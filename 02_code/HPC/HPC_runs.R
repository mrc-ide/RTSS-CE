# HPC set-up -------------------------------------------------------------------
library(didehpc)
setwd('Q:/GF-RTSS-CE')

options(didehpc.cluster = "fi--didemrchnb",
        didehpc.username = "htopazia")

source('./02_code/HPC/functions.R') # reference functions

# transfer the new malariasimulation folder manually to contexts or delete and reinstall using conan
# remotes::install_github('mrc-ide/malariasimulation@dev', force=T)

src <- conan::conan_sources("github::mrc-ide/malariasimulation@dev")

ctx <- context::context_save(path = "Q:/contexts",
                             sources = c('./02_code/HPC/functions.R'),
                             packages = c("dplyr", "malariasimulation"),
                             package_sources = src)

share <- didehpc::path_mapping('Home drive', "Q:", '//fi--san03.dide.ic.ac.uk/homes/htopazia', "M:")
config <- didehpc::didehpc_config(shares = share,
                                  use_rrq = FALSE,
                                  cores = 1,
                                  cluster = "fi--didemrchnb",
                                  parallel = FALSE)

# obj <- didehpc::queue_didehpc(ctx, config = config, provision = "upgrade")
obj <- didehpc::queue_didehpc(ctx, config = config)



# Set up your job --------------------------------------------------------------
library(tidyverse)
year <- 365

# population
population <- 200000

# seasonal profiles: c(g0, g[1], g[2], g[3], h[1], h[2], h[3])
# drawn from mlgts: https://github.com/mrc-ide/mlgts/tree/master/data
# g0 = a0, a = g, b = h
pfpr <- c(.10, .20, .40) # start at 10% (min rec for RTSS use)

  # FIRST
  seas_name <- 'highly seasonal'
  seasonality <- list(c(0.284596,-0.317878,-0.0017527,0.116455,-0.331361,0.293128,-0.0617547))
  s1 <- crossing(seasonality, seas_name, pfpr)

  # SECOND
  seas_name <- 'seasonal'
  seasonality <- list(c(0.285505,-0.325352,-0.0109352,0.0779865,-0.132815,0.104675,-0.013919))
  s2 <- crossing(seasonality, seas_name, pfpr)

  # THIRD
  seas_name <- 'perennial'
  seasonality <- list(c(0.2852770,-0.0248801,-0.0529426,-0.0168910,-0.0216681,-0.0242904,-0.0073646))
  s3 <- crossing(seasonality, seas_name, pfpr)

stable <- rbind(s1, s2, s3)

# vectors: list(arab_params, fun_params, gamb_params)
speciesprop <- data.frame(speciesprop = rbind(list(c(0.25, 0.25, 0.5))),
                          row.names = NULL)

# run time
warmup <- 6*year # needs to be multiple of 3 so that ITN dist. will line up with first timestep
sim_length <- 15*year

# interventions
ITN <- c('pyr', 'pbo')
ITNuse <- c(0,0.25, 0.50, 0.75)
ITNboost <- c(0,1)
resistance <- c(0, 0.4, 0.8)
IRS <-  c(0)
treatment <- c(0.30, 0.45, 0.60)
SMC <- c(0, 0.85)
RTSS <- c("none", "EPI", "SV") # leave out hybrid for now
RTSScov <- c(0, 0.85) # MVIP: dose 1 range 74-93%, dose 3 63-82%, dose 4 42-46%; first half of 2021
fifth <- c(0) # only one booster for now

interventions <-
  crossing(ITN, ITNuse, ITNboost, resistance, IRS, treatment, SMC, RTSS, RTSScov, fifth)

name <- "general"

# create combination of all runs and remove non-applicable scenarios
combo <- crossing(population, stable, warmup, sim_length, speciesprop, interventions) %>%
  mutate(name = paste0(name, "_", row_number())) %>%
  filter(!(RTSS == "none" & RTSScov > 0)) %>% # cannot set RTSS coverage when there is no RTSS
  filter(!(RTSScov == 0 & (RTSS == "EPI" | RTSS == 'SV' | RTSS == "hybrid"))) %>% # cannot set 0% coverage if RTSS is implemented
  filter(!(fifth == 1 & (RTSS == "EPI" | RTSS == 'none'))) %>% # no administration of fifth doses with EPI or no RTSS
  filter(!(seas_name == "perennial" & (RTSS == 'SV' | RTSS == "hybrid"))) %>% # do not introduce seasonal vaccination in perennial settings
  filter(!(SMC > 0 & seas_name == "perennial")) %>% # do not administer SMC in perennial settings
  filter(!(SMC == 0 & seas_name == "highly seasonal")) %>% # always introduce SMC in highly seasonal settings
  filter(!(ITNuse == 0 & resistance != 0)) %>% # do not introduce resistance when ITNuse==0
  filter(!(ITN == 'pbo' & resistance == 0)) %>% # only introduce PBO in areas that have resistance
  filter(!(ITN == 'pbo' & ITNboost == 1)) %>% # cannot have ITN boost + ITN PBO
  filter(!(ITN == 'pbo' & ITNuse == 0)) # cannot switch to PBO when ITNuse==0

# EIR / prev match from 'PfPR_EIR_match.R'
match <- readRDS("./02_code/HPC/EIRestimates.rds")

combo <- combo %>% left_join(match %>% select(-ITN)) %>%
  # put variables into the same order as function arguments
  select(population,        # simulation population
         seasonality,       # seasonal profile
         seas_name,         # name of seasonal profile
         starting_EIR,      # equilibrium EIR
         pfpr,              # corresponding PfPR
         warmup,            # warm-up period
         sim_length,        # length of simulation run
         speciesprop,       # proportion of each vector species
         ITN,               # ITN type - pyr, pbo
         ITNuse,            # ITN usage
         ITNboost,          # if ITN usage is boosted by 10%
         resistance,        # resistance level - none 0, med 0.4, high 0.8
         IRS,               # IRS coverage
         treatment,         # treatment coverage
         SMC,               # SMC coverage
         RTSS,              # RTS,S strategy
         RTSScov,           # RTS,S coverage
         fifth,             # status of 5th dose for SV or hybrid strategies
         name               # name of output file
  ) %>% as.data.frame()


# Run tasks --------------------------------------------------------------------
combo <- combo %>% mutate(f = paste0("./03_output/HPC/",combo$name,".rds")) %>%
  mutate(exist=case_when(file.exists(f) ~ 1, !file.exists(f) ~ 0)) %>%
  filter(exist==0) %>%
  select(-f, -exist)

t <- obj$enqueue_bulk(combo[1:2574,], runsimGF) # run 500 at a time [1:500,]
t$status()

beepr::beep(1)


# -#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# test run by hand
# x <- 4
# test <- runsimGF(population = combo[[x,1]],
#                  seasonality = combo[[x,2]],
#                  seas_name = combo[[x,3]],
#                  starting_EIR = combo[[x,4]],
#                  pfpr = combo[[x,5]],
#                  warmup = combo[[x,6]],
#                  sim_length = combo[[x,7]],
#                  speciesprop = combo[[x,8]],
#                  ITN = combo[[x,9]],
#                  ITNuse = combo[[x,10]],
#                  ITNboost = combo[[x,11]],
#                  resistance = combo[[x,12]],
#                  IRS = combo[[x,13]],
#                  treatment = combo[[x,14]],
#                  SMC = combo[[x,15]],
#                  RTSS = combo[[x,16]],
#                  RTSScov = combo[[x,17]],
#                  fifth = combo[[x,18]],
#                  name = combo[[x,19]])
# -#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
