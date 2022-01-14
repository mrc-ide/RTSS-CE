# Set-up ----------------------------------------------------------------------
library(didehpc)
setwd('M:/Hillary/GF-RTSS-CE')

options(didehpc.cluster = "fi--didemrchnb",
        didehpc.username = "htopazia")

# remotes::install_github('mrc-ide/malariasimulation@test/severe_demography', force=T)
src <- conan::conan_sources("github::mrc-ide/malariasimulation@dev")

ctx <- context::context_save(path = "M:/Hillary/contexts",
                             sources = c('./02_code/HPC/functions.R'),
                             packages = c("tidyverse", "malariasimulation"),
                             package_sources = src)

share <- didehpc::path_mapping("malaria", "M:", "//fi--didef3.dide.ic.ac.uk/malaria", "M:")
config <- didehpc::didehpc_config(shares = share,
                                  use_rrq = FALSE,
                                  cores = 1,
                                  cluster = "fi--didemrchnb",
                                  parallel = FALSE)

obj <- didehpc::queue_didehpc(ctx, config = config)


# Now set up your job ---------------------------------------------------------
library(tidyverse)
year <- 365

# population
population <- 10000

# seasonal profiles: c(g0, g[1], g[2], g[3], h[1], h[2], h[3])
# drawn from mlgts: https://github.com/mrc-ide/mlgts/tree/master/data
# g0 = a0, a = g, b = h
pfpr <- c(10,30,60) # start at 10% (min rec for RTSS use)

  # FIRST
  seas_name <- 'highly seasonal'
  seasonality <- list(c(0.284596,-0.317878,-0.0017527,0.116455,-0.331361,0.293128,-0.0617547))
  starting_EIR <- c(1.01,4.71,28.5)
  s1 <- crossing(seasonality, seas_name, starting_EIR) %>% cbind(pfpr)

  # SECOND
  seas_name <- 'seasonal'
  seasonality <- list(c(0.285505,-0.325352,-0.0109352,0.0779865,-0.132815,0.104675,-0.013919))
  starting_EIR <- c(1.11,5.41,34.0)
  s2 <- crossing(seasonality, seas_name, starting_EIR) %>% cbind(pfpr)

  # THIRD
  seas_name <- 'perennial'
  seasonality <- list(c(0.2852770,-0.0248801,-0.0529426,-0.0168910,-0.0216681,-0.0242904,-0.0073646))
  starting_EIR <- c(1.01,4.61,28.4)
  s3 <- crossing(seasonality, seas_name, starting_EIR) %>% cbind(pfpr)

stable <- rbind(s1, s2, s3)

# vectors: list(arab_params, fun_params, gamb_params)
speciesprop <- data.frame(speciesprop = rbind(list(c(0.25,0.25,0.5))),
                          row.names = NULL)

# run time
warmup <- 5*year
sim_length <- 20*year

# interventions
ITN <- c(0,0.25,0.50,0.75)
IRS <-  c(0)
treatment <- c(0.45)
SMC <- c(0,0.85)
RTSS <- c("none", "EPI", "SV") # leave out hybrid for now
RTSScov <- c(0,0.85) # MVIP: dose 1 range 74-93%, dose 3 63-82%, dose 4 42-46% first half of 2021
fifth <- c(0) # one booster for now

interventions <-
  crossing(ITN, IRS, treatment, SMC, RTSS, RTSScov, fifth) %>%
  filter(!(RTSS=="none" & RTSScov >0)) %>%                # can't have coverage when no RTSS
  filter(!(RTSScov==0 & (RTSS=="EPI" | RTSS=='SV' | RTSS=="hybrid"))) %>%  # can't have 0 coverage with RTSS
  filter(!(fifth==1 & (RTSS=="EPI" | RTSS=='none'))) %>%  # can't have fifth doses with EPI or none
  filter(!(SMC>0 & (seas_name=="seasonal" | seas_name=="perennial"))) # only administer SMC in highly seasonal settings

name <- "test"

# COMBINE into all combinations of runs
combo <- crossing(population, stable, warmup, sim_length, speciesprop, interventions) %>%
  mutate(name = paste0(name, "_", row_number())) %>%
  # put into correct order of function arguments
  select(population,        # simulation population
         seasonality,       # seasonal profile
         seas_name,         # name of seasonal profile
         starting_EIR,      # equilibrium EIR
         pfpr,              # corresponding PfPR
         warmup,            # warm-up period
         sim_length,        # length of simulation run
         speciesprop,       # proportion of each vector species
         ITN,               # ITN coverage
         IRS,               # IRS coverage
         treatment,         # treatment coverage
         SMC,               # SMC coverage
         RTSS,              # RTS,S strategy
         RTSScov,           # RTS,S coverage
         fifth,             # status of 5th dose for SV or hybrid strategies
         name               # name of output file
  ) %>% as.data.frame()


# Run tasks -------------------------------------------------------------------
t <- obj$enqueue_bulk(combo, runsimGF)
t$status()


# Still to add:
# insecticide resistance - ask Ellie - we want insecticide resistance to be an option. none / med / high
# new ITN types to try? - Ask Ellie if we would expect a big diff w/ PBO and IG2. If not a big difference, just choose one as a comparison example.

test <- runsimGF(combo[[1,1]],
                 combo[[1,2]],
                 combo[[1,3]],
                 combo[[1,4]],
                 combo[[1,5]],
                 combo[[1,6]],
                 combo[[1,7]],
                 combo[[1,8]],
                 combo[[1,9]],
                 combo[[1,10]],
                 combo[[1,11]],
                 combo[[1,12]],
                 combo[[1,13]],
                 combo[[1,14]],
                 combo[[1,15]],
                 combo[[1,16]])
