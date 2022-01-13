# Set-up ----------------------------------------------------------------------
library(didehpc)
setwd('M:/Hillary/rtss_malariasimulation')

options(didehpc.cluster = "fi--didemrchnb",
        didehpc.username = "htopazia")

# remotes::install_github('mrc-ide/malariasimulation@test/severe_demography', force=T)
src <- conan::conan_sources("github::mrc-ide/malariasimulation@feat/demography")

ctx <- context::context_save(root = "context",
                             sources = c('./02_code/functions.R'),
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

# population
population <- 100000

# seasonal profiles: c(g0, g[1], g[2], g[3], h[1], h[2], h[3])
pfpr <- c(5,10,15)

  # FIRST
  seas_name <- 'highly seasonal'
  seasonality <- list(c(0.284596,-0.317878,-0.0017527,0.116455,-0.331361,0.293128,-0.0617547))
  starting_EIR <- c(1.1,2.1,3.1)
  s1 <- crossing(seasonality, seas_name, starting_EIR) %>% cbind(pfpr)

  # SECOND
  seas_name <- 'seasonal'
  seasonality <- list(c(0.285505,-0.325352,-0.0109352,0.0779865,-0.132815,0.104675,-0.013919))
  starting_EIR <- c(1.2,2.2,3.2)
  s2 <- crossing(seasonality, seas_name, starting_EIR) %>% cbind(pfpr)

stable <- rbind(s1,s2)

# vectors: list(arab_params, fun_params, gamb_params)
speciesprop <- data.frame(speciesprop = rbind(list(c(0.25,0.25,0.5)),
                                              list(c(0.3,0.3,0.4))),
                          row.names = NULL)

# run time
warmup <- 1*year
sim_length <- 5*year

# interventions
ITN <- c(0,0.15,0.30,0.50,0.55,0.60,0.65,0.70,0.75) # Winskill et al. 2017
IRS <-  c(0,0.25,0.50,0.75,0.90) # Winskill et al. 2017
treatment <- c(0.60) # Winskill et al. 2017
SMC <- c(0,0.25,0.50,0.75,0.90) # Winskill et al. 2017
RTSS <- c("none", "EPI", "SV", "hybrid")
RTSScov <- c(0,0.25,0.50,0.75,0.90) # Winskill et al. 2017
fifth <- c(0,1)

interventions <-
  crossing(ITN, IRS, treatment, SMC, RTSS, RTSScov, fifth) %>%
  filter(!(RTSS=="none" & RTSScov >0)) %>%            # can't have coverage when no RTSS
  filter(!(fifth==1 & (RTSS=="EPI" | RTSS=='none')))  # can't have fifth doses with EPI or none

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
         ITN,               # ITN status
         IRS,               # IRS status
         treatment,         # treatment status
         SMC,               # SMC status
         RTSS,              # RTS,S strategy
         RTSScov,           # RTS,S coverage
         fifth,             # status of 5th dose for SV or hybrid strategies
         name               # name of output file
  ) %>% as.data.frame()


# Run tasks -------------------------------------------------------------------
t <- obj$enqueue_bulk(combo, runsimGF)
t$status()

test <- runsimGF(combotest[[1]],
                 combotest[[2]],
                 combotest[[3]],
                 combotest[[4]],
                 combotest[[5]],
                 combotest[[6]],
                 combotest[[7]],
                 combotest[[8]],
                 combotest[[9]],
                 combotest[[10]],
                 combotest[[11]],
                 combotest[[12]],
                 combotest[[13]],
                 combotest[[14]],
                 combotest[[15]],
                 combotest[[16]])
combotest <- combo[1,]; combotest

