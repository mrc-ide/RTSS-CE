# Specify HPC options ----------------------------------------------------------
library(didehpc)
setwd('M:/Hillary/GF-RTSS-CE')

options(didehpc.cluster = "fi--didemrchnb",
        didehpc.username = "htopazia")

# transfer the new malariasimulation folder manually to contexts or delete and reinstall using conan
# remotes::install_github('mrc-ide/malariasimulation@dev', force=T)
src <- conan::conan_sources(c("github::mrc-ide/malariasimulation@dev", "github::mrc-ide/cali"))

ctx <- context::context_save(path = "Q:/contexts",
                             sources = c('./02_code/HPC_draws/functions_draws.R'),
                             packages = c("dplyr", "malariasimulation", "cali"),
                             package_sources = src)

share <- didehpc::path_mapping("malaria", "M:", "//fi--didef3.dide.ic.ac.uk/malaria", "M:")
config <- didehpc::didehpc_config(shares = share,
                                  use_rrq = FALSE,
                                  cores = 1,
                                  cluster = "fi--didemrchnb", # fi--dideclusthn OR fi--didemrchnb
                                  parallel = FALSE)

obj <- didehpc::queue_didehpc(ctx, config = config) # check for latest v. of packages


# Set up your job --------------------------------------------------------------
# make all combinations of baseline scenarios
library(tidyverse)
library(malariasimulation)
year <- 365
month <- year/12

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

# baseline scenarios
ITN <- c('pyr')
ITNuse <- c(0, 0.25, 0.30, 0.50, 0.60, 0.75)
resistance <- c(0, 0.4, 0.8)
treatment <- c(0.30, 0.45, 0.60)
SMC <- c(0, 0.85)

interventions <-
  crossing(ITN, ITNuse, resistance, treatment, SMC)


# create combination of all runs and remove non-applicable scenarios
baseline <- crossing(population, stable, warmup, sim_length, speciesprop, interventions) %>%
  filter(
    (SMC == 0 & seas_name %in% c("perennial", "seasonal")) |
      (SMC > 0 & seas_name == 'highly seasonal')) %>% # SMC only in highly seasonal settings
  filter(!(ITNuse == 0 & resistance != 0)) %>% # do not introduce resistance when ITNuse==0
  filter(!(ITNuse %in% c(0.30, 0.60) & treatment %in% c(0.30, 0.60))) %>% # no 30% or 60% treatment in case studies
  filter(!(ITNuse %in% c(0.30, 0.60) & seas_name == 'seasonal')) # no SMC in case studies


dat1 <- baseline %>% filter(ITNuse %in% c(0, 0.25, 0.50, 0.75))
dat2 <- baseline %>% filter(ITNuse %in% c(0.30, 0.60))

baseline <- rbind(dat1, dat2)

# function to get parameters for baseline scenarios
getparams_baseline <- function(x){

  data <- baseline[x,]

  population = data$population
  seasonality = data$seasonality
  seas_name = data$seas_name
  pfpr = data$pfpr
  warmup = data$warmup
  sim_length = data$sim_length
  speciesprop = data$speciesprop
  ITN = data$ITN
  ITNuse = data$ITNuse
  resistance = data$resistance
  treatment = data$treatment
  SMC = data$SMC

  year <- 365
  month <- year/12

  # get starting parameters ----------
  params <- get_parameters(list(
    human_population = population,
    model_seasonality = TRUE,
    # rainfall fourier parameters
    g0 = unlist(seasonality)[1],
    g = unlist(seasonality)[2:4],
    h = unlist(seasonality)[5:7],
    individual_mosquitoes = FALSE))

  # outcome definitions ----------
  params$prevalence_rendering_min_ages = 2 * year
  params$prevalence_rendering_max_ages = 10 * year

  # demography ----------
  flat_demog <- read.table('./01_data/Flat_demog.txt') # from mlgts
  ages <- round(flat_demog$V3 * year) # top of age bracket
  deathrates <- flat_demog$V5 / 365 # age-specific death rates

  params <- set_demography(
    params,
    agegroups = ages,
    timesteps = 1,
    deathrates = matrix(deathrates, nrow = 1),
    birthrates = find_birthrates(population, ages, deathrates)
  )

  # vectors ----------
  params <- set_species(
    parameters = params,
    species = list(arab_params, fun_params, gamb_params),
    proportions = unlist(speciesprop))

  # proportion of bites taken in bed for each species
  # find values in S.I. of 10.1038/s41467-018-07357-w Table 3
  params$phi_bednets <- c(0.9, 0.9, 0.89) # Hogan et al. 2020
  # proportion of bites taken indoors for each species
  params$phi_indoors <- c(0.96, 0.98, 0.97) # Hogan et al. 2020

  # ITNs ----------
  # find values in S.I. of 10.1038/s41467-018-07357-w
  # or in Table S1.3 of Ellie's 2021 paper
  # same value for all species
  bednet_timesteps <- c(0)

  # no resistance
  dn0_1 <- 0.387 # pyr, 0 resistance

  # resistance
  dn0_2 <- case_when(ITN=='pyr' & resistance==0 ~ 0.387,
                     ITN=='pyr' & resistance==0.4 ~ 0.352,
                     ITN=='pyr' & resistance==0.8 ~ 0.270,
                     ITN=='pbo' & resistance==0 ~ 0.517,
                     ITN=='pbo' & resistance==0.4 ~ 0.494,
                     ITN=='pbo' & resistance==0.8 ~ 0.419)
  # no resistance
  rn_1 <- 0.563 # pyr, 0 resistance

  # resistance
  rn_2 <- case_when(ITN=='pyr' & resistance==0 ~ 0.563,
                    ITN=='pyr' & resistance==0.4 ~ 0.568,
                    ITN=='pyr' & resistance==0.8 ~ 0.626,
                    ITN=='pbo' & resistance==0 ~ 0.474,
                    ITN=='pbo' & resistance==0.4 ~ 0.493,
                    ITN=='pbo' & resistance==0.8 ~ 0.525)

  # no resistance
  gamman_1 <- 2.64 # pyr, 0 resistance

  # resistance
  gamman_2 <- case_when(ITN=='pyr' & resistance==0 ~ 2.64,
                        ITN=='pyr' & resistance==0.4 ~ 2.226,
                        ITN=='pyr' & resistance==0.8 ~ 1.616,
                        ITN=='pbo' & resistance==0 ~ 2.64,
                        ITN=='pbo' & resistance==0.4 ~ 2.160,
                        ITN=='pbo' & resistance==0.8 ~ 1.311)

  ITNuse1 = ITNuse
  ITNuse2 = ITNuse1

  npre <- ceiling(warmup/(3*year))      # number of distributions during warmup
  npost <- ceiling(sim_length/(3*year)) # number of distributions during sim_length

  params <- set_bednets(
    parameters = params,
    timesteps = c(seq(1, (warmup), 3*year),  # baseline coverage starts at 1
                  seq(warmup + 1, (warmup + sim_length), 3*year)), # intervention coverage starts at sim_length

    coverages = c(rep(ITNuse1, npre),         # set baseline coverage
                  rep(ITNuse2, npost)),    # set intervention coverage
    retention = 3 * year,
    dn0 = matrix(c(rep(dn0_1, npre), rep(dn0_2, npost),
                   rep(dn0_1, npre), rep(dn0_2, npost),
                   rep(dn0_1, npre), rep(dn0_2, npost)),
                 nrow=npre + npost, ncol=3),
    rn = matrix(c(rep(rn_1, npre), rep(rn_2, npost),
                  rep(rn_1, npre), rep(rn_2, npost),
                  rep(rn_1, npre), rep(rn_2, npost)),
                nrow=npre + npost, ncol=3),
    rnm = matrix(c(rep(.24, npre + npost),
                   rep(.24, npre + npost),
                   rep(.24, npre + npost)),
                 nrow=npre + npost, ncol=3),
    gamman = c(rep(gamman_1 * 365, npre), rep(gamman_2 * 365, npost))
  )

  bednet_timesteps <- params$bednet_timesteps - warmup


  # treatment ----------
  if (treatment > 0) {
    params <- set_drugs(
      parameters = params,
      list(AL_params, SP_AQ_params))

    params$drug_prophylaxis_scale <- c(10.6, 39.34)
    params$drug_prophylaxis_shape <- c(11.3, 3.40)

    params <- set_clinical_treatment(
      parameters = params,
      drug = 1,
      timesteps = c(1),
      coverages = c(treatment)
    )  }

  # SMC ----------
  smc_timesteps <- 0

  if (SMC > 0 & seas_name == 'highly seasonal') {
    peak <- peak_season_offset(params)
    first <- round(c(peak+c(-1,0,1,2)*month),0)
    firststeps <- sort(rep(first, (warmup+sim_length)/year))
    yearsteps <- rep(c(0, seq(year, (warmup+sim_length) - year, year)), length(first))
    timesteps <- yearsteps + firststeps

    params <- set_drugs(
      parameters = params,
      list(AL_params, SP_AQ_params))

    params$drug_prophylaxis_scale <- c(10.6, 39.34)
    params$drug_prophylaxis_shape <- c(11.3, 3.40)

    params <- set_smc(
      parameters = params,
      drug = 2,
      timesteps = sort(timesteps),
      coverages = rep(SMC, length(timesteps)),
      min_age = round(0.25*year),
      max_age = round(5*year))

    smc_timesteps <- params$smc_timesteps - warmup
  }

  parameters = data.frame(params = c(0), pfpr = c(0))
  parameters$params <- list(params)
  parameters$pfpr <- pfpr
  parameters$ID <- paste(data$pfpr, data$seas_name, data$ITNuse,
                         data$resistance, data$treatment, sep="_")
  parameters

}

# run intervention combinations through parameter function to assign draws
out <- map_dfr(seq(1, nrow(baseline), 1), getparams_baseline)

# save
saveRDS(out, './02_code/HPC_draws/baselinescenarios.rds')



# Run tasks -------------------------------------------------------------------
x = c(1:nrow(out)) # 306 baseline scenarios
y = c(1:50) # 50 parameter draws

# define all combinations of scenarios and draws

index <- crossing(x,y)
index$x = index$x

# remove ones that have already been run
index <- index %>% mutate(f = paste0("./03_output/PR_EIR/PRmatch_draws_", index$x, '_', index$y, ".rds")) %>%
  mutate(exist=case_when(file.exists(f) ~ 1, !file.exists(f) ~ 0)) %>%
  filter(exist==0) %>%
  select(-f, -exist)

# run a test with the first scenario
# index <- index %>% filter(x == 271)

# submit all remaining tasks
# t <- obj$enqueue_bulk(index, PRmatch)

# submit jobs, 100 as a time
sjob <- function(x, y){

  t <- obj$enqueue_bulk(index[x:y,], PRmatch)
  1

}

map2_dfr(seq(0, nrow(index) - 100, 100),
         seq(99, nrow(index), 100),
         sjob)




# Results ----------------------------------------------------------------------
# read in results
files <- list.files(path = "./03_output/PR_EIR", pattern = "PRmatch_draws_", full.names = TRUE)
dat_list <- lapply(files, function (x) readRDS(x))

# concatenate
match <-  do.call("rbind", dat_list) %>% as_tibble()

summary(match$starting_EIR)


# save EIR estimates
saveRDS(match, "./02_code/HPC_draws/EIRestimates.rds")



