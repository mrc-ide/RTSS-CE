
# Set interventions and run malariasimulation ----------------------------------

runsimGF <- function(population,        # simulation population
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
                     resistance,        # resistance level - 0, 0.4, 0.8
                     IRS,               # IRS coverage
                     treatment,         # treatment coverage
                     SMC,               # SMC coverage
                     RTSS,              # RTS,S strategy
                     RTSScov,           # RTS,S coverage
                     fifth,             # status of 5th dose for SV or hybrid
                     name               # name of output file
                     ){

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
    # incidence for every 5 year age group
  params$clinical_incidence_rendering_min_ages = c(0, 0.25, seq(5,100,5))*year
  params$clinical_incidence_rendering_max_ages = c(0.25, seq(5,100,5),200)*year
  params$severe_incidence_rendering_min_ages = c(0, 0.25, seq(5,100,5))*year
  params$severe_incidence_rendering_max_ages = c(0.25, seq(5,100,5),200)*year
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
  ITNuse2 = ITNuse1 + .10   # ITN boost by 10%

  if (ITNboost == 0) {      # if ITNs are not boosted, keep ITN use constant
    ITNuse2 = ITNuse1
  }

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


  # IRS ----------
  # find values in S.I. of 10.1038/s41467-018-07357-w Table 3
  if (IRS > 0) {
  params <- set_spraying(
    parameters = params,
    timesteps = seq(1, sim_length, year),
    coverages = rep(IRS, sim_length/year),
    ls_theta = matrix(c(rep(2.025, sim_length/year),
                        rep(2.025, sim_length/year),
                        rep(2.025, sim_length/year)),
                     nrow=sim_length/year, ncol=3),
    ls_gamma = matrix(c(rep(-0.009, sim_length/year),
                        rep(-0.009, sim_length/year),
                        rep(-0.009, sim_length/year)),
                      nrow=sim_length/year, ncol=3),
    ks_theta = matrix(c(rep(-2.222, sim_length/year),
                        rep(-2.222, sim_length/year),
                        rep(-2.222, sim_length/year)),
                      nrow=sim_length/year, ncol=3),
    ks_gamma = matrix(c(rep(0.008, sim_length/year),
                        rep(0.008, sim_length/year),
                        rep(0.008, sim_length/year)),
                      nrow=sim_length/year, ncol=3),
    ms_theta = matrix(c(rep(-1.232, sim_length/year),
                        rep(-1.232, sim_length/year),
                        rep(-1.232, sim_length/year)),
                      nrow=sim_length/year, ncol=3),
    ms_gamma = matrix(c(rep(-0.009, sim_length/year),
                        rep(-0.009, sim_length/year),
                        rep(-0.009, sim_length/year)),
                      nrow=sim_length/year, ncol=3)
  )  }

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

  if (SMC > 0 & seas_name == 'seasonal') {
    peak <- peak_season_offset(params)
    first <- round(warmup+c(peak+c(-2,-1,0,1,2)*month),0)
    firststeps <- sort(rep(first, sim_length/year))
    yearsteps <- rep(c(0, seq(year, sim_length - year, year)), length(first))
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



  # EPI ----------
  if (RTSS == "EPI") {
  params$rtss_doses <- round(c(0,1.5*month,3*month))
  boosters <- round(c(18*month))

  params <- set_rtss_epi(
    parameters = params,
    start = warmup,
    end = warmup + sim_length,
    coverage = RTSScov,
    age = round(6*month),
    min_wait = 0,
    boosters = boosters,
    booster_coverage = rep(.80, 1),
    seasonal_boosters = FALSE
  )  }

  # SV ----------
  rtss_mass_timesteps <- 0

  if (RTSS == "SV") {
  peak <- peak_season_offset(params)
  first <- round(warmup+(peak-month*3.5),0)
  timesteps <- c(first, first+seq(year, sim_length, year))
  params$rtss_doses <- round(c(0,1*month,2*month))

  boosters <- if(fifth==0) round(c(12*month+2*month)) else round(c(12*month+2*month, 24*month+2*month))

  params <- set_mass_rtss(
    parameters = params,
    timesteps = timesteps,
    coverages = rep(RTSScov,length(timesteps)),
    min_ages = round(5*month),
    max_ages = round(17*month),
    min_wait = 0,
    boosters = boosters,
    booster_coverage = rep(.80, length(boosters)))

  rtss_mass_timesteps <- params$rtss_mass_timesteps - warmup

  }

  # hybrid ----------
  if (RTSS == "hybrid") {
  params$rtss_doses <- round(c(0,1.5*month,3*month))

  peak <- peak_season_offset(params)
  first <- round(warmup+(peak-month*3.5),0)
  boosters <- if(fifth==0) round(c(first+3*month),0) else round(((first+3*month) + c(0, year)),0)

  params <- set_rtss_epi(
    parameters = params,
    start = warmup,
    end = warmup + sim_length,
    coverage = RTSScov,
    age = round(6*month),
    min_wait = 0,
    boosters = boosters,
    booster_coverage = rep(.80, length(boosters)),
    seasonal_boosters = TRUE)
  }

    # synergy ----------
  if (SMC > 0 & RTSS %in% c("EPI", "SV", "hybrid")) {

  params$rtss_beta <- 70.9
  params$rtss_alpha <- 0.868
  params$rtss_vmax <- 0.843
  params$rtss_cs_boost <- c(6.37008, 0.35)

  params$drug_prophylaxis_scale <- c(10.6, 45.76)
  params$drug_prophylaxis_shape <- c(11.3, 2.87)

  }

  # correlate interventions  ----------
  # correlations not working correctly with set.seed()
  # correlations <- get_correlation_parameters(params)

  # if (RTSScov == 0.77 & pfpr == 0.40) {
  #   correlations$inter_intervention_rho('rtss', 'bednets', 0.04)
  # }
  #
  # if (RTSScov == 0.72 & pfpr == 0.18) {
  #   correlations$inter_intervention_rho('rtss', 'bednets', 0.07)
  # }

  # EIR equilibrium ----------
  params <- set_equilibrium(params, as.numeric(starting_EIR))

  # run simulation ----------
  set.seed(123)

  output <- run_simulation(
    timesteps = warmup + sim_length,
    # correlations = correlations,
    parameters = params) %>%

    # add vars to output
    mutate(EIR = starting_EIR,
           warmup = warmup,
           sim_length = sim_length,
           population = population,
           pfpr = pfpr,
           timestep = timestep - warmup,
           seasonality = seas_name,
           speciesprop = paste(speciesprop, sep = ',', collapse = ''),
           ITN = ITN,
           ITNuse = ITNuse,
           ITNboost = ITNboost,
           resistance,
           IRS = IRS,
           treatment = treatment,
           SMC = SMC,
           RTSS = RTSS,
           RTSScov = RTSScov,
           fifth = fifth,
           bednet_timesteps = list(bednet_timesteps),
           smc_timesteps = list(smc_timesteps),
           rtss_mass_timesteps = list(rtss_mass_timesteps)) %>%
    ungroup() %>%
    filter(timestep > 0) %>% # remove warmup period

    # statistics by month
    mutate(year = ceiling(timestep/year),
           month = ceiling(timestep/month)) %>%

    # only necessary variables
    dplyr::select(EIR, warmup, sim_length, population, pfpr, month, year, seasonality, speciesprop,
                  ITN, ITNuse, ITNboost, resistance, IRS, treatment, SMC, RTSS, RTSScov, fifth,
                  starts_with("n_inc_severe"), starts_with("p_inc_severe"),
                  starts_with("n_rtss"),
                  starts_with("n_inc"), starts_with("p_inc"),
                  starts_with("n_detect"), starts_with("p_detect"),
                  starts_with("n_"), -n_bitten, n_treated, n_infections, bednet_timesteps,
                  smc_timesteps, rtss_mass_timesteps) %>%

    # take means of populations and sums of cases by month
    group_by(EIR, warmup, sim_length, population, pfpr, month, year, seasonality, speciesprop,
             ITN, ITNuse, ITNboost, resistance, IRS, treatment, SMC, RTSS, RTSScov, fifth,
             bednet_timesteps, smc_timesteps, rtss_mass_timesteps) %>%

    mutate_at(vars(n_0_91.25:n_36500_73000, n_730_3650,
                   n_detect_730_3650, p_detect_730_3650), mean, na.rm = TRUE) %>%
    mutate_at(vars(n_inc_severe_0_91.25:p_inc_clinical_36500_73000,
                   n_treated, n_infections), sum, na.rm = TRUE) %>%

    dplyr::select(n_0_91.25:n_36500_73000,
                  n_inc_severe_0_91.25:p_inc_clinical_36500_73000,
                  n_detect_730_3650, p_detect_730_3650,
                  n_730_3650,
                  n_treated, n_infections) %>%

    distinct()

 # save output ----------
  saveRDS(output, paste0('./03_output/HPC/', name,'.rds'))

}


# PfPR EIR match ---------------------------------------------------------------

PRmatch <- function(seasonality, seas_name, init_EIR, ITN, ITNuse, treatment, name){
  # parameters
  year <- 365
  month <- year/12
  human_population <- 100000

  params <- get_parameters(list(
    human_population = human_population,
    model_seasonality = TRUE,   # assign seasonality
    g0 = unlist(seasonality)[1],
    g = unlist(seasonality)[2:4],
    h = unlist(seasonality)[5:7],
    prevalence_rendering_min_ages = 2 * year,
    prevalence_rendering_max_ages = 10 * year,
    individual_mosquitoes = FALSE))

  flat_demog <- read.table('./01_data/Flat_demog.txt') # from mlgts
  ages <- round(flat_demog$V3 * year) # top of age bracket
  deathrates <- flat_demog$V5 / 365 # age-specific death rates
  params <- set_demography(
    params,
    agegroups = ages,
    timesteps = 1,
    deathrates = matrix(deathrates, nrow = 1),
    birthrates = find_birthrates(human_population, ages, deathrates)
  )

  params <- set_species(params, species = list(arab_params, fun_params, gamb_params),
                        proportions = c(0.25, 0.25, 0.5))
  params$phi_bednets <- c(0.9, 0.9, 0.89)
  params$phi_indoors <- c(0.96, 0.98, 0.97)

  params <- set_drugs(params, list(AL_params, SP_AQ_params))
  params <- set_clinical_treatment(params, 1, c(1), treatment)

  if (seas_name == 'highly seasonal') {
  peak <- peak_season_offset(params)
  first <- round(c(peak+c(-1,0,1,2)*month),0)
  firststeps <- sort(rep(first, (9*year)/year))
  yearsteps <- rep(c(0, seq(year, (9*year) - year, year)), length(first))
  timesteps <- yearsteps + firststeps

  params <- set_smc(
    parameters = params,
    drug = 2,
    timesteps = sort(timesteps),
    coverages = rep(.85, length(timesteps)),
    min_age = round(0.25*year),
    max_age = round(5*year))

  params$drug_prophylaxis_scale <- c(10.6, 39.34)
  params$drug_prophylaxis_shape <- c(11.3, 3.40)

  }

  # no resistance
  dn0_1 <- case_when(ITN=='pyr' ~ 0.387,
                     ITN=='pbo' ~ 0.517)

  # no resistance
  rn_1 <- case_when(ITN=='pyr' ~ 0.563,
                    ITN=='pbo' ~ 0.474)

  # no resistance
  gamman_1 <- case_when(ITN=='pyr'~ 2.64,
                        ITN=='pbo' ~ 2.64)

  params <- set_bednets(
    parameters = params,
    timesteps = seq(1, 9*year, 3*year),
    coverages = rep(ITNuse, 3),
    retention = 3 * year,
    dn0 = matrix(c(rep(dn0_1, 3),
                   rep(dn0_1, 3),
                   rep(dn0_1, 3)),
                 nrow=3, ncol=3),
    rn = matrix(c(rep(rn_1, 3),
                  rep(rn_1, 3),
                  rep(rn_1, 3)),
                nrow=3, ncol=3),
    rnm = matrix(c(rep(.24, 3),
                   rep(.24, 3),
                   rep(.24, 3)),
                 nrow=3, ncol=3),
    gamman = rep(gamman_1 * 365, 3))

  params <- set_equilibrium(params, as.numeric(init_EIR))

  set.seed(123)

  output <- run_simulation(
    timesteps = 9 * year,
    parameters = params,
    correlations = NULL)

  # output prev 2-10 values
  prev <-
      mean(
        output[
          output$timestep %in% seq(7 * year, 9 * year), # averaging a multiple of 3 for bednets
          'n_detect_730_3650'
        ] / output[
          output$timestep %in% seq(7 * year, 9 * year),
          'n_730_3650'
        ]
      )

  # create dataframe of initial EIR, output EIR, and prev 2-10 results
  EIR_prev <- cbind(seas_name, ITN, ITNuse, treatment, init_EIR, prev)

  saveRDS(EIR_prev, paste0('./03_output/PR_EIR/', name,'.rds'))

}

# PR match for DHS data, microscopy prev for 5-59 months
# changes from PRmatch include age group for prev calculations and smc coverage to 87%
PRmatch_u5 <- function(seasonality, seas_name, init_EIR, ITN, ITNuse, name){
  # parameters
  year <- 365
  month <- year/12
  human_population <- 100000

  params <- get_parameters(list(
    human_population = human_population,
    model_seasonality = TRUE,   # assign seasonality
    g0 = unlist(seasonality)[1],
    g = unlist(seasonality)[2:4],
    h = unlist(seasonality)[5:7],
    prevalence_rendering_min_ages = round(5 * month),
    prevalence_rendering_max_ages = 5 * year,
    individual_mosquitoes = FALSE))

  flat_demog <- read.table('./01_data/Flat_demog.txt') # from mlgts
  ages <- round(flat_demog$V3 * year) # top of age bracket
  deathrates <- flat_demog$V5 / 365 # age-specific death rates
  params <- set_demography(
    params,
    agegroups = ages,
    timesteps = 1,
    deathrates = matrix(deathrates, nrow = 1),
    birthrates = find_birthrates(human_population, ages, deathrates)
  )

  params <- set_species(params, species = list(arab_params, fun_params, gamb_params),
                        proportions = c(0.25, 0.25, 0.5))
  params$phi_bednets <- c(0.9, 0.9, 0.89)
  params$phi_indoors <- c(0.96, 0.98, 0.97)

  params <- set_drugs(params, list(AL_params, SP_AQ_params))
  params <- set_clinical_treatment(params, 1, c(1), c(0.45))

  if (seas_name == 'highly seasonal') {
    peak <- peak_season_offset(params)
    first <- round(c(peak+c(-1,0,1,2)*month),0)
    firststeps <- sort(rep(first, (9*year)/year))
    yearsteps <- rep(c(0, seq(year, (9*year) - year, year)), length(first))
    timesteps <- yearsteps + firststeps

    params <- set_smc(
      parameters = params,
      drug = 2,
      timesteps = sort(timesteps),
      coverages = rep(.87, length(timesteps)),
      min_age = round(0.25*year),
      max_age = round(5*year))

    params$drug_prophylaxis_scale <- c(10.6, 39.34)
    params$drug_prophylaxis_shape <- c(11.3, 3.40)

  }

  # no resistance
  dn0_1 <- case_when(ITN=='pyr' ~ 0.387,
                     ITN=='pbo' ~ 0.517)

  # no resistance
  rn_1 <- case_when(ITN=='pyr' ~ 0.563,
                    ITN=='pbo' ~ 0.474)

  # no resistance
  gamman_1 <- case_when(ITN=='pyr'~ 2.64,
                        ITN=='pbo' ~ 2.64)

  params <- set_bednets(
    parameters = params,
    timesteps = seq(1, 9*year, 3*year),
    coverages = rep(ITNuse, 3),
    retention = 3 * year,
    dn0 = matrix(c(rep(dn0_1, 3),
                   rep(dn0_1, 3),
                   rep(dn0_1, 3)),
                 nrow=3, ncol=3),
    rn = matrix(c(rep(rn_1, 3),
                  rep(rn_1, 3),
                  rep(rn_1, 3)),
                nrow=3, ncol=3),
    rnm = matrix(c(rep(.24, 3),
                   rep(.24, 3),
                   rep(.24, 3)),
                 nrow=3, ncol=3),
    gamman = rep(gamman_1 * 365, 3))

  params <- set_equilibrium(params, as.numeric(init_EIR))

  set.seed(123)
  output <- run_simulation(
    timesteps = 9 * year,
    parameters = params,
    correlations = NULL)

  # output prev 5-59 months values
  prev <-
    mean(
      output[
        output$timestep %in% seq(7 * year, 9 * year), # averaging a multiple of 3 for bednets
        'n_detect_152_1825'
      ] / output[
        output$timestep %in% seq(7 * year, 9 * year),
        'n_152_1825'
      ]
    )

  # create dataframe of initial EIR, output EIR, and prev 2-10 results
  EIR_prev <- cbind(seas_name, ITN, ITNuse, init_EIR, prev)

  saveRDS(EIR_prev, paste0('./03_output/PR_EIR/', name,'.rds'))

}
