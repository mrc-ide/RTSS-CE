
# set interventions and run malariasimulation ----------------------------------

runsimGF <- function(population,        # simulation population
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
  params$clinical_incidence_rendering_min_ages = seq(0,100,5)*year
  params$clinical_incidence_rendering_max_ages = c(seq(5,100,5),200)*year
  params$severe_incidence_rendering_min_ages = seq(0,100,5)*year
  params$severe_incidence_rendering_max_ages = c(seq(5,100,5),200)*year

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
  # same value for all species
  if (ITN > 0) {
  params <- set_bednets(
    parameters = params,
    timesteps = seq(1, sim_length, year),
    coverages = rep(ITN, sim_length/year),
    retention = 3 * year,
    dn0 = matrix(c(rep(.39, sim_length/year),
                   rep(.39, sim_length/year),
                   rep(.39, sim_length/year)),
                 nrow=sim_length/year, ncol=3),
    rn = matrix(c(rep(.56, sim_length/year),
                  rep(.56, sim_length/year),
                  rep(.56, sim_length/year)),
                nrow=sim_length/year, ncol=3),
    rnm = matrix(c(rep(.24, sim_length/year),
                   rep(.24, sim_length/year),
                   rep(.24, sim_length/year)),
                 nrow=sim_length/year, ncol=3),
    gamman = rep(2.64 * 365, sim_length/year)
  )  }

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

  params <- set_clinical_treatment(
    parameters = params,
    drug = 1,
    timesteps = c(1),
    coverages = c(treatment)
  )  }

  # SMC ----------
  if (SMC > 0) {
    peak <- peak_season_offset(params)
    first <- round(warmup+c(peak+c(-1.5,-0.5,0.5,1.5)*month),0) # check this to see where these rounds fall on the seasonal settings - if it doesn't cover the whole season in the low seasonality setting, we may need to add another round. Or have one SMC 4 rounds and one SMC 5 rounds and see if the additional round makes a difference.
    firststeps <- sort(rep(first, sim_length/year))
    yearsteps <- rep(c(0, seq(year, sim_length - year, year)), sim_length/year)
    timesteps <- yearsteps + firststeps

    params <- set_drugs(
      parameters = params,
      list(AL_params, SP_AQ_params))

    params <- set_smc(
      parameters = params,
      drug = 2,
      timesteps = sort(timesteps),
      coverages = rep(SMC, length(timesteps)),
      min_age = round(0.5*year),
      max_age = round(5*year))
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

  # EIR equilibrium ----------
  params <- set_equilibrium(params, as.numeric(starting_EIR))

  # run simulation ----------
  output <- run_simulation(
    timesteps = warmup + sim_length,
    parameters = params,
    correlations = NULL) %>%
    # add vars to output
    mutate(EIR = starting_EIR,
           warmup = warmup,
           sim_length = sim_length,
           pfpr = pfpr,
           timestep = timestep - warmup,
           seasonality = seas_name,
           speciesprop = paste(speciesprop, sep = ',', collapse = ''),
           ITN = ITN,
           IRS = IRS,
           treatment = treatment,
           SMC = SMC,
           RTSS = RTSS,
           RTSScov = RTSScov,
           fifth = fifth) %>%
    ungroup() %>%
    filter(timestep > 0) %>% # remove warmup period

    # statistics by month
    mutate(year = ceiling(timestep/year),
           month = ceiling(timestep/month)) %>%
    # only necessary variables
    dplyr::select(EIR, warmup, sim_length, pfpr, month, year, seasonality, speciesprop,
                  ITN, IRS, treatment, SMC, RTSS, RTSScov, fifth, starts_with("n_inc_severe"),
                  starts_with("n_rtss"), starts_with("p_inc"), starts_with("n_inc"), starts_with("n_detect"),
                  starts_with("p_detect"), starts_with("n_"), -n_bitten, -n_treated, -n_infections) %>%
    # take means of populations and sums of cases by month
    group_by(EIR, warmup, sim_length, pfpr, month, year, seasonality, speciesprop,
             ITN, IRS, treatment, SMC, RTSS, RTSScov, fifth) %>%
    mutate_at(vars(n_0_1825:n_36500_73000), mean, na.rm = TRUE) %>%
    mutate_at(vars(n_inc_severe_0_1825:n_inc_36500_73000), sum, na.rm = TRUE) %>%
    distinct()

 # save output ----------
  saveRDS(output, paste0('./03_output/HPC/', name,'.rds'))

}
