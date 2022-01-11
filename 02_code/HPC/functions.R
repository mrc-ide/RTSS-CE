# plotting --------------------------------------------------------------------
plotfun <- function(data,model){
  data %>%
    filter(model==model) %>%
    mutate(month=floor(timestep/(365/12))) %>%
    group_by(eir,month) %>%
    summarize(inc_month=sum(n_inc_0_36500, na.rm=T)/n_0_36500) %>%
    left_join(match, by="eir") %>%

    ggplot() +
    geom_line(aes(x=month, y=inc_month, color=as_factor(pfpr))) +
    labs(x="Timesteps (month)",
         y="Monthly \nclinical incidence \n(0-5 years)",
         color=expression(PfPR[2-10]),
         title=model) +
    scale_y_continuous(limits=c(0,1.7),breaks=c(0,0.5,1,1.5)) +
    theme_classic()
}

# no intervention -------------------------------------------------------------
runsim_none <- function(params, starting_EIR, warmup, sim_length, season){
  year <- 365
  month <- year/12

  params <- set_species(params, species = list(arab_params, fun_params, gamb_params),
                        proportions = c(0.25, 0.25, 0.5))
  params$phi_bednets <- c(0.9, 0.9, 0.89)
  params$phi_indoors <- c(0.96, 0.98, 0.97)

  params <- set_drugs(params, list(AL_params))
  params <- set_clinical_treatment(params, 1, c(1), c(0.45))
  params <- set_equilibrium(params, starting_EIR)
  output <- run_simulation(warmup + sim_length, params) %>%
    mutate(eir=starting_EIR,
           timestep=timestep-warmup,
           model=season) %>%
    filter(timestep > 0)

  saveRDS(output, paste0('./rds/HPC/none_',season,starting_EIR,'.rds'))
}

# EPI -------------------------------------------------------------------------
runsim_epi <- function(params, starting_EIR, warmup, sim_length, season){
  year <- 365
  month <- year/12

  params <- set_species(params, species = list(arab_params, fun_params, gamb_params),
                        proportions = c(0.25, 0.25, 0.5))
  params$phi_bednets <- c(0.9, 0.9, 0.89)
  params$phi_indoors <- c(0.96, 0.98, 0.97)

  params <- set_drugs(params, list(AL_params))
  params <- set_clinical_treatment(params, 1, c(1), c(0.45))

  params$rtss_doses <- round(c(0,1.5*month,3*month))
  boosters <- round(c(18*month))

  params <- set_rtss_epi(
    params,
    start = warmup,
    end = warmup + sim_length,
    coverage = 0.8,
    age = round(6*month),
    min_wait = 0,
    boosters = boosters,
    booster_coverage = rep(.80, 1),
    seasonal_boosters = FALSE)

  params <- set_equilibrium(params, starting_EIR)
  output <- run_simulation(warmup + sim_length, params) %>%
    mutate(eir=starting_EIR,
           timestep=timestep-warmup,
           model=season) %>%
    filter(timestep > 0)

  saveRDS(output, paste0('./rds/HPC/EPIalone_',season,starting_EIR,'.rds'))
}

# SV -------------------------------------------------------------------------
runsim_SV <- function(params, starting_EIR, warmup, sim_length, season, boosters, booster_coverage, rtss_cs_boost, name){
  year <- 365
  month <- year/12

  params <- set_species(params, species = list(arab_params, fun_params, gamb_params),
                        proportions = c(0.25, 0.25, 0.5))
  params$phi_bednets <- c(0.9, 0.9, 0.89)
  params$phi_indoors <- c(0.96, 0.98, 0.97)

  params <- set_drugs(params, list(AL_params))
  params <- set_clinical_treatment(params, 1, c(1), c(0.45))

  peak <- peak_season_offset(params)
  first <- ifelse(season=='high_seas', round(warmup+(peak-month*3.5),0), round(warmup+(peak-month*5.5),0))
  timesteps <- c(first, first+seq(1:14)*year)
  params$rtss_doses <- round(c(0,1*month,2*month))
  params$rtss_cs_boost <- c(rtss_cs_boost, 0.35)

  params <- set_mass_rtss(
    params,
    timesteps = timesteps,
    coverages = rep(0.8,length(timesteps)),
    min_ages = round(5*month),
    max_ages = round(17*month),
    min_wait = 0,
    boosters = as.vector(unlist(boosters)),
    booster_coverage = as.vector(unlist(booster_coverage)))

  params <- set_equilibrium(params, starting_EIR)
  output <- run_simulation(warmup + sim_length, params) %>%
    mutate(eir=starting_EIR,
           timestep=timestep-warmup,
           model=season) %>%
    filter(timestep > 0)

  saveRDS(output, paste0('./rds/HPC/',name,season,starting_EIR,'.rds'))
}

# SMC -------------------------------------------------------------------------
runsim_SMC <- function(params, starting_EIR, warmup, sim_length, season, shape, scale, name){
  year <- 365
  month <- year/12

  params <- set_species(params, species = list(arab_params, fun_params, gamb_params),
                        proportions = c(0.25, 0.25, 0.5))
  params$phi_bednets <- c(0.9, 0.9, 0.89)
  params$phi_indoors <- c(0.96, 0.98, 0.97)

  params <- set_drugs(params, list(AL_params, SP_AQ_params))
  params <- set_clinical_treatment(params, 1, c(1), c(0.45))

  peak <- peak_season_offset(params)
  first <- round(warmup+c(peak+c(-1.5,-0.5,0.5,1.5)*month),0)
  timesteps <- c(c(0,seq(1:14))*year + rep(first,15))

  params <- set_smc(
    params,
    drug = 2,
    timesteps = sort(timesteps),
    coverages = rep(0.75, length(timesteps)),
    min_age = round(0.25*year),
    max_age = round(5*year))

  params$drug_prophylaxis_shape <- c(11.3, shape)
  params$drug_prophylaxis_scale <- c(10.6, scale)

  params <- set_equilibrium(params, starting_EIR)
  output <- run_simulation(warmup + sim_length, params) %>%
    mutate(eir=starting_EIR,
           timestep=timestep-warmup,
           model=season) %>%
    filter(timestep > 0)

  saveRDS(output, paste0('./rds/HPC/',name,season,starting_EIR,'.rds'))
}

# SMC + EPI -------------------------------------------------------------------
runsim_SMCEPI <- function(params, starting_EIR, warmup, sim_length, season, shape, scale, name){
  year <- 365
  month <- year/12

  params <- set_species(params, species = list(arab_params, fun_params, gamb_params),
                        proportions = c(0.25, 0.25, 0.5))
  params$phi_bednets <- c(0.9, 0.9, 0.89)
  params$phi_indoors <- c(0.96, 0.98, 0.97)

  params <- set_drugs(params, list(AL_params, SP_AQ_params))
  params <- set_clinical_treatment(params, 1, c(1), c(0.45))

  peak <- peak_season_offset(params)
  first <- round(warmup+c(peak+c(-1.5,-0.5,0.5,1.5)*month),0)
  timesteps <- c(c(0,seq(1:14))*year + rep(first,15))

  params <- set_smc(
    params,
    drug = 2,
    timesteps = sort(timesteps),
    coverages = rep(0.75,length(timesteps)),
    min_age = round(0.25*year),
    max_age = round(5*year))

  params$drug_prophylaxis_shape <- c(11.3, shape)
  params$drug_prophylaxis_scale <- c(10.6, scale)

  params$rtss_doses <- round(c(0,1.5*month,3*month))
  boosters <- round(c(18*month))

  params <- set_rtss_epi(
    params,
    start = warmup,
    end = warmup + sim_length,
    coverage = 0.8,
    age = round(6*month),
    min_wait = 0,
    boosters = boosters,
    booster_coverage = rep(.80, 1),
    seasonal_boosters = FALSE)

  params <- set_equilibrium(params, starting_EIR)
  output <- run_simulation(warmup + sim_length, params) %>%
    mutate(eir=starting_EIR,
           timestep=timestep-warmup,
           model=season) %>%
    filter(timestep > 0)

  saveRDS(output, paste0('./rds/HPC/',name,season,starting_EIR,'.rds'))
}

# SMC + SV --------------------------------------------------------------------
runsim_SMCSV <- function(params, starting_EIR, warmup, sim_length, season, shape, scale, boosters, booster_coverage, rtss_cs_boost, name){
  year <- 365
  month <- year/12

  params <- set_species(params, species = list(arab_params, fun_params, gamb_params),
                        proportions = c(0.25, 0.25, 0.5))
  params$phi_bednets <- c(0.9, 0.9, 0.89)
  params$phi_indoors <- c(0.96, 0.98, 0.97)

  params <- set_drugs(params, list(AL_params, SP_AQ_params))
  params <- set_clinical_treatment(params, 1, c(1), c(0.45))

  params$drug_prophylaxis_shape <- c(11.3, shape)
  params$drug_prophylaxis_scale <- c(10.6, scale)

  peak <- peak_season_offset(params)
  first <- round(warmup+c(peak+c(-1.5,-0.5,0.5,1.5)*month),0)
  timesteps <- c(c(0,seq(1:14))*year + rep(first,15))

  params <- set_smc(
    params,
    drug = 2,
    timesteps = sort(timesteps),
    coverages = rep(0.75,length(timesteps)),
    min_age = round(0.25*year),
    max_age = round(5*year))

  peak <- peak_season_offset(params)
  first <- ifelse(season=='high_seas', round(warmup+(peak-month*3.5),0), round(warmup+(peak-month*5.5),0))
  timesteps <- c(first, first+seq(1:14)*year)
  params$rtss_doses <- round(c(0,1*month,2*month))
  params$rtss_cs_boost <- c(rtss_cs_boost, 0.35)

  params <- set_mass_rtss(
    params,
    timesteps = timesteps,
    coverages = rep(0.8,length(timesteps)),
    min_ages = round(5*month),
    max_ages = round(17*month),
    min_wait = 0,
    boosters = as.vector(unlist(boosters)),
    booster_coverage = as.vector(unlist(booster_coverage)))

  params$rtss_vmax <- if(name %in% c('SV4SMCsynergy_', 'SV5SMCsynergy_')) 0.911028 else .93
  params$rtss_alpha <- if(name %in% c('SV4SMCsynergy_', 'SV5SMCsynergy_')) 0.75303 else .74
  params$rtss_beta <- if(name %in% c('SV4SMCsynergy_', 'SV5SMCsynergy_')) 62.8525 else 99.4

  params <- set_equilibrium(params, starting_EIR)
  output <- run_simulation(warmup + sim_length, params) %>%
    mutate(eir=starting_EIR,
           timestep=timestep-warmup,
           model=season) %>%
    filter(timestep > 0)

  saveRDS(output, paste0('./rds/HPC/',name,season,starting_EIR,'.rds'))
}


# HYBRID --------------------------------------------------------------------
runsim_hybrid <- function(params, minwait, starting_EIR, warmup, sim_length, season, fifth){
  year <- 365
  month <- year/12

  params <- set_species(params, species = list(arab_params, fun_params, gamb_params),
                        proportions = c(0.25, 0.25, 0.5))
  params$phi_bednets <- c(0.9, 0.9, 0.89)
  params$phi_indoors <- c(0.96, 0.98, 0.97)

  params <- set_drugs(params, list(AL_params))
  params <- set_clinical_treatment(params, 1, c(1), c(0.45))

  params$rtss_doses <- round(c(0,1.5*month,3*month))

  peak <- peak_season_offset(params)
  first <- ifelse(season=='high_seas', round((peak-month*3.5),0), round((peak-month*5.5),0))
  boosters <- if(fifth==0) c(first+3*month) else ((first+3*month) + c(0, year))

  params <- set_rtss_epi(
    params,
    start = warmup,
    end = warmup + sim_length,
    coverage = 0.8,
    age = round(6*month),
    min_wait = minwait,
    boosters = boosters,
    booster_coverage = rep(.80, length(boosters)),
    seasonal_boosters = TRUE)

  params <- set_equilibrium(params, starting_EIR)
  output <- run_simulation(warmup + sim_length, params) %>%
    mutate(eir=starting_EIR,
           minwait=minwait,
           timestep=timestep-warmup,
           model=season,
           fifth=fifth) %>%
    filter(timestep > 0)

  saveRDS(output, paste0('./rds/HPC/hybrid_',minwait,'wait_',season,starting_EIR,'_',fifth,'.rds'))
}

# SEVERE age breakdown -------------------------------------------------------------
runsim_severe <- function(params, starting_EIR, warmup, sim_length, season){
  year <- 365
  month <- year/12

  params$severe_incidence_rendering_min_ages = seq(0,99,1)*year
  params$severe_incidence_rendering_max_ages = seq(1,100,1)*year

  params <- set_species(params, species = list(arab_params, fun_params, gamb_params),
                        proportions = c(0.25, 0.25, 0.5))
  params$phi_bednets <- c(0.9, 0.9, 0.89)
  params$phi_indoors <- c(0.96, 0.98, 0.97)

  params <- set_drugs(params, list(AL_params))
  params <- set_clinical_treatment(params, 1, c(1), c(0.45))
  params <- set_equilibrium(params, starting_EIR)
  output <- run_simulation(warmup + sim_length, params) %>%
    mutate(eir=starting_EIR,
           timestep=timestep-warmup,
           model=season) %>%
    filter(timestep > 0)

  saveRDS(output, paste0('./rds/severe/none_',season,starting_EIR,'.rds'))
}
