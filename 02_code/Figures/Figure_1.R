# Figure 1
# Scenarios, interventions, and timings

# set-up
source("./02_code/Figures/data_and_libraries.R")

library(malariasimulation)

year <- 365
month <- year / 12
warmup <- 6 * year
sim_length <- 1 * year

# create combination of all runs and remove non-applicable scenarios
seas_name <- c('highly seasonal', 'seasonal', 'perennial')

# EIR / prev match from 'PfPR_EIR_match.R'
match <- readRDS("./02_code/HPC_median/EIRestimates.rds") %>%
  filter(treatment == 0.45 & ITN == 'pyr' & ITNuse == 0 & pfpr == 0.4)

combo <- tibble(seas_name) %>%
  left_join(match %>%
            ungroup() %>%
            select(-ITN, -treatment, -ITNuse, -pfpr))


runsimGF <- function(seas_name,         # name of seasonal profile
                     starting_EIR      # equilibrium EIR
                     ){

  # define seasonality
  # seasonal profiles: c(g0, g[1], g[2], g[3], h[1], h[2], h[3])
  # drawn from mlgts: https://github.com/mrc-ide/mlgts/tree/master/data
  if(seas_name == 'highly seasonal'){
    seasonality <- list(c(0.284596,-0.317878,-0.0017527,0.116455,-0.331361,0.293128,-0))
  }

  if(seas_name == 'seasonal'){
    seasonality <- list(c(0.285505,-0.325352,-0.0109352,0.0779865,-0.132815,0.104675,-0.013919))
  }

  if(seas_name == 'perennial'){
    seasonality <- list(c(0.2852770,-0.0248801,-0.0529426,-0.0168910,-0.0216681,-0.0242904,-0.0073646))
  }

  # get starting parameters ----------
  params <- get_parameters(list(
    human_population = 20000,
    model_seasonality = TRUE,
    # rainfall fourier parameters
    g0 = unlist(seasonality)[1],
    g = unlist(seasonality)[2:4],
    h = unlist(seasonality)[5:7],
    individual_mosquitoes = FALSE))

  # outcome definitions ----------
  # incidence for every 5 year age group
  params$clinical_incidence_rendering_min_ages = 0*year
  params$clinical_incidence_rendering_max_ages = 5*year

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
    proportions = c(0.25, 0.25, 0.5))

  # proportion of bites taken in bed for each species
  # find values in S.I. of 10.1038/s41467-018-07357-w Table 3
  params$phi_bednets <- c(0.9, 0.9, 0.89) # Hogan et al. 2020
  # proportion of bites taken indoors for each species
  params$phi_indoors <- c(0.96, 0.98, 0.97) # Hogan et al. 2020

  # ITNs ----------
  # find values in S.I. of 10.1038/s41467-018-07357-w
  # same value for all species
  itn_timesteps <- 1

  # treatment ----------
  params <- set_drugs(
    parameters = params,
    list(AL_params, SP_AQ_params))

  params$drug_prophylaxis_scale <- c(10.6, 39.34)
  params$drug_prophylaxis_shape <- c(11.3, 3.40)

  params <- set_clinical_treatment(
    parameters = params,
    drug = 1,
    timesteps = c(1),
    coverages = 0.45)

  # SMC ----------
  smc_timesteps <- 0

  if (seas_name == 'seasonal') {
    peak <- peak_season_offset(params)
    first <- round(c(peak+c(-2,-1,0,1,2)*month),0)
    smc_timesteps <- sort(rep(first, sim_length/year))
  }

  if (seas_name == 'highly seasonal') {
    peak <- peak_season_offset(params)
    first <- round(c(peak+c(-1,0,1,2)*month),0)
    smc_timesteps <- sort(rep(first, sim_length/year))

    params <- set_drugs(
      parameters = params,
      list(AL_params, SP_AQ_params))

    params$drug_prophylaxis_scale <- c(10.6, 39.34)
    params$drug_prophylaxis_shape <- c(11.3, 3.40)

    params <- set_smc(
      parameters = params,
      drug = 2,
      timesteps = smc_timesteps,
      coverages = rep(0.85, length(smc_timesteps)),
      min_age = round(0.25*year),
      max_age = round(5*year))
  }

  if (seas_name == 'perennial') {
    smc_timesteps <- 0
  }

  # SV ----------
  rtss_mass_timesteps <- 0
  peak <- peak_season_offset(params)
  rtss_mass_timesteps <- round((peak-month*3.5) ,0)

  # EIR equilibrium ----------
  params <- set_equilibrium(params, as.numeric(starting_EIR))

  # run simulation ----------
  set.seed(123)

  output <- run_simulation(
    timesteps = warmup + sim_length,
    parameters = params)

  results <- output %>%

    # add vars to output
    mutate(timestep = timestep - warmup,
           seasonality = seas_name) %>%
    filter(timestep > 0) %>% # remove warmup period

    # statistics by month
    mutate(year = ceiling(timestep/year),
           month = ceiling(timestep/month),
           CI_0_5 = p_inc_clinical_0_1825 / n_0_1825) %>%

    group_by(seasonality, month) %>%
    summarize(CI_0_5 = sum(CI_0_5)) %>%

    # add intervention timesteps
    mutate(smc_timesteps = list(smc_timesteps),
           itn_timesteps = itn_timesteps,
           rtss_mass_timesteps = rtss_mass_timesteps)

  return(results)

}

# run function
output <- map2_dfr(combo$seas_name, combo$starting_EIR, runsimGF)

# copy dataset twice for pre- and post-intervention
output <- rbind(output, output %>% mutate(month = month-13))

# pull out intervention timings
SMCtime <- output %>% select(smc_timesteps, seasonality) %>%
  group_by(seasonality) %>%
  filter(seasonality != 'perennial') %>%
  mutate(t1 = unlist(smc_timesteps)[[1]],
         t2 = unlist(smc_timesteps)[[2]],
         t3 = unlist(smc_timesteps)[[3]],
         t4 = unlist(smc_timesteps)[[4]],
         t5 = unlist(smc_timesteps)[[5]]) %>%
  distinct() %>%
  pivot_longer(cols = t1:t5, names_to = "time", values_to = "month") %>%
  mutate(month = month / (365/12) + 1, intervention='SMC') %>% select(-smc_timesteps, -time)

SMCtime <- rbind(SMCtime, SMCtime %>% filter(seasonality == 'highly seasonal') %>% mutate(month = month-13)) %>%
  filter(!is.na(month))

RTSStime <- output %>% select(rtss_mass_timesteps, seasonality) %>%
  group_by(seasonality) %>%
  mutate(month = unlist(rtss_mass_timesteps)[[1]]) %>%
  distinct() %>%
  mutate(month = month / (365/12) + 1, intervention='RTS,S seasonal dose 3',
         month = ifelse(seasonality == 'perennial', 0, month)) %>% select(-rtss_mass_timesteps)

ITNtime <- output %>% select(bednet_timesteps, seasonality) %>%
  group_by(seasonality) %>%
  mutate(month = bednet_timesteps) %>%
  distinct() %>%
  mutate(month = month / (365/12) + 1, intervention='ITN *') %>% select(-bednet_timesteps)

ITNtime <- rbind(ITNtime, ITNtime %>% mutate(month = month-13))

interventions <- rbind(SMCtime, RTSStime, ITNtime)

my_text <- data_frame(seasonality = 'perennial',
                      lab = c('pre-intervention', 'post-intervention'),
                      x = c(-6, 7),
                      y = c(0.25,0.25))

# plot
ggplot(data = output) + # %>% filter(seasonality != 'perennial')
  geom_line(aes(x = month, y = CI_0_5), alpha = 0.8) +
  geom_rect(data = interventions %>% filter(intervention == 'RTS,S seasonal dose 3'), aes(xmin=1, xmax=12, ymin=0.01, ymax=0.03, fill = 'RTSS'), alpha = 0.1) +
  geom_rect(data = interventions %>% filter(seasonality=='highly seasonal' & intervention=='SMC' & month<0),
            aes(xmin = min(month, na.rm = T), xmax = max(month, na.rm = T)+1, ymin = 0, ymax = 0.3, fill = intervention), lty = 2, alpha = 0.02) +
  geom_rect(data = interventions %>% filter(seasonality=='highly seasonal' & intervention=='SMC' & month>0),
            aes(xmin = min(month, na.rm = T), xmax = max(month, na.rm = T)+1, ymin = 0, ymax = 0.3, fill = intervention), lty = 2, alpha = 0.02) +
  geom_rect(data = interventions %>% filter(seasonality=='seasonal' & intervention=='SMC'),
            aes(xmin = min(month, na.rm = T), xmax = max(month, na.rm = T)+1, ymin = 0, ymax = 0.3, fill = intervention), lty = 2, alpha = 0.02) +
  geom_vline(data = interventions, aes(xintercept = month, color = intervention), lty = 2) +
  geom_rect(aes(xmin = -0.9, xmax = 0.9, ymin = -1, ymax = 0.3), fill = 'white') +
  geom_vline(aes(xintercept=0)) +
  geom_text(data = my_text, aes(x = x,  y = y, label = lab)) +
  labs(x = 'month', y = 'monthly clinical incidence, 0-5 years', color = '', fill = '') +
  scale_x_continuous(breaks = seq(-12, 12, 1)) +
  facet_grid(factor(seasonality,
                    levels = c('perennial', 'seasonal', 'highly seasonal')) ~ .) +
  coord_cartesian(xlim = c(-12, 12), ylim = c(0, 0.3)) +
  scale_fill_manual(values = c('#088BBE', '#F6A1A5'),
                    labels = c('RTS,S age-based','SMC coverage')) +
  scale_color_manual(values = c('#1BB6AF','#088BBE','#F6A1A5')) +
  theme_classic()

ggsave('./03_output/figure1.pdf', width = 8, height = 3)


