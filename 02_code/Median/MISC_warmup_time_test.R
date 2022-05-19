# run malariasimulation to assess how long warmup time should be w/ ITNs
# run until PfPR stabilizes

library(malariasimulation)
library(tidyverse)

year <- 365
month <- year/12
human_population <- 50000
seasonality <- list(c(0.284596,-0.317878,-0.0017527,0.116455,-0.331361,0.293128,-0.0617547))
ITN = 'pyr'
ITNuse = 0.75
init_EIR = 50

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
params <- set_clinical_treatment(params, 1, c(1), c(0.45))


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

output <- run_simulation(
  timesteps = 9 * year,
  parameters = params,
  correlations = NULL)

beepr::beep(1) # notify when run is finished


# plot by day
ggplot(output) +
  geom_line(aes(x=timestep, y=n_detect_730_3650/n_730_3650))

# calculate prev over 3 year periods
output %>%
  mutate(yearblock = ceiling(timestep/(365*3))) %>%
  group_by(yearblock) %>%
  summarize(prev = mean(n_detect_730_3650/n_730_3650))

# appears to be stable between the second and third set of 3-year blocks
# 3 year warm-up is sufficient; run 6 to be safe
# tested EIR 5, EIR 50, and EIR 100
