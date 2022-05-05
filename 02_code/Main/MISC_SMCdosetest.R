
# SMC - is adding another dose better?

library(tidyverse)
library(malariasimulation)

year <- 365
month <- year/12

match <- readRDS("./02_code/HPC/EIRestimates.rds") %>% filter(ITNuse==0.5 & seas_name!='perennial'); match

# FIRST
seas_name <- 'highly seasonal'
seasonality <- list(c(0.284596,-0.317878,-0.0017527,0.116455,-0.331361,0.293128,-0.0617547))
starting_EIR <- c(10.2, 30.3, 167)
s1 <- crossing(seasonality, seas_name, starting_EIR)

# SECOND
seas_name <- 'seasonal'
seasonality <- list(c(0.285505,-0.325352,-0.0109352,0.0779865,-0.132815,0.104675,-0.013919))
starting_EIR <- c(4.42, 11.4, 43.6)
s2 <- crossing(seasonality, seas_name, starting_EIR)

fifth <- c(0,1,2)

stable <- rbind(s1, s2)
combo <- crossing(stable, fifth)

# FUNCTION ---------------------------------------------------------------------
runsim <- function(seasonality, seas_name, starting_EIR, fifth){

# run time
warmup <- 5*year
sim_length <- 5*year
population <- 50000

# get starting parameters ----------
params <- get_parameters(list(
  human_population = population,
  model_seasonality = TRUE,
  g0 = unlist(seasonality)[1],
  g = unlist(seasonality)[2:4],
  h = unlist(seasonality)[5:7],
  individual_mosquitoes = FALSE))

# outcome definitions ----------
params$clinical_incidence_rendering_min_ages = c(0,0)*year
params$clinical_incidence_rendering_max_ages = c(5,100)*year
params$severe_incidence_rendering_min_ages = c(0,0)*year
params$severe_incidence_rendering_max_ages = c(5,100)*year

# vectors ----------
params <- set_species(
  parameters = params,
  species = list(arab_params, fun_params, gamb_params),
  proportions = c(0.25,0.25,0.5))

# proportion of bites taken in bed for each species
# find values in S.I. of 10.1038/s41467-018-07357-w Table 3
params$phi_bednets <- c(0.9, 0.9, 0.89) # Hogan et al. 2020
# proportion of bites taken indoors for each species
params$phi_indoors <- c(0.96, 0.98, 0.97) # Hogan et al. 2020

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
  coverages = c(.45)
)

# ITNs ----------
params <- set_bednets(
  parameters = params,
  timesteps = seq(1, warmup + sim_length, year),
  coverages = rep(0.50, (warmup + sim_length)/year),
  retention = 1 * year,
  dn0 = matrix(c(rep(0.387, (warmup + sim_length)/year),
                 rep(0.387, (warmup + sim_length)/year),
                 rep(0.387, (warmup + sim_length)/year)),
               nrow=(warmup + sim_length)/year, ncol=3),
  rn = matrix(c(rep(0.563, (warmup + sim_length)/year),
                rep(0.563, (warmup + sim_length)/year),
                rep(0.563, (warmup + sim_length)/year)),
              nrow=(warmup + sim_length)/year, ncol=3),
  rnm = matrix(c(rep(.24, (warmup + sim_length)/year),
                 rep(.24,(warmup + sim_length)/year),
                 rep(.24, (warmup + sim_length)/year)),
               nrow=(warmup + sim_length)/year, ncol=3),
  gamman = rep(2.64 * 365, (warmup + sim_length)/year))

# SMC ----------
if (fifth==0) {
  peak <- peak_season_offset(params)
  first <- round(warmup+c(peak+c(-2,-1,0,1)*month),0)
  firststeps <- sort(rep(first, sim_length/year))
  yearsteps <- rep(c(0, seq(year, sim_length - year, year)), length(first))
  timesteps <- yearsteps + firststeps

  params <- set_smc(
    parameters = params,
    drug = 2,
    timesteps = sort(timesteps),
    coverages = rep(0.85, length(timesteps)),
    min_age = round(0.25*year),
    max_age = round(5*year))
}

if (fifth==1) {
  peak <- peak_season_offset(params)
  first <- round(warmup+c(peak+c(-1,0,1,2)*month),0)
  firststeps <- sort(rep(first, sim_length/year))
  yearsteps <- rep(c(0, seq(year, sim_length - year, year)), length(first))
  timesteps <- yearsteps + firststeps

  params <- set_smc(
    parameters = params,
    drug = 2,
    timesteps = sort(timesteps),
    coverages = rep(0.85, length(timesteps)),
    min_age = round(0.25*year),
    max_age = round(5*year))
}

if (fifth==2) {
  peak <- peak_season_offset(params)
  first <- round(warmup+c(peak+c(-2,-1,0,1,2)*month),0)
  firststeps <- sort(rep(first, sim_length/year))
  yearsteps <- rep(c(0, seq(year, sim_length - year, year)), length(first))
  timesteps <- yearsteps + firststeps

  params <- set_smc(
    parameters = params,
    drug = 2,
    timesteps = sort(timesteps),
    coverages = rep(0.85, length(timesteps)),
    min_age = round(0.25*year),
    max_age = round(5*year))
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
         timestep = timestep - warmup,
         seasonality = seas_name,
         fifth = fifth) %>%
  ungroup() %>%
  filter(timestep > 0)
}


runsim2 <- function(x){
  runsim(combo[[x,1]], combo[[x,2]], combo[[x,3]], combo[[x,4]])
}

# run
output <- map_dfr(c(1:nrow(combo)), runsim2)

# save
saveRDS(output,"./03_output/smc_test.rds")


# ANALYSIS ---------------------------------------------------------------------
output <- readRDS("./03_output/smc_test.rds")

# statistics by month
output2 <- output %>%
  mutate(year = ceiling(timestep/year),
         month = ceiling(timestep/month)) %>%
  # only necessary variables
  dplyr::select(EIR, month, year, seasonality, fifth, starts_with("n_inc_severe"),
                starts_with("p_inc"), starts_with("n_inc"), starts_with("n_")) %>%

  # take means of populations and sums of cases by month
  group_by(EIR, month, year, seasonality, fifth) %>%
  mutate_at(vars(n_0_1825:n_0_36500), mean, na.rm = TRUE) %>%
  mutate_at(vars(n_inc_severe_0_1825:n_inc_clinical_0_36500), sum, na.rm = TRUE) %>%
  dplyr::select(EIR, month, year, seasonality, fifth,
                n_0_1825:n_0_36500, n_inc_severe_0_1825:n_inc_clinical_0_36500) %>%
  distinct() %>%
  mutate(scenario = case_when(fifth==0 ~ '4 rounds normal',
                              fifth==1 ~ '4 rounds shifted',
                              fifth==2 ~ '5 rounds'),
         pfpr = case_when(EIR %in% c(10.2,4.42) ~ 'PfPR: 10',
                          EIR %in% c(30.3,11.4) ~ 'PfPR: 20',
                          EIR %in% c(167,43.6) ~ 'PfPR: 40'))

# plot
ggplot(data = output2) +
  geom_line(aes(x = month, y = p_inc_clinical_0_1825, color = scenario)) +
  labs(title = 'SMC roll-out',
       x = 'month',
       y = '# cases, 0-5 years',
       color = 'scenario',
       alpha = 'EIR') +
  facet_grid(rows=vars(pfpr), cols=vars(seasonality), scales='free') +
  theme_classic()

ggsave('./03_output/smc_test.pdf', width=6, height=6)


output3 <- output2 %>% group_by(pfpr, year, seasonality, fifth, scenario) %>%
  summarize(p_inc_clinical_0_1825 = sum(p_inc_clinical_0_1825, na.rm=T))

ggplot(data = output3) +
  geom_col(aes(x = scenario, y = p_inc_clinical_0_1825, fill = scenario)) +
  labs(title = 'SMC roll-out',
       x = 'scenario',
       y = '# cases, 0-5 years',
       color = 'scenario',
       alpha = 'EIR') +
  facet_grid(rows=vars(pfpr), cols=vars(seasonality), scales='free') +
  theme_classic() +
  theme(axis.text.x = element_blank())

ggsave('./03_output/smc_test2.pdf', width=6, height=6)

