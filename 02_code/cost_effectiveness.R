# Cost effectiveness -----------------------------------------------------------
library(tidyverse)

# pull in data from simulation runs (all interventions)
dat <- readRDS("./03_output/rtss_raw.rds")

# averted cases / severe-cases / deaths ----------------------------------------
# summarize over first 5 years
dat2 <- dat %>%
  filter(year <= 5) %>% # first 5 years
  group_by(file) %>%
  mutate_at(vars(n_0_1825:n_36500_73000), mean, na.rm = TRUE) %>% # mean of n in each age group
  mutate_at(vars(n_inc_severe_0_1825:dose4), sum, na.rm = TRUE) %>% # sum of cases and vax doses
  select(-month, -year) %>%
  distinct()

# calculate outputs by age
dat3 <- dat2 %>%
  dplyr::select(file:n_36500_73000, p_inc_clinical_0_1825:p_inc_severe_36500_73000, dose1:dose4) %>%
  pivot_longer(cols = c(n_0_1825:n_36500_73000, p_inc_clinical_0_1825:p_inc_clinical_36500_73000, p_inc_severe_0_1825:p_inc_severe_36500_73000), names_to = c('age'), values_to = c('value')) %>% # moving to long age groups
  mutate(n = ifelse(grepl('n_', age), value, NA),
         inc_clinical = ifelse(grepl('p_inc_clinical', age), value, NA), # creating var for inc_clinical
         inc_severe = ifelse(grepl('p_inc_severe', age), value, NA), # creating var for inc_severe
         age = gsub('n_', '', age), # creating var for age group
         age = gsub('p_inc_clinical_', '', age),
         age = gsub('p_inc_severe_', '', age),
         age = gsub('_', '-', age)) %>%
  group_by(file, age) %>%
  select(-value) %>%
  mutate_at(vars(n:inc_severe), sum, na.rm = TRUE) %>% # consolidate
  distinct() %>% ungroup()

# calculate outputs averted
baseline <- dat3 %>% # summarizing outputs when there is no intervention
  filter(ITNuse==0, RTSS=='none', ITN=='pyr', SMC==0) %>%
  select(file, EIR, pfpr, seasonality, age:inc_severe) %>%
  rename(n_baseline = n, inc_clinical_baseline = inc_clinical, inc_severe_baseline = inc_severe)

averted <- dat3 %>%
  filter(!(file %in% baseline$file)) %>% # taking out scenarios with no intervention
  left_join(baseline %>% select(-file), by=c('EIR', 'pfpr', 'seasonality', 'age')) %>% # adding baseline data
  mutate(case_avert = ((inc_clinical/n) - (inc_clinical_baseline/n_baseline)) * 10000, # cases averted per 100,000 people
         severe_avert = ((inc_severe/n) - (inc_severe_baseline/n_baseline)) * 10000) # severe cases averted per 100,000 people

# DALYs ------------------------------------------------------------------------
# DALYs = Years of life lost (YLL) + Years of live with disease (YLD)
# YLL = Deaths * remaining years of life
# YLD = cases and severe cases * disability weighting  * episode_length
# CE = $ per event (case, death DALY) averted

# Pete code: https://github.com/mrc-ide/gf/blob/69910e798a2ddce240c238d291bc36ea40661b90/R/epi.R#L89-L118
# Weights from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4772264/ {Gunda _et al_, 2016}
lifespan = 63,
episode_length = 0.01375
severe_episode_length = 0.04795
weight1 = 0.211     # Disability weight age group 1
weight2 = 0.195     # Disability weight age group 2
weight3 = 0.172     # Disability weight age group 3
severe_weight = 0.6 # Disability weight severe malaria

scaler = 0.215           #severe case to death modifier
treatment_scaler = 0.5   #treatment modifier

mortality_rate <- function(x, scaler = 0.215, treatment_scaler = 0.5){
  x %>%
    dplyr::mutate(mortality_rate = (1 - (treatment_scaler * .data$treatment_coverage)) * scaler * .data$sev)
}


dplyr::mutate(
  yll = .data$deaths * (lifespan - ((.data$age_lower + .data$age_upper) / 2)),
  yll_lower = .data$deaths_lower * (lifespan - ((.data$age_lower + .data$age_upper) / 2)),
  yll_upper = .data$deaths_upper * (lifespan - ((.data$age_lower + .data$age_upper) / 2)),

  yld = dplyr::case_when(
    .data$age_upper <= 5 ~ .data$cases * episode_length * weight1 + .data$severe_cases * severe_episode_length * severe_weight,
    .data$age_upper > 5 & .data$age_upper <= 15 ~ .data$cases * episode_length * weight2 + .data$severe_cases * severe_episode_length * severe_weight,
    .data$age_upper > 15 ~ .data$cases * episode_length * weight3 + .data$severe_cases * severe_episode_length * severe_weight),

  yld_lower = dplyr::case_when(
    .data$age_upper <= 5 ~ .data$cases_lower * episode_length * weight1 + .data$severe_cases * severe_episode_length * severe_weight,
    .data$age_upper > 5 & .data$age_upper <= 15 ~ .data$cases_lower * episode_length * weight2 + .data$severe_cases * severe_episode_length * severe_weight,
    .data$age_upper > 15 ~ .data$cases_lower * episode_length * weight3 + .data$severe_cases * severe_episode_length * severe_weight),

  yld_upper = dplyr::case_when(
    .data$age_upper <= 5 ~ .data$cases_upper * episode_length * weight1 + .data$severe_cases * severe_episode_length * severe_weight,
    .data$age_upper > 5 & .data$age_upper <= 15 ~ .data$cases_upper * episode_length * weight2 + .data$severe_cases * severe_episode_length * severe_weight,
    .data$age_upper > 15 ~ .data$cases_upper * episode_length * weight3 + .data$severe_cases * severe_episode_length * severe_weight))

daly = yyl + yld


test <- averted %>% filter(file==1)

#-costing data-----------------------------------------------------------------
cost_per_dose <- c(2.69,6.52,12.91) #cost per vaccine - including the cost of wastage etc $2,$5,$10 per dose
delivery_cost <- c(0.96,1.62,2.67)  #EPI delivery costs
tx_unit_cost  <- 1.47               #clinical treatment cost
severe_unit_cost <- 22.41           #severe treatment cost

cost_df <- expand_grid(cost_per_dose = cost_per_dose, delivery_cost = delivery_cost)



# notes for RTS,S call
costing of interventions
how to calculate mortality
life tables
getting upper and lower confidence estimates
