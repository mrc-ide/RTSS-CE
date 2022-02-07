# Cost effectiveness -----------------------------------------------------------
library(tidyverse)
# devtools::install_github("mrc-ide/netz@usage_to_npc") # preliminary version
library(netz)

# pull in data from simulation runs (all interventions)
dalyoutput <- readRDS("./03_output/rtss_long.rds") %>%
  separate(col = age, into = c("age_lower", "age_upper"), sep="-", remove = F) %>%
  mutate(age_lower = as.numeric(age_lower)/365,
         age_upper = as.numeric(age_upper)/365,
         inc = inc_clinical / n,
         sev = inc_severe / n,
         cases = inc_clinical,
         severe_cases = inc_severe)

# make sure age is in years for calculating DALYs
summary(dalyoutput$age_lower); summary(dalyoutput$age_upper)


# DALYs ------------------------------------------------------------------------
# DALYs = Years of life lost (YLL) + Years of live with disease (YLD)
# YLL = Deaths * remaining years of life
# YLD = cases and severe cases * disability weighting  * episode_length
# CE = $ per event (case, death DALY) averted

# Pete code: https://github.com/mrc-ide/gf/blob/69910e798a2ddce240c238d291bc36ea40661b90/R/epi.R
# Weights from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4772264/ {Gunda _et al_, 2016}

# mortality
mortality_rate <- function(x,
                           scaler = 0.215,           # severe case to death scaler
                           treatment_scaler = 0.5) { # treatment modifier
  x %>%
    dplyr::mutate(mortality_rate = (1 - (treatment_scaler * .data$treatment)) * scaler * .data$sev) %>% # mortality rate
    dplyr::mutate(deaths = .data$mortality_rate * .data$n)  # deaths
}


# case and death uncertainty
outcome_uncertainty <- function(x,
                                cases_cv = 0.227,   # case uncertainty SD scaler
                                deaths_cv = 0.265){ # death uncertainty SD scaler
  x %>%
    dplyr::mutate(cases_lower = round(pmax(0, stats::qnorm(0.025, .data$cases, .data$cases * cases_cv))),
                  cases_upper = round(stats::qnorm(0.975, .data$cases, .data$cases * cases_cv)),
                  deaths_lower = round(pmax(0, stats::qnorm(0.025, .data$deaths, .data$deaths * deaths_cv))),
                  deaths_upper = round(stats::qnorm(0.975, .data$deaths, .data$deaths * deaths_cv)))
}

# DALY components
daly_components <- function(x,
                            lifespan = 63,                   # average life expectancy
                            episode_length = 0.01375,        # average length of clinical episode
                            severe_episode_length = 0.04795, # average length of severe episode
                            weight1 = 0.211,      # Disability weight age group 1
                            weight2 = 0.195,      # Disability weight age group 2
                            weight3 = 0.172,      # Disability weight age group 3
                            severe_weight = 0.6){ # Disability weight severe malaria
  x %>%
    dplyr::mutate(yll = .data$deaths * (lifespan - ((.data$age_lower + .data$age_upper) / 2)),
                  yll_lower = .data$deaths_lower * (lifespan - ((.data$age_lower + .data$age_upper) / 2)),
                  yll_upper = .data$deaths_upper * (lifespan - ((.data$age_lower + .data$age_upper) / 2)),

                  yll = ifelse(yll < 0, 0, yll),                    # should be no negative yll from older age groups
                  yll_lower = ifelse(yll_lower < 0, 0, yll_lower),  # should be no negative yll from older age groups
                  yll_upper =  ifelse(yll_upper < 0, 0, yll_upper), # should be no negative yll from older age groups

                  yld = dplyr::case_when(.data$age_upper <= 5 ~ .data$cases * episode_length * weight1 + .data$severe_cases * severe_episode_length * severe_weight,
                                         .data$age_upper > 5 & .data$age_upper <= 15 ~ .data$cases * episode_length * weight2 + .data$severe_cases * severe_episode_length * severe_weight,
                                         .data$age_upper > 15 ~ .data$cases * episode_length * weight3 + .data$severe_cases * severe_episode_length * severe_weight),

                  yld_lower = dplyr::case_when(.data$age_upper <= 5 ~ .data$cases_lower * episode_length * weight1 + .data$severe_cases * severe_episode_length * severe_weight,
                                               .data$age_upper > 5 & .data$age_upper <= 15 ~ .data$cases_lower * episode_length * weight2 + .data$severe_cases * severe_episode_length * severe_weight,
                                               .data$age_upper > 15 ~ .data$cases_lower * episode_length * weight3 + .data$severe_cases * severe_episode_length * severe_weight),

                  yld_upper = dplyr::case_when(.data$age_upper <= 5 ~ .data$cases_upper * episode_length * weight1 + .data$severe_cases * severe_episode_length * severe_weight,
                                               .data$age_upper > 5 & .data$age_upper <= 15 ~ .data$cases_upper * episode_length * weight2 + .data$severe_cases * severe_episode_length * severe_weight,
                                               .data$age_upper > 15 ~ .data$cases_upper * episode_length * weight3 + .data$severe_cases * severe_episode_length * severe_weight)) %>%

    dplyr::mutate(daly = yll + yld)

  # Need to make lifespan more specific
  # Some age groups appear to go up to 200 years (if in days)?

}

# run functions
dalyoutput <- mortality_rate(dalyoutput)
dalyoutput <- outcome_uncertainty(dalyoutput)
dalyoutput <- daly_components(dalyoutput)

# check that there are no negative values
summary(dalyoutput$yll)
summary(dalyoutput[dalyoutput$yll==0,]$yll_lower)
summary(dalyoutput[dalyoutput$yll==0,]$yll_upper)
summary(dalyoutput$yld)

# consolidate across ages
dalyoutput <- dalyoutput %>%
  select(-inc, -sev, -mortality_rate) %>% # get rid of rate vars
  group_by(file) %>%                      # group to condense to one record per run
  mutate(n_182.5_1825 = ifelse(age=='182.5-1825', n, 0)) %>% # create variable for n ages 0.5-5 years to use in costing
  mutate_at(vars(n, n_182.5_1825, inc_clinical:daly), sum, na.rm=T) %>%      # condense outputs over all ages in population
  select(-age, -age_upper, -age_lower) %>%
  distinct()

# check that n_182.5_1825 var is created correctly
summary(dalyoutput$n_182.5_1825)

saveRDS(dalyoutput, './03_output/dalyoutput.rds')


# costing data------------------------------------------------------------------
# costs
# https://github.com/mrc-ide/gf/blob/master/data/unit_costs.rda
# https://github.com/mrc-ide/gf/blob/master/data/treatment_unit_costs.rda
load('./01_data/unit_costs.rda')

PYRcost <- unit_costs$cost_per_pyrethoid_net_delivered      # 5.13
PBOcost <- unit_costs$cost_per_pyrethroid_pbo_net_delivered # 7.05
TREATcost <- 1.47          # clinical treatment cost
SEVcost <- 22.41           # severe treatment cost
SMCcost <- unit_costs$cost_per_smc_dose_delivered           # 1.44
cost_per_dose <- c(2.69, 6.52, 12.91)
delivery_cost <- c(0.96, 1.62, 2.67)

# create combinations of dose cost and delivery cost
rtsscost_df <- expand_grid(cost_per_dose = cost_per_dose, delivery_cost = delivery_cost)

population <- 10000
sim_length <- 15*365

# Prepare to add costs to dataset
dalyoutput_cost <- dalyoutput %>%

  mutate(ITNuse2 = ifelse(ITNboost==1, ITNuse + .10, ITNuse), # account for booster coverage
         ITNcost = case_when(ITN=='pyr' ~ PYRcost,        # account for ITN type-specific cost
                             ITN=='pbo' ~ PBOcost)) %>%

  # count the number of interventions administered and the frequency of ITN dist
  mutate(bednet_distribution_frequency = as.numeric(lapply(lapply(bednet_timesteps,diff), unique)),
         bednet_timesteps = length(Filter(function(x) (x>0 & x<=15*365), bednet_timesteps)),
         smc_timesteps = length(Filter(function(x) (x>0 & x<=15*365), smc_timesteps)),
         rtss_mass_timesteps = length(Filter(function(x) (x>0 & x<=15*365), rtss_mass_timesteps))) %>%

  # merge in RTSS costing dataframe
  merge(rtsscost_df) %>%

  ungroup() %>% rowwise()


# check that the number of ITN distribution times is correct for ITNuse2 == 0.1
dalyoutput_cost %>% filter(ITNboost==1 & ITNuse==0) %>% group_by(bednet_timesteps) %>% summarize()


# read in netz package data to find the annual nets to distribute to give the simulated usage
nets_data <- netz::prepare_data()

# get nets to be distributed for each ITN usage
# assume maximum observed use rate and median bednet half life (across Africa)
nets_distributed <-
  convert_usage_to_annual_nets_distributed(
    target_usage = unique(dalyoutput_cost$ITNuse2),
    distribution_freq = unique(dalyoutput_cost$bednet_distribution_frequency)[
      !(is.na(unique(dalyoutput_cost$bednet_distribution_frequency)))],
    use_rate_data = max(nets_data$use_rate_by_country$use_rate),
    half_life_data = median(nets_data$half_life_data$half_life),
    extrapolate_npc = "loess",
    net_loss_function = net_loss_exp)

# Assumptions to be revised and discussed:
# Using median half life
# Using maximum use rate (with median, can only go up to usage of 81%)
# Extrapolating Loess curve according to curve trend
# Assuming exponential net loss
# check calculations excluding negative DALYS

# save output in case changes in package require changes in code:
saveRDS(nets_distributed, './03_output/net_usage_vs_nets_distributed.rds')

# create cost variables
dalyoutput_cost <- dalyoutput_cost %>%
  left_join(select(nets_distributed, target_use, annual_percapita_nets_distributed),
            by=c('ITNuse2' = 'target_use')) %>%
  mutate(annual_percapita_nets_distributed = ifelse(ITNuse2==0, 0,
                                                    annual_percapita_nets_distributed),
         cost_ITN = population * annual_percapita_nets_distributed * sim_length/365 *
           ITNcost,  # true net cost accounting for non-linear relationship
         cost_ITN_linear = population * ITNuse2 * bednet_timesteps * ITNcost,          # ITN linear
         cost_clinical = (cases-severe_cases) * treatment * TREATcost,                 # non-severe treatment
         cost_severe = severe_cases * treatment * SEVcost,                             # severe treatment
         cost_SMC = n_182.5_1825 * SMC * SMCcost * smc_timesteps,                      # SMC
         cost_vax = (dose1 + dose2 + dose3 + dose4) * (cost_per_dose + delivery_cost), # RTSS

         cost_total = cost_ITN + cost_clinical + cost_severe + cost_SMC + cost_vax)    # TOTAL

saveRDS(dalyoutput_cost, './03_output/dalyoutput_cost.rds')



# group data by scenarios ------------------------------------------------------

output <- dalyoutput_cost %>%
  filter(cost_per_dose==6.52 & delivery_cost==1.62) %>%
  mutate(ID = paste0(pfpr, seasonality, ITNuse, sep="_")) # create unique identifier

# there should be 36 baseline scenarios. 3 pfpr x 3 seasonality x 4 ITN usage
none <- output %>%
  filter(ITNboost==0 & ITN=='pyr' & resistance==0 & RTSS=='none' & # filter out interventions
           (SMC==0 | (seasonality=='highly seasonal'))) %>%
  rename(daly_baseline = daly,
         cost_total_baseline = cost_total) %>%
  select(file, ID, daly_baseline, cost_total_baseline)

base_IDs <- none$file

scenarios <- output %>% filter(!(file %in% base_IDs)) %>%
  left_join(none %>% select(-file), by=c('ID')) %>%
  mutate(CE = (cost_total - cost_total_baseline) / (daly_baseline - daly)) %>% # ICER
  mutate(intervention = case_when(

    ITN=='pyr' & ITNboost==0 & (SMC==0 | (seasonality=='highly seasonal')) & RTSS=='none' ~ 'none',

    ITN=='pbo' & ITNboost==0 & (SMC==0 | (seasonality=='highly seasonal')) & RTSS=='none' ~ 'ITN PBO',
    ITN=='pyr' & ITNboost==1 & (SMC==0 | (seasonality=='highly seasonal')) & RTSS=='none' ~ 'ITN 10% boost',
    ITN=='pyr' & ITNboost==0 & (SMC==0.85 & seasonality=='seasonal') & RTSS=='none' ~ 'SMC',
    ITN=='pyr' & ITNboost==0 & (SMC==0 | (seasonality=='highly seasonal')) & RTSS=='SV' ~ 'RTS,S SV',
    ITN=='pyr' & ITNboost==0 & (SMC==0 | (seasonality=='highly seasonal')) & RTSS=='EPI' ~ 'RTS,S EPI',

    ITN=='pbo' & ITNboost==0 & (SMC==0.85 & (seasonality=='seasonal')) & RTSS=='none' ~ 'ITN PBO + SMC',
    ITN=='pyr' & ITNboost==1 & (SMC==0.85 & (seasonality=='seasonal')) & RTSS=='none' ~ 'ITN 10% boost + SMC',

    ITN=='pbo' & ITNboost==0 & (SMC==0 | (seasonality=='highly seasonal')) & RTSS!='none' ~ 'ITN PBO + RTS,S',
    ITN=='pyr' & ITNboost==1 & (SMC==0 | (seasonality=='highly seasonal')) & RTSS!='none' ~ 'ITN 10% boost + RTS,S',

    ITN=='pyr' & ITNboost==0 & (SMC==0.85 & (seasonality=='seasonal')) & RTSS!='none' ~ 'RTS,S + SMC',
    ITN=='pyr' & ITNboost==1 & (SMC==0.85 & (seasonality=='seasonal')) & RTSS!='none' ~ 'ITN 10% boost + RTS,S + SMC',
    ITN=='pbo' & ITNboost==0 & (SMC==0.85 & (seasonality=='seasonal')) & RTSS!='none' ~ 'ITN PBO + RTS,S + SMC')) %>%

  mutate(intervention = factor(intervention, levels=c('none', 'ITN 10% boost', 'ITN PBO', 'RTS,S EPI', 'RTS,S SV', 'SMC', 'ITN 10% boost + RTS,S', 'ITN PBO + RTS,S', 'ITN 10% boost + SMC', 'ITN PBO + SMC', 'RTS,S + SMC', 'ITN 10% boost + RTS,S + SMC', 'ITN PBO + RTS,S + SMC')))

# check intervention
table(scenarios$intervention, useNA='always')


# inspect range of scenarios
table(scenarios$ID)
summary(scenarios$CE)
test <- scenarios %>% filter(CE<0)
table(scenarios$ID, scenarios$intervention) # three of scenarios with resistance, and one with 0 resistance (ITNuse=0)

saveRDS(scenarios, './03_output/scenarios.rds')


# data checks ##################################################################

# prev & daly
none <- dalyoutput_cost %>%
  filter(ITNboost==0 & RTSS=='none' & ITN=='pyr' & resistance==0 & (SMC==0 | (seasonality=='highly seasonal'))) %>%
  rename(daly_baseline = daly,
         cost_total_baseline = cost_total)

  ggplot(none, aes(x=pfpr, y=daly_baseline, color=factor(ITNuse))) +
  geom_point() +
  geom_smooth(method = 'lm', se=F) +
  facet_wrap(~seasonality) + theme_classic() +
  labs(x='PfPR', y='DALYs', color='ITN use',
       title = "Baseline scenarios, DALYs by PfPR and seasonality")

  ggplot(none, aes(x=pfpr, y=yll, color=factor(ITNuse))) +
    geom_point() +
    geom_smooth(method = 'lm', se=F) +
    facet_wrap(~seasonality) + theme_classic() +
    labs(x='PfPR', y='YLL', color='ITN use',
         title = "Baseline scenarios, YLLs by PfPR and seasonality")

  ggplot(none, aes(x=pfpr, y=yld, color=factor(ITNuse))) +
    geom_point() +
    geom_smooth(method = 'lm', se=F) +
    facet_wrap(~seasonality) + theme_classic() +
    labs(x='PfPR', y='YLD', color='ITN use',
         title = "Baseline scenarios, YLDs by PfPR and seasonality")

# RTSS doses are stable
table(dalyoutput_cost$RTSS, dalyoutput_cost$RTSScov, dalyoutput_cost$rtss_mass_timesteps)

# SMC doses are stable
table(dalyoutput_cost$seasonality, dalyoutput_cost$smc_timesteps)

# ITNs are stable
table(dalyoutput_cost$ITNuse, dalyoutput_cost$bednet_timesteps)

# no negative DALYs
summary(dalyoutput_cost$daly)

# DALYs averted dist
test <- scenarios %>% mutate(deltadaly = daly_baseline - daly)
summary(test$deltadaly)
test <- test %>% filter(deltadaly < 0)
table(test$resistance) # most negative DALYs are in the high resistance scenarios

# double-check relationship between cost_ITN_linear and cost_ITN:
# similar at low usage but divering at higher usage
ggplot(dalyoutput_cost) +
  geom_point(aes(x=cost_ITN_linear, y = cost_ITN, colour=ITNuse2)) +
  geom_abline(slope=1) +
  theme_classic()

################################################################################


