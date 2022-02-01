# Cost effectiveness -----------------------------------------------------------
library(tidyverse)

# pull in data from simulation runs (all interventions)
dalyoutput <- readRDS("C:/Users/htopazia/OneDrive - Imperial College London/Github/GF-RTSS-CE/03_output/rtss_long.rds") %>%
  separate(col = age, into = c("age_lower", "age_upper"), sep="-", remove = F) %>%
  mutate(age_lower = as.numeric(age_lower),
         age_upper = as.numeric(age_upper),
         inc = inc_clinical / n,
         sev = inc_severe / n,
         cases = inc_clinical,
         severe_cases = inc_severe)


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
}

# run functions
dalyoutput <- mortality_rate(dalyoutput)
dalyoutput <- outcome_uncertainty(dalyoutput)
dalyoutput <- daly_components(dalyoutput)

# consolidate across ages
dalyoutput <- dalyoutput %>%
  select(-inc, -sev, -mortality_rate) %>% # get rid of rate vars
  group_by(file) %>%                      # group to condense to one record per run
  mutate(n_182.5_1825 = n_182.5_1825) %>% # create variable for n ages 0.5-5 years to use in costing
  mutate_at(n, inc_clinical:daly)         # condense outputs over all ages in population

saveRDS(dalyoutput, './03_output/dalyoutput.rds')


# costing data------------------------------------------------------------------
# costs
# https://github.com/mrc-ide/gf/blob/master/data/unit_costs.rda
# https://github.com/mrc-ide/gf/blob/master/data/treatment_unit_costs.rda
load('C:/Users/htopazia/OneDrive - Imperial College London/Github/GF-RTSS-CE/01_data/unit_costs.rda')

PYRcost <- unit_costs$cost_per_pyrethoid_net_delivered
PBOcost <- unit_costs$cost_per_pyrethroid_pbo_net_delivered
TREATcost <- 1.47          # clinical treatment cost
SEVcost <- 22.41           # severe treatment cost
SMCcost <- unit_costs$cost_per_smc_dose_delivered
cost_per_dose <- c(2.69, 6.52, 12.91)
delivery_cost <- c(0.96, 1.62, 2.67)

# create combinations of dose cost and delivery cost
rtsscost_df <- expand_grid(cost_per_dose = cost_per_dose, delivery_cost = delivery_cost)

population <- 10000
sim_length <- 15*365

# add costs to dataset
dalyoutput_allages <- dalyoutput %>%

  mutate(ITNuse = ifelse(boost==1, ITNuse + .10, ITNuse), # account for booster coverage
         ITNcost = case_when(ITN=='pyr' ~ PYRcost,        # account for ITN type-specific cost
                             ITN=='pbo' ~ PBOcost)) %>%

  # count number of interventions administered
  mutate(bednet_timesteps = as_tibble(unlist(bednet_timesteps)) %>%
           filter(value>0 & value<=sim_length) %>% count() %>% as.numeric(),
         smc_timesteps = as_tibble(unlist(smc_timesteps)) %>%
           filter(value>0 & value<=sim_length) %>% count() %>% as.numeric()) %>%

  # merge in RTSS costing dataframe
  merge(rtsscost_df) %>%

  # create cost variables
  mutate(cost_novax = population * ITNuse * bednet_timesteps * ITNcost + # ITN
                      n_treated * TREATcost +                            # treatment
                      severe_cases * treatment * SEVcost +               # severe treatment
                      n_182.5_1825 * SMC * SMCcost * smc_timesteps,      # SMC

         cost_vax = (dose1 + dose2 + dose3 + dose4) * (cost_per_dose + delivery_cost),

         cost_total = cost_novax + cost_vax)


saveRDS(dalyoutput_allages, './03_output/dalyoutput_allages.rds')

