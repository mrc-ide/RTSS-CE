# Cost effectiveness -----------------------------------------------------------
library(tidyverse)
library(netz)
library(treasure)
# devtools::install_github("mrc-ide/netz@usage_to_npc") # preliminary version
# devtools::install_github("mrc-ide/treasure")


# pull in data from simulation runs (all interventions)
dalyoutput <- readRDS("./03_output/rtss_long.rds") %>%
  separate(col = age, into = c("age_lower", "age_upper"), sep="-", remove = F) %>%
  mutate(age_lower = as.numeric(age_lower)/365,
         age_upper = as.numeric(age_upper)/365,
         inc = inc_clinical / n,
         sev = inc_severe / n,
         cases = inc_clinical,
         severe_cases = inc_severe)

dalyoutput2 <- filter(dalyoutput, pfpr == 0.4, seasonality == "seasonal",
                     ITNuse == 0.5, treatment == 0.45, resistance == 0.4,
                     ITNboost == 0, SMC == 0, RTSScov == 0, ITN == "pyr",
                     RTSS == 'none')


# make sure age is in years for calculating DALYs
summary(dalyoutput$age_lower); summary(dalyoutput$age_upper)


# DALYs ------------------------------------------------------------------------
# DALYs = Years of life lost (YLL) + Years of life with disease (YLD)
# YLL = Deaths * remaining years of life
# YLD = cases and severe cases * disability weighting  * episode_length
# CE = $ per event (case, death DALY) averted

# Reference code: https://github.com/mrc-ide/gf/blob/69910e798a2ddce240c238d291bc36ea40661b90/R/epi.R
# Weights from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4772264/ {Gunda _et al_, 2016}

# mortality
mortality_rate <- function(x,
                           scaler = 0.215,           # severe case to death scaler
                           treatment_scaler = 0.5) { # treatment modifier
  x %>%
    dplyr::mutate(mortality_rate = scaler * .data$sev) %>% # mortality rate alternative (consistent with old ICL analysis)
    # dplyr::mutate(mortality_rate = (1 - (treatment_scaler * .data$treatment)) * scaler * .data$sev) %>% # mortality rate
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
                            lifespan = 64.49,                   # average life expectancy
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

    dplyr::mutate(daly = yll + yld,
                  daly_upper = yll_lower + yld_lower,
                  daly_lower = yll_upper + yld_upper)

}

# run functions
dalyoutput <- mortality_rate(dalyoutput)
dalyoutput <- outcome_uncertainty(dalyoutput)
dalyoutput <- daly_components(dalyoutput)

# check that there are no negative values
summary(dalyoutput$yll)
summary(dalyoutput$yld)

# consolidate across ages
dalyoutput <- dalyoutput %>%
  select(-inc, -sev, -mortality_rate) %>% # get rid of rate vars
  group_by(file) %>%                      # group to condense to one record per run

  # create vars for childhood cases for PCV comparison
  mutate(n_0_1825 = ifelse(age %in% c('0-91.25', '91.25-1825'), n, 0),
         n_91.25_1825 = ifelse(age=='91.25-1825', n, 0), # SMC denominator
         u5_cases = ifelse(age %in% c('0-91.25', '91.25-1825'), cases, 0),
         u5_severe = ifelse(age %in% c('0-91.25', '91.25-1825'), severe_cases, 0),
         u5_dalys = ifelse(age %in% c('0-91.25', '91.25-1825'), daly, 0)) %>%

  mutate_at(vars(n, n_0_1825, n_91.25_1825, u5_cases, u5_severe, u5_dalys, inc_clinical:daly_lower), sum, na.rm=T) %>%  # condense outputs over all ages in population
  select(-age, -age_upper, -age_lower) %>%
  distinct()


saveRDS(dalyoutput, './03_output/dalyoutput.rds')


# costing data------------------------------------------------------------------

# references in: https://mrc-ide.github.io/treasure/reference/index.html
# treatment costs from Penny et al. 2016

PYRcost <- 2.52 + 1.50                # pyrethroid $2.52 per net and $1.50 delivery cost
PBOcost <- 3.51 + 1.50                # PBO $3.51 per net and $1.50 delivery cost

SMCcost <- 0.9075                     # SMC $0.9075 per dose

cost_per_dose <- c(2.69, 6.52, 12.91) # RTS,S per dose $2, $5, $10 + consumables cost
delivery_cost <- 1.62                 # RTS,S delivery cost range c(0.96, 1.62, 2.67)

# create combinations of dose cost and delivery cost
rtsscost_df <- expand_grid(cost_per_dose = cost_per_dose, delivery_cost = delivery_cost)

RDT <- 0.46 + (0.46 * 0.15)           # RDT $0.46 unit cost + (unit cost * 15% delivery markup)
AL_adult <- 0.3 * 24                  # $7.2 clinical treatment cost ($0.3 * 24 doses)
AL_child <- 0.3 * 12                  # $3.6 clinical treatment cost ($0.3 * 12 doses)

outpatient <- 1.87                    # (Median WHO Choice cost for SSA)
inpatient <- 8.71                    # (Median WHO Choice cost for SSA, assuming average duration of stay of 3 days)

# clinical: RDT cost + Drug course cost + facility cost (outpatient)
TREATcost_adult <- RDT + AL_adult + outpatient
TREATcost_child <- RDT + AL_child + outpatient

# severe: RDT cost + Drug course cost + facility cost (inpatient)
SEVcost_adult <- RDT + AL_adult + inpatient
SEVcost_child <- RDT + AL_child + inpatient


# population and simulation length
population <- dalyoutput$population[[1]]
sim_length <- dalyoutput$sim_length[[1]]

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


# read in netz package data to find the annual nets to distribute to give the simulated usage
nets_data <- netz::prepare_data()
summary(nets_data$use_rate_by_country)
# min use_rate = .66, so can only go up to 66% usage
# max use_rate = .96


# get nets to be distributed for each ITN usage
ndist <- function(x) {

  convert_usage_to_annual_nets_distributed(
    target_usage = unique(dalyoutput_cost$ITNuse2),
    distribution_freq = unique(dalyoutput_cost$bednet_distribution_frequency)[
      !(is.na(unique(dalyoutput_cost$bednet_distribution_frequency)))],
    use_rate_data = x,
    half_life_data = median(nets_data$half_life_data$half_life),
    extrapolate_npc = "loess",
    net_loss_function = net_loss_map) %>%
  select(target_use, annual_percapita_nets_distributed)

}

# assume maximum observed use rate and median bednet half life (across Africa)
nets_distributed <- ndist(0.88)

# assume observed rate is the min in Africa
nets_distributed_min <- ndist(min(nets_data$use_rate_by_country$use_rate))
nets_distributed_min <- rename(nets_distributed_min, annual_percapita_nets_distmin = annual_percapita_nets_distributed)

# assume observed rate is the max in Africa
nets_distributed_max <- ndist(max(nets_data$use_rate_by_country$use_rate))
nets_distributed_max <- rename(nets_distributed_max, annual_percapita_nets_distmax = annual_percapita_nets_distributed)

nets_distributed <- full_join(nets_distributed, nets_distributed_min) %>% full_join(nets_distributed_max)


# Assumptions:
# Using median half life
# Using minimum use rate (88%) allowing to give 85% usage (this is between median and max)
# Extrapolating Loess curve according to curve trend
# Assuming net loss is like in MAP paper (smooth compact)

# save output in case changes in package require changes in code:
# saveRDS(nets_distributed, './03_output/net_usage_vs_nets_distributed.rds')


# create cost variables
# 77% of treatment costs are from the public sector (DHS, SSA)
dalyoutput_cost <- dalyoutput_cost %>%
  left_join(nets_distributed,
            by=c('ITNuse2' = 'target_use')) %>%
  mutate(annual_percapita_nets_distributed = ifelse(ITNuse2==0, 0,
                                                    annual_percapita_nets_distributed),
 # calculate cost of interventions
 # true net cost accounting for non-linear relationship
 cost_ITN = population * annual_percapita_nets_distributed * sim_length/365 * ITNcost,
 # true net cost MIN
 cost_ITNmin = population * annual_percapita_nets_distmin * sim_length/365 * ITNcost,
 # true net cost MAX
 cost_ITNmax = population * annual_percapita_nets_distmax * sim_length/365 * ITNcost,
 # ITN linear
 cost_ITN_linear = population * ITNuse2 * bednet_timesteps * ITNcost,
 # non-severe treatment
 cost_clinical = ((cases - severe_cases - u5_cases) * treatment * TREATcost_adult) * .77 + # adult
   ((u5_cases - u5_severe) * treatment * TREATcost_child) * .77, # child
 # severe treatment
 cost_severe = (severe_cases * treatment * SEVcost_adult) * .77 + # adult
   (u5_severe * treatment * SEVcost_child) * .77, # child
 # SMC
 cost_SMC = n_91.25_1825 * SMC * SMCcost * smc_timesteps,
 # RTSS
 cost_vax = (dose1 + dose2 + dose3 + dose4) * (cost_per_dose + delivery_cost),

 # TOTAL
 cost_total = cost_ITN + cost_clinical + cost_severe + cost_SMC + cost_vax,
 cost_total_ITNmin = cost_ITNmin + cost_clinical + cost_severe + cost_SMC + cost_vax,
 cost_total_ITNmax = cost_ITNmax + cost_clinical + cost_severe + cost_SMC + cost_vax,

 # cost just among children
 cost_total_u5 =
   # cost ITNs
   n_0_1825 * annual_percapita_nets_distributed * sim_length/365 * ITNcost +
   # cost clinical
   ((u5_cases-u5_severe) * treatment * TREATcost_child)*.77 +
   # cost severe
   (u5_severe * treatment * SEVcost_child)*.77 +
   # cost SMC and cost VAX are all among children
   cost_SMC + cost_vax)

saveRDS(dalyoutput_cost, './03_output/dalyoutput_cost.rds')



# assign scenarios -------------------------------------------------------------
output <- dalyoutput_cost %>%
  filter(cost_per_dose == 6.52) %>%
  mutate(ID = paste(pfpr, seasonality, ITNuse, resistance, treatment, sep="_")) # create unique identifier

# there should be 270 baseline scenarios. 3 pfpr x 3 seasonality x 4 ITN usage x 3 treatmet x 2.5 resistance (only pbo in resistance scenarios)
none <- output %>%
  filter(ITNboost==0 & ITN=='pyr' & RTSS=='none' & # filter out interventions
           (SMC==0 | (seasonality=='highly seasonal'))) %>%
  rename(daly_baseline = daly,
         cases_baseline = cases,
         severe_baseline = severe_cases,
         deaths_baseline = deaths,

         u5_dalys_baseline = u5_dalys,
         u5_cases_baseline = u5_cases,
         u5_severe_baseline = u5_severe,

         cost_total_baseline = cost_total,
         cost_total_ITNmin_baseline = cost_total_ITNmin,
         cost_total_ITNmax_baseline = cost_total_ITNmax,
         cost_total_u5_baseline = cost_total_u5,
         ) %>%
  select(file, ID, daly_baseline, cases_baseline, severe_baseline, deaths_baseline,
         u5_dalys_baseline, u5_cases_baseline, u5_severe_baseline,
         cost_total_baseline, cost_total_ITNmin_baseline, cost_total_ITNmax_baseline, cost_total_u5_baseline)

base_IDs <- none$file

scenarios <- output %>% filter(!(file %in% base_IDs)) %>%
  left_join(none %>% select(-file), by=c('ID')) %>%
  mutate(CE = (cost_total - cost_total_baseline) / (daly_baseline - daly),
         CE_u5 = (cost_total - cost_total_baseline) / (u5_dalys_baseline - u5_dalys),
         CE_ITNmin = (cost_total_ITNmin - cost_total_ITNmin_baseline) / (daly_baseline - daly),
         CE_ITNmax = (cost_total_ITNmax - cost_total_ITNmax_baseline) / (daly_baseline - daly),
         CE_case = (cost_total - cost_total_baseline) / (cases_baseline - cases),
         CE_u5_case = (cost_total - cost_total_baseline) / (u5_cases_baseline - u5_cases)
         ) %>% # ICER
  mutate(intervention = case_when(

    ITN=='pyr' & ITNboost==0 & (SMC==0 | (seasonality=='highly seasonal')) & RTSS=='none' ~ 'none',

    ITN=='pyr' & ITNboost==1 & (SMC==0 | (seasonality=='highly seasonal')) & RTSS=='none' ~ 'ITN 10% increase',
    ITN=='pbo' & ITNboost==0 & (SMC==0 | (seasonality=='highly seasonal')) & RTSS=='none' ~ 'ITN PBO',
    ITN=='pyr' & ITNboost==0 & (SMC==0 | (seasonality=='highly seasonal')) & RTSS=='EPI' ~ 'RTS,S age-based',
    ITN=='pyr' & ITNboost==0 & (SMC==0 | (seasonality=='highly seasonal')) & RTSS=='SV' ~ 'RTS,S seasonal',
    ITN=='pyr' & ITNboost==0 & (SMC==0.85 & seasonality=='seasonal') & RTSS=='none' ~ 'SMC',
    ITN=='pyr' & ITNboost==1 & (SMC==0 | (seasonality=='highly seasonal')) & RTSS!='none' ~ 'ITN 10% increase + RTS,S',
    ITN=='pbo' & ITNboost==0 & (SMC==0 | (seasonality=='highly seasonal')) & RTSS!='none' ~ 'ITN PBO + RTS,S',
    ITN=='pyr' & ITNboost==1 & (SMC==0.85 & (seasonality=='seasonal')) & RTSS=='none' ~ 'ITN 10% increase + SMC',
    ITN=='pbo' & ITNboost==0 & (SMC==0.85 & (seasonality=='seasonal')) & RTSS=='none' ~ 'ITN PBO + SMC',
    ITN=='pyr' & ITNboost==0 & (SMC==0.85 & (seasonality=='seasonal')) & RTSS!='none' ~ 'RTS,S + SMC',
    ITN=='pyr' & ITNboost==1 & (SMC==0.85 & (seasonality=='seasonal')) & RTSS!='none' ~ 'ITN 10% increase + RTS,S + SMC',
    ITN=='pbo' & ITNboost==0 & (SMC==0.85 & (seasonality=='seasonal')) & RTSS!='none' ~ 'ITN PBO + RTS,S + SMC')) %>%

  mutate(intervention_f = factor(intervention, levels=c('none', 'ITN 10% increase', 'ITN PBO', 'RTS,S age-based', 'RTS,S seasonal', 'SMC', 'ITN 10% increase + RTS,S', 'ITN PBO + RTS,S', 'ITN 10% increase + SMC', 'ITN PBO + SMC', 'RTS,S + SMC', 'ITN 10% increase + RTS,S + SMC', 'ITN PBO + RTS,S + SMC'))) %>%

  mutate(rank=as.numeric(intervention_f))

table(scenarios$intervention_f, scenarios$rank, useNA = 'always')


# inspect range of scenarios
table(scenarios$ID)

summary(scenarios$CE[scenarios$resistance==0])

test <- scenarios %>% filter(CE<0) %>% select(pfpr, seasonality, resistance, intervention, daly, daly_baseline, CE); table(test$resistance) # all negative CE scenarios have resistance except 2;

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
       title = "DALYs increase by PfPR and seasonality")

  ggplot(none, aes(x=pfpr, y=yll, color=factor(ITNuse))) +
    geom_point() +
    geom_smooth(method = 'lm', se=F) +
    facet_wrap(~seasonality) + theme_classic() +
    labs(x='PfPR', y='YLL', color='ITN use',
         title = "YLLs increase by PfPR and seasonality")

  ggplot(none, aes(x=pfpr, y=yld, color=factor(ITNuse))) +
    geom_point() +
    geom_smooth(method = 'lm', se=F) +
    facet_wrap(~seasonality) + theme_classic() +
    labs(x='PfPR', y='YLD', color='ITN use',
         title = "YLDs increase by PfPR and seasonality")

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
table(test$resistance) # all negative DALYs are in resistance scenarios

# double-check relationship between cost_ITN_linear and cost_ITN:
# similar at low usage but divering at higher usage
ggplot(dalyoutput_cost) +
  geom_point(aes(x=cost_ITN_linear, y = cost_ITN, colour=ITNuse2)) +
  facet_wrap(~ITN) +
  geom_abline(slope=1) +
  theme_classic()

#------------------------------------------------------------------------------#




# Case-study -------------------------------------------------------------------

# pull in data from simulation runs (all interventions)
dalyoutput <- readRDS("./03_output/rtss_long_casestudy.rds") %>%
  separate(col = age, into = c("age_lower", "age_upper"), sep="-", remove = F) %>%
  mutate(age_lower = as.numeric(age_lower)/365,
         age_upper = as.numeric(age_upper)/365,
         inc = inc_clinical / n,
         sev = inc_severe / n,
         cases = inc_clinical,
         severe_cases = inc_severe)

# run functions
dalyoutput <- mortality_rate(dalyoutput)
dalyoutput <- outcome_uncertainty(dalyoutput)
dalyoutput <- daly_components(dalyoutput)

# check that there are no negative values
summary(dalyoutput$yll)
summary(dalyoutput$yld)

# consolidate across ages
dalyoutput <- dalyoutput %>%
  select(-inc, -sev, -mortality_rate) %>% # get rid of rate vars
  group_by(file) %>%                      # group to condense to one record per run

  # create vars for childhood cases for PCV comparison
  mutate(n_0_1825 = ifelse(age %in% c('0-91.25', '91.25-1825'), n, 0),
         n_91.25_1825 = ifelse(age=='91.25-1825', n, 0), # SMC denominator
         u5_cases = ifelse(age %in% c('0-91.25', '91.25-1825'), cases, 0),
         u5_severe = ifelse(age %in% c('0-91.25', '91.25-1825'), severe_cases, 0),
         u5_deaths = ifelse(age %in% c('0-91.25', '91.25-1825'), deaths, 0),
         u5_dalys = ifelse(age %in% c('0-91.25', '91.25-1825'), daly, 0),

         u5_cases_lower = ifelse(age %in% c('0-91.25', '91.25-1825'), cases_lower, 0),
         u5_deaths_lower = ifelse(age %in% c('0-91.25', '91.25-1825'), deaths_lower, 0),
         u5_dalys_lower = ifelse(age %in% c('0-91.25', '91.25-1825'), daly_lower, 0),

         u5_cases_upper = ifelse(age %in% c('0-91.25', '91.25-1825'), cases_upper, 0),
         u5_deaths_upper = ifelse(age %in% c('0-91.25', '91.25-1825'), deaths_upper, 0),
         u5_dalys_upper = ifelse(age %in% c('0-91.25', '91.25-1825'), daly_upper, 0)
         ) %>%

  mutate(across(c(n, n_0_1825, n_91.25_1825, u5_cases:u5_dalys_upper, inc_clinical:daly_lower), ~sum(.x, na.rm=T))) %>%  # condense outputs over all ages in population
  select(-age, -age_upper, -age_lower) %>%
  distinct()

# costing data
# references in: https://mrc-ide.github.io/treasure/reference/index.html
# treatment costs from Penny et al. 2016

PYRcost <- 2.52 + 1.50                # pyrethroid $2.52 per net and $1.50 delivery cost
PBOcost <- 3.51 + 1.50                # PBO $3.51 per net and $1.50 delivery cost

SMCcost <- 0.9075                     # SMC $0.9075 per dose

cost_per_dose <- c(2.69, 6.52, 12.91) # RTS,S per dose $2, $5, $10 + consumables cost
delivery_cost <- 1.62                 # RTS,S delivery cost range c(0.96, 1.62, 2.67)

# create combinations of dose cost and delivery cost
rtsscost_df <- expand_grid(cost_per_dose = cost_per_dose, delivery_cost = delivery_cost)

RDT <- 0.46 + (0.46 * 0.15)           # RDT $0.46 unit cost + (unit cost * 15% delivery markup)
AL_adult <- 0.3 * 24                  # $7.2 clinical treatment cost ($0.3 * 24 doses)
AL_child <- 0.3 * 12                  # $3.6 clinical treatment cost ($0.3 * 12 doses)

outpatient <- 1.87                    # (Median WHO Choice cost for SSA)
inpatient <- 8.71                    # (Median WHO Choice cost for SSA, assuming average duration of stay of 3 days)

# clinical: RDT cost + Drug course cost + facility cost (outpatient)
TREATcost_adult <- RDT + AL_adult + outpatient
TREATcost_child <- RDT + AL_child + outpatient

# severe: RDT cost + Drug course cost + facility cost (inpatient)
SEVcost_adult <- RDT + AL_adult + inpatient
SEVcost_child <- RDT + AL_child + inpatient

population <- dalyoutput$population[[1]]
sim_length <- dalyoutput$sim_length[[1]]

# Prepare to add costs to dataset
dalyoutput_cost <- dalyoutput %>%

  mutate(ITNuse2 = ifelse(ITNboost==1, ITNuse + .10, ITNuse), # account for booster coverage
         ITNuse2 = round(ITNuse2, 2),
         ITNcost = case_when(ITN=='pyr' ~ PYRcost)) %>%

  # count the number of interventions administered and the frequency of ITN dist
  mutate(bednet_distribution_frequency = as.numeric(lapply(lapply(bednet_timesteps,diff), unique)),
         bednet_timesteps = length(Filter(function(x) (x>0 & x<=15*365), bednet_timesteps)),
         smc_timesteps = length(Filter(function(x) (x>0 & x<=15*365), smc_timesteps)),
         rtss_mass_timesteps = length(Filter(function(x) (x>0 & x<=15*365), rtss_mass_timesteps))) %>%

  # merge in RTSS costing dataframe
  merge(rtsscost_df) %>%
  ungroup() %>% rowwise()

# assume maximum observed use rate and median bednet half life (across Africa)
nets_distributed <- ndist(0.88)

# assume observed rate is the min in Africa
nets_distributed_min <- ndist(min(nets_data$use_rate_by_country$use_rate))
nets_distributed_min <- rename(nets_distributed_min, annual_percapita_nets_distmin = annual_percapita_nets_distributed)

# assume observed rate is the max in Africa
nets_distributed_max <- ndist(max(nets_data$use_rate_by_country$use_rate))
nets_distributed_max <- rename(nets_distributed_max, annual_percapita_nets_distmax = annual_percapita_nets_distributed)

nets_distributed <- full_join(nets_distributed, nets_distributed_min) %>% full_join(nets_distributed_max)


# create cost variables
# 77% of treatment costs are from the public sector (DHS, SSA)
dalyoutput_cost <- dalyoutput_cost %>%
  left_join(nets_distributed,
            by=c('ITNuse2' = 'target_use')) %>%
  mutate(annual_percapita_nets_distributed = ifelse(ITNuse2==0, 0,
                                                    annual_percapita_nets_distributed),
         # calculate cost of interventions
         # true net cost accounting for non-linear relationship
         cost_ITN = population * annual_percapita_nets_distributed * sim_length/365 * ITNcost,
         # true net cost MIN
         cost_ITNmin = population * annual_percapita_nets_distmin * sim_length/365 * ITNcost,
         # true net cost MAX
         cost_ITNmax = population * annual_percapita_nets_distmax * sim_length/365 * ITNcost,
         # ITN linear
         cost_ITN_linear = population * ITNuse2 * bednet_timesteps * ITNcost,
         # non-severe treatment
         cost_clinical = ((cases - severe_cases - u5_cases) * treatment * TREATcost_adult) * .77 +
           ((u5_cases - u5_severe) * treatment * TREATcost_child) * .77,
         # severe treatment
         cost_severe = (severe_cases * treatment * SEVcost_adult) * .77 +
           (u5_severe * treatment * SEVcost_child) * .77,
         # SMC
         cost_SMC = n_91.25_1825 * SMC * SMCcost * smc_timesteps,
         # RTSS
         cost_vax = (dose1 + dose2 + dose3 + dose4) * (cost_per_dose + delivery_cost),

         # TOTAL
         cost_total = cost_ITN + cost_clinical + cost_severe + cost_SMC + cost_vax,
         cost_total_ITNmin = cost_ITNmin + cost_clinical + cost_severe + cost_SMC + cost_vax,
         cost_total_ITNmax = cost_ITNmax + cost_clinical + cost_severe + cost_SMC + cost_vax,

         # cost just among children
         cost_total_u5 =
           # cost ITNs
           n_0_1825 * annual_percapita_nets_distributed * sim_length/365 * ITNcost +
           # cost clinical
           ((u5_cases-u5_severe) * treatment * TREATcost_child)*.77 +
           # cost severe
           (u5_severe * treatment * SEVcost_child)*.77 +
           # cost SMC and cost VAX are all among children
           cost_SMC + cost_vax)

# assign scenarios
output <- dalyoutput_cost %>%
  filter(cost_per_dose == 6.52) %>%
  mutate(ID = paste(pfpr, seasonality, ITNuse, resistance, sep="_")) # create unique identifier

saveRDS(output, './03_output/scenarios_casestudy.rds')


# repeat exercise by grouping rural / urban together for analysis
# categorize grouped scenarios

# gap in PfPR
output_pfpr <- output %>%
  filter(pfpr %in% c(0.40, 0.10) & RTSScov %in% c(0, 0.80) & ITNuse == 0.50)

# gap in ITN use
output_itn <- output %>%
  filter(RTSScov %in% c(0, 0.80) & ((pfpr == 0.20 & ITNuse == 0.60) | (pfpr == 0.10 & ITNuse == 0.30)))

# gap in RTSS coverage
output_rtss <- output %>%
  filter(ITNuse == 0.50 & ((pfpr == 0.20 & RTSScov %in% c(0, 0.50)) | (pfpr == 0.10 & RTSScov %in% c(0, 0.80))))


assign_scenarios <- function(output){
# assign scenarios
scenario1 <- output %>%
  filter(ITNboost==1 & RTSS=='none') %>% mutate(scenario=1)

scenario2 <- output %>%
  filter(RTSS=='EPI' & ITNboost==0) %>% mutate(scenario=2)

scenario3 <- output %>%
  filter((pfpr %in% c(0.20, 0.40) & ITNboost==1 & RTSS=='none') | (pfpr %in% c(0.10) & ITNboost==0 & RTSS=='none')) %>%
  mutate(scenario=3)

scenario4 <- output %>%
  filter((pfpr %in% c(0.20, 0.40) & ITNboost==0 & RTSS=='EPI') | (pfpr %in% c(0.10) & ITNboost==0 & RTSS=='none')) %>%
  mutate(scenario=4)

scenarios <- full_join(scenario1, scenario2) %>% full_join(scenario3) %>% full_join(scenario4) %>%
  mutate(scenario_f = factor(scenario,
                             levels=c(1,2,3,4),
                             labels=c('mass ITN 10% increase', 'mass age-based RTS,S', 'targeted ITN 10% increase', 'targeted age-based RTS,S')))

none <- output %>%
  filter((seasonality %in% c('perennial', 'highly seasonal') & ITNboost==0 & RTSS=='none') |
           (seasonality == 'seasonal' & ITNboost==0 & RTSS=='none' & SMC==0)) %>%
  mutate(daly_baseline = daly,
         daly_lower_baseline = daly_lower,
         daly_upper_baseline = daly_upper,

         cases_baseline = cases,
         cases_lower_baseline = cases_lower,
         cases_upper_baseline = cases_upper,

         severe_baseline = severe_cases,

         deaths_baseline = deaths,
         deaths_lower_baseline = deaths_lower,
         deaths_upper_baseline = deaths_upper,

         u5_dalys_baseline = u5_dalys,
         u5_dalys_lower_baseline = u5_dalys_lower,
         u5_dalys_upper_baseline = u5_dalys_upper,

         u5_cases_baseline = u5_cases,
         u5_cases_lower_baseline = u5_cases_lower,
         u5_cases_upper_baseline = u5_cases_upper,

         u5_severe_baseline = u5_severe,

         u5_deaths_baseline = u5_deaths,
         u5_deaths_lower_baseline = u5_deaths_lower,
         u5_deaths_upper_baseline = u5_deaths_upper,

         n_0_1825_baseline = n_0_1825,
         n_baseline = n,

         cost_total_baseline = cost_total
  ) %>%
  select(file, ID, daly_baseline:cost_total_baseline)

base_IDs <- none$file

# sum measures and calculate CE by scenario
scenarios2 <- scenarios %>% filter(!(file %in% base_IDs)) %>%
  left_join(none %>% dplyr::select(-file), by=c('ID')) %>%
  group_by(scenario, scenario_f, seasonality) %>%
  summarize(across(c(cases:u5_dalys_upper,
                     cost_total,
                     daly_baseline:cost_total_baseline), ~sum(.x, na.rm=T))) %>%

  mutate(CE_daly = (cost_total - cost_total_baseline) / (daly_baseline - daly),
        CE_daly_lower = (cost_total - cost_total_baseline) / (daly_lower_baseline - daly_lower),
        CE_daly_upper = (cost_total - cost_total_baseline) / (daly_upper_baseline - daly_upper),

        CE_daly_u5 = (cost_total - cost_total_baseline) / (u5_dalys_baseline - u5_dalys),
        CE_daly_u5_lower = (cost_total - cost_total_baseline) / (u5_dalys_lower_baseline - u5_dalys_lower),
        CE_daly_u5_upper = (cost_total - cost_total_baseline) / (u5_dalys_upper_baseline - u5_dalys_upper),

        CE_case = (cost_total - cost_total_baseline) / (cases_baseline - cases),
        CE_case_lower = (cost_total - cost_total_baseline) / (cases_lower_baseline - cases_lower),
        CE_case_upper = (cost_total - cost_total_baseline) / (cases_upper_baseline - cases_upper),

        CE_u5_case = (cost_total - cost_total_baseline) / (u5_cases_baseline - u5_cases),
        CE_u5_case_lower = (cost_total - cost_total_baseline) / (u5_cases_lower_baseline - u5_cases_lower),
        CE_u5_case_upper = (cost_total - cost_total_baseline) / (u5_cases_upper_baseline - u5_cases_upper),

        CE_death = (cost_total - cost_total_baseline) / (deaths_baseline - deaths),
        CE_death_lower = (cost_total - cost_total_baseline) / (deaths_lower_baseline - deaths_lower),
        CE_death_upper = (cost_total - cost_total_baseline) / (deaths_upper_baseline - deaths_upper),

        CE_u5_death = (cost_total - cost_total_baseline) / (u5_deaths_baseline - u5_deaths),
        CE_u5_death_lower = (cost_total - cost_total_baseline) / (u5_deaths_lower_baseline - u5_deaths_lower),
        CE_u5_death_upper = (cost_total - cost_total_baseline) / (u5_deaths_upper_baseline - u5_deaths_upper)
)

}

test1 <- assign_scenarios(output_pfpr) %>% mutate(scenario2 = 1)
test2 <- assign_scenarios(output_itn) %>% mutate(scenario2 = 2)
test3 <- assign_scenarios(output_rtss) %>% mutate(scenario2 = 3)


scenarios2 <- full_join(test1, test2) %>% full_join(test3) %>%
  mutate(scenario2_f = factor(scenario2,
                             levels=c(1,2,3),
                             labels=c('gap in PfPR', 'gap in ITN use', 'gap in vaccination')))

saveRDS(scenarios2, './03_output/scenarios2_casestudy.rds')

