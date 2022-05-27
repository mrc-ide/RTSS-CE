# Cost effectiveness -----------------------------------------------------------
library(tidyverse)
# devtools::install_github("mrc-ide/netz@usage_to_npc") # preliminary version

source("./02_code/Parameter draws/HPC_processing.R")  # process raw data
source("./02_code/Parameter draws/deaths_dalys.R")    # add deaths and dalys
source("./02_code/Parameter draws/add_costs.R")       # add costs
source("./02_code/Parameter draws/outcome_averted.R") # calculate case / daly averted

# read in a list of all malariasimulation outputs
files <- list.files(path = "Q:/GF-RTSS-CE/03_output/HPC/", pattern = "general_*", full.names = TRUE)

# create file index
index <- c(1:length(files))

# function for processing HPC data
cost_effectiveness <- function(x){ # input = index of file to process

  # process HPC output: condense by age group over simulation time ----
  output <- HPC_processing(x)

  # run mortality and DALY functions ----
  output <- mortality_rate(output)
  output <- outcome_uncertainty(output)
  output <- daly_components(output)

  # condense into one line per run ----
  output <- output %>%
    select(-inc, -sev, -mortality_rate) %>% # get rid of rate vars

    # create vars for childhood cases
    mutate(n_0_1825 = ifelse(age %in% c('0-91.25', '91.25-1825'), n, 0), # u5 denominator
           n_91.25_1825 = ifelse(age=='91.25-1825', n, 0), # SMC denominator
           u5_cases = ifelse(age %in% c('0-91.25', '91.25-1825'), cases, 0),
           u5_severe = ifelse(age %in% c('0-91.25', '91.25-1825'), severe_cases, 0),
           u5_dalys = ifelse(age %in% c('0-91.25', '91.25-1825'), daly, 0)) %>%

    mutate_at(vars(n, n_0_1825, n_91.25_1825,
                   u5_cases, u5_severe, u5_dalys,
                   inc_clinical, inc_severe,
                   cases, cases_lower, cases_upper, severe_cases,
                   deaths, deaths_lower, deaths_upper,
                   yll:daly_lower), sum, na.rm=T) %>%  # condense outputs over all ages in population
    select(-age, -age_upper, -age_lower) %>%
    distinct()

  # add costs (intervention-specific and total) ----
  output <- add_costs(output)

  print(x)

  return(output)

}

# run cost_effectiveness function
dalyoutput <- map_dfr(index, cost_effectiveness)

# save output
saveRDS(dalyoutput, './03_output/test_dalyoutput.rds')

# calculate cases / DALYs averted
output <- outcome_averted(dalyoutput)

# save output
saveRDS(output, './03_output/test_scenarios.rds')



# data checks ##################################################################

# check prevalence values

none <- output %>%
  filter(ITNboost==0 & ITN=='pyr' & RTSS=='none' & # filter out interventions
           (SMC==0 | (seasonality=='highly seasonal'))) %>%
  mutate(PR = n_detect_730_3650 / n_730_3650) %>%
  mutate(PRdiff = pfpr - PR)

summary(test$PRdiff)

# check that there are no negative DALY values
if(any(output$yll < 0) | any(output$yld < 0)){
  stop("DALYs must be greater than 0")
}

# RTSS doses are stable
table(output$RTSS, output$RTSScov, output$rtss_mass_timesteps)

# SMC doses are stable
table(output$seasonality, output$smc_timesteps)

# ITNs are stable
table(output$ITNuse, output$bednet_timesteps)

# DALYs averted dist
test <- output %>% mutate(deltadaly = daly_baseline - daly)
summary(test$deltadaly)
test <- test %>% filter(deltadaly < 0 & resistance == 0)
table(test$resistance) # all negative DALYs are in resistance scenarios
test$deaths; test$deaths_baseline
test$cases; test$cases_baseline
test$severe_cases; test$severe_baseline


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
load('./01_data/unit_costs.rda')

PYRcost <- 3.50        # $2.00 per net and $1.50 delivery cost
TREATcost <- 1.47      # clinical treatment cost
SEVcost <- 22.41       # severe treatment cost
SMCcost <- unit_costs$cost_per_smc_dose_delivered  # 1.44
cost_per_dose <- c(2.69, 6.52, 12.91)
delivery_cost <- c(0.96, 1.62, 2.67)

# create combinations of dose cost and delivery cost
rtsscost_df <- expand_grid(cost_per_dose = cost_per_dose, delivery_cost = delivery_cost)

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
  mutate(
         # calculate cost of interventions
         cost_ITN = population * annual_percapita_nets_distributed * sim_length/365 * ITNcost,  # true net cost accounting for non-linear relationship
         cost_clinical = ((cases-severe_cases) * treatment * TREATcost)*.77, # non-severe treatment
         cost_severe = (severe_cases * treatment * SEVcost)*.77,             # severe treatment
         cost_vax = (dose1 + dose2 + dose3 + dose4) * (cost_per_dose + delivery_cost), # RTSS
         cost_SMC = n_91.25_1825 * SMC * SMCcost * smc_timesteps,                      # SMC

         cost_total = cost_ITN + cost_clinical + cost_severe + cost_vax + cost_SMC)

# assign scenarios
output <- dalyoutput_cost %>%
  filter(cost_per_dose==6.52 & delivery_cost==1.62) %>%
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

