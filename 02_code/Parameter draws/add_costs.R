# input: data frame
# process: reads in dataframe, calculates cost of interventions and cost-effectiveness
# output: data frame with additional variables for cost and cost-effectiveness

add_costs <- function(x # dataframe to read in and process
                      ){
  # COST REFERENCES
  # https://github.com/mrc-ide/gf/blob/master/data/unit_costs.rda
  # https://github.com/mrc-ide/gf/blob/master/data/treatment_unit_costs.rda
  # load('./01_data/unit_costs.rda')
  # ITNs: https://www.thelancet.com/journals/lanplh/article/PIIS2542-5196(21)00296-5/fulltext

  PYRcost <- 3.50        # $2.00 per net and $1.50 delivery cost
  PBOcost <- 3.80        # $2.30 per net and $1.50 delivery cost
  TREATcost <- 1.47      # clinical treatment cost
  SEVcost <- 22.41       # severe treatment cost
  SMCcost <- 1.44        # 1.44 per dose
  cost_per_dose <- c(2.69, 6.52, 12.91) # cost per dose RTS,S
  delivery_cost <- 1.62  # delivery cost options RTS,S c(0.96, 1.62, 2.67)

  # create combinations of dose cost and delivery cost
  rtsscost_df <- expand_grid(cost_per_dose = cost_per_dose, delivery_cost = delivery_cost)

  population <- x$population[1]
  sim_length <- x$sim_length[1]

  # Prepare to add costs to dataset
  dalyoutput_cost <- x %>%

    mutate(ITNuse2 = ifelse(ITNboost == 1, ITNuse + .10, ITNuse), # account for booster coverage
           ITNcost = case_when(ITN == 'pyr' ~ PYRcost,            # account for ITN type-specific cost
                               ITN == 'pbo' ~ PBOcost)) %>%

    # count the number of interventions administered and the frequency of ITN dist
    mutate(bednet_distribution_frequency = as.numeric(lapply(lapply(bednet_timesteps,diff), unique)),
           bednet_timesteps = length(Filter(function(x) (x > 0 & x <= sim_length), bednet_timesteps)),
           smc_timesteps = length(Filter(function(x) (x > 0 & x <= sim_length), smc_timesteps)),
           rtss_mass_timesteps = length(Filter(function(x) (x > 0 & x<= sim_length), rtss_mass_timesteps))) %>%

    # merge in RTSS costing dataframe
    merge(rtsscost_df) %>%

    ungroup() %>% rowwise()


  # read in net distribution conversions, calculated through 'netz_dist.R'
  nets_distributed <- readRDS('./03_output/netz_data')

  # create cost variables
  # 77% of treatment costs are from the public sector (DHS, SSA)
  dalyoutput_cost <- dalyoutput_cost %>%
    left_join(nets_distributed,
              by=c('ITNuse2' = 'target_use')) %>%
    mutate(annual_percapita_nets_distributed = ifelse(ITNuse2==0, 0,
                                                      annual_percapita_nets_distributed),

           # calculate cost of interventions
           cost_ITN = population * annual_percapita_nets_distributed * sim_length/365 * ITNcost,  # true net cost accounting for non-linear relationship
           cost_ITNmin = population * annual_percapita_nets_distmin * sim_length/365 * ITNcost,  # true net cost MIN
           cost_ITNmax = population * annual_percapita_nets_distmax * sim_length/365 * ITNcost,  # true net cost MAX
           cost_ITN_linear = population * ITNuse2 * bednet_timesteps * ITNcost,          # ITN linear
           cost_clinical = ((cases-severe_cases) * treatment * TREATcost)*.77, # non-severe treatment
           cost_severe = (severe_cases * treatment * SEVcost)*.77,             # severe treatment
           cost_SMC = n_91.25_1825 * SMC * SMCcost * smc_timesteps,                      # SMC
           cost_vax = (dose1 + dose2 + dose3 + dose4) * (cost_per_dose + delivery_cost), # RTSS

           cost_total = cost_ITN + cost_clinical + cost_severe + cost_SMC + cost_vax, # TOTAL
           cost_total_ITNmin = cost_ITNmin + cost_clinical + cost_severe + cost_SMC + cost_vax,
           cost_total_ITNmax = cost_ITNmax + cost_clinical + cost_severe + cost_SMC + cost_vax,
           # cost just among children
           cost_total_u5 = n_0_1825 * annual_percapita_nets_distributed * sim_length/365 * ITNcost + # cost ITNs
             ((u5_cases-u5_severe) * treatment * TREATcost)*.77 + # cost clinical
             (u5_severe * treatment * SEVcost)*.77 +  # cost severe
             cost_SMC + cost_vax) # cost SMC and cost VAX are all among children

  # NOTE that values are NA for the target rates 0.75, 0.85 with the min SSA rate #

  return(dalyoutput_cost)

}
