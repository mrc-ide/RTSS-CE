# input: data frame
# process: reads in dataframe, calculates cost of interventions and cost-effectiveness
# output: data frame with additional variables for cost and cost-effectiveness

add_costs <- function(x # dataframe to read in and process
                      ){
  # COST REFERENCES
  # references in: https://mrc-ide.github.io/treasure/reference/index.html
  # treatment costs from Penny et al. 2016

  PYRcost <- 2.52 + 1.50                # pyrethroid $2.52 per net and $1.50 delivery cost
  PBOcost <- 3.51 + 1.50                # PBO $3.51 per net and $1.50 delivery cost

  SMCcost <- 0.9075                     # SMC $0.9075 per dose

  cost_per_dose <- c(2.69, 6.52, 12.91, 17.36) # RTS,S per dose $2, $5, $10 + consumables cost
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
  # NOTE that values are NA for the target rates 0.75, 0.85 with the min SSA rate #

  return(dalyoutput_cost)

}
