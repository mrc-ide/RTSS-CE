# input: data frame
# process: reads in dataframe, calculates cost-effectiveness of interventions compared to baseline
# output: data frame with additional variables for cost-effectiveness per case or DALY averted

outcome_averted <- function(x # dataframe to read in and process
                            ){

  output <- x

  # separate out baseline scenarios
  # there should be 270 baseline scenarios. 3 pfpr x 3 seasonality x 4 ITN usage x 3 treatment x 2.5 resistance (only pbo in resistance scenarios)
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
           cost_total_u5_baseline = cost_total_u5) %>%

    select(file, ID, scenario, drawID, cost_per_dose, daly_baseline, cases_baseline,
           severe_baseline, deaths_baseline, u5_dalys_baseline,
           u5_cases_baseline, u5_severe_baseline, cost_total_baseline,
           cost_total_ITNmin_baseline, cost_total_ITNmax_baseline,
           cost_total_u5_baseline)

  # separate out non baseline scenarios and merge
  base_IDs <- none$file

  scenarios <- output %>% filter(!(file %in% base_IDs)) %>%
    left_join(none %>% select(-file, -scenario), by=c('ID', 'drawID', 'cost_per_dose')) %>%

    # calculate cost-effectiveness
    mutate(CE = (cost_total - cost_total_baseline) / (daly_baseline - daly),
           CE_u5 = (cost_total - cost_total_baseline) / (u5_dalys_baseline - u5_dalys),
           CE_ITNmin = (cost_total_ITNmin - cost_total_ITNmin_baseline) / (daly_baseline - daly),
           CE_ITNmax = (cost_total_ITNmax - cost_total_ITNmax_baseline) / (daly_baseline - daly),
           CE_case = (cost_total - cost_total_baseline) / (cases_baseline - cases),
           CE_u5_case = (cost_total - cost_total_baseline) / (u5_cases_baseline - u5_cases)
    ) %>%

    # assign scenarios
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

  return(scenarios)

}
