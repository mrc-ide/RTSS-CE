# Tabular results

# set-up
source("./02_code/Figures/data_and_libraries.R")


# < the impact of RTSS on top of other interventions --------------
output <- scenarios %>%
  filter(cost_per_dose == 12.01 & delivery_cost == 1.62) %>%
  mutate(ID = paste(pfpr, seasonality, ITNuse, resistance, treatment, ITN, sep = "_")) %>%
  filter(ITNuse == 0.75)

none <- output %>%
  mutate(set = case_when(seasonality %in% c('perennial', 'highly seasonal') & intervention %in%
                           c('ITN 10% increase','ITN PBO') ~ 1,
                         seasonality == 'seasonal' & intervention %in%
                           c('ITN 10% increase + SMC','ITN PBO + SMC') ~ 1)) %>%
  filter(set == 1) %>%
  dplyr::select(file, ID, drawID, daly, cases, cost_total, u5_dalys, n_0_1825) %>%
  rename(daly_baseline = daly,
         cases_baseline = cases,
         cost_total_baseline = cost_total,
         u5_daly_baseline = u5_dalys) %>%
  dplyr::select(file, ID, drawID, daly_baseline, cases_baseline, cost_total_baseline, u5_daly_baseline)

base_IDs <- none$file

output2 <- output %>% filter(!(file %in% base_IDs)) %>%
  mutate(set = case_when(seasonality %in% c('perennial', 'highly seasonal') & intervention %in%
                           c('ITN 10% increase + RTS,S','ITN PBO + RTS,S') ~ 1,
                         seasonality=='seasonal' & intervention %in%
                           c('ITN 10% increase + RTS,S + SMC','ITN PBO + RTS,S + SMC') ~ 1)) %>%
  filter(set == 1) %>%
  dplyr::select(file, ID, drawID, pfpr, seasonality, intervention, daly, cases, cost_total, u5_dalys, dose3) %>%
  left_join(none %>% dplyr::select(-file), by=c('ID', 'drawID')) %>%
  mutate(CE = (cost_total - cost_total_baseline) / (daly_baseline - daly),
         deltadaly = daly_baseline - daly,
         deltacases = cases_baseline - cases,
         CE_u5 = (cost_total - cost_total_baseline) / (u5_daly_baseline - u5_dalys))


# adjusting for 15 year simulation period and 200,000 population arguments
summary(output2$deltadaly / (2*15)) # additional dalys averted per year in a population of 100,000 people
summary(output2$deltacases / (2*15)) # additional cases averted per year in a population of 100,000 people
summary(output2$deltacases / output2$dose3 * 100000) # additional cases averted per year per fully vaccinated child (dose3)
summary(output2$CE_u5) # additional cases averted per year in a population of 100,000 people


# < ICER table ----
# calculate change in dalys and cost
output <- scenarios %>%
  filter(cost_per_dose==6.52 & delivery_cost==1.62) %>% filter(resistance==0) %>%
  filter(intervention!='none') %>%
  mutate(deltadaly = daly_baseline - daly,
         deltacost = cost_total - cost_total_baseline,
         cost_daly_averted = (cost_total - cost_total_baseline) / deltadaly) %>%
  mutate(seasonality = factor(seasonality, levels = c('perennial', 'seasonal', 'highly seasonal'))) %>%
  group_by(ID, drawID, seasonality, intervention, intervention_f) %>%
  summarize(deltadaly = median(deltadaly),
            deltacost = median(deltacost),
            cost_daly_averted = median(cost_daly_averted))


final <- output %>%
  # filter out mixed strategies
  group_by(ID, drawID) %>% arrange(ID, drawID, deltacost) %>%
  filter(!(deltadaly < 0 & deltacost > 0)) %>%
  # filter out dominated strategies
  mutate(dominate = case_when(deltadaly < lag(deltadaly,n=12L) ~ 1,
                              deltadaly < lag(deltadaly,n=11L) ~ 1,
                              deltadaly < lag(deltadaly,n=10L) ~ 1,
                              deltadaly < lag(deltadaly,n=9L) ~ 1,
                              deltadaly < lag(deltadaly,n=8L) ~ 1,
                              deltadaly < lag(deltadaly,n=7L) ~ 1,
                              deltadaly < lag(deltadaly,n=6L) ~ 1,
                              deltadaly < lag(deltadaly,n=5L) ~ 1,
                              deltadaly < lag(deltadaly,n=4L) ~ 1,
                              deltadaly < lag(deltadaly,n=3L) ~ 1,
                              deltadaly < lag(deltadaly,n=2L) ~ 1,
                              deltadaly < lag(deltadaly,n=1L) ~ 1,
                              TRUE ~ 0)) %>%
  filter(dominate==0) %>%
  # marking extended dominated strategies
  mutate(ICER = (deltacost-lag(deltacost)) / (deltadaly-lag(deltadaly)),
         dominate = case_when(ICER > lead(ICER,n=12L) ~ 1,
                              ICER > lead(ICER,n=11L) ~ 1,
                              ICER > lead(ICER,n=10L) ~ 1,
                              ICER > lead(ICER,n=9L) ~ 1,
                              ICER > lead(ICER,n=8L) ~ 1,
                              ICER > lead(ICER,n=7L) ~ 1,
                              ICER > lead(ICER,n=6L) ~ 1,
                              ICER > lead(ICER,n=5L) ~ 1,
                              ICER > lead(ICER,n=4L) ~ 1,
                              ICER > lead(ICER,n=3L) ~ 1,
                              ICER > lead(ICER,n=2L) ~ 1,
                              ICER > lead(ICER,n=1L) ~ 1,
                              TRUE ~ 0)) %>%
  filter(dominate==0) %>%
  # loop through extendedly dominated strategies again
  mutate(ICER = (deltacost-lag(deltacost)) / (deltadaly-lag(deltadaly)),
         dominate = case_when(ICER > lead(ICER,n=12L) ~ 1,
                              ICER > lead(ICER,n=11L) ~ 1,
                              ICER > lead(ICER,n=10L) ~ 1,
                              ICER > lead(ICER,n=9L) ~ 1,
                              ICER > lead(ICER,n=8L) ~ 1,
                              ICER > lead(ICER,n=7L) ~ 1,
                              ICER > lead(ICER,n=6L) ~ 1,
                              ICER > lead(ICER,n=5L) ~ 1,
                              ICER > lead(ICER,n=4L) ~ 1,
                              ICER > lead(ICER,n=3L) ~ 1,
                              ICER > lead(ICER,n=2L) ~ 1,
                              ICER > lead(ICER,n=1L) ~ 1,
                              TRUE ~ 0)) %>%
  filter(dominate==0) %>%
  mutate(ICER = (deltacost-lag(deltacost)) / (deltadaly-lag(deltadaly)),
         dominate = case_when(ICER > lead(ICER,n=12L) ~ 1,
                              ICER > lead(ICER,n=11L) ~ 1,
                              ICER > lead(ICER,n=10L) ~ 1,
                              ICER > lead(ICER,n=9L) ~ 1,
                              ICER > lead(ICER,n=8L) ~ 1,
                              ICER > lead(ICER,n=7L) ~ 1,
                              ICER > lead(ICER,n=6L) ~ 1,
                              ICER > lead(ICER,n=5L) ~ 1,
                              ICER > lead(ICER,n=4L) ~ 1,
                              ICER > lead(ICER,n=3L) ~ 1,
                              ICER > lead(ICER,n=2L) ~ 1,
                              ICER > lead(ICER,n=1L) ~ 1,
                              TRUE ~ 0)) %>%
  filter(dominate==0) %>%
  mutate(ICER = (deltacost-lag(deltacost)) / (deltadaly-lag(deltadaly)),
         dominate = case_when(ICER > lead(ICER,n=12L) ~ 1,
                              ICER > lead(ICER,n=11L) ~ 1,
                              ICER > lead(ICER,n=10L) ~ 1,
                              ICER > lead(ICER,n=9L) ~ 1,
                              ICER > lead(ICER,n=8L) ~ 1,
                              ICER > lead(ICER,n=7L) ~ 1,
                              ICER > lead(ICER,n=6L) ~ 1,
                              ICER > lead(ICER,n=5L) ~ 1,
                              ICER > lead(ICER,n=4L) ~ 1,
                              ICER > lead(ICER,n=3L) ~ 1,
                              ICER > lead(ICER,n=2L) ~ 1,
                              ICER > lead(ICER,n=1L) ~ 1,
                              TRUE ~ 0)) %>%
  filter(dominate==0) %>%
  mutate(ICER = (deltacost-lag(deltacost)) / (deltadaly-lag(deltadaly)),
         dominate = case_when(ICER > lead(ICER,n=12L) ~ 1,
                              ICER > lead(ICER,n=11L) ~ 1,
                              ICER > lead(ICER,n=10L) ~ 1,
                              ICER > lead(ICER,n=9L) ~ 1,
                              ICER > lead(ICER,n=8L) ~ 1,
                              ICER > lead(ICER,n=7L) ~ 1,
                              ICER > lead(ICER,n=6L) ~ 1,
                              ICER > lead(ICER,n=5L) ~ 1,
                              ICER > lead(ICER,n=4L) ~ 1,
                              ICER > lead(ICER,n=3L) ~ 1,
                              ICER > lead(ICER,n=2L) ~ 1,
                              ICER > lead(ICER,n=1L) ~ 1,
                              TRUE ~ 0)) %>%
  filter(dominate==0) %>%
  mutate(ICER = ifelse(is.na(ICER), (deltacost) / (deltadaly), ICER),
         dominate = 0) %>%
  dplyr::select(ID, drawID, intervention, ICER, dominate)

merge <- output %>% left_join(final, by=c('ID', 'drawID', 'intervention')) %>%
  mutate(dominate = ifelse(is.na(dominate), 1, dominate))

# group results by intervention
merge %>% group_by(intervention_f) %>%
  summarize(t = n(),
            ndominate = n()-sum(dominate),
            p_ndominate = ndominate / t*100,
            ICER_m = median(ICER, na.rm=T),
            ICER_25 = quantile(ICER, prob=0.25, na.rm=T),
            ICER_75 = quantile(ICER, prob=0.75, na.rm=T))

# group results by seasonality and intervention
merge %>% group_by(seasonality, intervention_f) %>%
  summarize(t = n(),
            ndominate = n()-sum(dominate),
            p_ndominate = ndominate / t*100,
            ICER_m = median(ICER, na.rm=T),
            ICER_25 = quantile(ICER, prob=0.25, na.rm=T),
            ICER_75 = quantile(ICER, prob=0.75, na.rm=T)) %>%
  write.table("clipboard", sep="\t", row.names=FALSE, col.names=FALSE)


# ICER just among non-dominated strategies
merge %>% filter(dominate==0) %>% group_by(intervention) %>%
  summarize(ICER_m = median(ICER, na.rm=T),
            ICER_25 = quantile(ICER, prob=0.25, na.rm=T),
            ICER_75 = quantile(ICER, prob=0.75, na.rm=T))



# < ICERs among children ----
scenarios %>% ungroup() %>%
  #group_by(seasonality) %>%
  filter(intervention %in% c('RTS,S age-based', 'RTS,S seasonal')) %>%
  summarize(n = n(),
            median = median(CE_u5, na.rm=T),
            q25 = quantile(CE_u5, prob=0.25, na.rm=T),
            q75 = quantile(CE_u5, prob=0.75, na.rm=T))

scenarios %>% ungroup() %>% filter(resistance==0) %>%
  group_by(ITNuse) %>%
  summarize(n = n(),
            median = median(CE_u5, na.rm=T),
            q25 = quantile(CE_u5, prob=0.25, na.rm=T),
            q75 = quantile(CE_u5, prob=0.75, na.rm=T))

scenarios %>% ungroup() %>% filter(resistance==0) %>%
  group_by(ITNuse) %>%
  summarize(n = n(),
            median = median(CE, na.rm=T),
            q25 = quantile(CE, prob=0.25, na.rm=T),
            q75 = quantile(CE, prob=0.75, na.rm=T))

