# Figures & Tables -------------------------------------------------------------
# packages
library(tidyverse)
library(fuzzyjoin)
library(plotly)
library(kableExtra)
library(malariasimulation)
library(data.table)

library(patchwork)
library(scales)
library(LaCroixColoR)

# devtools::install_github('mrc-ide/malariasimulation@dev', force=TRUE)
# devtools::install_github('johannesbjork/LaCroixColoR')

# Look here: https://github.com/htopazian/rtss_malariasimulation
# find ideas in rmarkdowns for figures and insert code if helpful

dat <- readRDS("./03_output/rtss_raw.rds")
averted <- readRDS("./03_output/rtss_avert.rds")
dalyoutput <- readRDS('./03_output/dalyoutput.rds')
dalyoutput_cost <- readRDS('./03_output/dalyoutput_cost.rds')



# seasonality ------------------------------------------------------------------

# get starting parameters
seas_name <- 'highly seasonal'
seasonality <- list(c(0.284596,-0.317878,-0.0017527,0.116455,-0.331361,0.293128,-0.0617547))
s1 <- crossing(seasonality, seas_name)

seas_name <- 'seasonal'
seasonality <- list(c(0.285505,-0.325352,-0.0109352,0.0779865,-0.132815,0.104675,-0.013919))
s2 <- crossing(seasonality, seas_name)

stable <- rbind(s1, s2)

find_peak <- function(seasonality, seas_name){
  month <- 365/12

  params <- get_parameters(list(
    human_population = 1000,
    model_seasonality = TRUE,
    g0 = unlist(seasonality)[1],
    g = unlist(seasonality)[2:4],
    h = unlist(seasonality)[5:7],
    individual_mosquitoes = FALSE))

  peak <- peak_season_offset(params)
  smc <- crossing('smc', round(c(peak+c(-1.5,-0.5,0.5,1.5)*month),0))
  colnames(smc) <- c('intervention', 'month')

  first <- round((peak-month*3.5),0)
  sv <- crossing('sv', first + round(c(0,1*month,2*month)))
  colnames(sv) <- c('intervention', 'month')

  interventions <- rbind(smc, sv) %>%
    mutate(seasonality = seas_name,
           month = month/365*12,
           color = ifelse(intervention=='smc', 'red', 'blue'))

  return(interventions)
}

output <- map2_dfr(stable$seasonality, stable$seas_name, find_peak)

baseline <- dat %>%
  filter(ITNuse==0 &  RTSS=='none' &  ITN=='pyr' &
           seasonality != 'perennial' & ITNboost == 0 & (seasonality=="seasonal" & SMC==0 | seasonality=="highly seasonal" & SMC==0.85)) %>%
  mutate(inc_month=n_inc_clinical_0_1825/n_0_1825, na.rm=T) %>%
  filter(month %in% seq(1,12,1))

ggplot(data=baseline) +
  geom_line(data=baseline, aes(x=month, y=inc_month, color=as.factor(pfpr))) +
  geom_vline(data=output[output$intervention=='smc',], aes(xintercept=month, alpha='SMC'), lty=2, color='red') +
  geom_vline(data=output[output$intervention=='sv',], aes(xintercept=month, alpha='RTS,S'), lty=2, color='blue') +
  labs(title='distribution of SMC and seasonal RTS,S; ages 0-5 years',
       caption='highly seasonal has SMC implemented at baseline',
       x="Timesteps (month)",
       y="Monthly clinical incidence",
       color = 'PfPR') +
  scale_y_continuous(limits=c(0,0.8), breaks=seq(0,0.8,0.2)) +
  scale_x_continuous(limits=c(1,12), breaks=seq(1,12,1)) +
  facet_wrap('seasonality', ncol = 1) +
  theme_classic() +
  scale_alpha_manual(values = c(rep(1,2))) +
  guides(alpha = guide_legend(title = 'intervention',
                               override.aes = list(color = c('blue','red'))))

ggsave('./03_output/seasonality.pdf', width=6, height=6)


# seasonality II ---------------------------------------------------------------

output <- dat %>% filter(SMC == 0.85 & RTSS == 'SV' & ITN == 'pyr' & ITNuse == 0.5 & resistance == 0 & ITNboost == 0 )%>%
  filter(seasonality == 'highly seasonal' & pfpr == 0.4)

# look at one file
interventiontime <- output  %>%
  select(smc_timesteps, rtss_mass_timesteps, bednet_timesteps) %>% distinct()

SMCtime <- do.call(cbind.data.frame, interventiontime$smc_timesteps)
colnames(SMCtime) <- "smc"

RTSStime <- do.call(cbind.data.frame, interventiontime$rtss_mass_timesteps)
colnames(RTSStime) <- "rtss"

ITNtime <- do.call(cbind.data.frame, interventiontime$bednet_timesteps)
colnames(ITNtime) <- "itn"


ggplot(data = output) +
  geom_line(aes(x=month, y=p_inc_clinical_0_1825/n_0_1825), alpha = 0.8) +
  geom_vline(data = SMCtime, aes(xintercept=smc/(365/12), alpha='SMC'), lty=2, color='red') +
  geom_vline(data = RTSStime, aes(xintercept=rtss/(365/12), alpha='RTS,S dose 3'), lty=2, color='blue') +
  geom_vline(data = ITNtime, aes(xintercept=itn/(365/12), alpha='ITN'), lty=2, color='green') +
  labs(x='month', y='clinical incidence (month), 0-5 years', subtitle = "highly seasonal, PfPR 0.4") +
  scale_x_continuous(breaks = seq(0,20*12,12)) +
  coord_cartesian(xlim = c(0,20*12)) +
  scale_alpha_manual(values = c(rep(1,3))) +
  guides(alpha = guide_legend(title = 'intervention',
                              override.aes = list(color = c('green','blue','red')))) +
  theme_classic()

ggsave('./03_output/seasonalityII.pdf', width=10, height=4)

# seasonality III --------------------------------------------------------------
# one year zoomed
output <- dat %>% filter(SMC == 0.85 & RTSS == 'SV' & ITN == 'pyr' & ITNuse == 0.5 & resistance == 0 & ITNboost == 0 )

SMCtime <- output %>% select(smc_timesteps, seasonality, pfpr) %>%
  group_by(smc_timesteps, seasonality, pfpr) %>%
  mutate(t1 = unlist(smc_timesteps)[[21]],
         t2 = unlist(smc_timesteps)[[22]],
         t3 = unlist(smc_timesteps)[[23]],
         t4 = unlist(smc_timesteps)[[24]]) %>%
  distinct() %>%
  pivot_longer(cols=t1:t4, names_to = "time", values_to = "smc") %>%
  mutate(smc = case_when(seasonality=="seasonal" ~ smc - 5*365,
                         TRUE ~ smc),
         smc = smc / (365/12) + 1)

RTSStime <- output %>% select(rtss_mass_timesteps, seasonality, pfpr) %>%
  group_by(rtss_mass_timesteps, seasonality, pfpr) %>%
  mutate(rtss = unlist(rtss_mass_timesteps)[[1]]) %>%
  distinct() %>%
  mutate(rtss = rtss / (365/12) + 1)

ITNtime <- output %>% select(bednet_timesteps, seasonality, pfpr) %>%
  group_by(bednet_timesteps, seasonality, pfpr) %>%
  mutate(itn = unlist(bednet_timesteps)[[6]]) %>%
  distinct() %>%
  mutate(itn = itn / (365/12) + 1)


ggplot(data = output) +
  geom_line(aes(x=month, y=p_inc_clinical_0_1825/n_0_1825, color=as.factor(pfpr)), alpha = 0.8) +
  geom_vline(data = SMCtime, aes(xintercept=smc, alpha='SMC'), lty=2, color='red') +
  geom_vline(data = RTSStime, aes(xintercept=rtss, alpha='RTS,S dose 3'), lty=2, color='blue') +
  geom_vline(data = ITNtime, aes(xintercept=itn, alpha='ITN'), lty=2, color='green') +
  labs(x='month', y='clinical incidence (month), 0-5 years', color='PfPR') +
  scale_x_continuous(breaks = seq(1,12,1)) +
  facet_wrap(seasonality~pfpr) +
  coord_cartesian(xlim = c(0,12)) +
  scale_alpha_manual(values = c(rep(1,3))) +
  guides(alpha = guide_legend(title = 'intervention',
                              override.aes = list(color = c('green','blue','red')))) +
  theme_classic()

ggsave('./03_output/seasonalityIII.pdf', width=10, height=5)


# bed nets and resistance ------------------------------------------------------

ITNdata <- dat %>%
  filter(RTSS=='none' & SMC==0) %>%
  filter(ITNuse==.50) %>%
  mutate(inc_month=n_inc_clinical_0_1825/n_0_1825, na.rm=T) %>%
  filter(month %in% seq(1,12,1))

ggplot(data=ITNdata) +
  geom_line(aes(x=month, y=inc_month, color=as.factor(interaction(ITN, resistance)))) +
  labs(title='distribution of SMC and seasonal RTS,S, ages 0-5 years',
       x="Timesteps (month)",
       y="Monthly clinical incidence",
       color = 'PfPR') +
  scale_y_continuous(limits=c(0,0.5), breaks=seq(0,0.5,0.1)) +
  scale_x_continuous(limits=c(1,12), breaks=seq(1,12,1)) +
  facet_wrap(c('seasonality','pfpr'), ncol = 3, nrow = 3) +
  theme_classic() +
  scale_alpha_manual(values = c(rep(1,2))) +
  guides(alpha = guide_legend(title = 'intervention',
                              override.aes = list(color = c('blue','red'))))



# delta change -----------------------------------------------------------------
# cost_per_dose <- c(2.69, 6.52, 12.91)
# delivery_cost <- c(0.96, 1.62, 2.67)

output <- dalyoutput_cost %>%
  filter(cost_per_dose==6.52 & delivery_cost==1.62)

none <- output %>%
  filter(ITNboost==0 & RTSS=='none' & ITN=='pyr' & resistance==0 & (SMC==0 | (seasonality=='highly seasonal'))) %>%
  select(seasonality, pfpr, ITNuse, daly, cost_total) %>%
  rename(daly_baseline = daly,
         cost__total_baseline = cost_total)

output <- output %>% left_join(none, by=c('seasonality', 'pfpr', 'ITNuse')) %>%
  mutate(deltadaly = abs(daly_baseline - daly),
         deltacost = cost_total - cost__total_baseline,
         cost_daly_averted = (cost_total - cost__total_baseline) / deltadaly)

ITNboost <- output %>% filter(RTSS=='none') %>% filter(!(seasonality=='seasonal' & SMC==.85)) %>% filter(ITNboost==1) %>% filter(resistance==0)

ITNpbo <- output %>% filter(RTSS=='none') %>% filter(!(seasonality=='seasonal' & SMC==.85)) %>% filter(ITN=='pbo') %>% filter(resistance==0)

SMC <- output %>% filter(RTSS=='none') %>% filter(ITNboost==0 & ITN=='pyr') %>%
  filter(seasonality!='highly seasonal') %>% filter(resistance==0) %>% filter(SMC==.85)

RTSSSV <- output %>% filter(RTSS=='SV') %>% filter(ITNboost!=1 & ITN!='pbo') %>%
  filter(!(seasonality=='seasonal' & SMC==.85)) %>% filter(resistance==0)

RTSSEPI <- output %>% filter(RTSS=='EPI') %>% filter(ITNboost!=1 & ITN!='pbo') %>%
  filter(!(seasonality=='seasonal' & SMC==.85)) %>% filter(resistance==0)

ggplot(mapping=aes(x=deltacost, y=deltadaly)) +
  geom_point(data=ITNboost, aes(color='ITN boost (10%)', alpha=ITNuse, shape=factor(pfpr))) +
  geom_point(data=ITNpbo, aes(color='ITN PBO', alpha=ITNuse, shape=factor(pfpr))) +
  geom_point(data=SMC, aes(color='SMC alone', alpha=ITNuse, shape=factor(pfpr))) +
  geom_point(data=RTSSSV, aes(color='RTS,S SV alone ($5)', alpha=ITNuse, shape=factor(pfpr))) +
  geom_point(data=RTSSEPI, aes(color='RTS,S EPI alone ($5)', alpha=ITNuse, shape=factor(pfpr))) +
  geom_hline(yintercept = 0, lty=2, color="black") +
  geom_vline(xintercept = 0, lty=2, color="black") +
  facet_wrap(~seasonality) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x='change in cost (USD)',
       y='change in DALYs',
       color='intervention',
       title='intervention impact',
       subtitle='matched on seasonality, PfPR, and ITN usage',
       caption='Assuming resistance == 0, SMC implemented in seasonal settings')

ggsave('./03_output/impact_cloud.pdf', width=8, height=5)


# resistance ITNs --------------------------------------------------------------
output <- dalyoutput_cost %>%
  filter(RTSS=='none' & (SMC==0 | (seasonality=='highly seasonal'))) %>%
  filter(cost_per_dose==6.52 & delivery_cost==1.62)

none <- output %>%
  filter(ITNboost==0 & ITN=='pyr' & resistance==0) %>%
  select(seasonality, pfpr, ITNuse, daly, cost_total) %>%
  rename(daly_baseline = daly,
         cost__total_baseline = cost_total)

output <- output %>% left_join(none, by=c('seasonality', 'pfpr', 'ITNuse')) %>%
  mutate(deltadaly = abs(daly_baseline - daly),
         deltacost = cost_total - cost__total_baseline,
         cost_daly_averted = (cost_total - cost__total_baseline) / deltadaly)

ITNboost <- output %>% filter(ITNboost==1)

ITNpbo <- output %>% filter(ITN=='pbo')

ggplot(mapping=aes(x=deltacost, y=deltadaly)) +
  geom_point(data=ITNboost, aes(color='ITN boost', alpha=pfpr)) +
  geom_point(data=ITNpbo, aes(color='ITN PBO', alpha=pfpr)) +
  geom_hline(yintercept = 0, lty=2, color="black") +
  geom_vline(xintercept = 0, lty=2, color="black") +
  theme_classic() +
  facet_wrap(~resistance) +
  labs(x='change in cost (USD)',
       y='change in DALYs',
       color='intervention',
       alpha='PfPR',
       title='intervention impact',
       subtitle='matched on seasonality, PfPR, and ITN usage',
       caption='')

table(output$ITNboost, output$ITN) # check ITNboost and ITN type relationship
table(output$seasonality, output$pfpr)

test <- output %>% filter(!(ITN=='pyr' & ITNboost==0)) %>% filter(ITNuse>0) %>%
  mutate(ITN = ifelse(ITN=='pyr', "pyrethroid boost (10%)", "pyrethroid + PBO"),
         scenario = case_when(ITNuse==0.25 & resistance==0 ~ "ITN 0.25 - none",
                              ITNuse==0.25 & resistance==.4 ~ "ITN 0.25 - low",
                              ITNuse==0.25 & resistance==.8 ~ "ITN 0.25 - high",
                              ITNuse==0.50 & resistance==0 ~ "ITN 0.50 - none",
                              ITNuse==0.50 & resistance==.4 ~ "ITN 0.50 - low",
                              ITNuse==0.50 & resistance==.8 ~ "ITN 0.50 - high",
                              ITNuse==0.75 & resistance==0 ~ "ITN 0.75 - none",
                              ITNuse==0.75 & resistance==.4 ~ "ITN 0.75 - low",
                              ITNuse==0.75 & resistance==.8 ~ "ITN 0.75 - high"))

table(test$ITNuse, test$resistance, test$scenario)

ggplot(data = test, aes(x=factor(scenario), y=deltadaly, alpha=factor(resistance))) +
  geom_col(aes(fill=ITN), position=position_dodge()) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  facet_wrap(seasonality~pfpr) +
  labs(x='scenario: ITNuse - resistance',
       y='change in DALYs',
       fill='ITN',
       alpha = 'resistance',
       title='Insecticide resistance',
       subtitle='matched on seasonality, PfPR, and ITN usage',
       caption='faceted by seasonality and PfPR')

ggsave('./03_output/resistance.pdf', width=8, height=5)


# cost per DALY averted --------------------------------------------------------
output <- dalyoutput_cost %>%
  filter(cost_per_dose==6.52 & delivery_cost==1.62)

none <- output %>%
  filter(ITNboost==0 & RTSS=='none' & ITN=='pyr' & resistance==0 & (SMC==0 | (seasonality=='highly seasonal'))) %>%
  select(seasonality, pfpr, ITNuse, daly, cost_total) %>%
  rename(daly_baseline = daly,
         cost__total_baseline = cost_total)

output <- output %>% left_join(none, by=c('seasonality', 'pfpr', 'ITNuse')) %>%
  mutate(deltadaly = abs(daly_baseline - daly),
         deltacost = cost_total - cost__total_baseline,
         cost_daly_averted = (cost_total - cost__total_baseline) / deltadaly)

ggplot(mapping=aes(x=deltacost, y=deltadaly)) +
  geom_point(data=ITNboost, aes(color='ITN boost')) +
  geom_point(data=ITNpbo, aes(color='ITN PBO')) +
  geom_point(data=SMC, aes(color='SMC alone')) +
  geom_point(data=RTSS, aes(color='RTS,S alone')) +
  geom_hline(yintercept = 0, lty=2, color="black") +
  geom_vline(xintercept = 0, lty=2, color="black") +
  theme_classic() +
  labs(x='change in cost (USD)',
       y='change in DALYs',
       color='intervention',
       title='intervention impact',
       subtitle='matched on seasonality, PfPR, and ITN usage',
       caption='Assuming resistance == 0, SMC implemented in seasonal settings')

ggsave('./03_output/daly_averted.pdf', width=8, height=5)


# RTS,S cost -------------------------------------------------------------------
output <- dalyoutput_cost %>%
  filter(delivery_cost==1.62)

# cost_per_dose <- c(2.69, 6.52, 12.91)
# delivery_cost <- c(0.96, 1.62, 2.67)



# IDEAS ###########
- resistance - ITNs, upgrading to PBO and boosting. no SMC add ins or RTSS

- ITN change vs. SMC vs. RTSS in seasonal settings
- ITN change vs. RTSS in highly seasonal settings
- ITN change vs RTSS in perennial settings

- RTSS changing cost - bar charts

- looking at the effects of ITN boost, ITN change, SMC introduction, RTSS - scatter plot with change in CE and change in impact
##################





