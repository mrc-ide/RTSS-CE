# Figures & Tables -------------------------------------------------------------
# packages
library(tidyverse)
library(malariasimulation)
library(data.table)
library(patchwork)
library(LaCroixColoR)

# devtools::install_github('mrc-ide/malariasimulation@dev', force=TRUE)
# devtools::install_github('johannesbjork/LaCroixColoR')

# Look here: https://github.com/htopazian/rtss_malariasimulation
# find ideas in rmarkdowns for figures and insert code if helpful

dat <- readRDS("./03_output/rtss_raw.rds")
averted <- readRDS("./03_output/rtss_avert.rds")
dalyoutput <- readRDS('./03_output/dalyoutput.rds')
dalyoutput_cost <- readRDS('./03_output/dalyoutput_cost.rds')
scenarios <- readRDS('./03_output/scenarios.rds')


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
  geom_line(aes(x=month, y=n_inc_clinical_182.5_1825/n_182.5_1825), alpha = 0.8) +
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
output <- dat %>% filter(SMC == 0.85 & RTSS == 'SV' & ITN == 'pyr' & ITNuse == 0.5 & resistance == 0 & ITNboost == 0 & ITNuse == .5)

SMCtime <- output %>% select(smc_timesteps, seasonality, pfpr) %>%
  group_by(smc_timesteps, seasonality, pfpr) %>%
  mutate(t1 = ifelse(seasonality=='seasonal', unlist(smc_timesteps)[[1]], unlist(smc_timesteps)[[25]]),
         t2 = ifelse(seasonality=='seasonal', unlist(smc_timesteps)[[2]], unlist(smc_timesteps)[[26]]),
         t3 = ifelse(seasonality=='seasonal', unlist(smc_timesteps)[[3]], unlist(smc_timesteps)[[27]]),
         t4 = ifelse(seasonality=='seasonal', unlist(smc_timesteps)[[4]], unlist(smc_timesteps)[[28]]),
         t5 = ifelse(seasonality=='seasonal', unlist(smc_timesteps)[[5]], NA)) %>%
  distinct() %>%
  pivot_longer(cols=t1:t5, names_to = "time", values_to = "smc") %>%
  mutate(smc = smc / (365/12) + 1)

RTSStime <- output %>% select(rtss_mass_timesteps, seasonality, pfpr) %>%
  group_by(rtss_mass_timesteps, seasonality, pfpr) %>%
  mutate(rtss = unlist(rtss_mass_timesteps)[[1]]) %>%
  distinct() %>%
  mutate(rtss = rtss / (365/12) + 1)

ITNtime <- output %>% select(bednet_timesteps, seasonality, pfpr) %>%
  group_by(bednet_timesteps, seasonality, pfpr) %>%
  mutate(itn = unlist(bednet_timesteps)[[3]]) %>%
  distinct() %>%
  mutate(itn = itn / (365/12) + 1)

none <- dat %>%
  filter(ITNboost==0 & ITNuse==0.5 & RTSS=='none' & ITN=='pyr' & resistance==0 & (SMC==0 | (seasonality=='highly seasonal')))

ggplot(data = none %>% filter(seasonality != 'perennial')) +
  geom_line(aes(x=month, y=n_inc_clinical_182.5_1825/n_182.5_1825, color=as.factor(pfpr)), alpha = 0.8) +
  geom_vline(data = SMCtime, aes(xintercept=smc, alpha='SMC'), lty=2, color='red') +
  geom_vline(data = RTSStime, aes(xintercept=rtss, alpha='SV RTS,S dose 3'), lty=2, color='blue') +
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


# delta change -----------------------------------------------------------------
# cost_per_dose <- c(2.69, 6.52, 12.91)
# delivery_cost <- c(0.96, 1.62, 2.67)

output <- dalyoutput_cost %>%
  filter(cost_per_dose==6.52 & delivery_cost==1.62)

none <- output %>%
  filter(ITNboost==0 & RTSS=='none' & ITN=='pyr' & resistance==0 & (SMC==0 | (seasonality=='highly seasonal'))) %>%
  select(seasonality, pfpr, ITNuse, daly, cost_total) %>%
  rename(daly_baseline = daly,
         cost_total_baseline = cost_total)

output <- output %>% left_join(none, by=c('seasonality', 'pfpr', 'ITNuse')) %>%
  mutate(deltadaly = daly_baseline - daly,
         deltacost = cost_total - cost_total_baseline,
         cost_daly_averted = (cost_total - cost_total_baseline) / deltadaly)

ITNboost <- output %>% filter(RTSS=='none') %>% filter(!(seasonality=='seasonal' & SMC==.85)) %>% filter(ITNboost==1) %>% filter(resistance==0) %>% mutate(intervention='boost')

ITNpbo <- output %>% filter(RTSS=='none') %>% filter(!(seasonality=='seasonal' & SMC==.85)) %>% filter(ITN=='pbo') %>% filter(resistance==0) %>% mutate(intervention='pbo')

SMC <- output %>% filter(RTSS=='none') %>% filter(ITNboost==0 & ITN=='pyr') %>%
  filter(seasonality!='highly seasonal') %>% filter(resistance==0) %>% filter(SMC==.85) %>% mutate(intervention='smc')

RTSSSV <- output %>% filter(RTSS=='SV') %>% filter(ITNboost!=1 & ITN!='pbo') %>%
  filter(!(seasonality=='seasonal' & SMC==.85)) %>% filter(resistance==0) %>% mutate(intervention='sv')

RTSSEPI <- output %>% filter(RTSS=='EPI') %>% filter(ITNboost!=1 & ITN!='pbo') %>%
  filter(!(seasonality=='seasonal' & SMC==.85)) %>% filter(resistance==0) %>% mutate(intervention='epi')

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
       y='absolute change in DALYs',
       color='intervention',
       title='intervention impact',
       shape="PfPR",
       alpha="baseline ITN usage",
       subtitle='matched on seasonality, PfPR, and ITN usage',
       caption='Assuming resistance == 0, SMC implemented in seasonal settings')

ggsave('./03_output/impact_cloudI.pdf', width=8, height=5)

test <- rbind(ITNboost, ITNpbo, SMC, RTSSSV, RTSSEPI) %>%
  group_by(seasonality, pfpr, intervention) %>%
  mutate(ID = paste(seasonality, pfpr, ITNuse, sep="_"))

table(test$ID)
summary(test$deltadaly)
negvalues <- test %>% filter(deltadaly<0)

supp.labs <- c("PfPR: 0.1", "PfPR: 0.2", "PfPR: 0.4")
names(supp.labs) <- c(0.1, 0.2, 0.4)

ggplot(mapping=aes(x=deltadaly, y=deltacost)) +
  geom_line(data=test, aes(group=as.factor(ID)), color='lightgrey') +
  geom_point(data=ITNboost, aes(color='ITN boost (10%)', alpha=ITNuse), size=1.3) +
  geom_point(data=ITNpbo, aes(color='ITN PBO', alpha=ITNuse), size=1.3) +
  geom_point(data=SMC, aes(color='SMC alone', alpha=ITNuse), size=1.3) +
  geom_point(data=RTSSSV, aes(color='RTS,S SV alone ($5)', alpha=ITNuse), size=1.3) +
  geom_point(data=RTSSEPI, aes(color='RTS,S EPI alone ($5)', alpha=ITNuse), size=1.3) +
  geom_hline(yintercept = 0, lty=2, color="black") +
  geom_vline(xintercept = 0, lty=2, color="black") +
  facet_grid(seasonality~pfpr, labeller = labeller(pfpr=supp.labs)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(y='change in cost (USD)',
       x='change in DALYs averted',
       color='intervention',
       title='intervention impact',
       alpha="baseline ITN usage",
       subtitle='matched on seasonality, PfPR, and ITN usage',
       caption='Assuming resistance == 0, SMC implemented in seasonal settings')

ggsave('./03_output/impact_cloudII.pdf', width=10, height=7)



# resistance ITNs --------------------------------------------------------------
output <- dalyoutput_cost %>%
  filter(RTSS=='none' & (SMC==0 | (seasonality=='highly seasonal'))) %>%
  filter(cost_per_dose==6.52 & delivery_cost==1.62) %>%
  filter(ITNuse > 0)

none <- output %>%
  filter(ITNboost==0 & ITN=='pyr' & resistance==0) %>%
  select(seasonality, pfpr, ITNuse, daly, cost_total) %>%
  rename(daly_baseline = daly,
         cost__total_baseline = cost_total)

output <- output %>% left_join(none, by=c('seasonality', 'pfpr', 'ITNuse')) %>%
  mutate(deltadaly = daly_baseline - daly,
         deltacost = cost_total - cost__total_baseline,
         cost_daly_averted = (cost_total - cost__total_baseline) / deltadaly)

output <- output %>% group_by(pfpr, seasonality, ITN, ITNuse, resistance) %>%
  filter(!(ITN=='pyr'& ITNboost==0)) %>%
  mutate(ITNintervention = case_when(ITN=="pbo"~'pbo',
                                     ITNboost==1~'boost',
                                     TRUE~NA_character_)) %>%
  mutate(ID = paste0(pfpr, seasonality, ITNuse, resistance, sep="_"))

table(test$ITNintervention, useNA ='always')
table(test$ID)
table(test$ITNboost, test$ITN)

ITNboost <- output %>% filter(ITNboost==1)
table(ITNboost$ITN) # only pyr nets

ITNpbo <- output %>% filter(ITN=='pbo')
table(ITNpbo$ITN) # only pbo nets


# plot
ggplot(mapping=aes(x=deltadaly, y=deltacost)) +
  geom_line(data=output, aes(group=as.factor(ID)), alpha=.5) +
  geom_point(data=ITNboost, aes(color='ITN PYR 10% boost', alpha=pfpr), size=1.2) +
  geom_point(data=ITNpbo, aes(color='ITN PBO switch', alpha=pfpr), size=1.2) +
  geom_hline(yintercept = 0, lty=2, color="black") +
  geom_vline(xintercept = 0, lty=2, color="black") +
  # geom_smooth(data=ITNboost, method='lm', se=F, aes(color='ITN PYR 10% boost')) +
  # geom_smooth(data=ITNpbo, method='lm', se=F, aes(color='ITN PBO switch')) +
  theme_classic() +
  facet_grid(ITNuse~resistance,
             labeller = label_both) +
  labs(y='change in cost (USD)',
       x='change in DALYs averted',
       color='intervention',
       alpha='PfPR',
       title='Insecticide resistance',
       subtitle='matched on seasonality, PfPR, and ITN usage',
       caption='')

ggsave('./03_output/resistanceI.pdf', width=10, height=8)



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
  facet_grid(seasonality~pfpr) +
  scale_alpha_discrete(range = c(.4, 1)) +
  labs(x='scenario: ITNuse - resistance',
       y='change in DALYs averted',
       fill='ITN',
       alpha = 'resistance',
       title='Insecticide resistance',
       subtitle='matched on seasonality, PfPR, and ITN usage',
       caption='faceted by seasonality and PfPR')

ggsave('./03_output/resistanceII.pdf', width=8, height=5)


# density plot -----------------------------------------------------------------

scenarios_univariate <- scenarios %>% filter(intervention != 'mixed') %>%
  select(file, ID, intervention, CE) %>%
  arrange(ID, CE)

summary(scenarios_univariate$CE)

ggplot(scenarios_univariate) +
  #geom_histogram(aes(x=CE, fill=intervention), binwidth = 1) +
  geom_density(aes(x=CE, y=..count..,
                   fill=intervention, color=intervention, group=intervention), alpha=0.2) +
  scale_x_continuous(limits=c(-1000, 1000)) +
  labs(x=expression(paste(Delta," cost / ", Delta, " DALYs")),
       y='density (count)',
       title='density plot of intervention cost-effectiveness',
       caption='range in x = -8031 to 19809') +
  theme_classic()

ggsave('./03_output/interventionCE_density.pdf', width=8, height=5)


# per dose RTS,S cost -------------------------------------------------------------------
# set up cost vector
cost_per_dose2 <- seq(-50,100,.5) %>% as_tibble %>% rename(cost_per_dose2=value)

# pull out univariate scenarios
output <- scenarios %>% filter(intervention != 'mixed') %>%
  select(file:dose4, daly, ITNcost:intervention) %>% merge(cost_per_dose2) %>%
  mutate(cost_vax2 = (dose1 + dose2 + dose3 + dose4) * (cost_per_dose2),
         cost_total2 = cost_ITN + cost_clinical + cost_severe + cost_SMC + cost_vax2)

# checks
test <- output %>% filter(RTSS=='SV' & pfpr==0.2 & seasonality=='seasonal') %>%
  select(cost_per_dose, delivery_cost, cost_vax, cost_per_dose2, cost_vax2) %>%
  mutate(checkdosecost=(cost_per_dose+delivery_cost)/cost_per_dose2,
         checkvaxcost=cost_vax/cost_vax2)

table(test$checkdosecost, test$checkvaxcost) # fractions should be the same

output %>% group_by(ID) %>% summarize(n=n())

# rank order CE
test <- output %>% select(ID, pfpr, seasonality, ITNuse, resistance, cost_per_dose2, intervention, daly, daly_baseline, cost_total2, cost_total_baseline) %>%
  mutate(CE = (cost_total2 - cost_total_baseline) / (daly_baseline - daly)) %>%
  ungroup() %>%
  group_by(ID, seasonality, pfpr, ITNuse, resistance, cost_per_dose2) %>%
  arrange(ID, seasonality, pfpr, ITNuse, resistance, cost_per_dose2, CE) %>%
  slice_head(n = 1) %>%
  mutate(rank = ifelse(intervention %in% c('RTSS EPI', 'RTSS SV'), cost_per_dose2, NA)) %>%
  group_by(ID, seasonality, pfpr, ITNuse, resistance) %>%
  summarize(rank2 = max(rank, na.rm=T))

table(test$rank2)
test2 <- test %>% filter(pfpr==0.2 & seasonality=='seasonal' & resistance==0 & ITNuse==0.5)

ggplot(data=test) +
  geom_bar(aes(x=rank2, y=..count..,fill=as_factor(resistance), group=as_factor(resistance))) +
  theme_classic() +
  labs(x='RTS,S cost (USD)',
       y='count',
       fill='resistance',
       title='Distribution of the maximium cost per dose where \nRTS,S is the most cost efficient strategy',
       caption='14 scenarios where RTS,S is most CE >=$100 (all at resistance = 0.8)') +
  scale_x_continuous(breaks=seq(-30,10,10), limits=c(-30,10))

ggsave('./03_output/RTSS_price_dist.pdf', width=5, height=3)




# IDEAS ###########
# - resistance - ITNs, upgrading to PBO and boosting. no SMC add ins or RTSS
#
# - ITN change vs. SMC vs. RTSS in seasonal settings
# - ITN change vs. RTSS in highly seasonal settings
# - ITN change vs RTSS in perennial settings
#
# - RTSS changing cost - bar charts
#
# - looking at the effects of ITN boost, ITN change, SMC introduction, RTSS - scatter plot with change in CE and change in impact
##################





