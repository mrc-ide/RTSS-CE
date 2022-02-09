# Figures & Tables -------------------------------------------------------------
# packages
library(tidyverse)
library(malariasimulation)
library(data.table)
library(patchwork)
library(grid)
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

output <- dat %>% filter(SMC == 0.85 & RTSS == 'SV' & ITN == 'pyr' & ITNuse == 0 & resistance == 0 & ITNboost == 0 )%>%
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
output <- dat %>% filter(SMC == 0.85 & RTSS == 'SV' & ITN == 'pyr' & ITNuse == 0 & resistance == 0 & ITNboost == 0)

supp.labs <- c("PfPR: 0.1", "PfPR: 0.2", "PfPR: 0.4")
names(supp.labs) <- c(0.1, 0.2, 0.4)

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
  geom_vline(data = RTSStime, aes(xintercept=rtss, alpha='RTS,S SV dose 3'), lty=2, color='blue') +
  geom_vline(data = ITNtime, aes(xintercept=itn, alpha='ITN'), lty=2, color='green') +
  labs(x='month', y='clinical incidence (month), 0-5 years', color='PfPR') +
  scale_x_continuous(breaks = seq(1,12,1)) +
  facet_wrap(seasonality~pfpr, labeller = labeller(pfpr=supp.labs)) +
  coord_cartesian(xlim = c(0,12)) +
  scale_alpha_manual(values = c(rep(1,3))) +
  guides(alpha = guide_legend(title = 'intervention',
                              override.aes = list(color = c('green','blue','red')))) +
  theme_classic()

ggsave('./03_output/seasonalityIII.pdf', width=10, height=5)


# delta change -----------------------------------------------------------------
# cost_per_dose <- c(2.69, 6.52, 12.91)
# delivery_cost <- c(0.96, 1.62, 2.67)

output <- scenarios %>%
  filter(cost_per_dose==6.52 & delivery_cost==1.62) %>% filter(resistance==0) %>%
  filter(intervention!='none') %>%
  mutate(deltadaly = daly_baseline - daly,
         deltacost = cost_total - cost_total_baseline,
         cost_daly_averted = (cost_total - cost_total_baseline) / deltadaly)

ggplot(mapping=aes(x=deltadaly, y=deltacost)) +
  geom_point(data=output, aes(color=intervention, alpha=ITNuse, shape=factor(pfpr))) +
  geom_hline(yintercept = 0, lty=2, color="black") +
  geom_vline(xintercept = 0, lty=2, color="black") +
  facet_wrap(~seasonality) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(y='change in cost (USD)',
       x='absolute change in DALYs',
       color='intervention',
       title='intervention impact',
       shape="PfPR",
       alpha="baseline ITN usage",
       subtitle='matched on seasonality, PfPR, and ITN usage',
       caption='Assuming resistance == 0, SMC implemented in seasonal settings')

ggsave('./03_output/impact_cloudI.pdf', width=8, height=7)

table(output$ID)
summary(output$deltadaly)
negvalues <- output %>% filter(deltadaly<0) # three negative scenarios

supp.labs <- c("ITN use: 0", "ITN use: 0.25", "ITN use: 0.50", "ITN use: 0.75")
names(supp.labs) <- c(0,.25,.50,.75)


# simple strategies, by baseline ITNuse and seasonality
# removing dominated strategies
output %>%
  # filter out mixed strategies
  filter(intervention %in% c('ITN 10% boost','ITN PBO','SMC','RTS,S EPI','RTS,S SV')) %>%
  group_by(ID) %>% arrange(ID, deltadaly) %>%
  # filter out dominated strategies
  mutate(dominate = case_when(deltacost > lead(deltacost,n=3L) ~ 1,
                              deltacost > lead(deltacost,n=2L) ~ 1,
                              deltacost > lead(deltacost,n=1L) ~ 1,
                              TRUE ~ 0)) %>%
  filter(dominate==0) %>%
ggplot(mapping=aes(x=deltadaly, y=deltacost)) +
  geom_line(aes(group=as.factor(ID)), color='lightgrey') +
  geom_point(aes(color=intervention), size=1.3) +
  geom_hline(yintercept = 0, lty=2, color="black") +
  geom_vline(xintercept = 0, lty=2, color="black") +
  facet_grid(seasonality~ITNuse, labeller = labeller(ITNuse=supp.labs), scales = "free") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(y='change in cost (USD)',
       x='change in DALYs averted',
       color='intervention',
       caption='Assuming resistance == 0, SMC newly implemented in seasonal settings only') +
  theme(plot.caption.position = "plot")


ggsave('./03_output/impact_cloudIII.pdf', width=15, height=7)


# All interventions, by baseline ITNuse and seasonality
output <- output %>% mutate(intervention2 = factor(intervention, levels = c('ITN 10% boost','ITN PBO','RTS,S EPI','RTS,S SV','ITN 10% boost + RTS,S','ITN PBO + RTS,S','SMC','ITN 10% boost + SMC','ITN PBO + SMC', 'RTS,S + SMC', 'ITN 10% boost + RTS,S + SMC','ITN PBO + RTS,S + SMC')))

RColorBrewer::brewer.pal(12, "Paired")
colors <- c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C",'deeppink',"#CAB2D6","#6A3D9A",'black',"#FDBF6F","#FF7F00")

# removing dominated strategies
output %>%
  # filter out mixed strategies
  group_by(ID) %>% arrange(ID, deltadaly) %>%
  # filter out dominated strategies
  mutate(dominate = case_when(deltacost > lead(deltacost,n=3L) ~ 1,
                              deltacost > lead(deltacost,n=2L) ~ 1,
                              deltacost > lead(deltacost,n=1L) ~ 1,
                              TRUE ~ 0)) %>%
  filter(dominate==0) %>%
ggplot(mapping=aes(x=deltadaly, y=deltacost)) +
  geom_line(aes(group=as.factor(ID)), color='lightgrey') +
  geom_point(aes(color=intervention2), size=1.3) +
  geom_hline(yintercept = 0, lty=2, color="black") +
  geom_vline(xintercept = 0, lty=2, color="black") +
  facet_grid(seasonality~ITNuse, labeller = labeller(ITNuse=supp.labs), scales = "free") +
  theme_classic() +
  scale_color_manual(values = colors) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(y='change in cost (USD)',
       x='change in DALYs averted',
       color='intervention',
       caption='Assuming resistance == 0, SMC implemented in seasonal settings') +
  theme(plot.caption.position = "plot")

ggsave('./03_output/impact_cloudIV.pdf', width=15, height=7)


table(output$intervention)
test <- output %>% filter(ITNuse==.75, seasonality=='highly seasonal', resistance==0) %>% group_by(ID) %>% summarize(n=n())



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
scenarios_univariate <- scenarios %>% filter(intervention != 'none') %>%
  select(file, ID, intervention, intervention_f, CE) %>%
  arrange(ID, CE)

summary(scenarios_univariate$CE)

# all interventions
ggplot(data=scenarios_univariate) +
  geom_density(aes(x=CE, y=..count..,
                   fill=intervention_f, color=intervention_f, group=intervention_f), alpha=0.1) +
  scale_x_continuous(limits=c(-1000, 1000)) +
  labs(x=expression(paste(Delta," cost / ", Delta, " DALYs")),
       y='density (count)',
       color = 'intervention',
       fill = 'intervention',
       caption=paste0('range in x = ', round(min(scenarios_univariate$CE)), ' to ', round(max(scenarios_univariate$CE)))) +
  theme_classic() +
  theme(plot.caption.position = "plot")

ggsave('./03_output/interventionCE_density.pdf', width=8, height=5)

# univariate interventions
scenarios_univariate2 <-  scenarios_univariate %>%
  filter(intervention %in% c('ITN 10% boost','ITN PBO','SMC','RTS,S EPI','RTS,S SV'))

ggplot(scenarios_univariate2) +
  geom_density(aes(x=CE, y=..count..,
                   fill=intervention, color=intervention, group=intervention), alpha=0.3) +
  scale_x_continuous(limits=c(-1000, 1000)) +
  labs(x=expression(paste(Delta," cost / ", Delta, " DALYs")),
       y='density (count)',
       caption=paste0('range in x = ', round(min(scenarios_univariate2$CE)), ' to ', round(max(scenarios_univariate2$CE)))) +
  theme_classic() +
  theme(plot.caption.position = "plot")

ggsave('./03_output/interventionCE_density_univariate.pdf', width=8, height=5)



# box and whisker delta cost / delta daly --------------------------------------
summary(scenarios$CE)

text_high <- textGrob("Highest\nvalue", gp=gpar(fontsize=13, fontface="bold"))
text_low <- textGrob("Lowest\nvalue", gp=gpar(fontsize=13, fontface="bold"))

ggplot(scenarios %>% filter(intervention != 'none')) +
  geom_hline(yintercept = 0, lty=2, color='grey') +
  geom_vline(xintercept = 6.5, lty=2, color='grey') +
  geom_boxplot(aes(x=rank, y=CE,
                   fill=intervention_f, color=intervention_f, group=intervention), alpha=0.3) +
  scale_y_continuous(limits=c(-1000, 1000)) +
  labs(x='',
       y=expression(paste(Delta," cost / ", Delta, " DALYs")),
       fill = 'intervention',
       color = 'intervention',
       caption=paste0('range in cost / DALYs: ', round(min(scenarios$CE)), ' to ', round(max(scenarios$CE)))) +
  annotation_custom(textGrob("Univariate strategies"),xmin=2,xmax=6,ymin=-1200,ymax=-1200) +
  annotation_custom(textGrob("Mixed strategies"),xmin=7,xmax=13,ymin=-1200,ymax=-1200) +
  scale_x_continuous(breaks=c(0)) +
  coord_cartesian(clip="off") +
  theme_classic() +
  theme(plot.caption.position = "plot")

ggsave('./03_output/box_whisker_CE.pdf', width=10, height=5)




# per dose RTS,S cost ----------------------------------------------------------
# set up cost vector
cost_per_dose2 <- seq(-50,100,.5) %>% as_tibble %>% rename(cost_per_dose2=value)

# pull out univariate scenarios
output <- scenarios %>% filter(intervention %in% c('RTS,S SV', 'RTS,S EPI', 'SMC', 'none', 'ITN 10% boost', 'ITN PBO')) %>%
  select(file:dose4, daly, ITNcost:intervention) %>% merge(cost_per_dose2) %>%
  mutate(cost_vax2 = (dose1 + dose2 + dose3 + dose4) * (cost_per_dose2),
         cost_total2 = cost_ITN + cost_clinical + cost_severe + cost_SMC + cost_vax2)

# checks
test <- output %>% filter(RTSS=='SV' & pfpr==0.2 & seasonality=='seasonal') %>%
  select(cost_per_dose, delivery_cost, cost_vax, cost_per_dose2, cost_vax2) %>%
  mutate(checkdosecost=(cost_per_dose+delivery_cost)/cost_per_dose2,
         checkvaxcost=cost_vax/cost_vax2)

head(table(test$checkdosecost, test$checkvaxcost)) # fractions should be the same

output %>% group_by(ID) %>% summarize(n=n())

# rank order CE
test <- output %>% select(ID, pfpr, seasonality, ITNuse, resistance, cost_per_dose2, intervention, daly, daly_baseline, cost_total2, cost_total_baseline) %>%
  mutate(CE = (cost_total2 - cost_total_baseline) / (daly_baseline - daly)) %>%
  ungroup() %>%
  group_by(ID, seasonality, pfpr, ITNuse, resistance, cost_per_dose2) %>%
  arrange(ID, seasonality, pfpr, ITNuse, resistance, cost_per_dose2, CE) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(rank = ifelse(intervention %in% c('RTS,S EPI', 'RTS,S SV'), cost_per_dose2, NA)) %>%
  group_by(ID, seasonality, pfpr, ITNuse, resistance) %>%
  summarize(rank2 = max(rank, na.rm=T))

table(test$rank2)

ggplot(data=test) +
  geom_bar(aes(x=rank2, y=..count..,fill=as_factor(resistance), group=as_factor(resistance))) +
  theme_classic() +
  labs(x='RTS,S cost (USD) per dose',
       y='count',
       fill='resistance',
       caption='13 scenarios where RTS,S is most CE >=$100 (all at resistance = 0.8)\n1 scenario where RTS,S is most CE == -$40') +
  scale_x_continuous(breaks=seq(-10,10,1), limits=c(-10,10)) +
  theme(plot.caption.position = "plot")

ggsave('./03_output/RTSS_price_dist.pdf', width=7, height=4)


# per dose ITN cost -------------------------------------------------------------------
# set up cost vector
population <- 10000
sim_length <- 15*365

cost_per_dose2 <- seq(-5,200,.5) %>% as_tibble %>% rename(cost_per_dose2=value)

# pull out univariate scenarios
output <- scenarios %>% filter(intervention %in% c('RTS,S SV', 'RTS,S EPI', 'SMC', 'none', 'ITN 10% boost', 'ITN PBO')) %>%
  merge(cost_per_dose2) %>%
  mutate(cost_ITN2 = population * annual_percapita_nets_distributed * sim_length/365 * cost_per_dose2,
         cost_total2 = cost_ITN2 + cost_clinical + cost_severe + cost_SMC + cost_vax)

# checks
output %>% group_by(ID) %>% summarize(n=n())

# rank order CE
test <- output %>% select(ID, pfpr, seasonality, ITNuse, resistance, cost_per_dose2, intervention, daly, daly_baseline, cost_total2, cost_total_baseline) %>%
  mutate(CE = (cost_total2 - cost_total_baseline) / (daly_baseline - daly)) %>%
  ungroup() %>%
  group_by(ID, seasonality, pfpr, ITNuse, resistance, cost_per_dose2) %>%
  arrange(ID, seasonality, pfpr, ITNuse, resistance, cost_per_dose2, CE) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(rank = ifelse(intervention %in% c('ITN 10% boost', 'ITN PBO'), cost_per_dose2, NA)) %>%
  group_by(ID, seasonality, pfpr, ITNuse, resistance) %>%
  summarize(rank2 = max(rank, na.rm=T))

table(test$rank2)
test %>% filter(rank2 == -Inf)

ggplot(data=test) +
  geom_bar(aes(x=rank2, y=..count..,fill=as_factor(resistance), group=as_factor(resistance))) +
  theme_classic() +
  labs(x='ITN cost (USD) per dose',
       y='count',
       fill='resistance',
       caption='40 scenarios where ITNs are most CE >=$100 \n13 scenarios where ITNs are never the most CE (all resistance==0.8)') +
  scale_x_continuous(breaks=seq(-10,200,10), limits=c(-5,200)) +
  scale_y_continuous(limits=c(0,5)) +
  theme(plot.caption.position = "plot")


ggsave('./03_output/ITN_price_dist.pdf', width=7, height=4)



# per dose SMC cost -------------------------------------------------------------------
# set up cost vector
cost_per_dose2 <- seq(-50,100,.5) %>% as_tibble %>% rename(cost_per_dose2=value)

# pull out univariate scenarios
output <- scenarios %>% filter(intervention %in% c('RTS,S SV', 'RTS,S EPI', 'SMC')) %>%
  filter(seasonality=='seasonal') %>%
  merge(cost_per_dose2) %>%
  mutate(cost_SMC2 = n_182.5_1825 * SMC * cost_per_dose2 * smc_timesteps,
         cost_total2 = cost_ITN + cost_clinical + cost_severe + cost_SMC2 + cost_vax)

summary(output$cost_SMC); summary(output$cost_SMC2)


# checks
test <- output %>% filter(SMC==0.85 & pfpr==0.2 & seasonality=='seasonal') %>%
  select(cost_per_dose, cost_SMC, cost_per_dose2, cost_SMC2) %>%
  mutate(checkdosecost=1.44/cost_per_dose2,
         checkvaxcost=cost_SMC/cost_SMC2)

head(table(test$checkdosecost, test$checkvaxcost)) # fractions should be the same

output %>% group_by(ID) %>% summarize(n=n())

# rank order CE
test <- output %>% select(ID, pfpr, seasonality, ITNuse, resistance, cost_per_dose2, intervention, daly, daly_baseline, cost_total2, cost_total_baseline) %>%
  mutate(CE = (cost_total2 - cost_total_baseline) / (daly_baseline - daly)) %>%
  ungroup() %>%
  group_by(ID, seasonality, pfpr, ITNuse, resistance, cost_per_dose2) %>%
  arrange(ID, seasonality, pfpr, ITNuse, resistance, cost_per_dose2, CE) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(rank = ifelse(intervention == 'SMC', cost_per_dose2, NA)) %>%
  group_by(ID, seasonality, pfpr, ITNuse, resistance) %>%
  summarize(rank2 = max(rank, na.rm=T))

table(test$rank2)

ggplot(data=test) +
  geom_bar(aes(x=rank2, y=..count..,fill=as_factor(resistance), group=as_factor(resistance))) +
  theme_classic() +
  labs(x='SMC cost per dose (USD)',
       y='count',
       fill='resistance',
       caption='1 scenario where SMC is most CE == -$19 \ncomparisons solely between RTS,S and SMC strategies in seasonal settings') +
  scale_x_continuous(breaks=seq(-10,10,1), limits=c(-10,10)) +
  theme(plot.caption.position = "plot")

ggsave('./03_output/SMC_price_dist.pdf', width=7, height=4)



# ITN dist vs ITN use ----------------------------------------------------------
# relationship between cost_ITN_linear and cost_ITN:
summary(dalyoutput_cost$cost_ITN_linear); summary(dalyoutput_cost$cost_ITN)

A <- ggplot(dalyoutput_cost %>% filter(ITN=='pyr')) +
  geom_point(aes(x=cost_ITN_linear, y = cost_ITN, colour=ITNuse2)) +
  geom_abline(slope=1) +
  scale_x_continuous(limits=c(0, 265000)) +
  scale_y_continuous(limits=c(0, 408000)) +
  scale_color_gradient(limits=c(0,0.85)) +
  labs(x='Linear ITN cost', y='Netz ITN cost', color='ITN use', title='Pyrethroid') +
  theme_classic()

B <- ggplot(dalyoutput_cost %>% filter(ITN=='pbo')) +
  geom_point(aes(x=cost_ITN_linear, y = cost_ITN, colour=ITNuse2)) +
  geom_abline(slope=1) +
  scale_x_continuous(limits=c(0, 265000)) +
  scale_y_continuous(limits=c(0, 408000)) +
  scale_color_gradient(limits=c(0,0.85)) +
  labs(x='Linear ITN cost', y='Netz ITN cost', color='ITN use', title='Pyrethroid + PBO') +
  theme_classic()

A + B + plot_layout(guides = "collect", nrow=1) + plot_annotation(tag_levels = 'A')

ggsave('./03_output/ITN_netz.pdf', width=7, height=4)


# ICER table RTS,S doses -------------------------------------------------------

cost_per_dose2 <- c(2, 5, 10) %>% as_tibble %>% rename(cost_per_dose2=value)

# pull out univariate scenarios
scenarios %>%
  filter(intervention %in% c('RTS,S SV', 'RTS,S EPI')) %>%
  filter(resistance == 0 & seasonality == 'perennial' & ITNuse == 0.25) %>%
  merge(cost_per_dose2) %>%
  mutate(cost_vax2 = (dose1 + dose2 + dose3 + dose4) * (cost_per_dose2),
         cost_total2 = cost_ITN + cost_clinical + cost_severe + cost_SMC + cost_vax2) %>%
  # calculate ICER
  mutate(CE_daly = (cost_total2 - cost_total_baseline) / (daly_baseline - daly),
         CE_cinc = (cost_total2 - cost_total_baseline) / (cases_baseline - cases)) %>%

  group_by(cost_per_dose2) %>%
  summarize(median_d = median(CE_daly), min_d = min(CE_daly), max_d = max(CE_daly),
            median_c = median(CE_cinc), min_c = min(CE_cinc), max_c = max(CE_cinc))

# how many observations - 1 x season, 3 x pfpr, 2 x RTSS, 4 x ITN use
nrow(scenarios %>%
       filter(intervention %in% c('RTS,S SV', 'RTS,S EPI')) %>%
       filter(resistance == 0 & seasonality == 'perennial'))





