# Figures & Tables -------------------------------------------------------------
# packages
library(tidyverse)
library(data.table)
library(patchwork)
library(grid)
library(LaCroixColoR)

# devtools::install_github('mrc-ide/malariasimulation@dev', force=TRUE)
# devtools::install_github('johannesbjork/LaCroixColoR')

# load data
scenarios <- readRDS('./03_output/test_scenarios.rds')
dalyoutput_cost <- readRDS('./03_output/test_dalyoutput.rds')


# delta CE change --------------------------------------------------------------
# < facet by season -----------------
output <- scenarios %>%
  filter(cost_per_dose==6.52) %>% filter(resistance==0) %>%
  filter(intervention!='none') %>%

  mutate(deltadaly = daly_baseline - daly,
         deltacost = cost_total - cost_total_baseline,
         cost_daly_averted = deltacost / deltadaly) %>%

  group_by(ID, seasonality, intervention) %>%

  summarize(n = n(),
            deltadaly_m = mean(deltadaly, na.rm = T),
            deltadaly_l = deltadaly_m - 1.96*(sd(deltadaly, na.rm = T)/sqrt(50)),
            deltadaly_u = deltadaly_m + 1.96*(sd(deltadaly, na.rm = T)/sqrt(50)),

            deltacost_m = mean(deltacost, na.rm = T),
            deltacost_l = deltacost_m - 1.96*(sd(deltacost, na.rm = T)/sqrt(50)),
            deltacost_u = deltacost_m + 1.96*(sd(deltacost, na.rm = T)/sqrt(50))) %>%

  ungroup() %>%

  mutate(seasonality = factor(seasonality, levels = c('perennial', 'seasonal', 'highly seasonal')))

# creating factor variable for interventions, removing 'none'
output <- output %>% mutate(intervention2 = factor(intervention, levels = c('ITN 10% increase','ITN PBO','RTS,S age-based','RTS,S seasonal','ITN 10% increase + RTS,S','ITN PBO + RTS,S','SMC','ITN 10% increase + SMC','ITN PBO + SMC', 'RTS,S + SMC', 'ITN 10% increase + RTS,S + SMC','ITN PBO + RTS,S + SMC')))

RColorBrewer::brewer.pal(12, "Paired")
# colors <- c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C",'deeppink',"#CAB2D6","#6A3D9A",'black',"#FDBF6F","#FF7F00")
# ITN 10% boost #1F78B4
# RTSS EPI #B2DF8A
# RTSS SV #33A02C
# SMC #FB9A99
# ITN 10% boost + RTSS #A6CEE3
# ITN 10% boost + SMC  #E31A1C
# RTSS + SMC 'deeppink'
# ITN 10% boost + RTSS + SMC #6A3D9A

  colors <- c('#1F78B4',  '#B2DF8A', '#33A02C', '#FB9A99', '#A6CEE3', '#E31A1C',  'deeppink', '#6A3D9A')

A <- ggplot(data = output, mapping=aes(x=deltadaly_m, y=deltacost_m)) +
  geom_line(aes(group=as.factor(ID)), color='lightgrey', size=.5, alpha=0.2) +
  geom_errorbarh(aes(y=deltacost_m, xmin = deltadaly_l, xmax = deltadaly_u, color = intervention2), alpha=0.4) +
  # change in cost is minimal between parameter draws b/c most cost comes from interventions which are constant
  geom_point(aes(color=intervention2), size=2, alpha=0.6) +
  geom_hline(yintercept = 0, lty=2, color="black") +
  geom_vline(xintercept = 0, lty=2, color="black") +
  facet_grid(~seasonality, scales = "free") +
  theme_classic() +
  scale_color_manual(values = colors) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(title='All strategies',
       y='change in cost (USD)',
       x='change in DALYs averted',
       color='intervention')

# removing dominated strategies
output2 <- output %>%
  # filter out mixed strategies
  group_by(ID) %>% arrange(ID, deltacost_m) %>%
  filter(deltadaly_m >= 0 & deltacost_m > 0) %>%
  # filter out dominated strategies
  mutate(dominate = case_when(deltadaly_m < lag(deltadaly_m,n=12L) ~ 1,
                              deltadaly_m < lag(deltadaly_m,n=11L) ~ 1,
                              deltadaly_m < lag(deltadaly_m,n=10L) ~ 1,
                              deltadaly_m < lag(deltadaly_m,n=9L) ~ 1,
                              deltadaly_m < lag(deltadaly_m,n=8L) ~ 1,
                              deltadaly_m < lag(deltadaly_m,n=7L) ~ 1,
                              deltadaly_m < lag(deltadaly_m,n=6L) ~ 1,
                              deltadaly_m < lag(deltadaly_m,n=5L) ~ 1,
                              deltadaly_m < lag(deltadaly_m,n=4L) ~ 1,
                              deltadaly_m < lag(deltadaly_m,n=3L) ~ 1,
                              deltadaly_m < lag(deltadaly_m,n=2L) ~ 1,
                              deltadaly_m < lag(deltadaly_m,n=1L) ~ 1,
                              TRUE ~ 0)) %>%
  filter(dominate==0) %>%
  # marking extended dominated strategies
  mutate(ICER = (deltacost_m-lag(deltacost_m)) / (deltadaly_m-lag(deltadaly_m)),
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
  mutate(ICER = (deltacost_m-lag(deltacost_m)) / (deltadaly_m-lag(deltadaly_m)),
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
  mutate(ICER = (deltacost_m-lag(deltacost_m)) / (deltadaly_m-lag(deltadaly_m)),
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
  mutate(ICER = (deltacost_m-lag(deltacost_m)) / (deltadaly_m-lag(deltadaly_m)),
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
  mutate(ICER = (deltacost_m-lag(deltacost_m)) / (deltadaly_m-lag(deltadaly_m)),
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
  filter(dominate==0)

B <- ggplot(data = output2, mapping=aes(x=deltadaly_m, y=deltacost_m)) +
  geom_line(aes(group=as.factor(ID)), color='lightgrey', size=.5, alpha=0.2) +
  geom_errorbarh(aes(y=deltacost_m, xmin = deltadaly_l, xmax = deltadaly_u, color = intervention2), alpha=0.4) +
  # change in cost is minimal between parameter draws b/c most cost comes from interventions which are constant
  geom_point(aes(color=intervention2), size=2, alpha=0.6) +
  geom_hline(yintercept = 0, lty=2, color="black") +
  geom_vline(xintercept = 0, lty=2, color="black") +
  facet_grid(~seasonality, scales = "free") +
  theme_classic() +
  scale_color_manual(values = colors) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(title='Dominated strategies removed',
       y='change in cost (USD)',
       x='change in DALYs averted',
       color='intervention')

A + B + plot_layout(guides = "collect", nrow=2) + plot_annotation(tag_levels = 'A')

ggsave(paste0('./03_output/impact_cloud_all.pdf'), width=12, height=7)



# box and whisker delta cost / delta daly --------------------------------------
# inspect range of CE values
summary(scenarios$CE)

test <- scenarios %>% filter(CE < 0) %>% select(CE,intervention,file:fifth, cost_total, cost_total_baseline, daly, daly_baseline)

scenarios %>% group_by(intervention_f) %>%
  summarize(n=n(),
            median = round(median(CE, na.rm=T)),
            q25 = round(quantile(CE, p=0.25, na.rm=T)),
            q75 = round(quantile(CE, p=0.75, na.rm=T)))

# set up text for plot
text_high <- textGrob("Highest\nvalue", gp=gpar(fontsize=13, fontface="bold"))
text_low <- textGrob("Lowest\nvalue", gp=gpar(fontsize=13, fontface="bold"))

# set up order of interventions by median CE for plot
levels <- scenarios %>%
  mutate(group = case_when(intervention_f %in% c('ITN 10% increase', 'ITN PBO', 'SMC', 'RTS,S seasonal', 'RTS,S age-based') ~ 'uni',
                           TRUE ~ 'multi')) %>%
  group_by(intervention_f, group) %>%
  summarize(med = median(CE, na.rm=T) )%>% ungroup() %>%
  arrange(desc(group), med)

# plot of cost per DALY averted
scenarios %>%
  mutate(intervention_f = factor(intervention, levels=levels$intervention_f)) %>%
  mutate(rank=as.numeric(intervention_f)) %>%

ggplot(aes(x=rank, y=CE, fill=intervention_f, color=intervention_f, group=intervention)) +
  geom_hline(yintercept = 0, lty=2, color='grey') +
  geom_vline(xintercept = 5.5, lty=2, color='grey') +
  geom_boxplot(alpha=0.3) +
  coord_cartesian(ylim=c(-100, 500), clip="off") +
  labs(x='',
       y=expression(paste(Delta," cost / ", Delta, " DALYs")),
       fill = 'intervention',
       # caption=paste0('range in cost / DALYs: ', round(min(scenarios$CE, na.rm = T)), ' to ', round(max(scenarios$CE, na.rm=T))),
       color = 'intervention') +
  annotation_custom(textGrob("Univariate strategies"),xmin=1,xmax=5,ymin=-150,ymax=-150) +
  annotation_custom(textGrob("Mixed strategies"),xmin=6,xmax=12,ymin=-150,ymax=-150) +
  scale_x_continuous(breaks=c(0)) +
  theme_classic() +
  theme(plot.caption.position = "plot")

ggsave('./03_output/box_whisker_CE.pdf', width=10, height=5)


# plotting with cost per cases averted
# inspect range of CE case values
summary(scenarios$CE_case)

scenarios %>% group_by(intervention_f) %>%
  summarize(n=n(),
            median = round(median(CE_case)),
            q25 = round(quantile(CE_case, p=0.25)),
            q75 = round(quantile(CE_case, p=0.75)))

scenarios %>% filter(intervention != 'none') %>%
  mutate(intervention_f = factor(intervention, levels=levels$intervention_f)) %>%
  mutate(rank=as.numeric(intervention_f)) %>%

ggplot(aes(x=rank, y=CE_case, fill=intervention_f, color=intervention_f, group=intervention)) +
  geom_hline(yintercept = 0, lty=2, color='grey') +
  geom_vline(xintercept = 5.5, lty=2, color='grey') +
  geom_boxplot(alpha=0.3) +
  coord_cartesian(ylim=c(-1, 35), clip="off") +
  labs(x='',
       y=expression(paste(Delta," cost / ", Delta, " cases")),
       fill = 'intervention',
       color = 'intervention',
       caption=paste0('range in cost / cases: ', round(min(scenarios$CE_case, na.rm = T)), ' to ', round(max(scenarios$CE_case, na.rm=T)))) +
  annotation_custom(textGrob("Univariate strategies"),xmin=1,xmax=5,ymin=-4,ymax=-4) +
  annotation_custom(textGrob("Mixed strategies"),xmin=6,xmax=12,ymin=-4,ymax=-4) +
  scale_x_continuous(breaks=c(0)) +
  theme_classic() +
  theme(plot.caption.position = "plot")

ggsave('./03_output/box_whisker_CE_cases.pdf', width=10, height=5)



# per dose RTS,S cost ----------------------------------------------------------
# function to find the most common character value in a group
calculate_mode <- function(x) {
  uniqx <- unique(na.omit(x))
  uniqx[which.max(tabulate(match(x, uniqx)))]
}

output <- scenarios %>%
  filter(cost_per_dose==6.52 & delivery_cost==1.62) %>%
  # choose univariate scenarios
  filter(intervention %in% c('RTS,S seasonal', 'RTS,S age-based', 'SMC', 'none', 'ITN 10% increase', 'ITN PBO')) %>%
  # combine RTS,S strategies
  mutate(intervention = ifelse(intervention=='RTS,S age-based' | intervention=='RTS,S seasonal', 'RTS,S', intervention)) %>%
  arrange(ID) %>%
  # find minimum CE in each group. If intervention == RTS,S find the second or third lowest CE
  group_by(ID) %>%
  mutate(CE = ifelse(CE == min(CE) & intervention == 'RTS,S', 100000, CE),
         CE = ifelse(CE == min(CE) & intervention == 'RTS,S', 100000, CE),
         CEmin = min(CE, na.rm=T),
         interventionmin = ifelse(CE==CEmin, intervention, NA),
         interventionmin = calculate_mode(interventionmin)) %>%
  filter(intervention %in% c('RTS,S')) %>%
  rowwise() %>%
  # calculate needed cost of RTS,S to match the min CE intervention
  # CEmin = (cost_total - cost_total_baseline) / (daly_baseline - daly)
  mutate(delta_cost = CEmin * (daly_baseline - daly),
         cost_total = delta_cost + cost_total_baseline,
         cost_vax = cost_total - (cost_ITN + cost_clinical + cost_severe + cost_SMC),
         per_dose = cost_vax / (dose1 + dose2 + dose3 + dose4),
         costRTSS = per_dose - 1.62   # subtracting delivery cost
  ) %>%
  select(ID, seasonality, pfpr, ITNuse, ITN, RTSS, RTSScov, resistance, SMC, treatment, intervention, interventionmin, CE, CEmin, costRTSS) %>%
  group_by(ID) %>%
  arrange(ID,costRTSS)

  # all seasons
seasoncosts <- output %>% group_by(seasonality) %>%
  summarize(q25 = round(quantile(costRTSS, probs = 0.25, na.rm=T),2),
            med = round(median(costRTSS, na.rm=T),2),
            q75 = round(quantile(costRTSS, probs = 0.75, na.rm=T),2))

ggplot(output) +
  geom_boxplot(aes(x=factor(seasonality, levels=c('perennial','seasonal','highly seasonal')), y=costRTSS), fill = 'cornflower blue', color = 'cornflower blue', alpha = 0.4) +
  geom_text(data = seasoncosts, aes(x=seasonality, y=q25, label=q25), size = 3, nudge_y = .5, nudge_x = -.2) +
  geom_text(data = seasoncosts, aes(x=seasonality, y=med, label=med), size = 3, nudge_y = .5, nudge_x = -.2) +
  geom_text(data = seasoncosts, aes(x=seasonality, y=q75, label=q75), size = 3, nudge_y = .5, nudge_x = -.2) +
  geom_vline(xintercept = 0, lty = 2, color = 'grey') +
  theme_classic() +
  labs(y='RTS,S cost (USD) per dose',
       #caption=paste0('range = ', round(min(output$costRTSS)), ' to ', round(max(output$costRTSS)))
       x=''
       ) +
  coord_cartesian(ylim = c(-5,15)) +
  theme(plot.caption.position = "plot")

ggsave('./03_output/RTSS_price_dist_all.pdf', width=6, height=4)


# table
seasoncosts

# how many are above $5
output %>% filter(costRTSS >= 5) %>% group_by(SMC, ITNuse) %>%
  summarize(n=n()) %>%
  ungroup() %>%
  mutate(t=sum(n))

output %>% filter(costRTSS >= 5) %>% group_by(SMC, seasonality) %>%
  summarize(n=n()) %>%
  ungroup() %>%
  mutate(t=sum(n))


# univariate stacked histogram patterns ----------------------------------------
# look at color choices
RColorBrewer::display.brewer.all()
RColorBrewer::brewer.pal(5, "Paired")
colors <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99")

colors <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99")

output <- scenarios %>%
  filter(cost_per_dose==6.52 & delivery_cost==1.62) %>%
  filter(intervention %in% c('ITN 10% increase','ITN PBO','SMC','RTS,S age-based','RTS,S seasonal')) %>%
  mutate(seasonality = factor(seasonality, levels = c('perennial', 'seasonal', 'highly seasonal'))) %>%
  group_by(ID) %>% arrange(ID, CE) %>%
  slice(1L)

table(output$intervention, useNA = 'always'); (8+19)/270

# by pfpr
A <- ggplot(output) +
  geom_bar(aes(x=as.factor(pfpr), fill=intervention_f), position="fill", show.legend = F) +
  labs(x='PfPR', y='Proportion most \ncost-effective choice', fill='intervention') +
  scale_fill_manual(values=colors) +
  facet_grid(~seasonality) +
  theme_classic()

# by current ITN usage
B <- ggplot(output) +
  geom_bar(aes(x=factor(ITNuse), fill=intervention_f), position="fill", show.legend=F) +
  labs(x='ITN use', y='Proportion most \ncost-effective choice', fill='intervention') +
  scale_fill_manual(values=colors) +
  facet_grid(~seasonality) +
  theme_classic()

# by insecticide resistance
C <- ggplot(output) +
  geom_bar(aes(x=factor(resistance), fill=intervention_f), position="fill", show.legend=F) +
  labs(x='Resistance', y='Proportion most \ncost-effective choice', fill='intervention') +
  scale_fill_manual(values=colors) +
  facet_grid(~seasonality) +
  theme_classic()

# by ITN distribution efficiency
ITNefficient <- function(var, label) {
  scenarios %>%
    filter(cost_per_dose==6.52 & delivery_cost==1.62) %>%
    filter(intervention %in% c('ITN 10% increase','ITN PBO','SMC','RTS,S age-based','RTS,S seasonal')) %>%
    mutate(seasonality = factor(seasonality, levels = c('perennial', 'seasonal', 'highly seasonal'))) %>%
    group_by(ID) %>% arrange(ID, {{var}}) %>%
    slice(1L) %>% select(intervention, intervention_f, seasonality, pfpr, {{var}}) %>%
    mutate(model=label) %>%
    rename(CE = {{var}})
}

output2 <- ITNefficient(CE, 'standard') %>%
  full_join(ITNefficient(CE_ITNmin, 'less efficient')) %>%
  full_join(ITNefficient(CE_ITNmax, 'more efficient'))


D <- ggplot(output2) +
  geom_bar(aes(x=factor(model, levels=c('less efficient', 'standard', 'more efficient'), labels=c('less', 'standard', 'more')), fill=intervention_f), position="fill") +
  labs(x='ITN efficiency', y='Proportion most \ncost-effective choice', fill='intervention') +
  scale_fill_manual(values=colors) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  facet_grid(~seasonality) +
  theme_classic()

E <- ggplot(output) +
  geom_bar(aes(x=factor(treatment, labels=c('low', 'medium', 'high')), fill=intervention_f), position="fill") +
  labs(x='Treatment coverage', y='Proportion most \ncost-effective choice', fill='intervention') +
  scale_fill_manual(values=colors) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  facet_grid(~seasonality) +
  theme_classic()

(A + B + C + D + E) + plot_layout(guides = "collect", nrow=2) + plot_annotation(tag_levels = 'A')

ggsave(paste0('./03_output/univariate_quad_by_season.pdf'), width=17, height=5)

# pfpr
test <- output %>%
  mutate(intervention =
           case_when(intervention %in% c('RTS,S age-based', 'RTS,S seasonal') ~ 'RTS,S',
                     TRUE ~ intervention)) %>%
  group_by(seasonality, pfpr) %>%
  mutate(n=n()) %>%
  group_by(seasonality, intervention, pfpr) %>%
  summarize(n=mean(n), t=n(), p=t/n*100)

# ITN use
test <- output %>%
  mutate(intervention =
           case_when(intervention %in% c('RTS,S age-based', 'RTS,S seasonal') ~ 'RTS,S',
                     TRUE ~ intervention)) %>%
  group_by(seasonality, ITNuse) %>%
  mutate(n=n()) %>%
  group_by(seasonality, intervention, ITNuse) %>%
  summarize(n=mean(n), t=n(), p=t/n*100)

# ITN efficiency
test <- output2 %>%
  mutate(intervention =
           case_when(intervention %in% c('RTS,S age-based', 'RTS,S seasonal') ~ 'RTS,S',
                     TRUE ~ intervention)) %>%
  group_by(seasonality, model) %>%
  mutate(n=n()) %>%
  group_by(seasonality, intervention, model) %>%
  summarize(n=mean(n), t=n(), p=t/n*100)




# CE table RTS,S doses ---------------------------------------------------------
# cost_per_dose <- c(2.69, 6.52, 12.91)
# delivery_cost <- c(0.96, 1.62, 2.67)
# cost_per_dose + delivery_cost

cost_per_dose2 <- c(2.69+1.62, 6.52+1.62, 12.91+1.62) %>% as_tibble %>% rename(cost_per_dose2=value)

# pull out univariate scenarios
scenarios %>%
  filter(intervention %in% c('RTS,S SV', 'RTS,S EPI')) %>%
  filter(resistance == 0 & seasonality == 'perennial') %>%
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




# impact RTSS on top of other interventions ------------------------------------
output <- scenarios %>%
  mutate(ID = paste(pfpr, seasonality, ITNuse, resistance, treatment, ITN, sep="_")) %>%
  filter(ITNuse==0.75)

none <- output %>%
  mutate(set = case_when(seasonality %in% c('perennial', 'highly seasonal') & intervention %in%
                           c('ITN 10% boost','ITN PBO') ~ 1,
                         seasonality=='seasonal' & intervention %in%
                           c('ITN 10% boost + SMC','ITN PBO + SMC') ~ 1)) %>%
  filter(set == 1) %>%
  dplyr::select(file, ID, daly, cases, cost_total, u5_dalys, n_0_1825) %>%
  rename(daly_baseline = daly,
         cases_baseline = cases,
         cost_total_baseline = cost_total,
         u5_daly_baseline = u5_dalys) %>%
  dplyr::select(file, ID, daly_baseline, cases_baseline, cost_total_baseline, u5_daly_baseline)

base_IDs <- none$file

output2 <- output %>% filter(!(file %in% base_IDs)) %>%
  mutate(set = case_when(seasonality %in% c('perennial', 'highly seasonal') & intervention %in%
                          c('ITN 10% boost + RTS,S','ITN PBO + RTS,S') ~ 1,
                         seasonality=='seasonal' & intervention %in%
                          c('ITN 10% boost + RTS,S + SMC','ITN PBO + RTS,S + SMC') ~ 1)) %>%
  filter(set == 1) %>%
  dplyr::select(file, ID, pfpr, seasonality, intervention, daly, cases, cost_total, u5_dalys) %>%
  left_join(none %>% dplyr::select(-file), by=c('ID')) %>%
  mutate(CE = (cost_total - cost_total_baseline) / (daly_baseline - daly),
         deltadaly = daly_baseline - daly,
         deltacases = cases_baseline - cases,
         CE_u5 = (cost_total - cost_total_baseline) / (u5_daly_baseline - u5_dalys))

summary(output2$deltadaly / 2)
summary(output2$CE)
summary(output2[output2$seasonality=='perennial',]$CE)
summary(output2[output2$seasonality=='seasonal',]$CE)
summary(output2[output2$seasonality=='highly seasonal',]$CE)

summary(output2$deltadaly / (2*15)) # additional dalys averted per year in a population of 100,000 people
summary(output2$deltacases / (2*15)) # additional cases averted per year in a population of 100,000 people
summary(output2$CE_u5) # additional cases averted per year in a population of 100,000 people

ggplot(data = output2) +
  geom_boxplot(aes(x=pfpr, y=deltadaly / (2*15), fill=as.factor(pfpr), color=as.factor(pfpr), group=pfpr), alpha = 0.3, show.legend = F) +
  facet_grid(~factor(seasonality, levels = c('perennial', 'seasonal', 'highly seasonal'))) +
  labs(x = 'PfPR', y = 'DALYs averted', title = 'DALYS averted per year per 100,000 people') +
  theme_classic()

ggsave('./03_output/RTSS_additional_impact.pdf', width=6, height=3)



# ICER table -------------------------------------------------------------------

# calculate change in dalys and cost
output <- scenarios %>%
  filter(cost_per_dose==6.52 & delivery_cost==1.62) %>%
  filter(intervention!='none') %>%
  mutate(deltadaly = daly_baseline - daly,
         deltacost = cost_total - cost_total_baseline,
         cost_daly_averted = (cost_total - cost_total_baseline) / deltadaly)

# all strategies, by baseline ITNuse and seasonality
output <- output %>% mutate(intervention2 = factor(intervention, levels = c('ITN 10% boost','ITN PBO','RTS,S EPI','RTS,S SV','ITN 10% boost + RTS,S','ITN PBO + RTS,S','SMC','ITN 10% boost + SMC','ITN PBO + SMC', 'RTS,S + SMC', 'ITN 10% boost + RTS,S + SMC','ITN PBO + RTS,S + SMC')))

RColorBrewer::brewer.pal(12, "Paired")
# colors <- c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C",'deeppink',"#CAB2D6","#6A3D9A",'black',"#FDBF6F","#FF7F00")
colors <- c("#1F78B4","#B2DF8A","#33A02C","#FB9A99", "#A6CEE3",
            "#E31A1C",'deeppink',"#CAB2D6","#6A3D9A")

# removing dominated strategies
final <- output %>%
  filter(deltadaly >= 0 & deltacost > 0) %>%
  # filter out mixed strategies
  group_by(ID) %>% arrange(ID, deltacost) %>%
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
  dplyr::select(ID, intervention, ICER, dominate)

merge <- output %>% left_join(final, by=c('ID', 'intervention')) %>% mutate(dominate = ifelse(is.na(dominate),1,dominate)) %>%
  mutate(resistance = ifelse(resistance!=0, 1, 0))

output %>% group_by(intervention) %>% filter(resistance==0) %>%
  summarize(cost_daly_averted = median(cost_daly_averted))

merge %>% group_by(intervention) %>%
  summarize(t = n(),
            ndominate = n()-sum(dominate),
            p_ndominate = ndominate / t*100,
            ICER_m = median(ICER, na.rm=T),
            ICER_25 = quantile(ICER, prob=0.25, na.rm=T),
            ICER_75 = quantile(ICER, prob=0.75, na.rm=T))

merge %>% filter(resistance==0) %>% group_by(intervention) %>%
  summarize(t = n(),
            ndominate = n()-sum(dominate),
            p_ndominate = ndominate / t*100,
            ICER_m = median(ICER, na.rm=T),
            ICER_25 = quantile(ICER, prob=0.25, na.rm=T),
            ICER_75 = quantile(ICER, prob=0.75, na.rm=T))

# ICER just among non-dominated strategies
merge %>% filter(dominate==0) %>% group_by(intervention) %>%
  summarize(ICER_m = median(ICER, na.rm=T),
            ICER_25 = quantile(ICER, prob=0.25, na.rm=T),
            ICER_75 = quantile(ICER, prob=0.75, na.rm=T))


merge %>% filter(dominate==0) %>%
  mutate(RTSSalone = ifelse(intervention %in% c('RTS,S age-based', 'RTS,S seasonal'), 1, 0)) %>%
  group_by(seasonality, ID) %>%
  summarize(rtss = sum(RTSSalone, na.rm=T)) %>%
  mutate(rtss = ifelse(rtss>=1,1,0)) %>%
  group_by(seasonality) %>%
  summarize(n=n(), t=sum(rtss, na.rm=T), p=t/n*100, q=100-p)



# ICERs among children ---------------------------------------------------------
scenarios %>% ungroup() %>%
  #group_by(seasonality) %>%
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

# Supplement -------------------------------------------------------------------

# < ITN dist vs ITN use ----------------------------------------------------------
# relationship between cost_ITN_linear and cost_ITN:
# read in netz package data to find the annual nets to distribute to give the simulated usage
nets_data <- netz::prepare_data()

# get nets to be distributed for each ITN usage
ndist <- function(x) {

  convert_usage_to_annual_nets_distributed(
    target_usage = seq(0,0.85,0.05),
    distribution_freq = 1095, # 3 years
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
net <- c('pyrethroid', 'pyrethroid + PBO')
output <- crossing(nets_distributed, net)

output <- output %>%
  mutate(ITNcost = case_when(net=='pyrethroid' ~ 3.50,            # $2.00 per net and $1.50 delivery cost
                             net=='pyrethroid + PBO' ~ 3.80)) %>% # $2.30 per net and $1.50 delivery cost
  mutate(cost_ITN = annual_percapita_nets_distributed * ITNcost * 3,  # true net cost accounting for non-linear relationship
         cost_ITNmin = annual_percapita_nets_distmin * ITNcost * 3,  # true net cost MIN
         cost_ITNmax = annual_percapita_nets_distmax * ITNcost * 3,  # true net cost MAX
         cost_ITNmin = ifelse(is.na(cost_ITNmin), cost_ITN, cost_ITNmin),
         cost_ITN_linear = target_use * ITNcost)         # ITN linear

output2 <- output %>%
  filter((net=='pyrethroid' & round(target_use,2) %in% c(0.25, 0.35, 0.50, 0.60, 0.75, 0.85)) |
           (net=='pyrethroid + PBO' & round(target_use,2) %in% c(0.25, 0.5, 0.75))) %>%
  mutate(shape='sim')

# plot
ggplot(output, aes(x=cost_ITN_linear, y = cost_ITN, color=target_use)) +
  geom_point(data=output2, aes(x=cost_ITN_linear, y = cost_ITN, shape=shape), size=4, color='chartreuse1', fill='chartreuse1') +
  geom_point(data=output2, aes(x=cost_ITN_linear, y = cost_ITN, shape=shape), size=3.5, color='chartreuse1', fill='chartreuse1')+
  geom_abline(slope=1, size=1, lty=2, alpha=0.3) +
  geom_point(size=2) +
  geom_line(alpha=0.5, size=1) +
  geom_ribbon(aes(ymin=cost_ITNmin, ymax=cost_ITNmax), color=NA, fill='cornflower blue', alpha=0.4) +
  scale_color_gradient(limits=c(0,0.85)) +
  scale_shape_manual(values=c(1), labels=c('value used in simulation')) +
  facet_grid(~net) +
  labs(x='Linear ITN cost', y='Netz ITN cost', color='target ITN use', shape='') +
  theme_bw()

ggsave('./03_output/ITN_netz.pdf', width=8, height=4)


# < PfPR by country --------------------------------------------------------------
nets_data <- read_csv('./01_data/Intervention_ITN.csv')

PR <- read_csv('./01_data/00_PfPR_table_Global_admin0_2000-2019.csv') %>%
  filter(Year == 2019) %>% arrange(PfPR_median) %>%
  mutate(SSA = case_when(Name_0 %in% nets_data$Name ~ 1,
                         grepl('Cape|Comor|Ivoire|Lesotho|Maurit|Tome|Seyc|Tanz', Name_0) ~ 1)) %>%
  filter(SSA == 1) %>%
  mutate(Name_f = factor(Name_0, levels=Name_0))


lacroix_palettes$Pamplemousse

ggplot(PR) +
  geom_rect(xmin=1, xmax=49, ymin=0.09, ymax=.11, fill="#F6A1A5", alpha=0.009) +
  geom_rect(xmin=1, xmax=49, ymin=.19, ymax=.21, fill="#F8CD9C", alpha=0.01) +
  geom_rect(xmin=1, xmax=49, ymin=.39, ymax=.41, fill="#1BB6AF", alpha=0.005) +
  geom_pointrange(aes(x=Name_f, y = PfPR_median, ymin = PfPR_LCI, ymax = PfPR_UCI), color='#088BBE') +
  labs(x='Country', y=expression(italic(Pf)~PR[2-10]~', 2019')) +
  scale_y_continuous(breaks = c(0,0.1,0.2,0.3,0.4)) +
  scale_x_discrete(labels = scales::wrap_format(20)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave('./03_output/PfPR.pdf', width=10, height=4)


# < PfPR and mortality by country ------------------------------------------------
nets_data <- read_csv('./01_data/Intervention_ITN.csv')

PR <- read_csv('./01_data/00_PfPR_table_Global_admin0_2000-2019.csv') %>%
  filter(Year == 2019) %>% dplyr::select(ISO, Name_0, starts_with('PfPR'))

mortality <- read_csv('./01_data/00_Pf_mortality_rate_table_Global_admin0_2000-2019.csv') %>%
  filter(Year == 2019) %>% dplyr::select(ISO, Name_0, starts_with('mortality'))

PRmort <- PR %>% left_join(mortality) %>%
  arrange(PfPR_median) %>%
  mutate(SSA = case_when(Name_0 %in% nets_data$Name ~ 1,
                         grepl('Cape|Comor|Ivoire|Lesotho|Maurit|Tome|Seyc|Tanz', Name_0) ~ 1)) %>%
  filter(SSA == 1) %>%
  mutate(Name_f = factor(Name_0, levels=Name_0))

lacroix_palettes$Pamplemousse

ggplot(PRmort) +
  geom_rect(xmin=1, xmax=49, ymin=0.09, ymax=.11, fill="#F6A1A5", alpha=0.009) +
  geom_rect(xmin=1, xmax=49, ymin=.19, ymax=.21, fill="#F8CD9C", alpha=0.01) +
  geom_rect(xmin=1, xmax=49, ymin=.39, ymax=.41, fill="#1BB6AF", alpha=0.005) +
  geom_pointrange(aes(x=Name_f, y = PfPR_median, ymin = PfPR_LCI, ymax = PfPR_UCI, color='PfPR'), alpha=0.4) +
  geom_pointrange(aes(x=Name_f, y = mortality_rate_median*200, ymin = mortality_rate_LCI*200, ymax = mortality_rate_UCI*200, color='mortality'), alpha=0.4) +
  labs(x='Country', y=expression(italic(Pf)~PR[2-10]~', 2019'), color='') +
  scale_y_continuous(sec.axis = sec_axis(~ . / 200,
                                         name = expression(italic(Pf)~'mortality rate (all ages)')),
                     breaks = c(0,0.1,0.2,0.3,0.4)) +
  scale_x_discrete(labels = scales::wrap_format(20)) +
  scale_color_manual(values = c('#088BBE','blue')) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave('./03_output/PfPR_mortality.pdf', width=12, height=4)


# < ITN usage and countries -----------------------------------------------------
# https://malariaatlas.org/research-project/the-impact-of-malaria-control-on-plasmodium-falciparum-in-africa-2000-2015/
lacroix_palettes$Pamplemousse

# import ITN data from MAP
nets_data <- read_csv('./01_data/Intervention_ITN.csv') %>% arrange(`2015`) %>%
  mutate(Name_f = factor(Name, levels=Name))

ggplot(nets_data) +
  geom_rect(xmin=1, xmax=43, ymin=-0.01, ymax=.25, fill="#F6A1A5", alpha=0.009) +
  geom_rect(xmin=1, xmax=43, ymin=.25, ymax=.50, fill="#F8CD9C", alpha=0.01) +
  geom_rect(xmin=1, xmax=43, ymin=.50, ymax=.75, fill="#1BB6AF", alpha=0.005) +
  geom_rect(xmin=1, xmax=43, ymin=.75, ymax=.90, fill="#088BBE", alpha=0.006) +
  geom_point(aes(x=Name_f, y = `2015`)) +
  labs(x='Country', y='ITN use by country, 2015') +
  scale_x_discrete(labels = scales::wrap_format(20)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# OR try MAP data
SSA_ITN <- readRDS('./03_output/MAP_ITN.rds') %>% arrange(median) %>%
  filter(!is.na(median)) %>%
  mutate(name_0 = factor(name_0, levels=name_0))

ggplot(SSA_ITN) +
  geom_rect(xmin=1, xmax=47, ymin=-0.01, ymax=.25, fill="#F6A1A5", alpha=0.009) +
  geom_rect(xmin=1, xmax=47, ymin=.25, ymax=.50, fill="#F8CD9C", alpha=0.01) +
  geom_rect(xmin=1, xmax=47, ymin=.50, ymax=.75, fill="#1BB6AF", alpha=0.005) +
  geom_rect(xmin=1, xmax=47, ymin=.75, ymax=.90, fill="#088BBE", alpha=0.006) +
  geom_pointrange(aes(x=name_0, y=median, ymin=min, ymax=max), lwd=.3) +
  labs(x='Country', y='ITN use by country, 2019') +
  scale_x_discrete(labels = scales::wrap_format(20)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


ggsave('./03_output/ITN_usage.pdf', width=10, height=4)




# Case Study -------------------------------------------------------------------

# < point range outcomes -------------------------------------------------------
scenarios2 <- readRDS('./03_output/scenarios2_casestudy.rds') %>%
  filter(seasonality == 'highly seasonal' & scenario2_f == 'gap in PfPR')

plot_pointrange <- function(y, ymin, ymax, var) {

  y <- sym(y)
  ymin <- sym(ymin)
  ymax <- sym(ymax)

  scenarios2 %>%
    ggplot(aes(color=scenario_f, group=scenario_f)) +
    geom_pointrange(aes(x=scenario, y=!!y, ymin=!!ymin, ymax=!!ymax), alpha=0.7) +
    geom_hline(yintercept = 0, lty=2, color='grey') +
    scale_color_manual(values = c("#EA7580","#1BB6AF","#F6A1A5","#088BBE")) +
    scale_x_continuous(limits=c(0.5,4.4), breaks=c(1,2,3,4)) +
    # facet_grid(.~seasonality) +
    labs(x='',
         y=expr(paste(Delta," cost / ", Delta, !!var)),
         fill = '',
         color = '') +
    theme_classic()
}

lacroix_palettes$Pamplemousse

A <- plot_pointrange('CE_daly', 'CE_daly_lower', 'CE_daly_upper', " DALYs")
B <- plot_pointrange('CE_case', 'CE_case_lower', 'CE_case_upper', " cases")
C <- plot_pointrange('CE_death', 'CE_death_lower', 'CE_death_upper', " deaths")
D <- plot_pointrange('CE_daly_u5', 'CE_daly_u5_lower', 'CE_daly_u5_upper', " u5 DALYs")
E <- plot_pointrange('CE_u5_case', 'CE_u5_case_lower', 'CE_u5_case_upper', " u5 cases")
F <- plot_pointrange('CE_u5_death', 'CE_u5_death_lower', 'CE_u5_death_upper', " u5 deaths")

(A + B + C) / (D + E + F) + plot_layout(guides = "collect", nrow=2) + plot_annotation(tag_levels = 'A')

ggsave('./03_output/pointrange_casestudy.pdf', width=8, height=4)


# < incidence plot -------------------------------------------------------------
output <- readRDS('./03_output/scenarios_casestudy.rds')

# gap in PfPR
output_pfpr <- output %>%
  filter(pfpr %in% c(0.40, 0.10) & RTSScov %in% c(0, 0.80) & ITNuse == 0.50) %>%
  mutate(scenario2 = 1)

# gap in ITN use
output_itn <- output %>%
  filter(RTSScov %in% c(0, 0.80) & ((pfpr == 0.20 & ITNuse == 0.60) | (pfpr == 0.10 & ITNuse == 0.30))) %>%
  mutate(scenario2 = 2)

# gap in RTSS coverage
output_rtss <- output %>%
  filter(ITNuse == 0.50 & ((pfpr == 0.20 & RTSScov %in% c(0, 0.50)) | (pfpr == 0.10 & RTSScov %in% c(0, 0.80)))) %>%
  mutate(scenario2 = 3)

# combine
output <- full_join(output_pfpr, output_itn) %>% full_join(output_rtss) %>%
  mutate(scenario2_f = factor(scenario2,
                              levels=c(1,2,3),
                              labels=c('gap in PfPR', 'gap in ITN use', 'gap in vaccination')))


  # assign scenarios
  scenario1 <- output %>%
    filter(ITNboost==0 & RTSS=='none') %>% mutate(scenario=1)

  scenario2 <- output %>%
    filter(ITNboost==1 & RTSS=='none') %>% mutate(scenario=2)

  scenario3 <- output %>%
    filter(RTSS=='EPI' & ITNboost==0) %>% mutate(scenario=3)

  scenario4 <- output %>%
    filter((pfpr %in% c(0.20, 0.40) & ITNboost==1 & RTSS=='none') | (pfpr %in% c(0.10) & ITNboost==0 & RTSS=='none')) %>%
    mutate(scenario=4)

  scenario5 <- output %>%
    filter((pfpr %in% c(0.20, 0.40) & ITNboost==0 & RTSS=='EPI') | (pfpr %in% c(0.10) & ITNboost==0 & RTSS=='none')) %>%
    mutate(scenario=5)

  scenarios <- full_join(scenario1, scenario2) %>% full_join(scenario3) %>%
    full_join(scenario4) %>% full_join(scenario5) %>%
    mutate(scenario_f = factor(scenario,
                               levels=c(1,2,3,4,5),
                               labels=c('baseline','mass ITN 10% increase', 'mass age-based RTS,S', 'targeted ITN 10% increase', 'targeted age-based RTS,S'))) %>%
  mutate(
          clin_inc = cases / n,
          clin_inc_lower = cases_lower / n,
          clin_inc_upper = cases_upper / n,

          clin_inc_u5 = u5_cases / n_0_1825,
          clin_inc_u5_lower = u5_cases_lower / n_0_1825,
          clin_inc_u5_upper = u5_cases_upper / n_0_1825,

          sev_inc = severe_cases / n,
          sev_inc_u5 = u5_severe / n_0_1825) %>%

  mutate(residence = ifelse(pfpr %in% c(0.20, 0.40), 'rural', 'urban')) %>%
  dplyr::select(seasonality, scenario, scenario_f, scenario2, scenario2_f, residence, cost_total, clin_inc:sev_inc_u5) %>%
  pivot_longer(cols = clin_inc:sev_inc_u5, names_to = 'var', values_to = 'value') %>%
  mutate(estimate = ifelse(!(grepl('_lower',var) | grepl('_upper',var)), value, NA),
         lower = ifelse(grepl('_lower',var), value, NA),
         upper = ifelse(grepl('_upper',var), value, NA),
         var = case_when(grepl('clin_inc_u5',var) ~ 'clinical incidence <5s',
                         grepl('clin_inc',var) ~ 'clinical incidence total',
                         grepl('sev_inc_u5',var) ~ 'severe incidence <5s',
                         grepl('sev_inc',var) ~ 'severe incidence total')) %>%
  dplyr::select(-value) %>%
  group_by(seasonality, scenario2, scenario2_f, scenario, scenario_f, residence, var) %>%
  summarize(estimate = mean(estimate, na.rm=T),
         lower = mean(lower, na.rm=T),
         upper = mean(upper, na.rm=T),
         cost_total = mean(cost_total))


lacroix_palettes$Pamplemousse

scenarios %>%
  filter(var == 'clinical incidence total') %>%
  ggplot(aes(color=residence)) +
  geom_vline(aes(xintercept=0), lty=2, color='darkgrey', size=.8) +
  geom_pointrange(aes(y=(scenario)*-1, x=estimate, xmin=lower, xmax=upper), alpha=0.7,
                  position = position_dodge(width = .2)) +
  scale_color_manual(values = c("#088BBE","#172869")) +
  facet_grid(seasonality~scenario2_f, scales = 'free') +
  scale_x_continuous(limits = c(0, max(scenarios$upper)), breaks = scales::trans_breaks(identity, identity, n = 3)) +
  scale_y_continuous(breaks = seq(-5,-1,1),
                     labels = c('targeted RTS,S', 'targeted ITNs', 'universal RTS,S', 'universal ITNs', 'baseline')) +
  labs(x='',
       y='clinical incidence (per person)',
       fill = '',
       color = '') +
  theme_classic() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey96"))

ggsave('./03_output/inc_outcomes_casestudy_season.pdf', width=8, height=4)

scenarios %>%
  filter(var == 'clinical incidence total') %>%
  ggplot(aes(color=scenario_f, shape=residence)) +
  geom_vline(aes(xintercept=0), lty=2, color='darkgrey', size=.8) +
  geom_pointrange(aes(y=(scenario)*-1, x=estimate, xmin=lower, xmax=upper), alpha=0.7,
                  position = position_dodge(width = .2)) +
  scale_color_manual(values = c('black',"#EA7580","#1BB6AF","#F8CD9C","#088BBE")) +
  facet_grid(seasonality~scenario2_f, scales = 'free') +
  scale_x_continuous(limits = c(0, max(scenarios$upper)), breaks = scales::trans_breaks(identity, identity, n = 3)) +
  scale_shape_manual(values = c(19,17), guide = guide_legend(reverse = TRUE)) +
  scale_y_continuous(breaks = seq(-5,-1,1),
                     labels = NULL) +
  labs(y='',
       x='clinical incidence',
       fill = '',
       color = '') +
  theme_classic() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey96"))

ggsave('./03_output/inc_outcomes_casestudy_season2.pdf', width=8, height=4)



# < equality measures ----------------------------------------------------------
# change in daly / change in cost
readRDS('./03_output/scenarios2_casestudy.rds') %>% ungroup() %>%
  dplyr::select(seasonality, scenario2_f, scenario_f, CE_daly) %>%
  arrange(seasonality, scenario2_f, scenario_f) %>%
  pivot_wider(names_from = scenario_f, values_from = CE_daly) %>%
  write.table("clipboard", sep="\t")

lacroix_palettes$Pamplemousse

A <- readRDS('./03_output/scenarios2_casestudy.rds')  %>%
  ggplot(aes(color=scenario_f, group=scenario_f)) +
  geom_pointrange(aes(y=scenario*-1, x=CE_case, xmin=CE_case_lower, xmax=CE_case_upper), alpha=0.7) +
  geom_vline(xintercept = 0, lty=2, color='grey') +
  scale_color_manual(values = c("#EA7580","#1BB6AF","#F8CD9C","#088BBE")) +
  scale_y_continuous(limits=c(-4.4, -0.5), breaks=c(-4,-3,-2,-1), labels = NULL) +
  facet_grid(seasonality~scenario2_f, scales = 'free') +
  labs(y='',
       x=expr(paste(Delta," cost / ", Delta, " cases")),
       fill = '',
       color = '') +
  theme_classic() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey96"))


B <- scenarios %>%
  group_by(seasonality, scenario, scenario2_f, scenario_f, var) %>%
  arrange(seasonality, scenario2_f, scenario, scenario_f, var, residence) %>%
  filter(var %in% c('clinical incidence total')) %>%
  summarize(change_estimate = (estimate - lead(estimate))*2000000,
            change_estimate_lower = (lower - lead(lower))*2000000,
            change_estimate_upper = (upper - lead(upper))*2000000,
            sum_cost = cost_total + lead(cost_total)) %>%
  filter(!is.na(change_estimate)) %>%
  group_by(seasonality, scenario2_f, var) %>% arrange(var) %>%
  mutate(base_estimate = head(change_estimate, n=1),
         base_estimate_lower = head(change_estimate_lower, n=1),
         base_estimate_upper = head(change_estimate_upper, n=1),
         base_cost = head(sum_cost, n=1),
         per_change = (base_estimate - change_estimate),
         change_cost = (sum_cost - base_cost),
         equality =  change_cost / per_change,
         equality_lower = change_cost / (base_estimate_lower - change_estimate_lower),
         equality_upper = change_cost / (base_estimate_upper - change_estimate_upper)) %>%
  dplyr::select(seasonality, scenario, scenario2_f, scenario_f, equality, equality_lower, equality_upper) %>%
  filter(scenario_f != 'baseline') %>%
  filter(equality >= 0 ) %>%

  ggplot(aes(color=scenario_f)) +
  geom_vline(xintercept = 0, lty=2, color='grey') +
  geom_pointrange(aes(y=scenario*-1, x=equality, xmin=equality_lower, xmax=equality_upper), alpha=0.7) +
  scale_color_manual(values = c("#EA7580","#1BB6AF","#F8CD9C","#088BBE")) +
  facet_grid(seasonality~scenario2_f, scales = 'free') +
  scale_x_continuous(limits = c(min(scenarios$upper), max(scenarios$upper)), breaks = scales::trans_breaks(identity, identity, n = 3)) +
  scale_y_continuous(limits=c(-5.4, -1.5), breaks=c(-5,-4,-3,-2), labels = NULL) +
  labs(y='',
       x=expr(paste(Delta," cost / ", Delta, " difference in cases (urban vs. rural)")),
       fill = '',
       color = '') +
  theme_classic() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey96"))

(A / B) + plot_layout(guides = "collect", nrow=2) + plot_annotation(tag_levels = 'A')

ggsave('./03_output/pointrange_casestudy_season.pdf', width=9, height=6)


# % change in outcome / % change in cost
readRDS('./03_output/scenarios2_casestudy.rds') %>% ungroup() %>%
  dplyr::select(seasonality, scenario2_f, scenario_f, CE_daly) %>%
  arrange(seasonality, scenario2_f, scenario_f) %>%
  pivot_wider(names_from = scenario_f, values_from = CE_daly) %>%
  write.table("clipboard", sep="\t")

readRDS('./03_output/scenarios2_casestudy.rds') %>% ungroup() %>%
  dplyr::select(seasonality, scenario2_f, scenario_f, CE_case) %>%
  arrange(seasonality, scenario2_f, scenario_f) %>%
  pivot_wider(names_from = scenario_f, values_from = CE_case)


scenarios %>% group_by(seasonality, scenario2_f, scenario_f, var) %>%
  arrange(seasonality, scenario2_f, scenario_f, var, residence) %>%
  filter(var %in% c('clinical incidence total')) %>%
  summarize(change_estimate = (estimate - lead(estimate))*2000000,
            sum_cost = cost_total + lead(cost_total)) %>%
  filter(!is.na(change_estimate)) %>%
  group_by(seasonality, scenario2_f, var) %>% arrange(var) %>%
  mutate(base_estimate = head(change_estimate, n=1),
         base_cost = head(sum_cost, n=1),
         per_change = (base_estimate - change_estimate),
         change_cost = (sum_cost - base_cost),
         equality =  change_cost / per_change) %>%
  dplyr::select(seasonality, scenario2_f, scenario_f, equality) %>%
  pivot_wider(names_from = scenario_f, values_from = equality)




scenarios %>% group_by(seasonality, scenario2_f, scenario_f, var) %>%
  arrange(seasonality, scenario2_f, scenario_f, var, residence) %>%
  filter(var %in% c('clinical incidence total', 'severe incidence <5s')) %>%
  summarize(change_estimate = estimate - lead(estimate),
            sum_cost = cost_total + lead(cost_total)) %>%
  filter(!is.na(change_estimate)) %>%
  group_by(seasonality, scenario2_f, var) %>% arrange(var) %>%
  mutate(base_estimate = head(change_estimate, n=1),
         base_cost = head(sum_cost, n=1),
         per_change = (base_estimate - change_estimate) / base_estimate*100,
         # per_cost = (sum_cost - base_cost) / base_cost*100,
         change_cost = (sum_cost - base_cost) / 1000000,
         equality =  per_change / change_cost) %>%
  dplyr::select(seasonality, scenario2_f, scenario_f, equality) %>%
  pivot_wider(names_from = scenario_f, values_from = equality) %>%
  write.table("clipboard",sep="\t")




# < CE table -------------------------------------------------------------------
# dalys / cases / deaths
scenarios2 <- readRDS('./03_output/scenarios2_admin1.rds')

scenarios2 %>%
  dplyr::select(scenario, scenario_f, CE_daly, CE_case, CE_death, CE_daly_u5, CE_u5_case, CE_u5_death) %>%
  write.table("clipboard",sep="\t")


# equity
# continue using scenarios datasets from above
scenarios <- full_join(scenario0, scenario1) %>%
  full_join(scenario2) %>%
  full_join(scenario3) %>%
  full_join(scenario4) %>%
  mutate(scenario_f = factor(scenario, levels=c(1,2,3,4,5),
                             labels=c('baseline','mass ITN boost', 'mass age-based RTS,S', 'targeted ITN boost', 'targeted age-based RTS,S'))) %>%
  mutate(
    clin_inc = cases / n,
    clin_inc_lower = cases_lower / n,
    clin_inc_upper = cases_upper / n,

    clin_inc_u5 = u5_cases / n_0_1825,
    clin_inc_u5_lower = u5_cases_lower / n_0_1825,
    clin_inc_u5_upper = u5_cases_upper / n_0_1825,

    sev_inc = severe_cases / n,
    sev_inc_u5 = u5_severe / n_0_1825)

scenarios %>%
  mutate(name = paste0(pfpr,'_',scenario)) %>%
  dplyr::select(name, clin_inc, clin_inc_u5,
    sev_inc, sev_inc_u5) %>%
  pivot_longer(cols = clin_inc:sev_inc_u5, names_to = 'var', values_to = 'value') %>%
  pivot_wider(names_from = name, values_from = value) %>%
  write.table("clipboard",sep="\t")

scenarios %>%
  mutate(name = paste0(pfpr,'_',scenario)) %>%
  dplyr::select(name, cost_total) %>%
  pivot_wider(names_from = name, values_from = cost_total) %>%
  write.table("clipboard",sep="\t")

output <- scenarios %>% dplyr::select(pfpr, scenario, cost_total, clin_inc, clin_inc_u5,
                                      sev_inc, sev_inc_u5)

none <- output %>% filter(scenario==1) %>%
  rename(cost_total_baseline = cost_total,
         clin_inc_baseline = clin_inc,
         clin_inc_u5_baseline = clin_inc_u5,
         sev_inc_baseline = sev_inc,
         sev_inc_u5_baseline = sev_inc_u5)

output2 <- output %>% filter(scenario!=1) %>%
  full_join(none %>% dplyr::select(-scenario)) %>%
  mutate(CE_clin_inc = (cost_total - cost_total_baseline) / (clin_inc_baseline - clin_inc) / 100000,
         CE_clin_inc_u5 = (cost_total - cost_total_baseline) / (clin_inc_u5_baseline - clin_inc_u5) / 100000,
         CE_sev_inc = (cost_total - cost_total_baseline) / (sev_inc_baseline - sev_inc) / 100000,
         CE_sev_inc_u5 = (cost_total - cost_total_baseline) / (sev_inc_u5_baseline - sev_inc_u5) / 100000)

# ------------------------------------------------------------------------------



