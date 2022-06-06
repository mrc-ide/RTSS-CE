# Figures & Tables -------------------------------------------------------------
# packages
library(tidyverse)
library(data.table)
library(patchwork)
library(grid)
library(LaCroixColoR)

# color palette
# lacroix_palette(type = "paired")
# c("#C70E7B", "#FC6882","#007BC3","#54BCD1","#EF7C12","#F4B95A","#009F3F","#8FDA04",
#   "#AF6125","#F4E3C7","#B25D91","#EFC7E6", "#EF7C12","#F4B95A")
# RColorBrewer::display.brewer.pal(12, "Paired")

# assign scenarios colors
smc <- "#FC6882"
pbo <- "#007BC3"
itn <- "#54BCD1"
rtss_sv <- "#009F3F"
rtss_age <- "#54E356"
pbo_smc <- "#C70E7B"
itn_smc <- "#B25D91"
rtss_smc <- "#EFC7E6"
pbo_rtss_smc <- "#AF6125"
itn_rtss_smc <- "#F5DD42"
pbo_rtss <- "#EF7C12"
itn_rtss <- "#F4B95A"

# devtools::install_github('mrc-ide/malariasimulation@dev', force=TRUE)
# devtools::install_github('johannesbjork/LaCroixColoR')

# load data
scenarios <- readRDS('./03_output/scenarios_draws.rds')
dalyoutput_cost <- readRDS('./03_output/dalyoutput_draws.rds')

test <- dalyoutput_cost %>%
  group_by(ID, pfpr, seasonality, treatment, resistance, ITN, ITNuse, ITNboost, SMC, RTSS) %>%
  # filter(cost_per_dose==6.52 & delivery_cost==1.62) %>%
  summarize(n = n()) # 2,574 scenarios matches median results

test <- dalyoutput_cost %>%
  filter(ITNboost==0 & ITN=='pyr' & RTSS=='none' & # filter out interventions
           (SMC==0 | (seasonality=='highly seasonal'))) %>%
  filter(cost_per_dose==6.52 & delivery_cost==1.62) %>%
  summarize(n = n()); nrow(test) / 50 # 270 baseline scenarios matches median results

test <- scenarios %>%
  group_by(ID, pfpr, seasonality, treatment, resistance, ITN, ITNuse, ITNboost, SMC, RTSS) %>%
  summarize(n = n()) # 2,304 scenarios matches median results


# delta CE change --------------------------------------------------------------

# < by season ------------------------------

supp.labs <- c("ITN use: 0", "ITN use: 0.25", "ITN use: 0.50", "ITN use: 0.75")
names(supp.labs) <- c(0,.25,.50,.75)

deltaseason <- function(season) {

  output <- scenarios %>%
    filter(seasonality==season) %>%
    filter(cost_per_dose==6.52 & delivery_cost==1.62) %>% filter(resistance==0) %>%
    filter(intervention!='none') %>%
    mutate(deltadaly = daly_baseline - daly,
           deltacost = cost_total - cost_total_baseline,
           cost_daly_averted = (cost_total - cost_total_baseline) / deltadaly) %>%
    group_by(ID, seasonality, intervention, intervention_f) %>%
    summarize(deltadaly = median(deltadaly),
           deltacost = median(deltacost),
           cost_daly_averted = median(cost_daly_averted))

  # All interventions, by baseline ITNuse and seasonality
  output <- output %>% mutate(intervention2 = factor(intervention, levels = c('ITN 10% boost','ITN PBO','RTS,S EPI','RTS,S SV','ITN 10% boost + RTS,S','ITN PBO + RTS,S','SMC','ITN 10% boost + SMC','ITN PBO + SMC', 'RTS,S + SMC', 'ITN 10% boost + RTS,S + SMC','ITN PBO + RTS,S + SMC')))

  # assign colors
  if(season=='perennial'){
    colors <- c(itn, rtss_age, itn_rtss)
  }

  if(season=='highly seasonal'){
    colors <- c(itn, rtss_age, rtss_sv, itn_rtss)
  }

  if(season=='seasonal'){
    colors <- c(itn, rtss_age, rtss_sv, smc, itn_rtss, itn_smc, rtss_smc, itn_rtss_smc)
  }

  A <- ggplot(data = output, mapping=aes(x=deltadaly, y=deltacost)) +
    geom_line(aes(group=as.factor(ID)), color='lightgrey', size=.5) +
    geom_point(aes(color=intervention_f), size=2) +
    geom_hline(yintercept = 0, lty=2, color="black") +
    geom_vline(xintercept = 0, lty=2, color="black") +
    #facet_grid(~ITNuse, labeller = labeller(ITNuse=supp.labs), scales = "free") +
    theme_classic() +
    scale_color_manual(values = colors) +
    theme(axis.text.x = element_text(angle = 45, hjust=1)) +
    labs(title='All strategies',
         y='change in cost (USD)',
         x='change in DALYs averted',
         color='intervention')

  if(season=='perennial'){
    colors <- c(itn, rtss_age, itn_rtss)
  }

  if(season=='highly seasonal'){
    colors <- c(itn, rtss_age, rtss_sv, itn_rtss)
  }

  if(season=='seasonal'){
    colors <- c(itn, smc, itn_smc, rtss_smc, itn_rtss_smc)
  }

  # removing dominated strategies
  B <- output %>%
    # filter out mixed strategies
    group_by(ID) %>% arrange(ID, deltacost) %>%
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
    ggplot(mapping=aes(x=deltadaly, y=deltacost)) +
    geom_line(aes(group=as.factor(ID)), color='lightgrey', size=.5) +
    geom_point(aes(color=intervention_f), size=2, show.legend = F) +
    geom_hline(yintercept = 0, lty=2, color="black") +
    geom_vline(xintercept = 0, lty=2, color="black") +
    #facet_grid(~ITNuse, labeller = labeller(ITNuse=supp.labs), scales = "free") +
    theme_classic() +
    scale_color_manual(values = colors) +
    theme(axis.text.x = element_text(angle = 45, hjust=1)) +
    labs(title='Dominated strategies removed',
         y='change in cost (USD)',
         x='change in DALYs averted',
         color='intervention',
         caption='Assuming resistance == 0') +
    theme(plot.caption.position = "plot")

  A + B + plot_layout(guides = "collect", nrow=1) + plot_annotation(tag_levels = 'A')

  # ggsave(paste0('./03_output/impact_cloud_',season,'.pdf'), width=10, height=6)
  ggsave(paste0('./03_output/plots_draws/impact_cloud_',season,'.pdf'), width=12, height=5)

}

deltaseason('highly seasonal')
deltaseason('seasonal')
deltaseason('perennial')


# < facet by season -----------------
output <- scenarios %>%
  filter(cost_per_dose==6.52 & delivery_cost==1.62) %>% filter(resistance==0) %>%
  filter(intervention!='none') %>%
  mutate(deltadaly = daly_baseline - daly,
         deltacost = cost_total - cost_total_baseline,
         cost_daly_averted = (cost_total - cost_total_baseline) / deltadaly) %>%
  mutate(seasonality = factor(seasonality, levels = c('perennial', 'seasonal', 'highly seasonal'))) %>%
  group_by(ID, seasonality, intervention, intervention_f) %>%
  summarize(deltadaly = median(deltadaly),
            deltacost = median(deltacost),
            cost_daly_averted = median(cost_daly_averted))

# All interventions, by baseline ITNuse and seasonality
output <- output %>% mutate(intervention2 = factor(intervention, levels = c('ITN 10% increase','ITN PBO','RTS,S age-based','RTS,S seasonal','ITN 10% increase + RTS,S','ITN PBO + RTS,S','SMC','ITN 10% increase + SMC','ITN PBO + SMC', 'RTS,S + SMC', 'ITN 10% increase + RTS,S + SMC','ITN PBO + RTS,S + SMC')))


colors <- c(itn, rtss_age, rtss_sv, smc, itn_rtss, itn_smc, rtss_smc, itn_rtss_smc)


A <- ggplot(data = output, mapping=aes(x=deltadaly, y=deltacost)) +
  geom_line(aes(group=as.factor(ID)), color='lightgrey', size=.5, alpha=0.2) +
  geom_point(aes(color=intervention_f), size=2, alpha=0.6) +
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
B <- output %>%
  # filter out mixed strategies
  group_by(ID) %>% arrange(ID, deltacost) %>%
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
  ggplot(mapping=aes(x=deltadaly, y=deltacost)) +
  geom_line(aes(group=as.factor(ID)), color='lightgrey', size=.5, alpha=0.2) +
  geom_point(aes(color=intervention_f), size=2, alpha=0.6, show.legend = F) +
  geom_hline(yintercept = 0, lty=2, color="black") +
  geom_vline(xintercept = 0, lty=2, color="black") +
  facet_grid(~seasonality, scales = "free") +
  theme_classic() +
  scale_color_manual(values = colors) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(title='Dominated strategies removed',
       y='change in cost (USD)',
       x='change in DALYs averted',
       color='intervention',
       caption='Assuming resistance == 0') +
  theme(plot.caption.position = "plot")

A + B + plot_layout(guides = "collect", nrow=2) + plot_annotation(tag_levels = 'A')

ggsave(paste0('./03_output/plots_draws/impact_cloud_all.pdf'), width=12, height=7)


# box and whisker delta cost / delta daly --------------------------------------
tabledat <- scenarios %>%
  filter(cost_per_dose == 6.52 & delivery_cost == 1.62) %>%
  group_by(intervention_f) %>%
  summarize(n=n(),
            median = round(median(CE)),
            q25 = round(quantile(CE, p=0.25)),
            q75 = round(quantile(CE, p=0.75)),
            min = round(min(CE)),
            max = round(max(CE))) %>%
  arrange(median)

tabledat; sum(tabledat$n)


# set up order of interventions by median CE for plot
levels <- scenarios %>%
  mutate(group = case_when(intervention_f %in% c('ITN 10% increase', 'ITN PBO', 'SMC', 'RTS,S seasonal', 'RTS,S age-based') ~ 'uni',
                           TRUE ~ 'multi')) %>%
  group_by(intervention_f, group) %>%
  summarize(med = median(CE, na.rm=T) )%>% ungroup() %>%
  arrange(desc(group), med)

# plot of cost per DALY averted
box_plot <- function(cost_dose){
  scenarios %>%
    filter(cost_per_dose == cost_dose & delivery_cost == 1.62) %>%
    mutate(intervention_f = factor(intervention, levels=levels$intervention_f)) %>%
    mutate(rank=as.numeric(intervention_f)) %>%

    ggplot(aes(x=rank, y=CE, fill=intervention_f, color=intervention_f, group=intervention)) +
    geom_hline(yintercept = 0, lty=2, color='grey') +
    geom_vline(xintercept = 5.5, lty=2, color='grey') +
    geom_boxplot(alpha=0.3, outlier.alpha = 0.05,  outlier.size = 0.05) +
    coord_cartesian(ylim=c(-100, 500), clip="off") +
    labs(x='',
         y=expression(paste(Delta," cost / ", Delta, " DALYs")),
         fill = 'intervention',
         # caption=paste0('range in cost / DALYs: ', round(min(scenarios$CE, na.rm = T)), ' to ', round(max(scenarios$CE, na.rm=T))),
         color = 'intervention') +
    annotation_custom(textGrob("Univariate strategies"),xmin=1,xmax=5,ymin=-150,ymax=-150) +
    annotation_custom(textGrob("Mixed strategies"),xmin=6,xmax=12,ymin=-150,ymax=-150) +
    scale_fill_manual(values = c(smc, pbo, itn, rtss_sv, rtss_age, pbo_smc, itn_smc, rtss_smc, pbo_rtss_smc, itn_rtss_smc, pbo_rtss, itn_rtss)) +
    scale_color_manual(values = c(smc, pbo, itn, rtss_sv, rtss_age, pbo_smc, itn_smc, rtss_smc, pbo_rtss_smc, itn_rtss_smc, pbo_rtss, itn_rtss)) +
    scale_x_continuous(breaks=c(0)) +
    theme_classic() +
    theme(plot.caption.position = "plot")

  ggsave(paste0('./03_output/plots_draws/box_whisker_CE_', cost_dose, '.pdf'), width=10, height=5)

}

table(scenarios$cost_per_dose)

box_plot(2.69)
box_plot(6.52)
box_plot(12.91)


# plotting with cost per cases averted
tabledat <- scenarios %>%
  group_by(intervention_f) %>%
  summarize(n=n(),
            median = round(median(CE_case)),
            q25 = round(quantile(CE_case, p=0.25)),
            q75 = round(quantile(CE_case, p=0.75)),
            min = round(min(CE_case)),
            max = round(max(CE_case))) %>%
  arrange(median)

tabledat; sum(tabledat$n)

scenarios %>% filter(intervention != 'none') %>%
  filter(cost_per_dose == 6.52 & delivery_cost == 1.62) %>%
  mutate(intervention_f = factor(intervention, levels=levels$intervention_f)) %>%
  mutate(rank=as.numeric(intervention_f)) %>%

  ggplot(aes(x=rank, y=CE_case, fill=intervention_f, color=intervention_f, group=intervention)) +
  geom_hline(yintercept = 0, lty=2, color='grey') +
  geom_vline(xintercept = 5.5, lty=2, color='grey') +
  geom_boxplot(alpha=0.3, outlier.alpha = 0.05,  outlier.size = 0.05) +
  coord_cartesian(ylim=c(-1, 80), clip="off") +
  labs(x='',
       y=expression(paste(Delta," cost / ", Delta, " cases")),
       fill = 'intervention',
       color = 'intervention') +
  annotation_custom(textGrob("Univariate strategies"),xmin=1,xmax=5,ymin=-7,ymax=-7) +
  annotation_custom(textGrob("Mixed strategies"),xmin=6,xmax=12,ymin=-7,ymax=-7) +
  scale_fill_manual(values = c(smc, pbo, itn, rtss_sv, rtss_age, pbo_smc, itn_smc, rtss_smc, pbo_rtss_smc, itn_rtss_smc, pbo_rtss, itn_rtss)) +
  scale_color_manual(values = c(smc, pbo, itn, rtss_sv, rtss_age, pbo_smc, itn_smc, rtss_smc, pbo_rtss_smc, itn_rtss_smc, pbo_rtss, itn_rtss)) +
  scale_x_continuous(breaks=c(0)) +
  theme_classic()

ggsave('./03_output/plots_draws/box_whisker_CE_cases.pdf', width=10, height=5)



# per dose RTS,S cost ----------------------------------------------------------
# function to find the most common character value in a group
calculate_mode <- function(x) {
  uniqx <- unique(na.omit(x))
  uniqx[which.max(tabulate(match(x, uniqx)))]
}

output <- scenarios %>%
  filter(cost_per_dose == 6.52 & delivery_cost == 1.62) %>%
  # choose univariate scenarios
  filter(intervention %in% c('RTS,S seasonal', 'RTS,S age-based', 'SMC', 'none', 'ITN 10% increase', 'ITN PBO')) %>%
  # combine RTS,S strategies
  mutate(intervention = ifelse(intervention=='RTS,S age-based' | intervention=='RTS,S seasonal', 'RTS,S', intervention)) %>%
  arrange(ID, drawID) %>%
  group_by(ID, drawID) %>%
  # find minimum CE in each group. If intervention == RTS,S find the second or third lowest CE and move RTS,S to the top
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
  arrange(ID, costRTSS)


# < all seasons box and whisker ####
seasoncosts <- output %>% group_by(seasonality) %>%
  summarize(q25 = round(quantile(costRTSS, probs = 0.25, na.rm = T), 2),
            med = round(median(costRTSS, na.rm = T),2),
            q75 = round(quantile(costRTSS, probs = 0.75, na.rm = T), 2),
            min = round(min(costRTSS, na.rm = T), 2),
            max = round(max(costRTSS, na.rm = T), 2))

A <- ggplot(output) +
  geom_hline(yintercept = 0, color = 'light grey') +
  geom_boxplot(aes(x = factor(seasonality, levels=c('perennial','seasonal','highly seasonal')), y = costRTSS),
               fill = 'cornflower blue', color = 'cornflower blue', alpha = 0.4, outlier.alpha = 0.1,  outlier.size = 0.1) +
  geom_text(data = seasoncosts, aes(x=seasonality, y=q25, label=q25), size = 3, nudge_y = .55, nudge_x = -.22) +
  geom_text(data = seasoncosts, aes(x=seasonality, y=med, label=med), size = 3, nudge_y = .55, nudge_x = -.2) +
  geom_text(data = seasoncosts %>% filter(seasonality != 'seasonal'), aes(x=seasonality, y=q75, label=q75), size = 3, nudge_y = .55, nudge_x = -.2) +
  geom_text(data = seasoncosts %>% filter(seasonality == 'seasonal'), aes(x=seasonality, y=q75, label=q75), size = 3, nudge_y = 1.2, nudge_x = -.2) +
  geom_vline(xintercept = 0, lty = 2, color = 'grey') +
  theme_classic() +
  labs(y='RTS,S cost (USD) per dose',
       #caption=paste0('range = ', round(min(output$costRTSS)), ' to ', round(max(output$costRTSS)))
       x=''
  ) +
  coord_cartesian(ylim = c(-15,15)) +
  theme(plot.caption.position = "plot")

ggsave(A, './03_output/plots_draws/RTSS_price_dist_all.pdf', width=6, height=4)


# table
seasoncosts

output %>% ungroup() %>%
  group_by(interventionmin) %>%
  summarize(q25 = round(quantile(costRTSS, probs = 0.25, na.rm=T),2),
            med = round(median(costRTSS, na.rm=T),2),
            q75 = round(quantile(costRTSS, probs = 0.75, na.rm=T),2))

# how many are above $5 - SMC / ITNuse
output %>% filter(costRTSS >= 5) %>% group_by(SMC, ITNuse) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  mutate(t = sum(n), p = n/t * 100)


# how many are above $5 - pfpr
output %>% filter(costRTSS >= 5) %>% group_by(pfpr) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  mutate(t = sum(n), p = n/t * 100)


# < line plot ####

output2 <- output %>% ungroup() %>% arrange(costRTSS) %>%
  mutate(p = (nrow(output) - row_number() + 1) / nrow(output) * 100) %>%
  ungroup()

dose2 <- output2 %>% filter(costRTSS <= 2) %>% top_n(-1) %>% select(p) %>% as.numeric()
dose5 <- output2 %>% filter(costRTSS <= 5) %>% top_n(-1) %>% select(p) %>% as.numeric()
dose10 <- output2 %>% filter(costRTSS <= 10) %>% top_n(-1) %>% select(p) %>% as.numeric()
dose135 <-  output2 %>% filter(costRTSS <= 13.5) %>% top_n(-1) %>% select(p) %>% as.numeric()

table(output$drawID) # 450

output3 <- output %>% ungroup() %>% group_by(drawID) %>% arrange(costRTSS) %>%
  mutate(p = (450 - row_number() + 1) / 450 * 100) %>%
  ungroup()

segments <- data.frame(
  x = c(2, -20, 5, -20, 10, -20, 13.5, -20),
  xend = c(2, 2, 5, 5, 10, 10, 13.5, 13.5),
  y = c(-10, dose2, -10, dose5, -10, dose10, -10, dose135),
  yend = c(dose2, dose2, dose5, dose5, dose10, dose10, dose135, dose135))

points <- data.frame(
  x = c(2, 5, 10, 13.5),
  y = c(dose2, dose5, dose10, dose135)
)

B <- ggplot(output3) +
  geom_vline(xintercept = 0, color = 'grey') +
  geom_segment(data = segments,
               aes(x = x, xend = xend, y = y, yend = yend),
               lty = 3, color = 'grey') +
  geom_line(aes(x = costRTSS, y = p, group = drawID), color = '#619CFF', size = 1, alpha = 0.1) +
  geom_line(data = output2, aes(x = costRTSS, y = p, group = drawID), color = '#619CFF', size = 1) +
  geom_point(data = points, aes(x = x, y = y), color = 'blue', size = 2) +
  theme_classic() +
  labs(y='% of scenarios where \nRTS,S is most cost-effective',
       x='RTS,S cost per dose (USD)'
  ) +
  scale_x_continuous(breaks = c(-5, 0, 2, 5, 10, 13.5)) +
  coord_cartesian(xlim = c(0, 15), ylim = c(0, 100))

ggsave('./03_output/plots_draws/RTSS_price_dist_lineplot.pdf', B, width=6, height=4)

# combined plot
A + B + plot_annotation(tag_levels = 'A')

ggsave('./03_output/plots_draws/RTSS_price_dist_AB.pdf', width=10, height=4)


# univariate stacked histogram patterns ----------------------------------------
# look at color choices

# < by seasonality alone ------------------------------
univariateseason <- function(season) {

  # assign colors
  if(season=='perennial'){
    colors <- c(itn, pbo, rtss_age)
  }

  if(season=='highly seasonal'){
    colors <- c(itn, pbo, rtss_age, rtss_sv)
  }

  if(season=='seasonal'){
    colors <- c(itn, pbo, smc)
  }

  # by pfpr
  output <- scenarios %>%
    filter(seasonality==season) %>%
    filter(cost_per_dose==6.52 & delivery_cost==1.62) %>%
    filter(intervention %in% c('ITN 10% increase','ITN PBO','SMC','RTS,S age-based','RTS,S seasonal')) %>%
    group_by(ID, drawID) %>% arrange(ID, CE) %>%
    slice(1L)

  A <- ggplot(output) +
    geom_bar(aes(x=as.factor(pfpr), fill=intervention_f), position="fill", show.legend = F) +
    labs(x='PfPR', y='Proportion most \ncost-effective choice', fill='intervention') +
    scale_fill_manual(values=colors) +
    theme_classic()

  # by current ITN usage
  B <- ggplot(output) +
    geom_bar(aes(x=factor(ITNuse), fill=intervention_f), position="fill", show.legend=F) +
    labs(x='ITN use', y='Proportion most \ncost-effective choice', fill='intervention') +
    scale_fill_manual(values=colors) +
    theme_classic()

  # by insecticide resistance
  C <- ggplot(output) +
    geom_bar(aes(x=factor(resistance), fill=intervention_f), position="fill", show.legend=F) +
    labs(x='Resistance', y='Proportion most \ncost-effective choice', fill='intervention') +
    scale_fill_manual(values=colors) +
    theme_classic()

  # by ITN distribution efficiency
  ITNefficient <- function(var, label) {
    scenarios %>%
      filter(seasonality==season) %>%
      filter(cost_per_dose==6.52 & delivery_cost==1.62) %>%
      filter(intervention %in% c('ITN 10% increase','ITN PBO','SMC','RTS,S age-based','RTS,S seasonal')) %>%
      group_by(ID, drawID) %>% arrange(ID, drawID, {{var}}) %>%
      slice(1L) %>% select(intervention_f, seasonality, pfpr, {{var}}) %>%
      mutate(model=label) %>%
      rename(CE = {{var}})
  }

  output2 <- ITNefficient(CE, 'standard') %>%
    full_join(ITNefficient(CE_ITNmin, 'less efficient')) %>%
    full_join(ITNefficient(CE_ITNmax, 'more efficient'))


  if(season=='seasonal'){
    colors <- c(itn, pbo, rtss_age, rtss_sv, smc)
  }

  D <- ggplot(output2) +
    geom_bar(aes(x=factor(model, levels=c('less efficient', 'standard', 'more efficient')), fill=intervention_f), position="fill") +
    labs(x='ITN efficiency', y='Proportion most \ncost-effective choice', fill='intervention') +
    scale_fill_manual(values=colors) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
    theme_classic()


  if(season=='seasonal'){
    colors <- c(itn, pbo, smc)
  }

  E <- ggplot(output) +
    geom_bar(aes(x=factor(treatment, labels=c('low', 'medium', 'high')), fill=intervention_f), position="fill", show.legend = F) +
    labs(x='Treatment coverage', y='Proportion most \ncost-effective choice', fill='intervention') +
    scale_fill_manual(values=colors) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
    theme_classic()

  (A + B + C + D + E) + plot_layout(guides = "collect", nrow=3, ncol=2) + plot_annotation(tag_levels = 'A')

  ggsave(paste0('./03_output/plots_draws/univariate_quad_', season, '.pdf'), width=9, height=7)

}


univariateseason('highly seasonal')
univariateseason('seasonal')
univariateseason('perennial')


# < faceted by season ----------------------------------------------------------
colors <- c(itn, pbo, rtss_age, rtss_sv, smc)

output <- scenarios %>%
  filter(cost_per_dose==6.52 & delivery_cost==1.62) %>%
  filter(intervention %in% c('ITN 10% increase','ITN PBO','SMC','RTS,S age-based','RTS,S seasonal')) %>%
  mutate(seasonality = factor(seasonality, levels = c('perennial', 'seasonal', 'highly seasonal'))) %>%
  group_by(ID, drawID) %>% arrange(ID, drawID, CE) %>%
  slice(1L)


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
    group_by(ID, drawID) %>% arrange(ID, drawID, {{var}}) %>%
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
  geom_bar(aes(x=factor(treatment, labels=c('low', 'medium', 'high')), fill=intervention_f), position="fill", show.legend=F) +
  labs(x='Treatment coverage', y='Proportion most \ncost-effective choice', fill='intervention') +
  scale_fill_manual(values=colors) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  facet_grid(~seasonality) +
  theme_classic()

legend <- cowplot::get_legend(D)

Dm <- D + theme(legend.position = "none")


(A + B + C + Dm + E + legend) +
  plot_layout(nrow = 3) +
  plot_annotation(tag_levels = list(c('A', 'B', 'C', 'D', 'E', '')))

ggsave(paste0('./03_output/plots_draws/univariate_quad_by_season.pdf'), width=10, height=8)


prop_CE <- function(data, var){
  data %>%
    mutate(intervention = ifelse(grepl('RTS,S', intervention), 'RTS,S', intervention),
           intervention_f = factor(intervention, levels = c('ITN 10% increase', 'ITN PBO', 'RTS,S', 'SMC'))) %>%
    group_by(seasonality, {{var}}, intervention_f) %>%
    summarize(n = n()) %>%
    ungroup() %>%
    group_by(seasonality, {{var}}) %>%
    mutate(t = sum(n), p = n / t * 100) %>%
    filter(intervention_f == 'RTS,S')
}

print(prop_CE(output, sim_length), n = 30) # overall
print(prop_CE(output, pfpr), n = 30) # overall
print(prop_CE(output, ITNuse), n = 40)
print(prop_CE(output, resistance), n = 30)
print(prop_CE(output2, model), n = 40)
print(prop_CE(output, treatment), n = 30)

output %>%
  group_by(seasonality, intervention) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  group_by(seasonality) %>%
  mutate(t = sum(n), p = n / t * 100)


# TABLES -----------------------------------------------------------------------
# < CE table RTS,S doses ----
# cost_per_dose <- c(2.69, 6.52, 12.91)
# delivery_cost <- c(0.96, 1.62, 2.67)
# cost_per_dose + delivery_cost

cost_per_dose2 <- c(2.69+1.62, 6.52+1.62, 12.91+1.62) %>% as_tibble %>% rename(cost_per_dose2=value)

# pull out univariate scenarios
scenarios %>%
  filter(intervention %in% c('RTS,S seasonal', 'RTS,S age-based')) %>%
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
       filter(intervention %in% c('RTS,S seasonal', 'RTS,S age-based')) %>%
       filter(resistance == 0 & seasonality == 'perennial'))


# < impact RTSS on top of other interventions ----------------------------------
output <- scenarios %>%
  filter(cost_per_dose==6.52 & delivery_cost==1.62) %>%
  mutate(ID = paste(pfpr, seasonality, ITNuse, resistance, treatment, ITN, sep="_")) %>%
  filter(ITNuse==0.75)

none <- output %>%
  mutate(set = case_when(seasonality %in% c('perennial', 'highly seasonal') & intervention %in%
                           c('ITN 10% increase','ITN PBO') ~ 1,
                         seasonality=='seasonal' & intervention %in%
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
  dplyr::select(file, ID, drawID, pfpr, seasonality, intervention, daly, cases, cost_total, u5_dalys) %>%
  left_join(none %>% dplyr::select(-file), by=c('ID', 'drawID')) %>%
  mutate(CE = (cost_total - cost_total_baseline) / (daly_baseline - daly),
         deltadaly = daly_baseline - daly,
         deltacases = cases_baseline - cases,
         CE_u5 = (cost_total - cost_total_baseline) / (u5_daly_baseline - u5_dalys))


# adjusting for 15 year simulation period and 200,000 population arguments
summary(output2$deltadaly / (2*15)) # additional dalys averted per year in a population of 100,000 people
summary(output2$deltacases / (2*15)) # additional cases averted per year in a population of 100,000 people
summary(output2$CE_u5) # additional cases averted per year in a population of 100,000 people

ggplot(data = output2) +
  geom_boxplot(aes(x=pfpr, y=deltadaly / (2*15), alpha = 0.3, group = pfpr), show.legend = F, fill = '#619CFF', color = '#619CFF') +
  facet_grid(~factor(seasonality, levels = c('perennial', 'seasonal', 'highly seasonal'))) +
  labs(x = 'PfPR', y = 'DALYs averted', title = 'DALYS averted per year per 100,000 people') +
  theme_classic()

ggsave('./03_output/RTSS_additional_impact.pdf', width=6, height=3)



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




# CASE STUDY -------------------------------------------------------------------

output <- readRDS('./03_output/dalyoutput_draws_casestudy.rds') %>%
  filter(cost_per_dose==6.52 & delivery_cost==1.62)

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

  mutate(residence = ifelse(pfpr %in% c(0.20, 0.40), 'rural', 'urban')) %>%

  dplyr::select(drawID, seasonality, scenario, scenario_f, scenario2, scenario2_f, residence, cost_total, daly, cases, deaths) %>%

  arrange(seasonality, scenario2, scenario2_f, scenario, scenario_f, drawID, residence) %>%
  group_by(seasonality, scenario2, scenario2_f, scenario, scenario_f, drawID) %>%

  mutate(disparity = abs(daly - lag(daly))) %>%

  summarize(across(c(cost_total, cases, deaths, daly, disparity), sum, na.rm = T))


none <- scenarios %>% ungroup() %>%
  filter(scenario == 1) %>%
  mutate(cost_total_baseline = cost_total,
         cases_baseline = cases,
         deaths_baseline = deaths,
         daly_baseline = daly,
         disparity_baseline = disparity) %>%
  dplyr::select(seasonality, scenario2, scenario2_f, drawID, cost_total_baseline:disparity_baseline)


scenarios2 <- scenarios %>%
  left_join(none, by = c("seasonality", "scenario2", "scenario2_f", "drawID")) %>%
  mutate(cost_diff = cost_total - cost_total_baseline,
         cases_diff = cases_baseline - cases,
         deaths_diff = deaths_baseline - deaths,
         daly_diff = daly_baseline - daly,
         disparity_per = (disparity - disparity_baseline) / disparity_baseline * 100,
         CE = cost_diff / daly_diff) %>%
  ungroup() %>%
  # summarize over drawID
  group_by(seasonality, scenario2, scenario2_f, scenario, scenario_f) %>%
  summarize(n = 50,
            estimate_CE = mean(CE, na.rm = T),
            lower_CE = estimate_CE - 1.96 * (sd(CE, na.rm = T) / sqrt(n)),
            upper_CE = estimate_CE + 1.96 * (sd(CE, na.rm = T) / sqrt(n)),

            estimate_disparity = mean(disparity_per, na.rm = T),
            lower_dis = estimate_disparity - 1.96 * (sd(disparity_per, na.rm = T) / sqrt(n)),
            upper_dis = estimate_disparity + 1.96 * (sd(disparity_per, na.rm = T) / sqrt(n))) %>%
  # put negative cost effectiveness values up to 0
  mutate(estimate_CE = case_when(estimate_CE < 0 ~ 0,
                                 TRUE ~ estimate_CE),
         lower_CE = case_when(estimate_CE == 0 ~ 0,
                                 TRUE ~ lower_CE),
         upper_CE = case_when(estimate_CE == 0 ~ 0,
                                 TRUE ~ upper_CE))


# < primer ----
library(ggforce)

a <- data.frame(
  x = c(0.5, 0.5, 1.5, 1.5),
  y = c(1.5, 2.5, 2.5, 1.5)
)

b <- data.frame(
  x = c(1.5, 1.5, 2.5, 2.5),
  y = c(1.5, 2.5, 2.5, 1.5)
)

c <- data.frame(
  x = c(0.5, 0.5, 1.5, 1.5),
  y = c(0.5, 1.5, 1.5, 0.5)
)

d <- data.frame(
  x = c(1.5, 1.5, 2.5, 2.5),
  y = c(0.5, 1.5, 1.5, 0.5)
)


A <- ggplot() +
  geom_shape(data = a, aes(x = x, y = y), expand = unit(-.2, 'cm'), radius = unit(0.5, 'cm'), fill = 'grey', alpha = 0.9) +
  geom_shape(data = b, aes(x = x, y = y), expand = unit(-.2, 'cm'), radius = unit(0.5, 'cm'), fill = 'tomato1', alpha = .95) +
  geom_shape(data = c, aes(x = x, y = y), expand = unit(-.2, 'cm'), radius = unit(0.5, 'cm'), fill = 'chartreuse4', alpha = 0.8) +
  geom_shape(data = d, aes(x = x, y = y), expand = unit(-.2, 'cm'), radius = unit(0.5, 'cm'), fill = 'grey', alpha = 0.9) +
  annotate("text", x = 1, y = 2, label = 'CE \u274c \nEquity \u2714', color = 'white') +
  annotate("text", x = 2, y = 2, label = 'CE \u274c \nEquity \u274c', color = 'white') +
  annotate("text", x = 1, y = 1, label = 'CE \u2714 \nEquity \u2714', color = 'white') +
  annotate("text", x = 2, y = 1, label = 'CE \u2714 \nEquity \u274c', color = 'white') +
  labs(y = 'cost-effectiveness', x = 'disparity between urban and rural') +
  theme_classic() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank())

# adding asterisks to negative CE values
stars <- scenarios2 %>% filter(estimate_CE == 0) %>%
  select(estimate_CE, estimate_disparity) %>%
  mutate(x = estimate_disparity + 2, y = 4)

B <- ggplot(data = scenarios2 %>% filter(scenario > 1)) + # remove baseline
  geom_vline(xintercept = 0, lty=2, color='grey') +
  geom_pointrange(aes(x = estimate_disparity, y = estimate_CE, ymin = lower_CE, ymax = upper_CE, shape = scenario2_f, color = scenario_f)) +
  geom_pointrange(aes(x = estimate_disparity, xmin = lower_dis, xmax = upper_dis, y = estimate_CE, shape = scenario2_f, color = scenario_f)) +
  geom_point(data = stars, aes(x = x, y = y), shape = "*", size = 4, color = "black") +
  labs(y = expression(paste('cost-effectiveness (', Delta," cost / ", Delta, " DALYs)")),
       x = 'disparity between urban and rural (% change in DALYs)',
       shape = 'baseline scenario',
       color = 'intervention') +
  scale_color_manual(values = c("#C70E7B","#007BC3", "#FC6882","#54BCD1")) +
  scale_x_continuous(labels = scales::percent_format(scale = 1)) +
  coord_cartesian(xlim = c(-40, 40), ylim = c(0, 200), clip = 'on') +
  theme_classic()

(A + B) + plot_layout(guides = "collect", nrow=1) + plot_annotation(tag_levels = 'A')

ggsave(paste0('./03_output/plots_draws/case_study.png'), width=10, height=4)
