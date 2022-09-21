# Figure 3
# Change in cost per change in DALYs averted

# set-up
source("./02_code/Figures/data_and_libraries.R")


# < plot by season  -----------------

deltaseason <- function(season){

  output <- scenarios %>%
    filter(seasonality == season) %>%
    filter(cost_per_dose == 12.01 & delivery_cost == 1.62) %>%
    filter(resistance == 0) %>%
    filter(intervention != 'none') %>%

    mutate(deltadaly = daly_baseline - daly,
           deltacost = cost_total - cost_total_baseline,
           cost_daly_averted = deltacost / deltadaly) %>%

    group_by(ID, seasonality, intervention, intervention_f) %>%
    summarize(deltadaly = median(deltadaly),
              deltacost = median(deltacost),
              cost_daly_averted = median(cost_daly_averted))

  # all interventions, by baseline ITNuse and seasonality
  output <- output %>% mutate(intervention2 = factor(intervention, levels = c('ITN 10% boost','ITN PBO','RTS,S EPI','RTS,S SV','ITN 10% boost + RTS,S','ITN PBO + RTS,S','SMC','ITN 10% boost + SMC','ITN PBO + SMC', 'RTS,S + SMC', 'ITN 10% boost + RTS,S + SMC','ITN PBO + RTS,S + SMC')))

  # assign colors
  if(season == 'perennial'){
    colors <- c(itn, rtss_age, itn_rtss)
  }

  if(season == 'highly seasonal'){
    colors <- c(itn, rtss_age, rtss_sv, itn_rtss)
  }

  if(season == 'seasonal'){
    colors <- c(itn, rtss_age, rtss_sv, smc, itn_rtss, itn_smc, rtss_smc, itn_rtss_smc)
  }

  A <- ggplot(data = output, mapping = aes(x = deltadaly, y = deltacost)) +
    geom_line(aes(group = as.factor(ID)), color ='lightgrey', size = .5) +
    geom_point(aes(color = intervention_f), size = 2) +
    geom_hline(yintercept = 0, lty = 2, color = "black") +
    geom_vline(xintercept = 0, lty = 2, color = "black") +
    theme_classic() +
    scale_color_manual(values = colors) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title='All strategies',
         y='change in cost (USD)',
         x='change in DALYs averted',
         color='intervention')

  if(season=='perennial'){
    colors <- c(itn, rtss_age, itn_rtss)
  }

  if(season=='highly seasonal'){
    colors <- c(itn, rtss_sv, itn_rtss)
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
    filter(dominate == 0) %>%
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
    filter(dominate == 0) %>%
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
    filter(dominate == 0) %>%
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
    filter(dominate == 0) %>%
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
    filter(dominate == 0) %>%
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
    filter(dominate == 0) %>%

    ggplot(mapping = aes(x = deltadaly, y = deltacost)) +
    geom_line(aes(group = as.factor(ID)), color = 'lightgrey', size = .5) +
    geom_point(aes(color = intervention_f), size=2, show.legend = F) +
    geom_hline(yintercept = 0, lty = 2, color = "black") +
    geom_vline(xintercept = 0, lty = 2, color = "black") +
    theme_classic() +
    scale_color_manual(values = colors) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = 'Dominated strategies removed',
         y = 'change in cost (USD)',
         x = 'change in DALYs averted',
         color = 'intervention',
         caption = 'Assuming resistance == 0') +
    theme(plot.caption.position = "plot")

  A + B + plot_layout(guides = "collect", nrow = 1) +
    plot_annotation(tag_levels = 'A')

  ggsave(paste0('./03_output/figure3_', season, '.pdf'),
         width = 12, height = 5)

}


deltaseason('highly seasonal')
deltaseason('seasonal')
deltaseason('perennial')


# < facet by season -----------------

# calculate median change in DALYs and median change in cost
output <- scenarios %>%
  filter(cost_per_dose == 12.01 & delivery_cost == 1.62) %>%
  filter(resistance == 0) %>%
  filter(intervention != 'none') %>%

  mutate(deltadaly = daly_baseline - daly,
         deltacost = cost_total - cost_total_baseline,
         cost_daly_averted = deltacost / deltadaly,
         seasonality = factor(seasonality, levels = c('perennial', 'seasonal', 'highly seasonal'))) %>%

  group_by(ID, seasonality, intervention, intervention_f) %>%
  summarize(deltadaly = median(deltadaly),
            deltacost = median(deltacost),
            cost_daly_averted = median(cost_daly_averted))

# list color order for plots
colors <- c(itn, rtss_age, rtss_sv, smc, itn_rtss, itn_smc, rtss_smc, itn_rtss_smc)

# plot all interventions, by baseline ITNuse and seasonality
A <- ggplot(data = output, aes(x = deltadaly, y = deltacost)) +
  geom_line(aes(group = as.factor(ID)), color = 'lightgrey',
            size = .5, alpha = 0.2) +
  geom_point(aes(color = intervention_f), size = 2, alpha = 0.6) +
  geom_hline(yintercept = 0, lty = 2, color = "black") +
  geom_vline(xintercept = 0, lty = 2, color = "black") +
  facet_grid(~ seasonality, scales = "free") +
  theme_classic() +
  scale_color_manual(values = colors) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 14)) +
  labs(title = 'All strategies',
       y = 'change in cost (USD)',
       x = 'change in DALYs averted',
       color = 'intervention')

# removing dominated strategies
B <- output %>%
  # filter out mixed strategies
  group_by(ID) %>% arrange(ID, deltacost) %>%
  filter(!(deltadaly < 0 & deltacost > 0)) %>%
  # filter out dominated strategies
  mutate(dominate = case_when(deltadaly < lag(deltadaly, n=12L) ~ 1,
                              deltadaly < lag(deltadaly, n=11L) ~ 1,
                              deltadaly < lag(deltadaly, n=10L) ~ 1,
                              deltadaly < lag(deltadaly, n=9L) ~ 1,
                              deltadaly < lag(deltadaly, n=8L) ~ 1,
                              deltadaly < lag(deltadaly, n=7L) ~ 1,
                              deltadaly < lag(deltadaly, n=6L) ~ 1,
                              deltadaly < lag(deltadaly, n=5L) ~ 1,
                              deltadaly < lag(deltadaly, n=4L) ~ 1,
                              deltadaly < lag(deltadaly, n=3L) ~ 1,
                              deltadaly < lag(deltadaly, n=2L) ~ 1,
                              deltadaly < lag(deltadaly, n=1L) ~ 1,
                              TRUE ~ 0)) %>%
  filter(dominate == 0) %>%
  # marking extended dominated strategies
  mutate(ICER = (deltacost - lag(deltacost)) / (deltadaly - lag(deltadaly)),
         dominate = case_when(ICER > lead(ICER, n=12L) ~ 1,
                              ICER > lead(ICER, n=11L) ~ 1,
                              ICER > lead(ICER, n=10L) ~ 1,
                              ICER > lead(ICER, n=9L) ~ 1,
                              ICER > lead(ICER, n=8L) ~ 1,
                              ICER > lead(ICER, n=7L) ~ 1,
                              ICER > lead(ICER, n=6L) ~ 1,
                              ICER > lead(ICER, n=5L) ~ 1,
                              ICER > lead(ICER, n=4L) ~ 1,
                              ICER > lead(ICER, n=3L) ~ 1,
                              ICER > lead(ICER, n=2L) ~ 1,
                              ICER > lead(ICER, n=1L) ~ 1,
                              TRUE ~ 0)) %>%
  filter(dominate == 0) %>%
  # loop through extendedly dominated strategies again
  mutate(ICER = (deltacost - lag(deltacost)) / (deltadaly - lag(deltadaly)),
         dominate = case_when(ICER > lead(ICER, n=12L) ~ 1,
                              ICER > lead(ICER, n=11L) ~ 1,
                              ICER > lead(ICER, n=10L) ~ 1,
                              ICER > lead(ICER, n=9L) ~ 1,
                              ICER > lead(ICER, n=8L) ~ 1,
                              ICER > lead(ICER, n=7L) ~ 1,
                              ICER > lead(ICER, n=6L) ~ 1,
                              ICER > lead(ICER, n=5L) ~ 1,
                              ICER > lead(ICER, n=4L) ~ 1,
                              ICER > lead(ICER, n=3L) ~ 1,
                              ICER > lead(ICER, n=2L) ~ 1,
                              ICER > lead(ICER, n=1L) ~ 1,
                              TRUE ~ 0)) %>%
  filter(dominate == 0) %>%
  mutate(ICER = (deltacost - lag(deltacost)) / (deltadaly - lag(deltadaly)),
         dominate = case_when(ICER > lead(ICER, n=12L) ~ 1,
                              ICER > lead(ICER, n=11L) ~ 1,
                              ICER > lead(ICER, n=10L) ~ 1,
                              ICER > lead(ICER, n=9L) ~ 1,
                              ICER > lead(ICER, n=8L) ~ 1,
                              ICER > lead(ICER, n=7L) ~ 1,
                              ICER > lead(ICER, n=6L) ~ 1,
                              ICER > lead(ICER, n=5L) ~ 1,
                              ICER > lead(ICER, n=4L) ~ 1,
                              ICER > lead(ICER, n=3L) ~ 1,
                              ICER > lead(ICER, n=2L) ~ 1,
                              ICER > lead(ICER, n=1L) ~ 1,
                              TRUE ~ 0)) %>%
  filter(dominate == 0) %>%
  mutate(ICER = (deltacost - lag(deltacost)) / (deltadaly - lag(deltadaly)),
         dominate = case_when(ICER > lead(ICER, n=12L) ~ 1,
                              ICER > lead(ICER, n=11L) ~ 1,
                              ICER > lead(ICER, n=10L) ~ 1,
                              ICER > lead(ICER, n=9L) ~ 1,
                              ICER > lead(ICER, n=8L) ~ 1,
                              ICER > lead(ICER, n=7L) ~ 1,
                              ICER > lead(ICER, n=6L) ~ 1,
                              ICER > lead(ICER, n=5L) ~ 1,
                              ICER > lead(ICER, n=4L) ~ 1,
                              ICER > lead(ICER, n=3L) ~ 1,
                              ICER > lead(ICER, n=2L) ~ 1,
                              ICER > lead(ICER, n=1L) ~ 1,
                              TRUE ~ 0)) %>%
  filter(dominate == 0) %>%
  mutate(ICER = (deltacost - lag(deltacost)) / (deltadaly - lag(deltadaly)),
         dominate = case_when(ICER > lead(ICER, n=12L) ~ 1,
                              ICER > lead(ICER, n=11L) ~ 1,
                              ICER > lead(ICER, n=10L) ~ 1,
                              ICER > lead(ICER, n=9L) ~ 1,
                              ICER > lead(ICER, n=8L) ~ 1,
                              ICER > lead(ICER, n=7L) ~ 1,
                              ICER > lead(ICER, n=6L) ~ 1,
                              ICER > lead(ICER, n=5L) ~ 1,
                              ICER > lead(ICER, n=4L) ~ 1,
                              ICER > lead(ICER, n=3L) ~ 1,
                              ICER > lead(ICER, n=2L) ~ 1,
                              ICER > lead(ICER, n=1L) ~ 1,
                              TRUE ~ 0)) %>%
  filter(dominate == 0) %>%

  ggplot(mapping = aes(x = deltadaly, y = deltacost)) +
  geom_line(aes(group = as.factor(ID)), color = 'lightgrey',
            size = .5, alpha = 0.2) +
  geom_point(aes(color = intervention_f), size = 2, alpha = 0.6, show.legend = F) +
  geom_hline(yintercept = 0, lty = 2, color = "black") +
  geom_vline(xintercept = 0, lty = 2, color = "black") +
  facet_grid(~ seasonality, scales = "free") +
  theme_classic() +
  scale_color_manual(values = colors) +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        text = element_text(size = 14)) +
  labs(title = 'Dominated strategies removed',
       y = 'change in cost (USD)',
       x = 'change in DALYs averted',
       color = 'intervention') +
  theme(plot.caption.position = "plot")

A + B + plot_layout(guides = "collect", nrow=2) + plot_annotation(tag_levels = 'A')

ggsave(paste0('./03_output/figure3_all.pdf'), width = 12, height = 7)

