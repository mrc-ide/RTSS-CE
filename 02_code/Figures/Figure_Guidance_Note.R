# Figure to support Technical Guidance Note
# Cost-effectiveness ratios (cost per DALY averted) of univariate and mixed strategies
# Proportion of times an intervention was the most cost-efficient out of all scenarios

# set-up
source("./02_code/Figures/data_and_libraries.R")

# A ----
# set up order of interventions by median CE for plot
levels <- scenarios %>%
  mutate(group = case_when(intervention_f %in% c('ITN 10% increase', 'ITN PBO', 'SMC', 'RTS,S seasonal', 'RTS,S age-based') ~ 'uni',
                           TRUE ~ 'multi')) %>%
  filter(group == 'uni') %>%
  group_by(intervention_f, group) %>%
  summarize(med = median(CE, na.rm = T) )%>% ungroup() %>%
  arrange(desc(group), med)


# list color order for plots
colors <- c(smc, pbo, itn, rtss_sv, rtss_age)

# plot of cost per DALY averted
A <- scenarios %>%
  filter(intervention_f %in% c('ITN 10% increase', 'ITN PBO', 'SMC', 'RTS,S seasonal', 'RTS,S age-based')) %>%
  filter(cost_per_dose == 12.01 & delivery_cost == 1.62) %>%
  mutate(intervention_f = factor(intervention, levels = levels$intervention_f)) %>%
  mutate(rank = as.numeric(intervention_f)) %>%

  ggplot(aes(x = rank, y = CE, fill = intervention_f,
             color = intervention_f, group = intervention)) +
  geom_hline(yintercept = 0, lty = 2, color = 'grey') +
  geom_boxplot(alpha = 0.3, outlier.alpha = 0.05,  outlier.size = 0.05, coef = 500) +
  coord_cartesian(ylim = c(-50, 500), clip = "off") +
  labs(x = '',
       y = expression(paste(Delta," cost / ", Delta, " DALYs")),
       fill = 'intervention',
       color = 'intervention') +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  scale_x_continuous(breaks = c(0)) +
  scale_y_continuous(breaks = seq(0, 500, 50)) +
  theme_classic() +
  theme(plot.caption.position = "plot",
        text = element_text(size = 14))


# B and C ----
output <- scenarios %>%
  filter(cost_per_dose == 12.01 & delivery_cost == 1.62) %>%
  filter(intervention %in% c('ITN 10% increase','ITN PBO','SMC','RTS,S age-based','RTS,S seasonal')) %>%
  mutate(seasonality = factor(seasonality, levels = c('perennial', 'seasonal', 'highly seasonal'))) %>%
  group_by(ID, drawID) %>% arrange(ID, drawID, CE) %>%
  slice(1L)

colors <- c(itn, pbo, rtss_age, rtss_sv, smc)


# by pfpr
B <- ggplot(output) +
  geom_bar(aes(x = as.factor(pfpr), fill = intervention_f), position = "fill") +
  labs(x = 'PfPR', y = 'Proportion most \ncost-effective choice', fill = 'intervention') +
  scale_fill_manual(values = colors) +
  facet_grid(~ seasonality) +
  theme_classic() +
  theme(text = element_text(size = 12))

# by current ITN usage
C <- ggplot(output) +
  geom_bar(aes(x = factor(ITNuse), fill = intervention_f), position = "fill") +
  labs(x = 'ITN use', y = 'Proportion most \ncost-effective choice', fill = 'intervention') +
  scale_fill_manual(values = colors) +
  facet_grid(~ seasonality) +
  theme_classic() +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(size = 9))


# final figure
A + (B / C) +
  plot_layout(guides = "collect", nrow = 1, widths = c(1, 1.3)) +
  plot_annotation(tag_levels = list(c('A')))

ggsave(paste0('./03_output/figure_guidancenote.pdf'), width = 10, height = 4.5)
