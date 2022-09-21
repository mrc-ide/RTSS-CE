# Figure 2
# Cost-effectiveness ratios (cost per DALY averted) of univariate and mixed strategies

# set-up
source("./02_code/Figures/data_and_libraries.R")

# print stats
tabledat <- scenarios %>%
  filter(cost_per_dose == 12.01 & delivery_cost == 1.62) %>%
  group_by(intervention_f) %>%
  summarize(n = n(),
            median = round(median(CE)),
            q25 = round(quantile(CE, p = 0.25)),
            q75 = round(quantile(CE, p = 0.75)),
            min = round(min(CE)),
            max = round(max(CE))) %>%
  arrange(median)

tabledat; sum(tabledat$n)


# set up order of interventions by median CE for plot
levels <- scenarios %>%
  mutate(group = case_when(intervention_f %in% c('ITN 10% increase', 'ITN PBO', 'SMC', 'RTS,S seasonal', 'RTS,S age-based') ~ 'uni',
                           TRUE ~ 'multi')) %>%
  group_by(intervention_f, group) %>%
  summarize(med = median(CE, na.rm = T) )%>% ungroup() %>%
  arrange(desc(group), med)


# list color order for plots
colors <- c(smc, pbo, itn, rtss_sv, rtss_age, pbo_smc, itn_smc, rtss_smc, pbo_rtss_smc, itn_rtss_smc, pbo_rtss, itn_rtss)


# plot of cost per DALY averted
box_plot <- function(cost_dose){

  scenarios %>%
    filter(cost_per_dose == cost_dose & delivery_cost == 1.62) %>%
    mutate(intervention_f = factor(intervention, levels = levels$intervention_f)) %>%
    mutate(rank = as.numeric(intervention_f)) %>%

    ggplot(aes(x = rank, y = CE, fill = intervention_f,
               color = intervention_f, group = intervention)) +
    geom_hline(yintercept = 0, lty = 2, color = 'grey') +
    geom_vline(xintercept = 5.5, lty = 2, color = 'grey') +
    geom_boxplot(alpha = 0.3, outlier.alpha = 0.05,  outlier.size = 0.05, coef = 500) +
    coord_cartesian(ylim = c(-50, 500), clip = "off") +
    labs(x = '',
         y = expression(paste(Delta," cost / ", Delta, " DALYs")),
         fill = 'intervention',
         color = 'intervention') +
    annotation_custom(textGrob("Univariate strategies"),
                      xmin = 1, xmax = 5,ymin = -100, ymax = -100) +
    annotation_custom(textGrob("Mixed strategies"),
                      xmin = 6, xmax = 12, ymin = -100, ymax = -100) +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    scale_x_continuous(breaks = c(0)) +
    scale_y_continuous(breaks = seq(0, 500, 50)) +
    theme_classic() +
    theme(plot.caption.position = "plot",
          text = element_text(size = 14))

  ggsave(paste0('./03_output/figure2_', cost_dose, '.pdf'), width = 9, height = 4)

}

table(scenarios$cost_per_dose)

box_plot(2.69)
box_plot(6.52)
box_plot(12.01)
box_plot(12.91)
box_plot(17.36)


scenarios %>%
  filter(cost_per_dose == 12.01 & delivery_cost == 1.62) %>%
  mutate(intervention_f = factor(intervention, levels = levels$intervention_f)) %>%
  mutate(rank = as.numeric(intervention_f)) %>%

  ggplot(aes(x = rank, y = CE_case, fill = intervention_f,
             color = intervention_f, group = intervention)) +
  geom_hline(yintercept = 0, lty = 2, color = 'grey') +
  geom_vline(xintercept = 5.5, lty = 2, color = 'grey') +
  geom_boxplot(alpha = 0.3, outlier.alpha = 0.05,  outlier.size = 0.05, coef = 500) +
  coord_cartesian(ylim = c(-50, 300), clip = "off") +
  labs(x = '',
       y = expression(paste(Delta," cost / ", Delta, " cases")),
       fill = 'intervention',
       color = 'intervention') +
  annotation_custom(textGrob("Univariate strategies"),
                    xmin = 1, xmax = 5,ymin = -80, ymax = -80) +
  annotation_custom(textGrob("Mixed strategies"),
                    xmin = 6, xmax = 12, ymin = -80, ymax = -80) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  scale_x_continuous(breaks = c(0)) +
  scale_y_continuous(breaks = seq(0, 300, 50)) +
  theme_classic() +
  theme(plot.caption.position = "plot",
        text = element_text(size = 14))

ggsave(paste0('./03_output/figure2_cases.pdf'), width = 9, height = 4)

