# Figure S4
# Cost-effectiveness ratios (cost per case averted) of univariate and mixed strategies

# set-up
source("./02_code/Figures/data_and_libraries.R")


# print stats
tabledat <- scenarios %>%
  # remove scenarios with a negative impact on DALYs
  mutate(casediff = cases_baseline - cases) %>%
  filter(casediff >= 0) %>%
  filter(cost_per_dose == 12.01 & delivery_cost == 1.62) %>%
  group_by(intervention_f) %>%
  summarize(n = n(),
            median = round(median(CE_case)),
            q25 = round(quantile(CE_case, p=0.25)),
            q75 = round(quantile(CE_case, p=0.75)),
            min = round(min(CE_case)),
            max = round(max(CE_case))) %>%
  arrange(median)

tabledat; sum(tabledat$n)

# < cost per case averted in cost-effective settings
scenarios %>%
  filter(cost_per_dose == 12.01 & delivery_cost == 1.62) %>%
  filter(pfpr == 0.40 & ITNuse == 0.75) %>%
  filter(seasonality == 'highly seasonal' | seasonality == 'perennial' | (seasonality == 'seasonal' & SMC == 0.85)) %>%
  group_by(intervention_f) %>%
  summarize(n = n(),
            median = round(median(CE_case)),
            q25 = round(quantile(CE_case, p = 0.25)),
            q75 = round(quantile(CE_case, p = 0.75)),
            min = round(min(CE_case)),
            max = round(max(CE_case))) %>%
  arrange(median)


scenarios %>% filter(intervention != 'none') %>%
  # remove scenarios with a negative impact on DALYs
  mutate(casediff = cases_baseline - cases) %>%
  filter(casediff >= 0) %>%
  filter(cost_per_dose == 12.01 & delivery_cost == 1.62) %>%
  mutate(intervention_f = factor(intervention, levels=levels$intervention_f)) %>%
  mutate(rank=as.numeric(intervention_f)) %>%

  ggplot(aes(x=rank, y=CE_case, fill=intervention_f, color=intervention_f, group=intervention)) +
  geom_hline(yintercept = 0, lty=2, color='grey') +
  geom_vline(xintercept = 5.5, lty=2, color='grey') +
  geom_boxplot(alpha=0.3, outlier.alpha = 0.05,  outlier.size = 0.05, coef = 500) +
  coord_cartesian(ylim=c(-1, 80), clip="off") +
  labs(x='',
       y=expression(paste(Delta," cost / ", Delta, " cases")),
       fill = 'intervention',
       color = 'intervention') +
  annotation_custom(textGrob("Univariate strategies"),xmin=1,xmax=5,ymin=-8,ymax=-8) +
  annotation_custom(textGrob("Mixed strategies"),xmin=6,xmax=12,ymin=-8,ymax=-8) +
  scale_fill_manual(values = c(smc, pbo, itn, rtss_sv, rtss_age, pbo_smc, itn_smc, rtss_smc, pbo_rtss_smc, itn_rtss_smc, pbo_rtss, itn_rtss)) +
  scale_color_manual(values = c(smc, pbo, itn, rtss_sv, rtss_age, pbo_smc, itn_smc, rtss_smc, pbo_rtss_smc, itn_rtss_smc, pbo_rtss, itn_rtss)) +
  scale_x_continuous(breaks=c(0)) +
  theme_classic() +
  theme(plot.caption.position = "plot",
        text = element_text(size = 14))

ggsave('./03_output/figureS4.pdf', width = 9, height = 4)


