# Figure 5
# Cost-per-dose

# set-up
source("./02_code/Figures/data_and_libraries.R")


# function to find the most common character value in a group
calculate_mode <- function(x) {

  uniqx <- unique(na.omit(x))
  uniqx[which.max(tabulate(match(x, uniqx)))]

}

output <- scenarios %>%
  filter(cost_per_dose == 12.01 & delivery_cost == 1.62) %>%
  # choose univariate scenarios
  filter(intervention %in% c('RTS,S seasonal', 'RTS,S age-based', 'SMC', 'none', 'ITN 10% increase', 'ITN PBO')) %>%
  # combine RTS,S strategies
  mutate(intervention = ifelse(intervention == 'RTS,S age-based' | intervention == 'RTS,S seasonal', 'RTS,S', intervention)) %>%
  arrange(ID, drawID) %>%
  group_by(ID, drawID) %>%
  # find minimum CE in each group. If intervention == RTS,S find the second or third lowest CE and move RTS,S to the top
  mutate(CE = ifelse(CE == min(CE) & intervention == 'RTS,S', 100000, CE),
         CE = ifelse(CE == min(CE) & intervention == 'RTS,S', 100000, CE),
         CEmin = min(CE, na.rm = T),
         interventionmin = ifelse(CE == CEmin, intervention, NA),
         interventionmin = calculate_mode(interventionmin)) %>%
  filter(intervention %in% c('RTS,S')) %>%
  rowwise() %>%
  # calculate needed cost of RTS,S to match the min CE intervention
  # CEmin = (cost_total - cost_total_baseline) / (daly_baseline - daly)
  mutate(delta_cost = CEmin * (daly_baseline - daly),
         cost_total = delta_cost + cost_total_baseline,
         cost_vax = cost_total - (cost_ITN + cost_clinical + cost_severe + cost_SMC),
         per_dose = cost_vax / (dose1 + dose2 + dose3 + dose4),
         costRTSScon = per_dose - 1.62,   # subtracting delivery cost
         # cost = dosecost + 0.111*(dosecost) + 0.165*(dosecost + 0.111*dosecost) + 0.14
         # (cost - 0.14) / 1.294315 = dose
         costRTSS = (costRTSScon - 0.14)  / 1.294315  # subtracting consumables cost
  ) %>%
  select(ID, seasonality, pfpr, ITNuse, ITN, RTSS, RTSScov, resistance, SMC, treatment, intervention, interventionmin, CE, CEmin, costRTSS) %>%
  group_by(ID) %>%
  arrange(ID, costRTSS)


# < all seasons box and whisker ----
seasoncosts <- output %>% group_by(seasonality) %>%
  summarize(q25 = round(quantile(costRTSS, probs = 0.25, na.rm = T), 2),
            med = round(median(costRTSS, na.rm = T),2),
            q75 = round(quantile(costRTSS, probs = 0.75, na.rm = T), 2),
            min = round(min(costRTSS, na.rm = T), 2),
            max = round(max(costRTSS, na.rm = T), 2))

A <- ggplot(output) +
  geom_hline(yintercept = 0, color = 'light grey') +
  geom_boxplot(aes(x = factor(seasonality, levels = c('perennial','seasonal','highly seasonal')), y = costRTSS),
               fill = 'cornflower blue', color = 'cornflower blue', alpha = 0.4, outlier.alpha = 0.1,  outlier.size = 0.1, coef = 100) +
  geom_text(data = seasoncosts, aes(x = seasonality, y = q25, label = q25), size = 3, nudge_y = -0.55, nudge_x = -.223) +
  geom_text(data = seasoncosts, aes(x = seasonality, y = med, label = med), size = 3, nudge_y = 0.1, nudge_x = -.2) +
  geom_text(data = seasoncosts %>% filter(seasonality != 'seasonal'),
            aes(x = seasonality, y = q75, label = q75), size = 3, nudge_y = .55, nudge_x = -.2) +
  geom_text(data = seasoncosts %>% filter(seasonality == 'seasonal'),
            aes(x= seasonality, y = q75, label = format(round(q75, 2), nsmall = 2)), size = 3, nudge_y = 1.2, nudge_x = -.2) +
  geom_vline(xintercept = 0, lty = 2, color = 'grey') +
  theme_classic() +
  labs(y='RTS,S cost (USD) per dose', x='') +
  coord_cartesian(ylim = c(-5, 20)) +
  theme(plot.caption.position = "plot",
        text = element_text(size = 14))


# print stats
seasoncosts

output %>% ungroup() %>%
  group_by(interventionmin) %>%
  summarize(q25 = round(quantile(costRTSS, probs = 0.25, na.rm = T),2),
            med = round(median(costRTSS, na.rm = T),2),
            q75 = round(quantile(costRTSS, probs = 0.75, na.rm = T),2))

# how many are above $9.3 - stratified by SMC / ITNuse
output %>% filter(costRTSS >= 9.3) %>% group_by(SMC, ITNuse) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  mutate(t = sum(n), p = n/t * 100)


# how many are above $9.3 - stratified by PfPR
output %>% filter(costRTSS >= 9.3) %>% group_by(pfpr) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  mutate(t = sum(n), p = n/t * 100)

# median cost when ITNs and SMC are not maximized
test <- output %>% filter(ITNuse == 0.75 | SMC == 0.85)
summary(test$costRTSS)

test <- output %>% filter(!(ITNuse == 0.75 | SMC == 0.85))
summary(test$costRTSS)



# < line plot ----
# calculate proportion of scenarios at each cost-per-dose value
output2 <- output %>% ungroup() %>% arrange(costRTSS) %>%
  mutate(p = (nrow(output) - row_number() + 1) / nrow(output) * 100) %>%
  ungroup()

dose2 <- output2 %>% filter(costRTSS <= 2) %>% top_n(-1) %>% select(p) %>% as.numeric()
dose5 <- output2 %>% filter(costRTSS <= 5) %>% top_n(-1) %>% select(p) %>% as.numeric()
dose93 <- output2 %>% filter(costRTSS <= 9.3) %>% top_n(-1) %>% select(p) %>% as.numeric()
dose10 <- output2 %>% filter(costRTSS <= 10) %>% top_n(-1) %>% select(p) %>% as.numeric()
dose135 <-  output2 %>% filter(costRTSS <= 13.5) %>% top_n(-1) %>% select(p) %>% as.numeric()

table(output$drawID) # 450 total

output3 <- output %>% ungroup() %>% group_by(drawID) %>% arrange(costRTSS) %>%
  mutate(p = (n() - row_number() + 1) / n() * 100) %>%
  ungroup()

# create median lines and points for plot
segments <- data.frame(
  x = c(2, -20, 5, -20, 9.3, -10, 13.5, -20),
  xend = c(2, 2, 5, 5, 9.3, 9.3, 13.5, 13.5),
  y = c(-10, dose2, -10, dose5, -10, dose93, -10,  dose135),
  yend = c(dose2, dose2, dose5, dose5, dose93, dose93, dose135, dose135))


points <- data.frame(
  x = c(2, 5, 9.3, 13.5),
  y = c(dose2, dose5, dose93, dose135)
)

B <- ggplot(output3) +
  geom_vline(xintercept = 0, color = 'grey') +
  geom_segment(data = segments,
               aes(x = x, xend = xend, y = y, yend = yend),
               lty = 3, color = 'grey') +
  geom_line(aes(x = costRTSS, y = p, group = drawID), color = 'cornflower blue', size = 1, alpha = 0.1) +
  geom_line(data = output2, aes(x = costRTSS, y = p, group = drawID), color = 'cornflower blue', size = 1) +
  geom_point(data = points, aes(x = x, y = y), color = 'blue', size = 2) +
  theme_classic() +
  labs(y='% of scenarios where \nRTS,S is most cost-effective',
       x='RTS,S cost per dose (USD)'
  ) +
  scale_x_continuous(breaks = c(-5, 0, 2, 5, 9.3, 13.5)) +
  coord_cartesian(xlim = c(0, 15), ylim = c(0, 100)) +
  theme(text = element_text(size = 14))


# combined plot
A + B + plot_annotation(tag_levels = 'A')

   # coord_cartesian(xlim = c(0, 12), ylim = c(0, 100))) +

ggsave('./03_output/figure5.pdf', width = 10, height = 4)


# print stats
100 - dose2
100 - dose5
100 - dose93
100 - dose10
100 - dose135

output2 %>% filter(p == 50) %>% select(costRTSS) %>% as.numeric()

summary(output2$costRTSS)

