# Figure S5
# Cost-per-dose

# set-up
source("./02_code/Figures/data_and_libraries.R")


# function to find the most common character value in a group
calculate_mode <- function(x) {

  uniqx <- unique(na.omit(x))
  uniqx[which.max(tabulate(match(x, uniqx)))]

}


compare_vax <- function(comparison){

  if(comparison == 'SMC'){
    scenarios2 <- scenarios %>% filter(seasonality == 'seasonal')
  }


  if(comparison == 'ITN PBO'){
    scenarios2 <- scenarios %>% filter(resistance > 0)
  }

  if(comparison == 'ITN 10% increase'){
    scenarios2 <- scenarios
  }

  output <- scenarios2 %>%
    filter(cost_per_dose == 12.01 & delivery_cost == 1.62) %>%
    # combine RTS,S strategies
    mutate(intervention = ifelse(intervention == 'RTS,S age-based' | intervention == 'RTS,S seasonal', 'RTS,S', intervention)) %>%
    filter(intervention %in% c('RTS,S', comparison)) %>%
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
    arrange(ID, costRTSS) %>%
    mutate(comparison = comparison)

  return(output)

}


output <- map_dfr(c('SMC', 'ITN PBO', 'ITN 10% increase'), compare_vax)

output %>% group_by(comparison) %>% summarize(median = median(costRTSS),
                                              min = min(costRTSS),
                                              max = max(costRTSS))


# < line plot ----
# calculate proportion of scenarios at each cost-per-dose value
output2 <- output %>% ungroup() %>% group_by(comparison) %>% arrange(costRTSS) %>%
  mutate(p = (n() - row_number() + 1) / n() * 100) %>%
  ungroup() %>%
  select(comparison, costRTSS, p) %>%
  add_row(comparison = 'SMC', costRTSS = 1000, p = 0)



points_select <- function(x){ # x cost RTS,S dose

  output2 %>% filter(costRTSS <= x) %>% group_by(comparison) %>% top_n(-1) %>% select(p) %>% mutate(cost = x)

}

points <- map_dfr(c(2, 5, 9.3, 13.5), points_select) %>% filter(p > 0.1)

table(output$drawID) # 900 total

output3 <- output %>% ungroup() %>% group_by(drawID, comparison) %>% arrange(costRTSS) %>%
  mutate(p = (n() - row_number() + 1) / n() * 100) %>%
  ungroup()

# create median lines and points for plot
segments <- data.frame(
  x = c(2, -20, 5, -20, 9.3, -10, 10, -20, 13.5, -20),
  xend = c(2, 2, 5, 5, 9.3, 9.3, 10, 10, 13.5, 13.5),
  y = c(-10, dose2, -10, dose5, -10, dose93, -10,  dose10, -10, dose135),
  yend = c(dose2, dose2, dose5, dose5, dose93, dose93, dose10, dose10, dose135, dose135))

segments <- data.frame(
  x = c(2, 5, 9.3, 13.5,
        rep(-10, length(points$cost))),
  xend = c(2, 5, 9.3, 13.5,
           points$cost),
  y = c(-10, -10, -10, -10,
        points$p),
  yend = c(max(points[points$cost == 2, ]$p), max(points[points$cost == 5, ]$p), max(points[points$cost == 9.3, ]$p), max(points[points$cost == 13.5, ]$p),
           points$p))


ggplot(output3) +
  geom_vline(xintercept = 0, color = 'grey') +
  geom_segment(data = segments %>% filter(xend != 13.5),
               aes(x = x, xend = xend, y = y, yend = yend),
               lty = 3, color = 'grey') +
  geom_line(aes(x = costRTSS, y = p, group = interaction(comparison, drawID), color = comparison), size = 1, alpha = 0.05) +
  geom_line(data = output2, aes(x = costRTSS, y = p, group = comparison, color = comparison), size = 1) +
  geom_point(data = points, aes(x = cost, y = p, color = comparison), size = 2) +
  scale_color_manual(values = c(itn, pbo, smc)) +
  theme_classic() +
  labs(y = '% of scenarios where RTS,S is implemented \nprior to an alternative intervention',
       x = 'RTS,S cost per dose (USD)',
       color = 'comparison \nintervention'
  ) +
  scale_x_continuous(breaks = c(-5, 0, 2, 5, 9.3)) +
  coord_cartesian(xlim = c(0, 12), ylim = c(0, 100)) +
  theme(text = element_text(size = 14))


ggsave('./03_output/figureS5.pdf', width = 8, height = 5)


