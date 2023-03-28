# Figure 6
# Measuring cost-effectiveness and equity

# set-up
source("./02_code/Figures/data_and_libraries.R")


# read in equity analysis results
output <- readRDS("./03_output/dalyoutput_draws_casestudy.rds") |>
  ungroup() |>
  filter(cost_per_dose == 12.01 & delivery_cost == 1.62)


# gap in PfPR
output_pfpr <- output |>
  filter(pfpr %in% c(0.40, 0.10) & RTSScov %in% c(0, 0.80) & ITNuse == 0.50) |>
  mutate(scenario2 = 1)

# gap in ITN use
output_itn <- output |>
  filter(RTSScov %in% c(0, 0.80) & ((pfpr == 0.20 & ITNuse == 0.60) | (pfpr == 0.10 & ITNuse == 0.30))) |>
  mutate(scenario2 = 2)

# gap in RTSS coverage
output_rtss <- output |>
  filter(ITNuse == 0.50 & ((pfpr == 0.20 & RTSScov %in% c(0, 0.50)) | (pfpr == 0.10 & RTSScov %in% c(0, 0.80)))) |>
  mutate(scenario2 = 3)

# combine
output <- full_join(output_pfpr, output_itn) |> full_join(output_rtss) |>
  mutate(scenario2_f = factor(scenario2,
                              levels=c(1,2,3),
                              labels=c("gap in PfPR", "gap in ITN use", "gap in vaccination")))


# assign scenarios
scenario1 <- output |>
  filter(ITNboost==0 & RTSS=="none") |> mutate(scenario=1)

scenario2 <- output |>
  filter(ITNboost==1 & RTSS=="none") |> mutate(scenario=2)

scenario3 <- output |>
  filter(RTSS=="EPI" & ITNboost==0) |> mutate(scenario=3)

scenario4 <- output |>
  filter((pfpr %in% c(0.20, 0.40) & ITNboost==1 & RTSS=="none") | (pfpr %in% c(0.10) & ITNboost==0 & RTSS=="none")) |>
  mutate(scenario=4)

scenario5 <- output |>
  filter((pfpr %in% c(0.20, 0.40) & ITNboost==0 & RTSS=="EPI") | (pfpr %in% c(0.10) & ITNboost==0 & RTSS=="none")) |>
  mutate(scenario=5)

scenariosb <- full_join(scenario1, scenario2) |> full_join(scenario3) |>
  full_join(scenario4) |> full_join(scenario5) |>
  mutate(scenario_f = factor(scenario,
                             levels=c(1,2,3,4,5),
                             labels=c("baseline","mass ITN 10% increase", "mass age-based RTS,S", "targeted ITN 10% increase", "targeted age-based RTS,S"))) |>
  mutate(residence = ifelse(pfpr %in% c(0.20, 0.40), "rural", "urban")) |>
  dplyr::select(drawID, seasonality, scenario, scenario_f, scenario2, scenario2_f, residence, cost_total, daly, cases, deaths)


# make plot for cases, deaths, or DALYs
plot_equity <- function(var){ # var = cases, deaths, or daly

  # calculate disparity between urban and rural
  scenarios <- scenariosb |> arrange(seasonality, scenario2, scenario2_f, scenario, scenario_f, drawID, residence) |>
    group_by(seasonality, scenario2, scenario2_f, scenario, scenario_f, drawID) |>
    mutate(disparity = abs({{var}} - lag({{var}}))) |>
    summarize(across(c(cost_total, cases, deaths, daly, disparity), sum, na.rm = T))

  # separate out and label baseline values
  none <- scenarios |> ungroup() |>
    filter(scenario == 1) |>
    mutate(cost_total_baseline = cost_total,
           outcome_baseline = {{var}},
           disparity_baseline = disparity) |>
    dplyr::select(seasonality, scenario2, scenario2_f, drawID, cost_total_baseline:disparity_baseline)

  # merge baseline and non-baseline scenarios and calculate difference
  scenarios2 <- scenarios |>
    left_join(none, by = c("seasonality", "scenario2", "scenario2_f", "drawID")) |>
    mutate(cost_diff = cost_total - cost_total_baseline,
           outcome_diff = outcome_baseline - {{var}},
           disparity_per = (disparity - disparity_baseline) / disparity_baseline * 100,
           CE_outcome = cost_diff / outcome_diff) |>
    filter(outcome_diff >= 0) |> # take out scenarios where DALYs increase
    ungroup() |>
    # summarize over drawID
    group_by(seasonality, scenario2, scenario2_f, scenario, scenario_f) |>
    filter(scenario_f != "baseline") |>
    summarize(n = 50,
              estimate_CE = median(CE_outcome, na.rm = T),
              lower_CE = quantile(CE_outcome, 0.025),
              upper_CE = quantile(CE_outcome, 0.975),

              estimate_disparity = median(disparity_per, na.rm = T),
              lower_dis = quantile(disparity_per, 0.025),
              upper_dis = quantile(disparity_per, 0.975)) |>
    # put negative cost effectiveness values up to 0
    mutate(estimate_CE = case_when(estimate_CE < 0 ~ 0,
                                   TRUE ~ estimate_CE),
           lower_CE = case_when(estimate_CE == 0 ~ 0,
                                TRUE ~ lower_CE),
           upper_CE = case_when(estimate_CE == 0 ~ 0,
                                TRUE ~ upper_CE))

  # adding asterisks to negative CE values
  stars <- scenarios2 |> filter(estimate_CE == 0) |>
    select(estimate_CE, estimate_disparity) |>
    mutate(x = estimate_disparity + 2, y = 4)

  if(deparse(substitute(var)) == "daly"){
    name <- "DALYs"
    xlim <- c(-50, 20)
    ylim <- c(0, 300)
  }
  if(deparse(substitute(var)) == "cases"){
    name <- "cases"
    xlim <- c(-65, 10)
    ylim <- c(0, 110)
  }
  if(deparse(substitute(var)) == "deaths"){
    name <- "deaths"
    xlim <- c(min(scenarios2$lower_dis), max(scenarios2$upper_dis))
    ylim <- c(min(scenarios2$lower_CE), max(scenarios2$upper_CE))
  }

  plot <- ggplot(data = scenarios2 |> filter(scenario > 1)) + # remove baseline
    geom_vline(xintercept = 0, lty=2, color="grey") +
    geom_pointrange(aes(x = estimate_disparity, y = estimate_CE, ymin = lower_CE, ymax = upper_CE, shape = scenario2_f, color = scenario_f)) +
    geom_pointrange(aes(x = estimate_disparity, xmin = lower_dis, xmax = upper_dis, y = estimate_CE, shape = scenario2_f, color = scenario_f)) +
    geom_point(data = stars, aes(x = x, y = y), shape = "*", size = 4, color = "black") +
    labs(y = expr(paste("cost-effectiveness (", Delta," cost / ", Delta, " ", !!name, ")")),
         x = paste0("disparity between urban and rural (% change in ", name, ")"),
         shape = "baseline scenario",
         color = "intervention") +
    scale_color_manual(values = c("#C70E7B","#007BC3", "#FC6882","#54BCD1")) +
    scale_x_continuous(labels = scales::percent_format(scale = 1)) +
    coord_cartesian(xlim = xlim, ylim = ylim, clip = "off") +
    theme_classic() +
    theme(text = element_text(size = 14))

  return(plot)

}


# plot cases, deaths, or DALYs
B <- plot_equity(daly)
C <- plot_equity(cases)
D <- plot_equity(deaths)


P1 <- B + geom_hline(yintercept = 0, lty=2, color="grey") +
  labs(y = "more cost-effective",
       x = "more equitable",
       title = "DALYs",
       shape = "baseline scenario",
       color = "intervention") +
  geom_segment(aes(x = 20, xend = -50, y = -40, yend = -40),
               arrow=arrow(length=unit(0.2,"cm")), color = "cornflowerblue", size = 1) +
  geom_segment(aes(x = -57, xend = -57, y = 300, yend = 20),
               arrow=arrow(length=unit(0.2,"cm")), color = "cornflowerblue", size = 1) +
  scale_x_continuous(breaks = c(0)) +
  scale_y_continuous(breaks = c(0))

P2 <- C + geom_hline(yintercept = 0, lty=2, color="grey") +
  labs(y = "more cost-effective",
       x = "more equitable",
       title = "cases",
       shape = "baseline scenario",
       color = "intervention") +
  geom_segment(aes(x = 10, xend = -65, y = -15, yend = -15),
               arrow=arrow(length=unit(0.2,"cm")), color = "cornflowerblue", size = 1) +
  geom_segment(aes(x = -72.5, xend = -72.5, y = 110, yend = 7),
               arrow=arrow(length=unit(0.2,"cm")), color = "cornflowerblue", size = 1) +
  scale_x_continuous(breaks = c(0)) +
  scale_y_continuous(breaks = c(0))


(P1 + P2) + plot_layout(guides = "collect", nrow=1) + plot_annotation(tag_levels = "A")

ggsave(paste0("./03_output/figure6.png"), width = 10, height = 4)

