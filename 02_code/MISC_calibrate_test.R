# testing calibration function for PfPR 2-10 from warmup years 6-9
library(tidyverse)
library(malariasimulation)
library(cali)

year <- 365
p <- get_parameters(list(human_population = 1000, individual_mosquitoes = FALSE))
p$timesteps <- 9 * year # simulation run time = 9 years
target <- 0.3

summary_mean_pfpr_2_10_6y9y <- function(x){

  x$year <- ceiling(x$timestep / year)
  x <- x |> filter(year >= 7)
  prev_2_10 <- mean(x$n_detect_730_3650 / x$n_730_3650)
  return(prev_2_10)

}

# run calibration model
set.seed(123)
out <- cali::calibrate(parameters = p,
                       target = target,
                       summary_function = summary_mean_pfpr_2_10_6y9y,
                       tolerance = 0.02,
                       low = 0.001,
                       high = 500)


# plot calibrated model to check match to target
parameters <- set_equilibrium(p, init_EIR = out)
raw <- run_simulation(parameters$timesteps + 100, parameters = parameters)
summary_mean_pfpr_2_10_6y9y(raw)
pfpr <- raw$n_detect_730_3650 / raw$n_730_3650

pd <- data.frame(time = 1:(parameters$timesteps + 100), pfpr = pfpr)

ggplot() +
  geom_hline(yintercept = target, col = "dodgerblue", lty = 2) +
  geom_line(data = pd, aes(x = time, y = pfpr), col = "deeppink", size = 1) +
  ylim(0, 1) +
  theme_bw()
