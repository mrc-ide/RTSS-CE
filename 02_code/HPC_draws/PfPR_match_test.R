# calibrate PfPRs to EIRs

# devtools::install_github('https://github.com/mrc-ide/cali')

library(cali)
library(malariasimulation)
library(ggplot2)

year <- 365

# define target: PfPR2-10 value
target <- c(0.3)
# time point at which to match target
target_tt <- seq(3*year, 6*year, 100)

parameters <- malariasimulation::get_parameters(list(

human_population = 1000,
model_seasonality = TRUE,

g0 = 0.284596,
g = c(-0.317878,-0.0017527,0.116455),
h = c(-0.331361,0.293128,-0.0617547)
))

params <- set_bednets(
  parameters = params,
  timesteps = c(1*year, 4*year),
  coverages = c(0.8, 0.8),
  retention = 3 * year,
  dn0 = matrix(c(0.387, 0.387), nrow=2, ncol=1),
  rn = matrix(c(0.563, 0.563), nrow=2, ncol=1),
  rnm = matrix(c(.24, .24), nrow=2, ncol=1),
  gamman = c(2.64 * 365, 2.64 * 365))


set.seed(123)
out <- calibrate(parameters = parameters,
                 target = target,
                 target_tt = target_tt,
                 summary_function = summary_pfpr_2_10,
                 tolerance = 0.02,
                 interval = c(1, 6))

parameters <- malariasimulation::set_equilibrium(parameters, init_EIR = out$root)
raw <- malariasimulation::run_simulation(timesteps = max(target_tt) + 100, parameters = parameters)
pfpr <- summary_pfpr_2_10(raw)

pd <- data.frame(time = 1:(max(target_tt) + 100), pfpr = pfpr)

pd %>% filter(time > 4*year & time < 6*year) %>% summarize(m = mean(pfpr))

ggplot() +
  geom_line(data = pd, aes(x = time, y = pfpr), col = "deeppink", size = 1) +
  geom_point(aes(x = target_tt, y = target), col = "dodgerblue", size = 4) +
  ylim(0, 1) +
  theme_bw()
