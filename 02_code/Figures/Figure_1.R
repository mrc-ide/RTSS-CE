# Figure 1
# Scenarios, interventions, and timings

# set-up
source("./02_code/Figures/data_and_libraries.R")
library(malariasimulation)

year <- 365
month <- year / 12

# seasonal profiles: c(g0, g[1], g[2], g[3], h[1], h[2], h[3])
# drawn from mlgts: https://github.com/mrc-ide/mlgts/tree/master/data
# g0 = a0, a = g, b = h
seas_name <- 'highly seasonal'
seasonality <- list(c(0.284596,-0.317878,-0.0017527,0.116455,-0.331361,0.293128,-0.0617547))
s1 <- tibble(seasonality, seas_name)

seas_name <- 'seasonal'
seasonality <- list(c(0.285505,-0.325352,-0.0109352,0.0779865,-0.132815,0.104675,-0.013919))
s2 <- tibble(seasonality, seas_name)

seas_name <- 'perennial'
seasonality <- list(c(0.2852770,-0.0248801,-0.0529426,-0.0168910,-0.0216681,-0.0242904,-0.0073646))
s3 <- tibble(seasonality, seas_name)

stable <- bind_rows(s1, s2, s3)

# determine intervention timings
add_dose <- function(x){

  data <- stable[x, ]

  params <- get_parameters(list(
    human_population = 1000,
    model_seasonality = TRUE,
    # rainfall fourier parameters
    g0 = unlist(data$seasonality)[1],
    g = unlist(data$seasonality)[2:4],
    h = unlist(data$seasonality)[5:7],
    individual_mosquitoes = FALSE))

  # assign RTS,S
  data$peak <- peak_season_offset(params) + year
  if(data$seas_name == "highly seasonal"){
    data$dose1 <- round((data$peak - month * 3.5), 0)}
  if(data$seas_name == "seasonal"){
    data$dose1 <- round((data$peak - month * 5.5), 0)}
  if(data$seas_name == "perennial"){
    data$dose1 <- NA}

  data$dose2 <- round(data$dose1 + 1*month)
  data$dose3 <- round(data$dose1 + 2*month)
  data$boost <- data$dose3

  data$EPI <- list(c(year, year * 2, 0, 0.03))

  # assign SMC
  data$smc <- list(rep(NA, 4)) # placeholder

  if(data$seas_name == "highly seasonal"){
    data$smc <- list(c(round(c(data$peak - year + c(-1, 0, 1, 2) * month), 0),
                     round(c(data$peak + c(-1, 0, 1, 2) * month), 0)))
  }
  if(data$seas_name == "seasonal"){
    data$smc <- list(c(round(c(data$peak + c(-2, -1, 0, 1, 2) * month), 0)))
  }

  # assign ITN
  data$ITN <- "ITN *"
  data$time3 <- 1

  return(data)
}

interventions <- map_dfr(c(1:3), add_dose) |>
  rename("seas" = "seasonality", "seasonality" = "seas_name") |>
  rowwise() |>
  mutate(round1 = unlist(smc)[1],
         round2 = unlist(smc)[2],
         round3 = unlist(smc)[3],
         round4 = unlist(smc)[4],
         round5 = unlist(smc)[5],
         round6 = unlist(smc)[6],
         round7 = unlist(smc)[7],
         round8 = unlist(smc)[8]) |>
  pivot_longer(dose1:boost, names_to = "vaccine", values_to = "time1") |>
  pivot_longer(round1:round8, names_to = "SMC", values_to = "time2"); interventions


# determine CI at baseline with no interventions
run_model <- function(x){

  data <- stable[x, ]

  params <- get_parameters(list(
    human_population = 10000,
    clinical_incidence_rendering_min_ages = c(0, 0),
    clinical_incidence_rendering_max_ages = c(5 * year, 200 * year),
    model_seasonality = TRUE,
    # rainfall fourier parameters
    g0 = unlist(data$seasonality)[1],
    g = unlist(data$seasonality)[2:4],
    h = unlist(data$seasonality)[5:7],
    individual_mosquitoes = FALSE))

  params <- set_equilibrium(params, init_EIR = 30)

  print(paste("running model", x))

  output <- run_simulation(6 * 365, params)

  output <- as_tibble(output)
  output$seasonality <- data$seas_name

  output <- output |>
    mutate(timestep = timestep - 4 * 365) |>
    filter(timestep > 0)

  return(output)
}

CI <- map_dfr(c(1:3), run_model)


# clean the output CI data with no interventions
CI_plot <- CI |>
  mutate(month = ceiling(timestep / month)) |>
  mutate(month = case_when(month >= 1 & month <= 12 ~ month - 13,
                           month >= 13 ~ month - 12)) |>
  mutate(c_inc = p_inc_clinical_0_73000 / n_0_73000,
         c_inc_u5 = p_inc_clinical_0_1825 / n_0_1825)


# plot CI and dose timings, seasonal settings only
# text for plot
my_text <- tibble(seasonality = "perennial",
                  lab = c("pre-intervention", "post-intervention"),
                  x = c(year * 0.5, year * 1.5),
                  y = c(0.02, 0.02))
# age-based rectangle
EPI <- tibble(xmin = year, xmax = year * 2, ymin = 0, ymax = 0.003, seasonality = c("perennial", "seasonal", "highly seasonal"))
# SMC coverage rectangle
interventions |> select(seasonality, SMC, time2) |> filter(!is.na(time2)) |> distinct()
SMC <- tibble(xmin = c(533, 201, 201 + year), xmax = c(624 + month, 292 + month, 292 + month + year),
              ymin = rep(0, 3), ymax = rep(Inf, 3), seasonality = c("seasonal", "highly seasonal", "highly seasonal"))

ggplot() +
  stat_smooth(data = CI_plot, aes(x = timestep, y = c_inc_u5), span = .2, se = F,
              color = "grey", fill = "grey", alpha = 0.5, geom = "area") +
  geom_rect(data = EPI,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = "RTS,S age-based"), alpha = 0.1) +
  geom_rect(data = SMC,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = "SMC coverage"), alpha = 0.1) +
  geom_vline(data = interventions, aes(xintercept = time1, color = "RTS,S seasonal"), lty = 2) +
  geom_vline(data = interventions, aes(xintercept = time2, color = "SMC"), lty = 2) +
  geom_vline(data = interventions, aes(xintercept = time3, color = "ITN*"), lty = 2) +
  geom_vline(xintercept = year + 1) +
  geom_text(data = my_text, aes(x = x,  y = y, label = lab)) +
  facet_grid(factor(seasonality,
                    levels = c("perennial", "seasonal", "highly seasonal")) ~ .) +
  coord_cartesian(xlim = c(1, 2 * year - 43), ylim = c(0, max(CI_plot$c_inc_u5))) +
  scale_color_manual(values = c('#1BB6AF','#088BBE','#F6A1A5')) +
  scale_fill_manual(values = c("#088BBE", "#F6A1A5")) +
  labs(x = "month", y = "clinical incidence, 0-5 years", color = "", fill = "") +
  scale_x_continuous(breaks = seq(1, 2 * year + month, month), labels = seq(-12, 12, 1)) +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_line(color = "#F5F5F5", size = 0.2))

# save
ggsave('./03_output/figure1v2.pdf', width = 8, height = 3.5)


