library(malariasimulation)
# remotes::install_github('mrc-ide/malariasimulation@test/severe_demography', force=T)
library(tidyverse)
library(mgcv)

# SET SEASONALITY --------------------------------------------------------------

# parameters
year <- 365
human_population <- 100000

params <- get_parameters(list(
  human_population = human_population,
  # average_age = 8453.323,     # to match flat_demog
  model_seasonality = TRUE,   # assign seasonality
  # g0 = 0.285505,                              # low
  # g = c(-0.325352, -0.0109352, 0.0779865),    # low
  # h = c(-0.132815, 0.104675, -0.013919),      # low
  g0 = 0.284596,                              # high
  g = c(-0.317878, -0.0017527, 0.116455),     # high
  h = c(-0.331361, 0.293128, -0.0617547),     # high
  prevalence_rendering_min_ages = 2 * year,
  prevalence_rendering_max_ages = 10 * year,
  individual_mosquitoes = FALSE,
  # individual_mosquitoes = TRUE,
  severe_enabled = T))

# set flat demography to match Haley's:
# https://github.com/ht1212/seasonal_use_case/blob/main/Part_1/2_model_function.R#L44
flat_demog <- read.table('.//01_data/Flat_demog.txt') # from mlgts

ages <- round(flat_demog$V3 * year) # top of age bracket

deathrates <- flat_demog$V5 / 365 # age-specific death rates

params <- set_demography(
  params,
  agegroups = ages,
  timesteps = 1,
  deathrates = matrix(deathrates, nrow = 1),
  birthrates = find_birthrates(human_population, ages, deathrates)
)

# set species / drugs / treatment parameters
params <- set_species(params, species = list(arab_params, fun_params, gamb_params),
                      proportions = c(0.25, 0.25, 0.5))
params <- set_drugs(params, list(AL_params))
params <- set_clinical_treatment(params, 1, c(1), c(0.45))


# RUN MODEL --------------------------------------------------------------------

# loop over malariasimulation runs
init_EIR <- c(0.01, 0.1, seq(1,9,1), seq(10, 100, by=5)) # set EIR values

outputs <- lapply(
  init_EIR,
  function(init) {
    p_i <- set_equilibrium(params, init)
    run_simulation(5 * year, p_i)
  }
)

# output EIR values
EIR <- lapply(
  outputs,
  function(output) {
    mean(
      rowSums(
        output[
          output$timestep %in% seq(4 * 365, 5 * 365, 1),
          grepl('EIR_', names(output))
        ] / human_population * year
      )
    )
  }
)

# output prev 2-10 values
prev <- lapply(
  outputs,
  function(output) {
    mean(
      output[
        output$timestep %in% seq(4 * 365, 5 * 365),
        'n_detect_730_3650'
      ] / output[
        output$timestep %in% seq(4 * 365, 5 * 365),
        'n_730_3650'
      ]
    )
  }
)

# create dataframe of initial EIR, output EIR, and prev 2-10 results
EIR_prev <- cbind(init_EIR, unlist(EIR), unlist(prev)) %>%
  as_tibble() %>% rename(EIR = V2, prev = V3)


# save data for later use
# saveRDS(EIR_prev, 'C:/Users/htopazia/OneDrive - Imperial College London/Github/rtss_malariasimulation/rds/EIR_prev_lowF.rds')    # low, individual_mosquitoes=F
# saveRDS(EIR_prev, 'C:/Users/htopazia/OneDrive - Imperial College London/Github/rtss_malariasimulation/rds/EIR_prev_highF.rds')   # high, individual_mosquitoes=F



# PLOT RESULTS -----------------------------------------------------------------
EIR_prev1 <- readRDS('C:/Users/htopazia/OneDrive - Imperial College London/Github/rtss_malariasimulation/rds/EIR_prev_highF.rds') %>% mutate(season = 'high seasonality')

EIR_prev2 <- readRDS('C:/Users/htopazia/OneDrive - Imperial College London/Github/rtss_malariasimulation/rds/EIR_prev_lowF.rds') %>% mutate(season = 'low seasonality')

EIR_prev <- rbind(EIR_prev1, EIR_prev2)

# calculate EIR / prev 2-10 relationship from malariEquilibrium
eir <- seq(from = .1, to = 100, by=.1)
eq_params <- malariaEquilibrium::load_parameter_set("Jamie_parameters.rds")

prev <- vapply( # calculate prevalence between 2:10 for a bunch of EIRs
  eir,
  function(eir) {
    eq <- malariaEquilibrium::human_equilibrium(
      eir,
      ft=0,
      p=eq_params,
      age=0:100
    )
    sum(eq$states[2:10, 'pos_M']) / sum(eq$states[2:10, 'prop'])
  },
  numeric(1)
)

prevmatch <- as_tibble(cbind(eir,prev))

colors <- c("malariaEquilibrium" = "blue", "malariasimulation" = "red")


p <- ggplot(data=EIR_prev) +
  geom_point(aes(x=init_EIR, y=prev), color='black') +
  stat_smooth(aes(x=init_EIR, y=prev, color='malariasimulation'), method = 'gam', n=1000) +
  geom_line(data=prevmatch, aes(x=eir, y=prev, color = 'malariaEquilibrium')) +
  labs(x='initial EIR', y=expression(paste(italic(Pf),"PR"[2-10])), color='model') +
  facet_wrap('season', nrow=1) +
  scale_color_manual(values = colors) +
  theme_classic()

p


# MATCH EIR AND PfPR -----------------------------------------------------------

# extract points from stat_smooth
p2 <- ggplot_build(p)
p2 <- p2$data[[2]]

p2 <- p2 %>% as_tibble() %>% select(x,y,PANEL) %>% rename(init_EIR=x, pred=y, season=PANEL) %>%
  mutate(season=ifelse(season==1,'high seasonality','low seasonality'))


# Pre-intervention baseline PfPR2-10 starting at values 10, 25, 35, 55
PfPR <- as_tibble_col(c(.10, .25, .35, .55), column_name="pfpr")

# match via points
PfPR %>%
  fuzzyjoin::difference_left_join(EIR_prev, by=c("pfpr"="prev"),
                                  max_dist=1, distance_col="dist") %>%
  group_by(pfpr, season) %>% slice_min(dist)


# match via stat_smooth predictions
PfPR %>%
  fuzzyjoin::difference_left_join(p2, by=c("pfpr"="pred"),
                                  max_dist=1, distance_col="dist") %>%
  group_by(pfpr, season) %>% slice_min(dist)

#

