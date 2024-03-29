# HPC set-up -------------------------------------------------------------------
library(didehpc)
setwd("M:/Hillary/RTSS-CE")

# to edit HPC username and password below
# usethis::edit_r_environ

sources <- c("./02_code/HPC_draws/functions_draws.R",               # parameter draw runs
             "./02_code/HPC_draws/Processing/HPC_processing.R",     # one line per age group
             "./02_code/HPC_draws/Processing/deaths_dalys.R",       # calc deaths & dalys
             "./02_code/HPC_draws/Processing/add_costs.R",          # add cost estimates
             "./02_code/HPC_draws/Processing/outcome_averted.R",    # cases & dalys averted
             "./02_code/HPC_draws/Processing/cost_effectiveness.R") # calc CE

ctx <- context::context_save(path = "M:/Hillary/RTSS-CE/contexts/",
                             sources = sources,
                             packages = c("dplyr", "malariasimulation", "purrr", "tidyr"))

share <- didehpc::path_mapping("malaria", "M:", "//fi--didenas1/malaria", "M:")

config <- didehpc::didehpc_config(credentials = list(
                                  username = Sys.getenv("DIDE_USERNAME"),
                                  password = Sys.getenv("DIDE_PASSWORD")),
                                  shares = share,
                                  use_rrq = FALSE,
                                  cores = 1,
                                  cluster = "fi--didemrchnb", # fi--dideclusthn OR fi--didemrchnb
                                  template = "GeneralNodes", # "GeneralNodes", "12Core", "16Core", "12and16Core", "20Core", "24Core", "32Core"
                                  parallel = FALSE)

# transfer the new malariasimulation folder manually to contexts or delete and re-install using conan
# obj <- didehpc::queue_didehpc(ctx, config = config, provision = "later")
# obj$install_packages("mrc-ide/cali")
# obj$install_packages("mrc-ide/individual")
# obj$install_packages("mrc-ide/malariaEquilibrium")
# obj$install_packages("mrc-ide/malariasimulation")

obj <- didehpc::queue_didehpc(ctx, config = config)



# Set up your job --------------------------------------------------------------
library(tidyverse)
year <- 365

# population
population <- 200000

# seasonal profiles: c(g0, g[1], g[2], g[3], h[1], h[2], h[3])
# drawn from mlgts: https://github.com/mrc-ide/mlgts/tree/master/data
# g0 = a0, a = g, b = h
pfpr <- c(.10, .20, .40) # start at 10% (min rec for RTSS use)

  # FIRST
  seas_name <- "highly seasonal"
  seasonality <- list(c(0.284596,-0.317878,-0.0017527,0.116455,-0.331361,0.293128,-0.0617547))
  s1 <- crossing(seasonality, seas_name, pfpr)

  # SECOND
  seas_name <- "seasonal"
  seasonality <- list(c(0.285505,-0.325352,-0.0109352,0.0779865,-0.132815,0.104675,-0.013919))
  s2 <- crossing(seasonality, seas_name, pfpr)

  # THIRD
  seas_name <- "perennial"
  seasonality <- list(c(0.2852770,-0.0248801,-0.0529426,-0.0168910,-0.0216681,-0.0242904,-0.0073646))
  s3 <- crossing(seasonality, seas_name, pfpr)

stable <- rbind(s1, s2, s3)

# vectors: list(arab_params, fun_params, gamb_params)
speciesprop <- data.frame(speciesprop = rbind(list(c(0.25, 0.25, 0.5))),
                          row.names = NULL)

# run time
warmup <- 21*year # needs to be multiple of 3 so that ITN dist. will line up with first timestep
sim_length <- 15*year

# interventions
ITN <- c("pyr", "pbo")
ITNuse <- c(0, 0.141, 0.339, 0.641) # values 0, 0.2, 0.4, 0.6 calculated using MISC_ITN_usage_distribution.R
ITNboost <- c(0, 1)
resistance <- c(0, 0.4, 0.8)
IRS <-  c(0)
treatment <- c(0.30, 0.45, 0.60)
SMC <- c(0, 0.85)
RTSS <- c("none", "EPI", "SV") # leave out hybrid for now
RTSScov <- c(0, 0.85) # MVIP: dose 1 range 74-93%, dose 3 63-82%, dose 4 42-46%; first half of 2021
fifth <- c(0) # only one booster for now

interventions <-
  crossing(ITN, ITNuse, ITNboost, resistance, IRS, treatment, SMC, RTSS, RTSScov, fifth)

drawID <- c(1:50)

# create combination of all runs and remove non-applicable scenarios
combo <- crossing(population, stable, warmup, sim_length, speciesprop, interventions, drawID) |>
  mutate(ID = paste(pfpr, seas_name, ITNuse,
        resistance, treatment, sep="_")) |>
  filter(!(RTSS == "none" & RTSScov > 0)) |> # cannot set RTSS coverage when there is no RTSS
  filter(!(RTSScov == 0 & (RTSS == "EPI" | RTSS == "SV" | RTSS == "hybrid"))) |> # cannot set 0% coverage if RTSS is implemented
  filter(!(fifth == 1 & (RTSS == "EPI" | RTSS == "none"))) |> # no administration of fifth doses with EPI or no RTSS
  filter(!(seas_name == "perennial" & (RTSS == "SV" | RTSS == "hybrid"))) |> # do not introduce seasonal vaccination in perennial settings
  filter(!(SMC > 0 & seas_name == "perennial")) |> # do not administer SMC in perennial settings
  filter(!(SMC == 0 & seas_name == "highly seasonal")) |> # always introduce SMC in highly seasonal settings
  filter(!(ITNuse == 0 & resistance != 0)) |> # do not introduce resistance when ITNuse==0
  filter(!(ITN == "pbo" & resistance == 0)) |> # only introduce PBO in areas that have resistance
  filter(!(ITN == "pbo" & ITNboost == 1)) |> # cannot have ITN boost + ITN PBO
  filter(!(ITN == "pbo" & ITNuse == 0)) # cannot switch to PBO when ITNuse==0

# EIR / prev match from "PfPR_EIR_match.R"
match <- readRDS("./02_code/HPC_draws/EIRestimates.rds")

combo <- combo |> left_join(match, by = c("drawID", "ID")) |>
  # put variables into the same order as function arguments
  select(population,        # simulation population
         seasonality,       # seasonal profile
         seas_name,         # name of seasonal profile
         starting_EIR,      # equilibrium EIR
         pfpr,              # corresponding PfPR
         warmup,            # warm-up period
         sim_length,        # length of simulation run
         speciesprop,       # proportion of each vector species
         ITN,               # ITN type - pyr, pbo
         ITNuse,            # ITN usage
         ITNboost,          # if ITN usage is boosted by 10%
         resistance,        # resistance level - none 0, med 0.4, high 0.8
         IRS,               # IRS coverage
         treatment,         # treatment coverage
         SMC,               # SMC coverage
         RTSS,              # RTS,S strategy
         RTSScov,           # RTS,S coverage
         fifth,             # status of 5th dose for SV or hybrid strategies
         ID,                # name of output file
         drawID             # parameter draw no.
  ) |> as.data.frame()


saveRDS(combo, "./02_code/HPC_draws/scenarios.rds")


# Run tasks --------------------------------------------------------------------
x = c(1:nrow(combo)) # 128,700 runs

index <- crossing(x)

# remove ones that have already been run
index <- index |> mutate(f = paste0("./03_output/HPC/general_", index$x, ".rds")) |>
  mutate(exist = case_when(file.exists(f) ~ 1, !file.exists(f) ~ 0)) |>
  filter(exist == 0) |>
  select(-f, -exist)

# run a test with the first scenario
# index <- index |> filter(x == 1)

# submit all remaining tasks
# t <- obj$enqueue_bulk(index, runsimGF)
# t$status()

# submit jobs, 100 as a time
sjob <- function(x, y){

  t <- obj$enqueue_bulk(index[x:y,], runsimGF)
  return(1)

}

map2_dfr(seq(0, nrow(index) - 100, 100),
         seq(99, nrow(index), 100),
         sjob)



# Processing -------------------------------------------------------------------
combo <- readRDS("./02_code/HPC_draws/scenarios.rds")

index <- seq(1, nrow(combo), 1)

x <- seq(1, length(index) - 1000, 1000)
y <- seq(1000, length(index), 1000)

data <- tibble(x, y)
data <- rbind(data, c(128001, length(index))) # add final row

data <- data |> mutate(f = paste0("./03_output/HPC_processing/run_", x, "_", y, ".rds")) |>
  mutate(exist = case_when(file.exists(f) ~ 1, !file.exists(f) ~ 0)) |>
  filter(exist == 0) |>
  select(-f, -exist)

t <- obj$enqueue_bulk(data, cost_effectiveness)
t$status()


# read in all parameter draw runs and process
files <- list.files(path = "./03_output/HPC_processing/", pattern = "run*", full.names = TRUE)
dat_list <- lapply(files, function (x) readRDS(x))
dalyoutput <- data.table::rbindlist(dat_list, fill = TRUE, idcol = "identifier")

# take out areas where ITN use == 0 and where resistance is high and ITN 10% increase intervention is used
dalyoutput <- dalyoutput |>
  filter(!(resistance == 0.8 & ITNboost == 1)) |>
  filter(ITNuse != 0)

# save output
saveRDS(dalyoutput, "./03_output/dalyoutput_draws.rds")

# calculate cases / DALYs averted
source("./02_code/HPC_draws/Processing/outcome_averted.R")
output <- outcome_averted(dalyoutput)

# save output
saveRDS(output, "./03_output/scenarios_draws.rds")



# Case study -------------------------------------------------------------------
year <- 365
population <- 200000
pfpr <- c(0.1, 0.2, 0.4)

# FIRST
seas_name <- "highly seasonal"
seasonality <- list(c(0.284596,-0.317878,-0.0017527,0.116455,-0.331361,0.293128,-0.0617547))
s1 <- crossing(seasonality, seas_name, pfpr)

# SECOND
seas_name <- "seasonal"
seasonality <- list(c(0.285505,-0.325352,-0.0109352,0.0779865,-0.132815,0.104675,-0.013919))
s2 <- crossing(seasonality, seas_name, pfpr)

# THIRD
seas_name <- "perennial"
seasonality <- list(c(0.2852770,-0.0248801,-0.0529426,-0.0168910,-0.0216681,-0.0242904,-0.0073646))

s3 <- crossing(seasonality, seas_name, pfpr)

stable <- rbind(s1, s2, s3)
speciesprop <- data.frame(speciesprop = rbind(list(c(0.25, 0.25, 0.5))),
                          row.names = NULL)

warmup <- 21*year
sim_length <- 15*year

# interventions
ITN <- c("pyr")
ITNuse <-  c(0.231, 0.473, 0.641) # 0.30, 0.50, 0.60
ITNboost <- c(0,1)
resistance <- c(0)
IRS <-  c(0)
treatment <- c(0.45)
SMC <- c(0, 0.85)
RTSS <- c("none", "EPI")
RTSScov <- c(0, 0.5, 0.8)
fifth <- c(0)

interventions <-
  crossing(ITN, ITNuse, ITNboost, resistance, IRS, treatment, SMC, RTSS, RTSScov, fifth)

drawID <- c(1:50)

name <- "casestudy"

# create combination of all runs and remove non-applicable scenarios
combo <- crossing(population, stable, warmup, sim_length, speciesprop, interventions, drawID) |>
  mutate(ID = paste(pfpr, seas_name, ITNuse,
                    resistance, treatment, sep="_")) |>
  mutate(name = paste0(name, "_", row_number())) |>
  filter(!(RTSS == "none" & RTSScov > 0)) |> # cannot set RTSS coverage when there is no RTSS
  filter(!(RTSScov == 0 & RTSS == "EPI")) |> # cannot set 0% coverage if RTSS is implemented
  filter(!(SMC > 0 & seas_name == "perennial")) |> # do not administer SMC in perennial settings
  filter(!(SMC == 0 & seas_name == "highly seasonal")) |> # always introduce SMC in highly seasonal settings
  filter(!(seas_name == "seasonal"))


# EIR / prev match from "PfPR_EIR_match.R"
match <- readRDS("./02_code/HPC_draws/EIRestimates.rds")

combo <- combo |> left_join(match, by = c("drawID", "ID"), multiple = "first") |>
  # put variables into the same order as function arguments
  select(population,        # simulation population
         seasonality,       # seasonal profile
         seas_name,         # name of seasonal profile
         starting_EIR,      # equilibrium EIR
         pfpr,              # corresponding PfPR
         warmup,            # warm-up period
         sim_length,        # length of simulation run
         speciesprop,       # proportion of each vector species
         ITN,               # ITN type - pyr, pbo
         ITNuse,            # ITN usage
         ITNboost,          # if ITN usage is boosted by 10%
         resistance,        # resistance level - none 0, med 0.4, high 0.8
         IRS,               # IRS coverage
         treatment,         # treatment coverage
         SMC,               # SMC coverage
         RTSS,              # RTS,S strategy
         RTSScov,           # RTS,S coverage
         fifth,             # status of 5th dose for SV or hybrid strategies
         ID,                # name of output file
         drawID             # parameter draw no.
  ) |> as.data.frame()


saveRDS(combo, "./02_code/HPC_draws/scenarios_casestudy.rds")


# < Run tasks ------------------------------------------------------------------
x = c(1:nrow(combo)) # 5,400 runs

index <- crossing(x)

# remove ones that have already been run
index <- index |> mutate(f = paste0("./03_output/HPC/casestudy_", index$x, ".rds")) |>
  mutate(exist = case_when(file.exists(f) ~ 1, !file.exists(f) ~ 0)) |>
  filter(exist == 0) |>
  select(-f, -exist)

# run a test with the first scenario
# index <- index |> filter(x == 1)

# submit all remaining tasks
# t <- obj$enqueue_bulk(index, runsimGF_casestudy)
# t$status()

# submit jobs, 100 as a time
sjob <- function(x, y){

  t <- obj$enqueue_bulk(index[x:y,], runsimGF_casestudy)
  return(1)

}

map2_dfr(seq(0, nrow(index) - 100, 100),
         seq(99, nrow(index), 100),
         sjob)



# < Processing -----------------------------------------------------------------
# read in all parameter draw runs and process

# tell R where to find new functions to process case study results
sources <- c("./02_code/HPC_draws/functions_draws.R",               # parameter draw runs
             "./02_code/HPC_draws/Processing/deaths_dalys.R",       # calc deaths & dalys
             "./02_code/HPC_draws/Processing/add_costs.R",          # add cost estimates
             "./02_code/HPC_draws/Processing/outcome_averted.R",    # cases & dalys averted
             "./02_code/HPC_draws/Processing_casestudy/HPC_processing.R",     # one line per age group
             "./02_code/HPC_draws/Processing_casestudy/cost_effectiveness.R") # calc CE

ctx <- context::context_save(path = "M:/Hillary/RTSS-CE/contexts/",
                             sources = sources,
                             packages = c("dplyr", "malariasimulation", "purrr", "tidyr"))

obj <- didehpc::queue_didehpc(ctx, config = config)


# read in all parameter draw runs and process
t <- obj$enqueue(cost_effectiveness(1))
t$status()

table(combo$seas_name, combo$ITNuse, useNA = "always")
table(combo$seas_name, combo$RTSScov, useNA = "always")
table(combo$seas_name, combo$pfpr, useNA = "always")
table(output$seasonality, output$ITNuse, useNA = "always")
table(output$seasonality, output$RTSScov, useNA = "always")
table(output$seasonality, output$pfpr, useNA = "always")


