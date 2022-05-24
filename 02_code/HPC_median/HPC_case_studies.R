# HPC set-up -------------------------------------------------------------------
library(didehpc)
setwd('Q:/GF-RTSS-CE')

options(didehpc.cluster = "fi--didemrchnb",
        didehpc.username = "htopazia")

source('./02_code/HPC/functions.R') # reference functions

# transfer the new malariasimulation folder manually to contexts or delete and reinstall using conan
# remotes::install_github('mrc-ide/malariasimulation@dev', force=T)

src <- conan::conan_sources("github::mrc-ide/malariasimulation@dev")

ctx <- context::context_save(path = "Q:/contexts",
                             sources = c('./02_code/HPC/functions.R'),
                             packages = c("dplyr", "malariasimulation"),
                             package_sources = src)

share <- didehpc::path_mapping('Home drive', "Q:", '//fi--san03.dide.ic.ac.uk/homes/htopazia', "M:")
config <- didehpc::didehpc_config(shares = share,
                                  use_rrq = FALSE,
                                  cores = 1,
                                  cluster = "fi--didemrchnb", # fi--dideclusthn OR fi--didemrchnb
                                  parallel = FALSE)

# obj <- didehpc::queue_didehpc(ctx, config = config, provision = "upgrade")
obj <- didehpc::queue_didehpc(ctx, config = config)


# Run tasks --------------------------------------------------------------------
library(tidyverse)
year <- 365
population <- 200000
pfpr <- c(0.1, 0.2, 0.4)

# FIRST
seas_name <- 'highly seasonal'
seasonality <- list(c(0.284596,-0.317878,-0.0017527,0.116455,-0.331361,0.293128,-0.0617547))
s1 <- crossing(seasonality, seas_name, pfpr)

# SECOND
seas_name <- 'seasonal'
seasonality <- list(c(0.285505,-0.325352,-0.0109352,0.0779865,-0.132815,0.104675,-0.013919))
s2 <- crossing(seasonality, seas_name, pfpr)

# THIRD
seas_name <- 'perennial'
seasonality <- list(c(0.2852770,-0.0248801,-0.0529426,-0.0168910,-0.0216681,-0.0242904,-0.0073646))

s3 <- crossing(seasonality, seas_name, pfpr)

stable <- rbind(s1, s2, s3)
speciesprop <- data.frame(speciesprop = rbind(list(c(0.25, 0.25, 0.5))),
                          row.names = NULL)

warmup <- 6*year
sim_length <- 15*year

# interventions
ITN <- c('pyr')
ITNuse <-  c(0.30, 0.50, 0.60)
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

name <- "casestudy"

# create combination of all runs and remove non-applicable scenarios
combo <- crossing(population, stable, warmup, sim_length, speciesprop, interventions) %>%
  mutate(name = paste0(name, "_", row_number())) %>%
  filter(!(RTSS == "none" & RTSScov > 0)) %>% # cannot set RTSS coverage when there is no RTSS
  filter(!(RTSScov == 0 & RTSS == "EPI")) %>% # cannot set 0% coverage if RTSS is implemented
  filter(!(SMC > 0 & seas_name == "perennial")) %>% # do not administer SMC in perennial settings
  filter(!(SMC == 0 & seas_name == "highly seasonal")) %>% # always introduce SMC in highly seasonal settings
  filter(!(seas_name == 'seasonal'))

# EIR / prev match from 'PfPR_EIR_match.R'
match <- readRDS("./02_code/HPC/EIRestimates.rds")

combo <- combo %>% left_join(match %>% select(-ITN)) %>%
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
         name               # name of output file
  ) %>% as.data.frame()

# Run tasks
combo <- combo %>% mutate(f = paste0("./03_output/HPC/",combo$name,".rds")) %>%
  mutate(exist=case_when(file.exists(f) ~ 1, !file.exists(f) ~ 0)) %>%
  filter(exist==0) %>%
  select(-f, -exist)

t <- obj$enqueue_bulk(combo, runsimGF)
t$status()

beepr::beep(1)

# -#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# test run by hand
# x <- 4
# test <- runsimGF(population = combo[[x,1]],
#                  seasonality = combo[[x,2]],
#                  seas_name = combo[[x,3]],
#                  starting_EIR = combo[[x,4]],
#                  pfpr = combo[[x,5]],
#                  warmup = combo[[x,6]],
#                  sim_length = combo[[x,7]],
#                  speciesprop = combo[[x,8]],
#                  ITN = combo[[x,9]],
#                  ITNuse = combo[[x,10]],
#                  ITNboost = combo[[x,11]],
#                  resistance = combo[[x,12]],
#                  IRS = combo[[x,13]],
#                  treatment = combo[[x,14]],
#                  SMC = combo[[x,15]],
#                  RTSS = combo[[x,16]],
#                  RTSScov = combo[[x,17]],
#                  fifth = combo[[x,18]],
#                  name = combo[[x,19]])
# -#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# DHS processing ---------------------------------------------------------------

require(data.table)
library(tidyverse)

data.dir <- 'C:/Users/htopazia/OneDrive - Imperial College London/Github/GF-RTSS-CE/'

# pull all .rds files from HPC output folder and combine
files <- list.files(path = "Q:/GF-RTSS-CE/03_output/HPC_median//", pattern = "casestudy_*", full.names = TRUE)
dat_list <- lapply(files, function (x) readRDS(x))
dat <- rbindlist(dat_list, fill = TRUE, idcol="file")

# add vaccine doses
dat <- dat %>% rowwise() %>%
  mutate(dose1 = case_when(RTSS=="none" ~ 0,
                           #RTSS=="SV" ~ n_rtss_mass_dose_1,
                           RTSS=="EPI" ~ n_rtss_epi_dose_1),
         dose2 = case_when(RTSS=="none" ~ 0,
                           #RTSS=="SV" ~ n_rtss_mass_dose_2,
                           RTSS=="EPI" ~ n_rtss_epi_dose_2),
         dose3 = case_when(RTSS=="none" ~ 0,
                           #RTSS=="SV" ~ n_rtss_mass_dose_3,
                           RTSS=="EPI" ~ n_rtss_epi_dose_3),
         dose4 = case_when(RTSS=="none" ~ 0,
                           #RTSS=="SV" ~ n_rtss_mass_booster_1,
                           RTSS=="EPI" ~ n_rtss_epi_booster_1)) %>% ungroup()

# check prevalence - that the EIR used in the simulation results in the matching prevalence value
test <- dat %>%
  filter((seasonality != 'highly seasonal' & RTSS=='none' & ITN=='pyr' & SMC==0) |
           (seasonality=="highly seasonal" & RTSS=='none' & ITN=='pyr')) %>%
  filter(year %in% c(1,2,3)) %>%
  group_by(seasonality, ITNuse, EIR, pfpr) %>%
  summarize(prev = mean(n_detect_730_3650/n_730_3650))


# summarize data over the first 15 years (mult of 3 for ITNs)
dat2 <- dat %>%
  filter(year <= 15) %>% # first 15 years
  group_by(file) %>%
  mutate_at(vars(n_0_91.25:n_36500_73000), mean, na.rm = TRUE) %>%   # mean of n in each age group
  mutate_at(vars(n_inc_severe_0_91.25:dose4), sum, na.rm = TRUE) %>% # sum of cases and vax doses
  select(-month, -year) %>%
  distinct()

# calculate outputs by age
dat3 <- dat2 %>%
  dplyr::select(file:n_36500_73000, n_inc_clinical_0_91.25:n_inc_clinical_36500_73000,
                n_inc_severe_0_91.25:n_inc_severe_36500_73000, n_treated, n_infections, dose1:dose4) %>%
  # moving to long age groups
  pivot_longer(cols = c(n_0_91.25:n_36500_73000,
                        n_inc_clinical_0_91.25:n_inc_clinical_36500_73000,
                        n_inc_severe_0_91.25:n_inc_severe_36500_73000),
               names_to = c('age'), values_to = c('value')) %>%
  mutate(n = ifelse(grepl('n_[[:digit:]]', age), value, NA),             # creating var for age group
         inc_clinical = ifelse(grepl('n_inc_clinical', age), value, NA), # creating var for inc_clinical
         inc_severe = ifelse(grepl('n_inc_severe', age), value, NA),     # creating var for inc_severe
         age = gsub('n_inc_clinical_', '', age),                         # combining age vars
         age = gsub('n_inc_severe_', '', age),
         age = gsub('n_', '', age),
         age = gsub('_', '-', age)) %>%
  group_by(file, age) %>%
  select(-value) %>%
  mutate_at(vars(inc_clinical:n), sum, na.rm = TRUE) %>% # consolidate
  distinct() %>% ungroup()

# save
saveRDS(dat3, paste0(data.dir, "03_output/rtss_long_casestudy.rds"))

