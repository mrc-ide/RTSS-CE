# Process HPC results
require(data.table)
library(tidyverse)

data.dir <- 'C:/Users/htopazia/OneDrive - Imperial College London/Github/GF-RTSS-CE/'


# raw data ---------------------------------------------------------------------
# pull all .rds files from HPC output folder and combine
files <- list.files(path = "Q:/GF-RTSS-CE/03_output/HPC_median/", pattern = "general_*", full.names = TRUE)
dat_list <- lapply(files, function (x) readRDS(x))
dat <- rbindlist(dat_list, fill = TRUE, idcol="file")

# add vaccine doses
dat <- dat %>% rowwise() %>%
  mutate(dose1 = case_when(RTSS=="none" ~ 0,
                           RTSS=="EPI" ~ n_rtss_epi_dose_1,
                           RTSS=="SV" ~ n_rtss_mass_dose_1),
         dose2 = case_when(RTSS=="none" ~ 0,
                           RTSS=="EPI" ~ n_rtss_epi_dose_2,
                           RTSS=="SV" ~ n_rtss_mass_dose_2),
         dose3 = case_when(RTSS=="none" ~ 0,
                           RTSS=="EPI" ~ n_rtss_epi_dose_3,
                           RTSS=="SV" ~ n_rtss_mass_dose_3),
         dose4 = case_when(RTSS=="none" ~ 0,
                           RTSS=="EPI" ~ n_rtss_epi_booster_1,
                           RTSS=="SV" ~ n_rtss_mass_booster_1)) %>% ungroup()

# save
saveRDS(dat, paste0(data.dir, "03_output/rtss_raw.rds"))

# check prevalence - that the EIR used in the simulation results in the matching prevalence value
test <- dat %>%
  filter((RTSS=='none' & ITN=='pyr' & resistance==0) &
           (SMC==0 | (seasonality=="highly seasonal" & SMC==0.85))) %>%
  filter(year %in% c(1,2,3)) %>%
  group_by(seasonality, ITNuse, treatment, EIR, pfpr) %>%
  summarize(prev = mean(n_detect_730_3650/n_730_3650))


# long format by age group -----------------------------------------------------
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
  mutate_at(vars(n, inc_clinical, inc_severe), sum, na.rm = TRUE) %>% # consolidate
  distinct() %>% ungroup()

# save
saveRDS(dat3, paste0(data.dir, "03_output/rtss_long.rds"))


