# Process HPC results ----------------------------------------------------------
# read in files
require(data.table)
library(tidyverse)

# check to make sure you are pulling files from the correct folder to process.
files <- list.files(path = "M:/Hillary/GF-RTSS-CE/03_output/HPC/", pattern = "*.rds", full.names = TRUE)
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

# preview
View(dat[1:1000,c('n_rtss_epi_dose_1', 'n_rtss_mass_booster_1', 'dose1', 'dose4')])

# save
saveRDS(dat,"C:/Users/htopazia/OneDrive - Imperial College London/Github/GF-RTSS-CE//03_output/rtss_raw.rds")

