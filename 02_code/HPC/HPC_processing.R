# Process HPC results ----------------------------------------------------------
# read in files
require(data.table)
library(tidyverse)

# check to make sure you are pulling files from the correct folder to process.
files <- list.files(path = "M:/Hillary/rtss_malariasimulation/rds/HPC", pattern = "*.rds", full.names = TRUE)
dat_list <- lapply(files, function (x) data.table(readRDS(x)))

dat <- rbindlist(dat_list, fill = TRUE, idcol="file") %>%
  mutate(file = files[file],
         file = gsub("M:/Hillary/rtss_malariasimulation/rds/HPC/","",file),
         intervention = case_when(grepl('EPIalone',file)~"EPI",
                                  grepl('SV4_',file)~"SV 4-dose",
                                  grepl('SV5_',file)~"SV 5-dose",
                                  grepl('SV4updated_',file)~"SV 4-dose - updated booster",
                                  grepl('SV5updated_',file)~"SV 5-dose - updated booster",
                                  grepl('SMCEPI_',file)~"EPI + SMC",
                                  grepl('SV4SMC_',file)~"SV 4-dose + SMC",
                                  grepl('SV5SMC_',file)~"SV 5-dose + SMC",
                                  grepl('SV4updatedSMC',file)~"SV 4-dose - updated booster + SMC",
                                  grepl('SV5updatedSMC',file)~"SV 5-dose - updated booster + SMC",
                                  grepl('SV4SMCsynergy',file)~"SV 4-dose synergy + SMC",
                                  grepl('SV5SMCsynergy',file)~"SV 5-dose synergy + SMC",
                                  grepl('SMCalone',file)~"SMC alone",
                                  grepl('none',file)~"none"),
         season = case_when(model=='high_seas'~"highly_seasonal",
                            model=='low_seas'~"seasonal"))

saveRDS(dat,"C:/Users/htopazia/OneDrive - Imperial College London/Github/rtss_malariasimulation/rds/rtss_smc_raw.rds")
#saveRDS(dat,"C:/Users/htopazia/OneDrive - Imperial College London/Github/rtss_malariasimulation/rds/severe_test.rds")

summary(dat$n_0_36500)

dat2 <- dat %>%
  mutate(cases = (n_inc_clinical_0_36500/n_0_36500) * human_population,
         n = n_0_36500,
         cases_p = (p_inc_clinical_0_36500/n_0_36500) * human_population,
         severe = (n_inc_severe_0_36500/n_0_36500) * human_population,
         severe_p = (p_inc_severe_0_36500/n_0_36500) * human_population,
         deaths = 0.215 * (n_inc_severe_0_36500/n_0_36500) * human_population,
         deaths_p = 0.215 * (p_inc_severe_0_36500/n_0_36500) * human_population) %>%
  rowwise() %>%
  mutate(dose1 = sum(n_rtss_epi_dose_1,n_rtss_mass_dose_1,na.rm=T),
         dose2 = sum(n_rtss_epi_dose_2,n_rtss_mass_dose_2,na.rm=T),
         dose3 = sum(n_rtss_epi_dose_3,n_rtss_mass_dose_3,na.rm=T),
         dose4 = sum(n_rtss_epi_booster_1,n_rtss_mass_booster_1,na.rm=T),
         dose5 = sum(n_rtss_mass_booster_2,na.rm=T),
         dosecomplete = dose3) %>%
  ungroup()

saveRDS(dat2,"C:/Users/htopazia/OneDrive - Imperial College London/Github/rtss_malariasimulation/rds/rtss_smc_all.rds")

none <- dat2 %>% filter(intervention == 'none') %>%
  rename(base_case = cases,
         base_n = n_0_36500,
         base_death = deaths) %>%
  select(timestep, base_case, base_n, base_death, eir, season)

dat3 <- dat2 %>% filter(intervention != 'none') %>% left_join(none) %>%
  group_by(eir, season, intervention) %>%
  summarize(base_case = sum(base_case, na.rm=T),
            base_death = sum(base_death, na.rm=T),
            dosecomplete = sum(dosecomplete, na.rm=T),
            cases = sum(cases, na.rm=T),
            pop = mean(n_0_36500, na.rm=T),
            severe = sum(severe, na.rm=T),
            deaths = sum(deaths, na.rm=T),
            cases_averted = base_case - cases,
            cases_averted_per_100000_fvp = (base_case - cases)/dosecomplete * human_population,
            deaths_averted = base_death - deaths,
            deaths_averted_per_100000_fvp = (base_death - deaths)/dosecomplete * human_population)

saveRDS(dat3,"C:/Users/htopazia/OneDrive - Imperial College London/Github/rtss_malariasimulation/rds/rtss_smc_averted.rds")
