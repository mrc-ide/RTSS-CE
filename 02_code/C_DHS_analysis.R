# country case-study -----------------------------------------------------------
library(rdhs)
library(malariaAtlas)
library(raster)
library(sf)
library(tidyverse)

# MAP data ---------------------------------------------------------------------
PR <- read_csv('./01_data/00_PfPR_table_Global_admin0_2000-2019.csv') %>%
  filter(Year == 2019) %>% arrange(PfPR_median) %>%
  mutate(SSA = case_when(Name_0 %in% nets_data$Name ~ 1,
                         grepl('Cape|Comor|Ivoire|Lesotho|Maurit|Tome|Seyc|Tanz', Name_0) ~ 1)) %>%
  filter(SSA == 1) %>%
  mutate(Name_f = factor(Name_0, levels=Name_0))


# get MAP shapefiles for SSA admin1
SSA_shp0 <- malariaAtlas::getShp(ISO = PR$ISO, admin_level = "admin0")
SSA_shp1 <- malariaAtlas::getShp(ISO = PR$ISO, admin_level = "admin1")

# get MAP raster for ITN use
SSA_ITN_2019 <- malariaAtlas::getRaster(surface = 'Insecticide treated bednet (ITN) use version 2020',
                                        shp = SSA_shp0,
                                        year = 2019)


v0 <- raster::extract(SSA_ITN_2019, SSA_shp0, weights=TRUE, fun=mean, na.rm=T)
v1 <- raster::extract(SSA_ITN_2019, SSA_shp1, weights=TRUE, fun=mean, na.rm=T)


# turn spatial object to dataframe and add raster values
SSA_ITN0 <- sf::st_drop_geometry(sf::st_as_sf(SSA_shp0))
SSA_ITN0$mean <- v0

SSA_ITN1 <- sf::st_drop_geometry(sf::st_as_sf(SSA_shp1))
SSA_ITN1$mean <- v1

# summarize data
SSA_ITN0 <- SSA_ITN0 %>% group_by(name_0) %>%
  summarize(n = n(), median = median(mean, na.rm=T))

SSA_ITN1 <- SSA_ITN1 %>% group_by(name_0) %>%
  summarize(min = min(mean, na.rm = T), max = max(mean, na.rm=T))

SSA_ITN <- SSA_ITN0 %>% full_join(SSA_ITN1)

saveRDS(SSA_ITN, './03_output/MAP_ITN.rds')




# DHS data ---------------------------------------------------------------------
set_rdhs_config(email = "h.topazian@imperial.ac.uk",
                project = "Strategic introduction of the RTS,S/AS01 malaria
                vaccine relative to scale-up of existing interventions")

# select country and year
survs <- dhs_surveys(countryIds = c("NG"), # c("GH","KE", "MW", "NG", "TZ", "CD"),
                     surveyYear = c(2018))

# select children's recode and individual recode
datasets_kr <- dhs_datasets(
  surveyIds = survs$SurveyId,
  fileFormat = "FL",
  fileType = "KR")

datasets_pr <- dhs_datasets(
  surveyIds = survs$SurveyId,
  fileFormat = "FL",
  fileType = "PR")

# download datsets and variables
downloads_kr <- get_datasets(datasets_kr$FileName)
vars_kr <- c("caseid", "v001", "v002", "v003", "midx", "v005", "b5",
             "b8", "v459", "v008", "h1", "h7", "h9", 'h10', "h22", "h32z")
questions_kr <- search_variables(datasets_kr$FileName, variables = vars_kr)

downloads_pr <- get_datasets(datasets_pr$FileName)
vars_pr <- c("hhid", "hv001", "hv002", "hvidx", "hv005", "hv103", "hml12", "hc1", "hml35")
questions_pr <- search_variables(datasets_pr$FileName, variables = vars_pr)


# print questions (note: assumes questions/codes are the same across all surveys, but this is not always the case)
print(questions_kr[1:length(vars_kr), 1:3])
print(questions_pr[1:length(vars_pr), 1:3])

# extract data
extract_kr <- extract_dhs(questions_kr, add_geo = TRUE)
extract_bound_kr <-
  rbind_labelled(
    extract_kr$NGKR7BFL
  )

extract_pr <- extract_dhs(questions_pr, add_geo = TRUE)
extract_bound_pr <-
  rbind_labelled(
    extract_pr$NGPR7BFL
  )

# save
saveRDS(extract_bound_kr,  file = "./03_output/DHS_nets_dtp3_KR.RDS")
saveRDS(extract_bound_pr,  file = "./03_output/DHS_nets_dtp3_PR.RDS")



# read in
extract_bound_kr <- readRDS("./03_output/DHS_nets_dtp3_KR.RDS")
extract_bound_pf <- readRDS("./03_output/DHS_nets_dtp3_PR.RDS")


# First filter to children aged 1 and 2 (and alive) in Children's recode, then link person's recode. Then keep those who slept in house last night.

dat_start <- extract_bound_kr %>%
  filter(b5 == 1, b8 %in% c(1,2)) %>%
  left_join(extract_bound_pr, by = c("v001" = "hv001", "v002" = "hv002", "v003" = "hvidx",
                                     "CLUSTER", "ALT_DEM", "LATNUM", "LONGNUM", "ADM1NAME",
                                     "DHSREGNA", "SurveyId")) %>%
  filter(hv103 == 1, hml12 != 9, v459 != 9)


# select either "admin1", "region", or "country". Only need country

level <- 'admin1'

if (level == "admin1"){
  grouping <- c("SurveyId", "DHSREGNA", "ADM1NAME")
} else if (level == "region") {
  grouping <- c("SurveyId", "DHSREGNA")
} else if (level == "country") {
  grouping <- c("SurveyId")
} else {
  stop("Invalid level specified")
}

# Summarise at selected level
dat <- dat_start %>%
  group_by_at(grouping) %>%
  summarise(num_children = sum(v005/1e6),
            date_vacc_survey_cmc = mean(v008))

dat_DTP_card <- dat_start %>%
  group_by_at(grouping) %>%
  filter((h7 %in% c(1, 3)) | (h1 == 1 & h7 == 2)) %>%
  summarise(num_DTP3_card = sum(v005 / 1e6))

dat_DTP_report <- dat_start %>%
  group_by_at(grouping) %>%
  filter(h7 == 2 & h1 != 1) %>%
  summarise(num_DTP3_report = sum(v005 / 1e6))

dat_measles_card <- dat_start %>%
  group_by_at(grouping) %>%
  filter((h9 %in% c(1, 3)) | (h1 == 1 & h9 == 2)) %>%
  summarise(num_measles_card = sum(v005 / 1e6))

dat_measles_report <- dat_start %>%
  group_by_at(grouping) %>%
  filter(h9 == 2 & h1 != 1) %>%
  summarise(num_measles_report = sum(v005 / 1e6))

dat_any_net <- dat_start %>%
  group_by_at(grouping) %>%
  # check if have net and if child slept under any net
  filter(v459 == 1, hml12 %in% c(1,2,3)) %>%
  summarise(slept_any_net = sum(v005/1e6))

dat_itns <- dat_start %>%
  group_by_at(grouping) %>%
  # check if have net and if child slept under treated net
  filter(v459 == 1, hml12 %in% c(1,2)) %>%
  summarise(slept_itn = sum(v005/1e6))

dat_access <- dat_start %>%
  group_by_at(grouping) %>%
  # remove missing for having a net
  filter(v459 %in% c(1)) %>%
  summarise(access_itn = sum(v005/1e6))

dat_vaccine_yes_itn_yes  <- dat_start %>%
  group_by_at(grouping) %>%
  # filter yes for vaccine and yes for nets
  filter(h7 %in% c(1, 2, 3) & hml12 %in% c(1,2)) %>%
  summarise(vac_y_itn_y = sum(v005/1e6))

dat_vaccine_no_itn_yes  <- dat_start %>%
  group_by_at(grouping) %>%
  # filter no for vaccine and yes for nets
  filter(h7 == 0 & hml12 %in% c(1,2)) %>%
  summarise(vac_n_itn_y = sum(v005/1e6))

dat_vaccine_yes_itn_no <- dat_start %>%
  group_by_at(grouping) %>%
  # filter yes for vaccine and no for nets
  filter(h7 %in% c(1, 2, 3) & hml12 %in% c(0,3)) %>%
  summarise(vac_y_itn_n = sum(v005/1e6))

dat_vaccine_no_itn_no <- dat_start %>%
  group_by_at(grouping) %>%
  # filter no for vaccine and no for net
  filter(h7 == 0 & hml12 %in% c(0,3)) %>%
  summarise(vac_n_itn_n = sum(v005/1e6))

dat_all <- left_join(dat, dat_DTP_card, by = grouping) %>%
  left_join(dat_DTP_report, by = grouping) %>%
  left_join(dat_measles_card, by = grouping) %>%
  left_join(dat_measles_report, by = grouping) %>%
  left_join(dat_any_net, by = grouping) %>%
  left_join(dat_itns, by = grouping) %>%
  left_join(dat_access, by = grouping) %>%
  left_join(dat_vaccine_yes_itn_yes, by = grouping) %>%
  left_join(dat_vaccine_no_itn_yes, by = grouping) %>%
  left_join(dat_vaccine_yes_itn_no, by = grouping) %>%
  left_join(dat_vaccine_no_itn_no, by = grouping) %>%
  mutate(prop_DTP3_vacc = (num_DTP3_card + num_DTP3_report)/num_children,
         prop_measles_vacc = (num_measles_card + num_measles_report)/num_children,
         prop_any_net = slept_any_net/num_children,
         prop_itn = slept_itn/num_children,
         prop_access = access_itn/num_children,
         prop_vac_y_int_y = vac_y_itn_y/num_children,
         prop_vac_n_int_y = vac_n_itn_y/num_children,
         prop_vac_y_int_n = vac_y_itn_n/num_children,
         prop_vac_n_int_n = vac_n_itn_n/num_children) %>%
  dplyr::select(-c(num_children, num_DTP3_card, num_DTP3_report, num_measles_card, num_measles_report, num_children, slept_any_net, slept_itn, num_children, vac_y_itn_y, vac_n_itn_y, vac_y_itn_n, vac_n_itn_n))

# add the countrycode
dat_all$countrycode <- substr(dat_all$SurveyId, 1, 2)

# remove SurveyId
dat_all <- dat_all %>% dplyr::select(-c(SurveyId))

dat_all <- dat_all %>%
  dplyr::select(c(countrycode, date_net_survey_cmc, date_treat_survey_cmc,
                  date_vacc_survey_cmc, everything()))

write_csv(dat_all, paste0("./03_output/combined_DHS_data_", level, ".csv"))

