# country case-study -----------------------------------------------------------
library(rdhs)
library(malariaAtlas)
# devtools::install_github("ropensci/rdhs")
library(raster)
library(sf)
library(tidyverse)
library(netz)
library(survey)
library(broom)


# MAP data ---------------------------------------------------------------------
# ITNs #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
# get MAP data from ITNs and from PfPR tables
nets_data <- read_csv('./01_data/Intervention_ITN.csv') %>% arrange(`2015`) %>%
  mutate(Name_f = factor(Name, levels=Name))

PR <- read_csv('./01_data/00_PfPR_table_Global_admin0_2000-2019.csv') %>%
  filter(Year == 2019) %>% arrange(PfPR_median) %>%
  # filter out only countries in SSA
  mutate(SSA = case_when(Name_0 %in% nets_data$Name ~ 1,
                         grepl('Cape|Comor|Ivoire|Lesotho|Maurit|Tome|Seyc|Tanz', Name_0) ~ 1)) %>% # add SSA countries missing ITN info
  filter(SSA == 1) %>%
  mutate(Name_f = factor(Name_0, levels=Name_0))


# get MAP shapefiles for SSA admin1
SSA_shp0 <- malariaAtlas::getShp(ISO = PR$ISO, admin_level = "admin0")
SSA_shp1 <- malariaAtlas::getShp(ISO = PR$ISO, admin_level = "admin1")

# get MAP raster for ITN use
SSA_ITN_2019 <- malariaAtlas::getRaster(surface = 'Insecticide treated bednet (ITN) use version 2020',
                                        shp = SSA_shp0,
                                        year = 2019)

# extract ITN information from rasters at admin0 and admin1 levels
v0 <- raster::extract(SSA_ITN_2019, SSA_shp0, weights=TRUE, fun=mean, na.rm=T)
v1 <- raster::extract(SSA_ITN_2019, SSA_shp1, weights=TRUE, fun=mean, na.rm=T)


# turn spatial object into a dataframe and add raster values
SSA_ITN0 <- sf::st_drop_geometry(sf::st_as_sf(SSA_shp0))
SSA_ITN0$mean <- v0

SSA_ITN1 <- sf::st_drop_geometry(sf::st_as_sf(SSA_shp1))
SSA_ITN1$mean <- v1

# summarize data
# find country medians and min and max admin1 values
SSA_ITN0 <- SSA_ITN0 %>% group_by(name_0) %>%
  summarize(n = n(),
            median = median(mean, na.rm=T))

SSA_ITN1 <- SSA_ITN1 %>% group_by(name_0) %>%
  summarize(min = min(mean, na.rm = T),
            max = max(mean, na.rm=T))

# combine country values and min / max admin values into one dataframe
SSA_ITN <- SSA_ITN0 %>% full_join(SSA_ITN1)

# save for ITN plots
saveRDS(SSA_ITN, './03_output/MAP_ITN.rds')


# TRAVEL TIMES #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
# download shapefile and travel times raster from MAP
shp <- malariaAtlas::getShp(ISO = c("NGA","GHA"), admin_level = "admin1")
plot(shp)

traveltime <- malariaAtlas::getRaster(
  surface = 'Global travel time to healthcare map with access to motorized transport',
  shp = shp)
malariaAtlas::autoplot_MAPraster(traveltime)


# PfPR #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
# get MAP PfPR for SSA admin1
shpNGA <- malariaAtlas::getShp(ISO = c("NGA"), admin_level = "admin1")

shpGHA <- malariaAtlas::getShp(ISO = c("GHA"), admin_level = "admin1")

PfPR_NGA <- malariaAtlas::getRaster(
  surface = 'Plasmodium falciparum PR2 - 10 version 2020',
  year = 2018,
  shp = shpNGA)
malariaAtlas::autoplot_MAPraster(PfPR_NGA)

PfPR_GHA <- malariaAtlas::getRaster(
  surface = 'Plasmodium falciparum PR2 - 10 version 2020',
  year = 2014,
  shp = shpGHA)
malariaAtlas::autoplot_MAPraster(PfPR_GHA)

PfPR2 <- do.call(merge, list(PfPR_NGA, PfPR_GHA))
plot(PfPR2)


# DHS data ---------------------------------------------------------------------
set_rdhs_config(email = "h.topazian@imperial.ac.uk",
                project = "Strategic introduction of the RTS,S/AS01 malaria
                vaccine relative to scale-up of existing interventions")
# select option 2 to not allow rdhs to store outside of temp directory

# select country and year
survs <- dhs_surveys(surveyCharacteristicIds = 39, # maalria parasitemia
                     surveyType = c('DHS'),
                     surveyYear = seq(2000,2020,1))
                     # countryIds = c("NG","GH"), # c("GH","KE", "MW", "NG", "TZ", "CD"),
                     # surveyYear = c(2018, 2014))

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
vars_kr <- c("caseid", "v001", "v002", "v003", "midx", "v005", 'b4', "b5", 'b16',
             "b8", "v459", "v008", "v190", "h1", "h7", "h9", 'h10', "h22", "h32z")
questions_kr <- search_variables(datasets_kr$FileName, variables = vars_kr)

downloads_pr <- get_datasets(datasets_pr$FileName)
vars_pr <- c("hhid", "hv001", "hv002", "hvidx", "hv005", "hv023", 'hv024', "hv025", "hv103", "hml12", "hc1", "hml32", "hml35", 'sh212i')
questions_pr <- search_variables(datasets_pr$FileName, variables = vars_pr)


# print questions (note: assumes questions/codes are the same across all surveys, but this is not always the case)
print(questions_kr[1:length(vars_kr), 1:3])
print(questions_pr[1:length(vars_pr), 1:3])

# extract data
extract_kr <- extract_dhs(questions_kr, add_geo = TRUE)
extract_bound_kr <-
  full_join(
    extract_kr$NGKR7BFL,
    extract_kr$GHKR72FL,
    by=c("caseid", "v001", "v002", "v003", "midx", "v005", "b4", "b5", "b8", 'b16', "v459", "v008", "v190", "h1", "h7", "h9", "h10", "h22", "h32z", "CLUSTER", "ALT_DEM", "LATNUM", "LONGNUM", "ADM1NAME", "DHSREGNA", "SurveyId")
  )

extract_pr <- extract_dhs(questions_pr, add_geo = TRUE)
extract_bound_pr <-
  full_join(
    extract_pr$NGPR7BFL,
    extract_pr$GHPR72FL,
    by = c("hhid", "hv001", "hv002", "hvidx", "hv005", "hv023", "hv024", "hv025", "hv103", "hml12", "hc1", "hml32", "hml35", "CLUSTER", "ALT_DEM", "LATNUM", "LONGNUM", "ADM1NAME", "DHSREGNA", "SurveyId")
  )

# save
saveRDS(extract_bound_kr,  file = "./03_output/DHS_nets_dtp3_KR.rds")
saveRDS(extract_bound_pr,  file = "./03_output/DHS_nets_dtp3_PR.rds")

# read in
extract_bound_kr <- readRDS("./03_output/DHS_nets_dtp3_KR.rds")
extract_bound_pr <- readRDS("./03_output/DHS_nets_dtp3_PR.rds")


# first filter to children aged 1 and 2 (and alive) in children's recode, then link person's recode. Then keep those who slept in house last night.
# link KR and PR datasets by v001 and v002 per: https://dhsprogram.com/data/Merging-datasets.cfm
dat_start <- extract_bound_kr %>%
  filter(b5 == 1) %>% # filtering to alive children
  left_join(extract_bound_pr, by = c("v001" = "hv001", "v002" = "hv002", "b16" = "hvidx",
                                     "CLUSTER", "ALT_DEM", "LATNUM", "LONGNUM", "ADM1NAME", "SurveyId",
                                     "DHSREGNA")) %>%
  mutate(wt = v005/1e6)

# summarize data by admin1
# select either "admin1", "region", or "country"
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
  # sum women's individual sample weight and take the mean interview date
  summarise(num_children = sum(v005 / 1e6),
            date_vacc_survey_cmc = mean(v008))

dat_DTP_card <- dat_start %>%
  group_by_at(grouping) %>%
  # received dpt3 (on vax card or reported by mother and health card is seen)
  filter((h7 %in% c(1, 3)) | (h1 == 1 & h7 == 2)) %>%
  summarise(num_DTP3_card = sum(v005 / 1e6))

dat_DTP_report <- dat_start %>%
  group_by_at(grouping) %>%
  # received dpt3 (reported by mother but health card not seen)
  filter(h7 == 2 & h1 != 1) %>%
  summarise(num_DTP3_report = sum(v005 / 1e6))

dat_measles_card <- dat_start %>%
  group_by_at(grouping) %>%
  # received measles1 (on vax card or reported by mother and health card is seen)
  filter((h9 %in% c(1, 3)) | (h1 == 1 & h9 == 2)) %>%
  summarise(num_measles_card = sum(v005 / 1e6))

dat_measles_report <- dat_start %>%
  group_by_at(grouping) %>%
  # received measles1 (reported by mother but health card not seen)
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

# join all vars into single dataset
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

# add coordinates
dat_all <- rbind(
  (shpNGA %>% st_as_sf() %>%
   dplyr::select(name_1) %>%
   mutate(name_1 = tolower(name_1))),
  (shpGHA %>% st_as_sf() %>%
   dplyr::select(name_1) %>%
   mutate(name_1 = tolower(name_1)))) %>%
  left_join(dat_all %>% mutate(ADM1NAME = case_when(ADM1NAME=='FCT ABUJA'~'abuja',
                                                    ADM1NAME=='NASARAWA'~'nassarawa',
                                                    TRUE ~ ADM1NAME),
                               ADM1NAME = tolower(ADM1NAME)), by = c('name_1'='ADM1NAME'))

# save to .rds
saveRDS(dat_all, paste0("./03_output/combined_DHS_data_", level, ".rds"))

# add travel times and PfPR to individual level data
# transforming to UTM 31N, to allow for creation of buffers in meters
dhs_map <- st_as_sf(dat_start, coords = c("LONGNUM", "LATNUM"), crs = 4326) # all DHS coords in WGS 1984
dhs_map <- st_transform(dhs_map %>% dplyr::select(SurveyId, CLUSTER) %>% dplyr::distinct(.keep_all=T), crs = 32631)
plot(dhs_map)

# extract travel times around cluster points, with a 1km buffer
travel <- raster::extract(traveltime, dhs_map, buffer=1000, fun=mean)

# extract PfPR around cluster points, with a 1km buffer
PfPRbuff <- raster::extract(PfPR2, dhs_map, buffer=1000, fun=mean)

dhs_map <- cbind(dhs_map, travel, PfPRbuff)
plot(dhs_map$travel); plot(dhs_map$PfPRbuff)

ggplot(dhs_map %>% filter(!is.na(travel))) +
  geom_sf(aes(color=travel)) +
  scale_color_distiller(palette = 'Spectral')

ggplot(dhs_map %>% filter(!is.na(PfPRbuff))) +
  geom_sf(aes(color=PfPRbuff)) +
  scale_color_distiller(palette = 'Spectral')


# inspect individual level DHS data
# look for associations with DPT3, ITN, travel time
# make variables for DPT3 and ITN status
dat_start2 <- dat_start %>%
  mutate(dpt3 = case_when(h7 %in% c(1, 2, 3) ~ 1,
                          h7 == 0 ~ 0,
                          TRUE ~ NA_real_),
         itn = case_when(hml12 %in% c(1,2) ~ 1,
                         hml12 %in% c(0,3) ~ 0,
                         TRUE ~ NA_real_),
         itn_access = case_when(v459 == 1 ~ 1,
                                v459 == 0 ~ 0,
                                TRUE ~ NA_real_),
         vac_y_itn_n = case_when(h7 %in% c(1, 2, 3) & hml12 %in% c(0,3) ~ 1,
                                 is.na(h7) | is.na(hml12) ~ NA_real_,
                                 TRUE ~ 0),
         microscopy = case_when(hml32 == 0 ~ 0,
                                hml32 == 1 ~ 1,
                                TRUE ~ NA_real_),
         mRDT = case_when(hml35 == 0 ~ 0,
                          hml35 == 1 ~ 1,
                          TRUE ~ NA_real_),
         ACT = case_when(sh212i == 0 ~ 0,
                         sh212i == 1 ~ 1,
                         TRUE ~ NA_real_),
         treatment= case_when(h32z == 0 ~ 0,
                              h32z == 1 ~ 1,
                              TRUE ~ NA_real_),
         wealth = case_when(v190 == 1 ~ 'poorest',
                            v190 == 2 ~ 'poorer',
                            v190 == 3 ~ 'middle',
                            v190 == 4 ~ 'richer',
                            v190 == 5 ~ 'richest'),
         residence = case_when(hv025 ==1 ~ 'urban',
                               hv025 == 2 ~ 'rural'),
         sex = case_when(b4 == 1 ~ 'male',
                         b4 == 2 ~ 'female')) %>%
  left_join(dhs_map, by = c('SurveyId', 'CLUSTER'))

# save to .rds
saveRDS(dat_start2, "./03_output/individual_DHS_data.rds")


# models -----------------------------------------------------------------------
# choose 1 option for NGA or GHA
dat_start <- readRDS("./03_output/individual_DHS_data.rds") %>%
  filter(SurveyId == 'GH2014DHS') %>% # NG2018DHS
  mutate(travel60 = travel/60) %>%
  filter(!is.na(PfPRbuff))


# plot clusters to ensure correct admin1 or country was selected
ggplot() +
  geom_sf(data = st_as_sf(shp)) +
  geom_sf(data = st_as_sf(dat_start))

# create survey design object
DHSdesign <- survey::svydesign(id = ~CLUSTER, strata = ~hv023, weights = ~wt, data = dat_start %>% filter(!is.na(hv023)))

# group data by cluster
plotd <- dat_start %>% group_by(CLUSTER) %>%
  summarize(dpt3 = mean(dpt3, na.rm=T),
            itn = mean(itn, na.rm=T),
            wealth = mean(v190, na.rm=T)) %>%
  left_join(dhs_map, by='CLUSTER')

ggplot(data = plotd, aes(x=travel, y=dpt3)) +
  geom_point() +
  geom_smooth(method='lm') +
  coord_cartesian(ylim = c(0,1)) + theme_classic()

ggplot(data = plotd, aes(x=travel, y=itn)) +
  geom_point() +
  geom_smooth(method='lm') +
  coord_cartesian(ylim = c(0,1)) + theme_classic()

ggplot(data = plotd, aes(x=travel, y=wealth)) +
  geom_point() +
  geom_smooth(method='lm') +
  coord_cartesian(ylim = c(1,5)) + theme_classic()

ggplot(data = plotd, aes(x=PfPRbuff, y=wealth)) +
  geom_point() +
  geom_smooth(method='lm') +
  coord_cartesian(ylim = c(1,5)) + theme_classic()

ggplot(data = plotd, aes(x=PfPRbuff, y=itn)) +
  geom_point() +
  geom_smooth(method='lm') +
  theme_classic()

ggplot(data = plotd, aes(x=PfPRbuff, y=travel)) +
  geom_point() +
  geom_smooth(method='lm') +
  theme_classic()

m <- svyglm(dpt3 ~ travel60, DHSdesign, family=binomial("logit"))
tidy(m, conf.int = T, conf.level = 0.95, exponentiate = T)

m <- svyglm(itn ~ travel60 + PfPRbuff, DHSdesign, family=binomial("logit"))
tidy(m, conf.int = T, conf.level = 0.95, exponentiate = T)

m <- svyglm(itn_access ~ b4, DHSdesign, family=binomial("logit"))
tidy(m, conf.int = T, conf.level = 0.95, exponentiate = T)

m <- svyglm(dpt3 ~ b4, DHSdesign, family=binomial("logit"))
tidy(m, conf.int = T, conf.level = 0.95, exponentiate = T)

m <- svyglm(itn_access ~ v190, DHSdesign, family=binomial("logit"))
tidy(m, conf.int = T, conf.level = 0.95, exponentiate = T)

m <- svyglm(dpt3 ~ factor(v190), DHSdesign, family=binomial("logit"))
tidy(m, conf.int = T, conf.level = 0.95, exponentiate = T)


# create function for printing out weighted means
output_admin1 <- function(admin1) {

  dat <- dat_start %>% filter(!is.na(hv023)) %>% filter(ADM1NAME == admin1)

  DHSdesign <- survey::svydesign(id = ~CLUSTER, strata = ~hv023, weights = ~wt, data = dat)

  survmean_factor <- function(var1, var2){

    m <- svyby(as.formula(paste0('~', var1)), as.formula(paste0('~', var2)), DHSdesign, svymean, na.rm=T, survey.lonely.psu="adjust")
    m <- as.data.frame(m)
    rownames(m) <- c()
    m %>% rename(value = var1, factor = var2) %>% mutate(var = var1)

  }

  t <- map2_dfr(rep(c('itn_access','itn','dpt3','vac_y_itn_n', 'microscopy', 'mRDT', 'ACT', 'travel60', 'treatment'), 3),
         c(rep('wealth', 9), rep('sex', 9), rep('residence', 9)),
         survmean_factor) %>% mutate(admin1 = admin1)

  t

}

# run function
vars <- unique(dat_start$ADM1NAME)
d <- map_dfr(vars, output_admin1)

# copy estimates to clipboard by factor variable
d %>% select(-se) %>%
  filter(factor %in% c('urban', 'rural')) %>%
  pivot_wider(names_from = factor, values_from = value) %>%
  arrange(var, admin1, rural, urban) %>% write.table("clipboard",sep="\t")

d %>% select(-se) %>%
  filter(factor %in% c('male', 'female')) %>%
  pivot_wider(names_from = factor, values_from = value) %>%
  arrange(var, admin1, male, female) %>% write.table("clipboard",sep="\t")

d %>% select(-se) %>%
  filter(factor %in% c('poorest', 'poorer', 'middle', 'richer', 'richest')) %>%
  pivot_wider(names_from = factor, values_from = value) %>%
  arrange(var, admin1) %>% write.table("clipboard",sep="\t")


# find correlation between ITN use and DPT3 vaccination
dat <- dat_start %>% filter(!is.na(hv023)) %>% filter(ADM1NAME == 'Eastern')

# unweighted
cor(dat$itn, dat$dpt3, method = 'pearson', use = 'complete.obs')

# weighted
DHSdesign <- survey::svydesign(id = ~CLUSTER, strata = ~hv023, weights = ~wt,
                               data = dat %>%
                                 filter(!is.na(itn) & !is.na(dpt3) & residence == 'rural'))
jtools::svycor(~dpt3 + itn, design = DHSdesign)

DHSdesign <- survey::svydesign(id = ~CLUSTER, strata = ~hv023, weights = ~wt,
                               data = dat %>%
                                 filter(!is.na(itn) & !is.na(dpt3) & residence == 'urban'))
jtools::svycor(~dpt3 + itn, design = DHSdesign)
