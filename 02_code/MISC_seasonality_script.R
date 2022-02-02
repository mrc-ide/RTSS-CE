# Set-up -----------------------------------------------------------------------
library(umbrella) # remotes::install_github('mrc-ide/umbrella@gee', force=T)
library(rgee) # remotes::install_github('r-spatial/rgee', force=T)
library(sf)
library(dplyr)
library(ggplot2)
library(reticulate)

setwd('M:/Hillary/rtss_malariasimulation')
# if earthengine is already installed, skip to line 34, ee_Initialize()

# Set python version (needs to be <=3.9!) and install earthengine --------------
py_discover_config() # check version
import('sys')$executable
reticulate::use_python('C:/Users/htopazia/Anaconda3/python.exe', required=T)
Sys.setenv(RETICULATE_PYTHON="C:/Users/htopazia/Anaconda3/envs/rgee/python.exe")
options(reticulate.conda_binary = "C:/Users/htopazia/Anaconda3/envs/rgee/python.exe")

# Use Anaconda prompt to create a new environment
# conda create --name rgee python=3.9
# conda activate rgee
# conda install -c conda-forge earthengine-api

# install earthengine
rgee::ee_install_set_pyenv(py_path = "C:/Users/htopazia/Anaconda3/envs/rgee/python.exe",
                           py_env = 'rgee')
# ::ee_clean_pyenv() # to clear and start over

# check
ee_check() # may encounter an error but keep going
# usethis::edit_r_environ() # to check environment

# sign into google earth account
ee_Initialize(user='htopazian@gmail.com', drive=T) # Users/htopazian_earthengine


# Start analysis ---------------------------------------------------------------
# import spatial data
admin0 <- readRDS("M:/Eradication/rds/admin0.RDS") %>% filter(Country %in% c('Mali', 'Burkina Faso', 'Senegal', 'Guinea', 'Mauritania', 'Algeria', 'Niger', "Côte d'Ivoire", "CÃ´te d'Ivoire", 'Ghana', 'Togo', 'Benin', 'Nigeria', 'Sierra Leone', 'Morocco'))
admin1 <- readRDS("M:/Eradication/rds/admin1.RDS") %>% filter(admin1 %in% c('Hauts-Bassins', 'Sikasso'))
countries <- filter(admin0, Country %in% c('Mali', 'Burkina Faso'))
admin2 <- readRDS("M:/Eradication/rds/admin2.RDS") %>%
  filter(admin2 %in% c('Tuy', 'Bougouni'))

# plot
ggplot() +
  geom_sf(data = admin0, fill="cornsilk2",color="cornsilk3") +
  geom_sf(data = countries, fill="cornsilk", color="tan4", size=0.5) +
  geom_sf(data = admin1, aes(fill = Country)) +
  geom_sf(data = admin2, color='black', fill=NA) +
  labs(fill = 'country sites', x = "", y = "") +
  theme_bw() +
  scale_x_continuous(limits=c(-12, 4)) +
  scale_y_continuous(limits=c(9.5, 25)) +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank())

# specify time-range
start_date <- "2015-01-01"
end_date <- "2016-12-31"

# extract
daily_rain <- pull_daily_rainfall(sf = admin2, start_date = start_date, end_date = end_date)

# process
seasonality <- process(gee_output = daily_rain, ISO, Country, admin1, admin2)

# inspect
seasonality %>%
  select(admin1, coefficients) %>%
  tidyr::unnest(cols = c(coefficients))

pd <- seasonality %>%
  select(admin1, profile) %>%
  tidyr::unnest(cols = c(profile)) %>%
  mutate(day = t * 365)

ggplot(pd, aes(x = day, y = y, col = admin1)) +
  geom_line() +
  xlab("Day") +
  ylab("Precipitation") +
  theme_bw()
