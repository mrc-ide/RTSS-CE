# Figures & Tables -------------------------------------------------------------

# packages
library(tidyverse)
library(data.table)
library(patchwork)
library(grid)
library(LaCroixColoR) # devtools::install_github("johannesbjork/LaCroixColoR")

# color palette
# lacroix_palette(type = "paired")

# assign scenarios to particular colors
smc          <- "#FC6882"
pbo          <- "#007BC3"
itn          <- "#54BCD1"
rtss_sv      <- "#009F3F"
rtss_age     <- "#54E356"
pbo_smc      <- "#C70E7B"
itn_smc      <- "#B25D91"
rtss_smc     <- "#EFC7E6"
pbo_rtss_smc <- "#AF6125"
itn_rtss_smc <- "#F5DD42"
pbo_rtss     <- "#EF7C12"
itn_rtss     <- "#F4B95A"

# load data
scenarios <- readRDS("./03_output/scenarios_draws.rds") |>
  # remove scenarios with a negative impact on DALYs
  mutate(dalydiff = daly_baseline - daly) |>
  filter(dalydiff >= 0)


# the data read in above is the main dataset for the paper
# it captures uncertainty from parameter draw runs


