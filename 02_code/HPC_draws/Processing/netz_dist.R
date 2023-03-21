# translate ITN usage to ITN annual distribution

# load netz package data
library(netz)
library(tidyverse)

# input: target ITN useage and observed use rate
# process: loads
# output: data frame with the median / min / max distribution rates required for each target ITN usage

netz_dist <- function(x  # observed use rate
                  ) {

  output <- expand_grid(target_use = c(0, 0.10, 0.20, 0.25, 0.30, 0.35, 0.40, 0.50, 0.60, 0.70)) |>
    mutate(half_life = median(get_halflife_data()$half_life),
           usage_rate = x) |>
    mutate(distribution_freq = 365 * 3) |>
    mutate(access = usage_to_access(target_use, usage_rate),
           npc = access_to_crop(access, type = "loess"),
           annual_percapita_nets_distributed = crop_to_distribution(npc, distribution_freq = distribution_freq,
                                        net_loss_function = net_loss_map,
                                        half_life = half_life)) |>
    select(target_use, annual_percapita_nets_distributed)

  # Assumptions:
  # Using median half life
  # Using minimum use rate (88%) allowing to give 85% usage (this is between median and max)
  # Extrapolating Loess curve according to curve trend
  # Assuming net loss is like in MAP paper (smooth compact)

  # rename var for min SSA use rate
  if(x == min(get_usage_rate_data()$usage_rate)){
    output <- rename(output, annual_percapita_nets_distmin = annual_percapita_nets_distributed)
  }

  # rename var for max SSA use rate
  if(x == max(get_usage_rate_data()$usage_rate)){
    output <- rename(output, annual_percapita_nets_distmax = annual_percapita_nets_distributed)
  }

  return(output)

}

target_use <- c(0.90, # assume maximum observed use rate and median bednet half life (across Africa)
                min(get_usage_rate_data()$usage_rate), # assume observed rate is the min in Africa
                max(get_usage_rate_data()$usage_rate)) # assume observed rate is the max in Africa

# run function
output <- map_dfr(target_use, netz_dist) |>
  group_by(target_use) |>
  mutate(across(1:3, mean, na.rm = T)) |>
  distinct()

output[1, 2:4] <- 0 # set NAs in first row with 0 usage to 0 nets distributed
output[10, 3] <- 0.618 # set last value of min observed rate to increase linearly re: target usage error

# save
saveRDS(output, "./03_output/netz_data.rds")

# NOTE that values are NA for the target rates 0.75, 0.85 with the min SSA rate #

