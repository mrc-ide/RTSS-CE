# translate ITN usage to ITN annual distribution

# load netz package data
library(netz)
nets_data <- netz::prepare_data()


# input: target ITN useage and observed use rate
# process: loads
# output: data frame with the median / min / max distribution rates required for each target ITN usage

netz_dist <- function(x  # observed use rate
                  ) {

  output <- convert_usage_to_annual_nets_distributed(
    target_usage = c(0, 0.10, 0.25, 0.35, 0.50, 0.60, 0.75, 0.85),
    distribution_freq = 3 * 365, # every three years
    use_rate_data = x,
    half_life_data = median(nets_data$half_life_data$half_life),
    extrapolate_npc = "loess",
    net_loss_function = net_loss_map) %>%
    select(target_use, annual_percapita_nets_distributed)

  # Assumptions:
  # Using median half life
  # Using minimum use rate (88%) allowing to give 85% usage (this is between median and max)
  # Extrapolating Loess curve according to curve trend
  # Assuming net loss is like in MAP paper (smooth compact)

  # rename var for min SSA use rate
  if(x == min(nets_data$use_rate_by_country$use_rate)){
    output <- rename(output, annual_percapita_nets_distmin = annual_percapita_nets_distributed)
  }

  # rename var for max SSA use rate
  if(x == max(nets_data$use_rate_by_country$use_rate)){
    output <- rename(output, annual_percapita_nets_distmax = annual_percapita_nets_distributed)
  }

  return(output)

}

# set target_use
target_use <- c(0.88, # assume maximum observed use rate and median bednet half life (across Africa)
                min(nets_data$use_rate_by_country$use_rate), # assume observed rate is the min in Africa
                max(nets_data$use_rate_by_country$use_rate)) # assume observed rate is the max in Africa

# run function
output <- map_dfr(target_use, netz_dist) %>%
  group_by(target_use) %>%
  mutate(across(1:3, mean, na.rm = T)) %>%
  distinct()

# save
saveRDS(output, './03_output/netz_data')

# NOTE that values are NA for the target rates 0.75, 0.85 with the min SSA rate #

