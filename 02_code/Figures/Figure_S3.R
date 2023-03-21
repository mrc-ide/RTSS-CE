# Figure S3
# Insecticide-treated net (ITN) distribution cost curves vs linear distribution cost curves

# set-up
source("./02_code/Figures/data_and_libraries.R")


# read in netz package data to find the annual nets to distribute to give the simulated usage
library(netz)

# input: target ITN useage and observed use rate
# process: loads
# output: data frame with the median / min / max distribution rates required for each target ITN usage

netz_dist <- function(x  # observed use rate
) {

  output <- expand_grid(target_use = seq(0,0.85,0.05)) |>
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
nets_distributed <- map_dfr(target_use, netz_dist) |>
  group_by(target_use) |>
  mutate(across(1:3, mean, na.rm = T)) |>
  distinct()

nets_distributed[1, 2:4] <- 0 # set NAs in first row with 0 usage to 0 nets distributed

net <- c("pyrethroid", "pyrethroid + PBO")
output <- crossing(nets_distributed, net)

output <- output |>
  mutate(ITNcost = case_when(net == "pyrethroid" ~ 3.50,            # $2.00 per net and $1.50 delivery cost
                             net == "pyrethroid + PBO" ~ 3.80)) |> # $2.30 per net and $1.50 delivery cost
  mutate(cost_ITN = annual_percapita_nets_distributed * ITNcost * 3,  # true net cost accounting for non-linear relationship
         cost_ITNmin = annual_percapita_nets_distmin * ITNcost * 3,  # true net cost MIN
         cost_ITNmax = annual_percapita_nets_distmax * ITNcost * 3,  # true net cost MAX
         cost_ITNmin = ifelse(is.na(cost_ITNmin), cost_ITN, cost_ITNmin),
         cost_ITNmax = ifelse(is.na(cost_ITNmax), cost_ITN, cost_ITNmax),
         cost_ITN_linear = target_use * ITNcost)         # ITN linear

output2 <- output |>
  filter((net == "pyrethroid" & round(target_use,2) %in% c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7)) |
           (net == "pyrethroid + PBO" & round(target_use,2) %in% c(0.2, 0.4, 0.6))) |>
  mutate(shape="sim")

# plot
ggplot(output, aes(x = cost_ITN_linear, y = cost_ITN, color = target_use)) +
  geom_point(data = output2, aes(x = cost_ITN_linear, y = cost_ITN, shape = shape), size = 4, color = "chartreuse1", fill = "chartreuse1") +
  geom_point(data = output2, aes(x = cost_ITN_linear, y = cost_ITN, shape = shape), size = 3.5, color = "chartreuse1", fill = "chartreuse1")+
  geom_abline(slope = 1, size = 1, lty = 2, alpha = 0.3) +
  geom_point(size = 2) +
  geom_line(data = output |> filter(target_use > 0), alpha = 0.5, size = 1) +
  geom_ribbon(data = output |> filter(target_use > 0), aes(ymin = cost_ITNmin, ymax = cost_ITNmax), color = NA, fill = "cornflower blue", alpha = 0.4) +
  scale_color_gradient(limits = c(0, 0.85)) +
  scale_shape_manual(values = c(1), labels = c("value used in simulation")) +
  facet_grid(~ net) +
  labs(x = "Linear ITN cost", y = "Netz ITN cost", color = "target ITN use", shape = "") +
  scale_y_continuous(limits = c(0, 6)) +
  theme_bw() +
  theme(text = element_text(size = 14))

ggsave("./03_output/figureS3.pdf", width = 8, height = 4)

