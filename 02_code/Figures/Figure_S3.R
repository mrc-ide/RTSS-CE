# Figure S3
# Insecticide-treated net (ITN) distribution cost curves vs linear distribution cost curves

# set-up
source("./02_code/Figures/data_and_libraries.R")


# read in netz package data to find the annual nets to distribute to give the simulated usage
devtools::install_github('https://github.com/mrc-ide/netz/tree/f65bc686243dc0f6f210e8d459de8e11358f246f', force = T)
nets_data <- netz::prepare_data()

# get nets to be distributed for each ITN usage
ndist <- function(x) {

  convert_usage_to_annual_nets_distributed(
    target_usage = seq(0,0.85,0.05),
    distribution_freq = 1095, # 3 years
    use_rate_data = x,
    half_life_data = median(nets_data$half_life_data$half_life),
    extrapolate_npc = "loess",
    net_loss_function = net_loss_map) %>%
    select(target_use, annual_percapita_nets_distributed)

}

# assume maximum observed use rate and median bednet half life (across Africa)
nets_distributed <- ndist(0.88)

# assume observed rate is the min in Africa
nets_distributed_min <- ndist(min(nets_data$use_rate_by_country$use_rate))
nets_distributed_min <- rename(nets_distributed_min, annual_percapita_nets_distmin = annual_percapita_nets_distributed)

# assume observed rate is the max in Africa
nets_distributed_max <- ndist(max(nets_data$use_rate_by_country$use_rate))
nets_distributed_max <- rename(nets_distributed_max, annual_percapita_nets_distmax = annual_percapita_nets_distributed)

nets_distributed <- full_join(nets_distributed, nets_distributed_min) %>% full_join(nets_distributed_max)

net <- c('pyrethroid', 'pyrethroid + PBO')
output <- crossing(nets_distributed, net)

output <- output %>%
  mutate(ITNcost = case_when(net=='pyrethroid' ~ 3.50,            # $2.00 per net and $1.50 delivery cost
                             net=='pyrethroid + PBO' ~ 3.80)) %>% # $2.30 per net and $1.50 delivery cost
  mutate(cost_ITN = annual_percapita_nets_distributed * ITNcost * 3,  # true net cost accounting for non-linear relationship
         cost_ITNmin = annual_percapita_nets_distmin * ITNcost * 3,  # true net cost MIN
         cost_ITNmax = annual_percapita_nets_distmax * ITNcost * 3,  # true net cost MAX
         cost_ITNmin = ifelse(is.na(cost_ITNmin), cost_ITN, cost_ITNmin),
         cost_ITN_linear = target_use * ITNcost)         # ITN linear

output2 <- output %>%
  filter((net=='pyrethroid' & round(target_use,2) %in% c(0.25, 0.35, 0.50, 0.60, 0.75, 0.85)) |
           (net=='pyrethroid + PBO' & round(target_use,2) %in% c(0.25, 0.5, 0.75))) %>%
  mutate(shape='sim')

# plot
ggplot(output, aes(x=cost_ITN_linear, y = cost_ITN, color=target_use)) +
  geom_point(data=output2, aes(x=cost_ITN_linear, y = cost_ITN, shape=shape), size=4, color='chartreuse1', fill='chartreuse1') +
  geom_point(data=output2, aes(x=cost_ITN_linear, y = cost_ITN, shape=shape), size=3.5, color='chartreuse1', fill='chartreuse1')+
  geom_abline(slope=1, size=1, lty=2, alpha=0.3) +
  geom_point(size=2) +
  geom_line(alpha=0.5, size=1) +
  geom_ribbon(aes(ymin=cost_ITNmin, ymax=cost_ITNmax), color=NA, fill='cornflower blue', alpha=0.4) +
  scale_color_gradient(limits=c(0,0.85)) +
  scale_shape_manual(values=c(1), labels=c('value used in simulation')) +
  facet_grid(~net) +
  labs(x='Linear ITN cost', y='Netz ITN cost', color='target ITN use', shape='') +
  theme_bw() +
  theme(text = element_text(size = 14))

ggsave('./03_output/figureS3.pdf', width = 8, height = 4)

