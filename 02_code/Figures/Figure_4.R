# Figure 4
# Proportion of times an intervention was the most cost-efficient out of all scenarios

# set-up
source("./02_code/Figures/data_and_libraries.R")


# < by seasonality alone ------------------------------
univariateseason <- function(season) {

  # assign colors
  if(season == 'perennial'){
    colors <- c(itn, pbo, rtss_age)
  }

  if(season == 'highly seasonal'){
    colors <- c(itn, pbo, rtss_age, rtss_sv)
  }

  if(season == 'seasonal'){
    colors <- c(itn, pbo, smc)
  }

  # by pfpr
  output <- scenarios %>%
    filter(seasonality == season) %>%
    filter(cost_per_dose == 6.52 & delivery_cost == 1.62) %>%
    filter(intervention %in% c('ITN 10% increase','ITN PBO','SMC','RTS,S age-based','RTS,S seasonal')) %>%
    group_by(ID, drawID) %>% arrange(ID, CE) %>%
    slice(1L)

  A <- ggplot(output) +
    geom_bar(aes(x = as.factor(pfpr), fill = intervention_f), position = "fill", show.legend = F) +
    labs(x = 'PfPR', y = 'Proportion most \ncost-effective choice', fill = 'intervention') +
    scale_fill_manual(values =colors) +
    theme_classic()

  # by current ITN usage
  B <- ggplot(output) +
    geom_bar(aes(x = factor(ITNuse), fill = intervention_f), position = "fill", show.legend = F) +
    labs(x = 'ITN use', y = 'Proportion most \ncost-effective choice', fill='intervention') +
    scale_fill_manual(values = colors) +
    theme_classic()

  # by insecticide resistance
  C <- ggplot(output) +
    geom_bar(aes(x = factor(resistance), fill = intervention_f), position = "fill", show.legend = F) +
    labs(x = 'Resistance', y = 'Proportion most \ncost-effective choice', fill = 'intervention') +
    scale_fill_manual(values = colors) +
    theme_classic()

  # by ITN distribution efficiency
  ITNefficient <- function(var, label) {
    scenarios %>%
      filter(seasonality == season) %>%
      filter(cost_per_dose == 6.52 & delivery_cost == 1.62) %>%
      filter(intervention %in% c('ITN 10% increase','ITN PBO','SMC','RTS,S age-based','RTS,S seasonal')) %>%
      group_by(ID, drawID) %>% arrange(ID, drawID, {{var}}) %>%
      slice(1L) %>% select(intervention_f, seasonality, pfpr, {{var}}) %>%
      mutate(model = label) %>%
      rename(CE = {{var}})
  }

  output2 <- ITNefficient(CE, 'standard') %>%
    full_join(ITNefficient(CE_ITNmin, 'less efficient')) %>%
    full_join(ITNefficient(CE_ITNmax, 'more efficient'))


  if(season == 'seasonal'){
    colors <- c(itn, pbo, rtss_age, rtss_sv, smc)
  }

  D <- ggplot(output2) +
    geom_bar(aes(x = factor(model, levels = c('less efficient', 'standard', 'more efficient')), fill = intervention_f), position = "fill") +
    labs(x = 'ITN efficiency', y = 'Proportion most \ncost-effective choice', fill = 'intervention') +
    scale_fill_manual(values = colors) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
    theme_classic()


  if(season == 'seasonal'){
    colors <- c(itn, pbo, smc)
  }

  E <- ggplot(output) +
    geom_bar(aes(x = factor(treatment, labels = c('low', 'medium', 'high')), fill = intervention_f), position = "fill", show.legend = F) +
    labs(x = 'Treatment coverage', y = 'Proportion most \ncost-effective choice', fill = 'intervention') +
    scale_fill_manual(values = colors) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
    theme_classic()

  (A + B + C + D + E) + plot_layout(guides = "collect", nrow = 3, ncol = 2) + plot_annotation(tag_levels = 'A')

  ggsave(paste0('./03_output/figure4_', season, '.pdf'), width = 9, height = 7)

}


univariateseason('highly seasonal')
univariateseason('seasonal')
univariateseason('perennial')


# < faceted by season ----------------------------------------------------------
stacked_plot <- function(cost_dose){

  colors <- c(itn, pbo, rtss_age, rtss_sv, smc)

  output <- scenarios %>%
    filter(cost_per_dose == cost_dose & delivery_cost == 1.62) %>%
    filter(intervention %in% c('ITN 10% increase','ITN PBO','SMC','RTS,S age-based','RTS,S seasonal')) %>%
    mutate(seasonality = factor(seasonality, levels = c('perennial', 'seasonal', 'highly seasonal'))) %>%
    group_by(ID, drawID) %>% arrange(ID, drawID, CE) %>%
    slice(1L)


  # by pfpr
  A <- ggplot(output) +
    geom_bar(aes(x = as.factor(pfpr), fill=intervention_f), position = "fill", show.legend = F) +
    labs(x = 'PfPR', y = 'Proportion most \ncost-effective choice', fill = 'intervention') +
    scale_fill_manual(values = colors) +
    facet_grid(~ seasonality) +
    theme_classic() +
    theme(text = element_text(size = 14))

  # by current ITN usage
  B <- ggplot(output) +
    geom_bar(aes(x = factor(ITNuse), fill = intervention_f), position = "fill", show.legend = F) +
    labs(x = 'ITN use', y = 'Proportion most \ncost-effective choice', fill = 'intervention') +
    scale_fill_manual(values = colors) +
    facet_grid(~ seasonality) +
    theme_classic() +
    theme(text = element_text(size = 14))

  # by insecticide resistance
  C <- ggplot(output) +
    geom_bar(aes(x = factor(resistance), fill = intervention_f), position = "fill", show.legend = F) +
    labs(x = 'Resistance', y='Proportion most \ncost-effective choice', fill = 'intervention') +
    scale_fill_manual(values = colors) +
    facet_grid(~ seasonality) +
    theme_classic() +
    theme(text = element_text(size = 14))

  # by ITN distribution efficiency
  ITNefficient <- function(var, label) {
    scenarios %>%
      filter(cost_per_dose == cost_dose & delivery_cost == 1.62) %>%
      filter(intervention %in% c('ITN 10% increase','ITN PBO','SMC','RTS,S age-based','RTS,S seasonal')) %>%
      mutate(seasonality = factor(seasonality, levels = c('perennial', 'seasonal', 'highly seasonal'))) %>%
      group_by(ID, drawID) %>% arrange(ID, drawID, {{var}}) %>%
      slice(1L) %>% select(intervention, intervention_f, seasonality, pfpr, {{var}}) %>%
      mutate(model = label) %>%
      rename(CE = {{var}})
  }

  output2 <- ITNefficient(CE, 'standard') %>%
    full_join(ITNefficient(CE_ITNmin, 'less efficient')) %>%
    full_join(ITNefficient(CE_ITNmax, 'more efficient'))


  D <- ggplot(output2) +
    geom_bar(aes(x = factor(model, levels = c('less efficient', 'standard', 'more efficient'),
                            labels = c('less', 'standard', 'more')), fill = intervention_f), position = "fill") +
    labs(x = 'ITN efficiency', y = 'Proportion most \ncost-effective choice', fill = 'intervention') +
    scale_fill_manual(values = colors) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
    facet_grid(~ seasonality) +
    theme_classic() +
    theme(text = element_text(size = 14),
          axis.text.x = element_text(size = 8))

  E <- ggplot(output) +
    geom_bar(aes(x = factor(treatment, labels = c('low', 'medium', 'high')), fill = intervention_f), position = "fill", show.legend = F) +
    labs(x = 'Treatment coverage', y = 'Proportion most \ncost-effective choice', fill = 'intervention') +
    scale_fill_manual(values = colors) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
    facet_grid(~ seasonality) +
    theme_classic() +
    theme(text = element_text(size = 14),
          axis.text.x = element_text(size = 9))

  legend <- cowplot::get_legend(D)

  Dm <- D + theme(legend.position = "none")


  (A + B + C + Dm + E + legend) +
    plot_layout(nrow = 3) +
    plot_annotation(tag_levels = list(c('A', 'B', 'C', 'D', 'E', '')))

  ggsave(paste0('./03_output/figure4_', cost_dose, '.pdf'), width = 10, height = 8)

}

stacked_plot(2.69)
stacked_plot(6.52)
stacked_plot(12.91)
stacked_plot(17.36)


# print stats
output <- scenarios %>%
  filter(cost_per_dose == 6.52 & delivery_cost == 1.62) %>%
  filter(intervention %in% c('ITN 10% increase','ITN PBO','SMC','RTS,S age-based','RTS,S seasonal')) %>%
  mutate(seasonality = factor(seasonality, levels = c('perennial', 'seasonal', 'highly seasonal'))) %>%
  group_by(ID, drawID) %>% arrange(ID, drawID, CE) %>%
  slice(1L)

# by ITN distribution efficiency
ITNefficient <- function(var, label) {
  scenarios %>%
    filter(cost_per_dose == 6.52 & delivery_cost == 1.62) %>%
    filter(intervention %in% c('ITN 10% increase','ITN PBO','SMC','RTS,S age-based','RTS,S seasonal')) %>%
    mutate(seasonality = factor(seasonality, levels = c('perennial', 'seasonal', 'highly seasonal'))) %>%
    group_by(ID, drawID) %>% arrange(ID, drawID, {{var}}) %>%
    slice(1L) %>% select(intervention, intervention_f, seasonality, pfpr, {{var}}) %>%
    mutate(model = label) %>%
    rename(CE = {{var}})
}

output2 <- ITNefficient(CE, 'standard') %>%
  full_join(ITNefficient(CE_ITNmin, 'less efficient')) %>%
  full_join(ITNefficient(CE_ITNmax, 'more efficient'))

prop_CE <- function(data, var){
  data %>%
    mutate(intervention = ifelse(grepl('RTS,S', intervention), 'RTS,S', intervention),
           intervention_f = factor(intervention, levels = c('ITN 10% increase', 'ITN PBO', 'RTS,S', 'SMC'))) %>%
    group_by(seasonality, {{var}}, intervention_f) %>%
    summarize(n = n()) %>%
    ungroup() %>%
    group_by(seasonality, {{var}}) %>%
    mutate(t = sum(n), p = n / t * 100) %>%
    filter(intervention_f == 'RTS,S')
}

print(prop_CE(output, sim_length), n = 30) # overall
print(prop_CE(output, pfpr), n = 30) # overall
print(prop_CE(output, ITNuse), n = 40)
print(prop_CE(output, resistance), n = 30)
print(prop_CE(output2, model), n = 40)
print(prop_CE(output, treatment), n = 30)

output %>%
  group_by(seasonality, intervention) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  group_by(seasonality) %>%
  mutate(t = sum(n), p = n / t * 100)

