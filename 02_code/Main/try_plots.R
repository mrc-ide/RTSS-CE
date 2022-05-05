dalyoutput_cost <- readRDS('./03_output/dalyoutput_cost.rds')

# Groups are defined by:
# pfpr, seasonality, ITNuse, resistance, SMC (depends on seasonality), RTSS cost per dose and delivery cost

# scenarios are: ITNboost, ITN=pbo, RTSS (RTSScov)
# baseline is:
# ITNboost=0
# RTSS = none
# RTSS_cov=0 (this is the same as RTSS=none)
# ITN=pyr
# Also group by:
# 3 options for dose cost and delivery cost

# scenarios are:
# ITNboost=1
# RTSS = EPI or SV
# RTSScov = 0.85
# ITN=pbo

# SMC = always 0 in perennial, always 0.85 in highly seasonal
# 0 or 0.85 in seasonal setting

View(group_by(dalyoutput_cost, pfpr, seasonality, ITNuse, resistance, cost_per_dose, delivery_cost) %>%
  summarise(n()))

# All have 9 entries, except seasonal has 18 because of SMC

View(filter(dalyoutput_cost, pfpr==0.2, seasonality=="perennial", ITNuse==0, resistance==0,
            cost_per_dose==2.69, delivery_cost==0.96))
# 3 pyr without boost (baseline), 3 pyr with boost (scenario), 3 pbo without boost (scenario)
# 3 options within each is the RTSS

dalyoutput_cost2 <- as_tibble(dalyoutput_cost) %>%
  group_by(pfpr, seasonality, ITNuse, resistance, SMC, cost_per_dose, delivery_cost) %>%
  mutate(ID=cur_group_id(),
         scenario=ifelse(test=(ITNboost==0 & RTSS == "none" & ITN== "pyr"),yes="baseline",
                          no=paste0("scenario",ITNboost,"_",ITN,"_",RTSS))) %>%
  ungroup() %>% select(ID, scenario, daly, cost_total) %>%
  pivot_wider(names_from = scenario,
              values_from = c(daly, cost_total))
