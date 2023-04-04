# Determine input ITN usage to reach target average population usage values

library(netz)

# test various dist values until values match target ITN values
model_usage <- netz::population_usage(distribution = rep(0.858, 10),
                                      timesteps = 30 * 365,
                                      distribution_timesteps = seq(1, 30, 3) * 365,
                                      half_life = 5 * 365)

mean(tail(model_usage, 15 * 365))
plot(model_usage, t = "l")

# ITN values
# 0 = 0
# 0.10 = 0.066
# 0.20 = 0.141
# 0.25 = 0.184
# 0.30 = 0.231
# 0.35 = 0.283
# 0.40 = 0.339
# 0.50 = 0.473
# 0.60 = 0.641
# 0.65 = 0.742
# 0.70 = 0.858
# 0.75 = 0.994
# 0.85 = 1 (will result in 0.752)
