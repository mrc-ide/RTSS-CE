# Set-up parameter draws for uncertainty runs
require(data.table)
library(tidyverse)

data.dir <- 'C:/Users/htopazia/OneDrive - Imperial College London/Github/MalariaLaunchR/inst/model_files/parameters/'


# get draws --------------------------------------------------------------------
# pull all parameter draw .txt files from MalariaLaunchR folder and combine
files <- list.files(path = data.dir, pattern = "*.txt", full.names = TRUE)
dat_list <- lapply(files, function (x) read.delim(x, header = F, col.names = c('var', 'value')))
dat <- rbindlist(dat_list, fill = TRUE, idcol="file")

dat <- dat %>% group_by(file) %>% pivot_wider(names_from = 'var', values_from = 'value')
head(dat)
str(dat)

saveRDS(dat, './02_code/HPC/parameter_draws.rds')


params <- malariasimulation::get_parameters(list(
  # IA0 = 0 for all
  ib0 = d$IB0,
  ic0 = d$IC0,
  id0 = d$ID0,
  iv0 = d$IV0,
  pcm = d$P_IC_M,
  pvm = d$P_IV_M,
  a0  = d$a0, # constant
  rb = d$ab0, # constant
  d$ac0,
  d$ad0,
  d$alpha_pcr,
  av = d$av0,
  d$beta_pcr
  b0 = d$bh,
  b1 = d$bmin, # constant
  cd = d$cD,
  ct = d$cT, # FIX**** mismatching parameter
  cu = d$cU,
  d$clin_immunity
  d$da
  d$db
  rc = d$dc, # constant
  d$dd
  d$det_immunity
  rm = d$dm,
  d1 = d$dmin,
  dd = d$dur_D, # constant
  de = d$dur_E, # constant
  d$dur_P
  dt = d$dur_T # constant
  du = d$dur_U,
  d$dur_min
  rva = d$dv,
  rvm = d$dvm,
  d$fb0
  d$fc0
  fd0 = d$fd0,
  fv0 = d$fv0,
  gamma1 = d$gamma_inf,
  d$gammab
  d$gammac
  gammad = d$gammad,
  gammav = d$gammav,
  d$inf_immunity
  d$ka
  kb = d$kb,
  kc = d$kc,
  kd = d$kd,
  kv = d$kv,
  delay_gam = d$latgam, # constant
  dem = d$latmosq, # constant
  d$par_A2
  d$par_immunity
  d$par_u
  phi0 = d$phi0
  phi1 = d$phi1
  d$tau_v
  rho = d$rho, # constant
  d$sev_immunity
  sigma_squared = d$sigma2, # constant
  d$separate_sev
  theta0 = d$theta0,
  theta1 = d$theta1,
  d$two_A
  ub = d$ub,
  uc = d$uc,
  ud = d$ud,
  uv = d$uv,
  d$wa
  d$dur_A # FIX**** mismatching parameter
  average_age = 1 / d$eta # constant

))

summary(d$theta0); table(params$theta0)
