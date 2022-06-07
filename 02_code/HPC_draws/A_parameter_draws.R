# set-up parameter draws for uncertainty runs
require(data.table)
library(tidyverse)

data.dir <- 'C:/Users/htopazia/OneDrive - Imperial College London/Github/MalariaLaunchR/inst/model_files/parameters/'


# get draws --------------------------------------------------------------------
# pull all parameter draw .txt files from MalariaLaunchR folder and combine
files <- list.files(path = data.dir, pattern = "*.txt", full.names = TRUE)
dat_list <- lapply(files, function (x) read.delim(x, header = F, col.names = c('var', 'value')))

dat <- rbindlist(dat_list, fill = TRUE, idcol="file")

d <- dat %>% group_by(file) %>% pivot_wider(names_from = 'var', values_from = 'value')
head(d)
str(d)

saveRDS(d, './02_code/HPC/parameter_draws.rds')

params <- malariasimulation::get_parameters(list())
summary(d$dur_A); table(params$da)

# model wiki https://github.com/mrc-ide/Malaria_model/wiki/model_parms
params <- malariasimulation::get_parameters(list(
  # https://github.com/mrc-ide/malariasimulation/blob/aace0c5f212285d22cff12b519cbc62daf79624b/R/compatibility.R
  dd = d$dur_D,
  dt = d$dur_T,
  da = d$dur_A, # value 195 from the old model is the right one! New says 200
  du = d$dur_U,
  del	#N/A
  dl	#N/A
  dpl	#N/A
  mup	#N/A
  mum	#N/A
  sigma_squared = d$sigma2,
  rm = d$dm,
  rvm = d$dvm,
  rb = d$db,
  rc = d$dc,
  rva = d$dv,
  rid =	d$dd,
  b0 = d$bh,
  b1 = d$bmin,
  ib0 = d$IB0,
  kb = d$kb,
  ub = d$ub,
  uc = d$uc,
  uv = d$uv,
  ud = d$ud,
  cd = d$cD,
  gamma1 = d$gamma_inf,
  cu = d$cU,
  ct = d$cT, # values in malariasim are more precise
  a0 = d$a0,
  rho = d$rho,
  phi0 = d$phi0,
  phi1 = d$phi1,
  ic0 = d$IC0,
  kc = d$kc,
  theta0 = d$theta0,
  theta1 = d$theta1,
  kv = d$kv,
  fv0 = d$fv0,
  av = d$av0,
  gammav = d$gammav,
  iv0 = d$IV0,
  de = d$dur_E,
  delay_gam = d$latgam,
  dem = d$latmosq,
  fd0 = d$fd0,
  ad = d$ad0,
  gammad = d$gammad,
  d1 = d$dmin,
  id0 = d$ID0,
  kd = d$kd,
  average_age = round( 1/ d$eta),
  pcm = d$P_IC_M,
  pvm = d$P_IV_M
))

