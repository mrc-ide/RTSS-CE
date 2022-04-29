# PfPR EIR match ---------------------------------------------------------------

PRmatch <- function(x, y){

  # read in selected scenario
  data <- readRDS('./02_code/HPC_draws/baselinescenarios.rds')[x,]

  # choose a parameter set from baseline scenarios
  p <- data$params

  # choose a parameter draw
  d <- readRDS('./02_code/HPC/parameter_draws.rds')[y,]

  # over-write malariasimulation parameters to match the parameter draw
  p$human_population = 10000
  p$dd = d$dur_D
  p$dt = d$dur_T
  p$da = d$dur_A # value 195 from the old model is the right one! New says 200
  p$du = d$dur_U
  p$sigma_squared = d$sigma2
  p$rm = d$dm
  p$rvm = d$dvm
  p$rb = d$db
  p$rc = d$dc
  p$rva = d$dv
  p$rid =	d$dd
  p$b0 = d$bh
  p$b1 = d$bmin
  p$ib0 = d$IB0
  p$kb = d$kb
  p$ub = d$ub
  p$uc = d$uc
  p$uv = d$uv
  p$ud = d$ud
  p$cd = d$cD
  p$gamma1 = d$gamma_inf
  p$cu = d$cU
  p$ct = d$cT # values in malariasim are more precise
  p$a0 = d$a0
  p$rho = d$rho
  p$phi0 = d$phi0
  p$phi1 = d$phi1
  p$ic0 = d$IC0
  p$kc = d$kc
  p$theta0 = d$theta0
  p$theta1 = d$theta1
  p$kv = d$kv
  p$fv0 = d$fv0
  p$av = d$av0
  p$gammav = d$gammav
  p$iv0 = d$IV0
  p$de = d$dur_E
  p$delay_gam = d$latgam
  p$dem = d$latmosq
  p$fd0 = d$fd0
  p$ad = d$ad0
  p$gammad = d$gammad
  p$d1 = d$dmin
  p$id0 = d$ID0
  p$kd = d$kd
  p$average_age = round(1 / d$eta)
  p$pcm = d$P_IC_M
  p$pvm = d$P_IV_M

  # calibration ref: https://mrc-ide.github.io/cali/articles/Basic_calibration.html
  # define target: PfPR2-10 value
  target <- data$pfpr

  # time points at which to match target = years 4 to 6
  year <- 365
  target_tt <- seq(4*year, 6*year, 100)

  # run calibration model
  set.seed(123)
  out <- calibrate(parameters = p,
                   target = target,
                   target_tt = target_tt,
                   summary_function = summary_pfpr_2_10,
                   tolerance = 0.02,
                   interval = c(1, 300))

  # store init_EIR results as an .rds file to be read in later
  PR <- data.frame(scenarioID = x, drawID = y)
  PR$init_EIR <- out$root
  PR$ID <- data$ID

  saveRDS(PR, paste0('./03_output/PR_EIR/PRmatch_', x , '_', y, '.rds')

}

