# Figures & Tables -------------------------------------------------------------
# packages
library(tidyverse)
library(fuzzyjoin)
library(plotly)
library(kableExtra)
library(malariasimulation)
library(data.table)
library(RcppRoll)

library(patchwork)
library(scales)
library(LaCroixColoR)

# devtools::install_github('mrc-ide/malariasimulation@dev', force=TRUE)
# devtools::install_github('johannesbjork/LaCroixColoR')

# Look here: https://github.com/htopazian/rtss_malariasimulation
# find ideas in rmarkdowns for figures and insert code if helpful

dat <- readRDS("./03_output/rtss_raw.rds")
averted <- readRDS("./03_output/rtss_avert.rds")


# seasonality ------------------------------------------------------------------

# get starting parameters
seas_name <- 'highly seasonal'
seasonality <- list(c(0.284596,-0.317878,-0.0017527,0.116455,-0.331361,0.293128,-0.0617547))
s1 <- crossing(seasonality, seas_name)

seas_name <- 'seasonal'
seasonality <- list(c(0.285505,-0.325352,-0.0109352,0.0779865,-0.132815,0.104675,-0.013919))
s2 <- crossing(seasonality, seas_name)

seas_name <- 'perennial'
seasonality <- list(c(0.2852770,-0.0248801,-0.0529426,-0.0168910,-0.0216681,-0.0242904,-0.0073646))
s3 <- crossing(seasonality, seas_name)

stable <- rbind(s1, s2, s3)

find_peak <- function(seasonality, seas_name){
  month <- 365/12

  params <- get_parameters(list(
    human_population = 1000,
    model_seasonality = TRUE,
    g0 = unlist(seasonality)[1],
    g = unlist(seasonality)[2:4],
    h = unlist(seasonality)[5:7],
    individual_mosquitoes = FALSE))

  peak <- peak_season_offset(params)
  smc <- crossing('smc', round(c(peak+c(-1.5,-0.5,0.5,1.5)*month),0))
  colnames(smc) <- c('intervention', 'month')

  first <- round((peak-month*3.5),0)
  sv <- crossing('sv', first + round(c(0,1*month,2*month)))
  colnames(sv) <- c('intervention', 'month')

  interventions <- rbind(smc, sv) %>%
    mutate(seasonality = seas_name,
           month = month/365*12,
           color = ifelse(intervention=='smc', 'red', 'blue'))

  return(interventions)
}

output <- map2_dfr(stable$seasonality, stable$seas_name, find_peak)

baseline <- dat %>%
  filter(ITNuse==0, RTSS=='none', ITN=='pyr') %>%
  mutate(inc_month=n_inc_clinical_0_1825/n_0_1825, na.rm=T) %>%
  filter(month %in% seq(1,12,1))


ggplot(data=baseline) +
  geom_line(data=baseline, aes(x=month, y=inc_month, color=as.factor(pfpr))) +
  geom_vline(data=output[output$intervention=='smc',], aes(xintercept=month, alpha='SMC'), lty=2, color='red') +
  geom_vline(data=output[output$intervention=='sv',], aes(xintercept=month, alpha='RTS,S'), lty=2, color='blue') +
  labs(title='distribution of SMC and seasonal RTS,S; ages 0-5 years',
       x="Timesteps (month)",
       y="Monthly clinical incidence",
       color = 'PfPR') +
  scale_y_continuous(limits=c(0,0.8), breaks=seq(0,0.8,0.2)) +
  scale_x_continuous(limits=c(1,12), breaks=seq(1,12,1)) +
  facet_wrap('seasonality', ncol = 1) +
  theme_classic() +
  scale_alpha_manual(values = c(rep(1,2))) +
  guides(alpha = guide_legend(title = 'intervention',
                               override.aes = list(color = c('blue','red'))))

ggsave('./03_output/seasonality.pdf', width=6, height=6)


# bed nets and resistance ------------------------------------------------------
ITNdata <- dat %>%
  filter(RTSS=='none' & SMC==0) %>%
  filter(ITNuse==.50) %>%
  mutate(inc_month=n_inc_clinical_0_1825/n_0_1825, na.rm=T) %>%
  filter(month %in% seq(1,12,1))

ggplot(data=ITNdata) +
  geom_line(aes(x=month, y=inc_month, color=as.factor(interaction(ITN, resistance)))) +
  labs(title='distribution of SMC and seasonal RTS,S, ages 0-5 years',
       x="Timesteps (month)",
       y="Monthly clinical incidence",
       color = 'PfPR') +
  scale_y_continuous(limits=c(0,0.5), breaks=seq(0,0.5,0.1)) +
  scale_x_continuous(limits=c(1,12), breaks=seq(1,12,1)) +
  facet_wrap(c('seasonality','pfpr'), ncol = 3, nrow = 3) +
  theme_classic() +
  scale_alpha_manual(values = c(rep(1,2))) +
  guides(alpha = guide_legend(title = 'intervention',
                              override.aes = list(color = c('blue','red'))))


# cases averted ----------------------------------------------------------------
avertdat <- averted %>%

ggplot(data = avertdat) +
  geom_bar(aes(x = , y = , fill = )) +
  labs(title = "cases averted, aggregated years 1-5",
       xlab = "scenario",
       ylab = "cases")



