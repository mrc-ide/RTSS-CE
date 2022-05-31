# Specify HPC options ----------------------------------------------------------
library(didehpc)
setwd('Q:/GF-RTSS-CE')

options(didehpc.cluster = "fi--didemrchnb",
        didehpc.username = "htopazia")

source('./02_code/HPC/functions.R')

# transfer the new malariasimulation folder manually to contexts or delete and reinstall using conan
# remotes::install_github('mrc-ide/malariasimulation@dev', force=T)
src <- conan::conan_sources("github::mrc-ide/malariasimulation@dev")

ctx <- context::context_save(path = "Q:/contexts",
                             sources = c('./02_code/HPC/functions.R'),
                             packages = c("dplyr", "malariasimulation"),
                             package_sources = src)

share <- didehpc::path_mapping('Home drive', "Q:", '//fi--san03.dide.ic.ac.uk/homes/htopazia', "M:")
config <- didehpc::didehpc_config(shares = share,
                                  use_rrq = FALSE,
                                  cores = 1,
                                  cluster = "fi--didemrchnb", # fi--dideclusthn OR fi--didemrchnb
                                  parallel = FALSE)

obj <- didehpc::queue_didehpc(ctx, config = config) # check for latest v. of packages


# Set up your job --------------------------------------------------------------

library(tidyverse)
year <- 365

# FIRST
seas_name <- 'highly seasonal'
seasonality <- list(c(0.284596,-0.317878,-0.0017527,0.116455,-0.331361,0.293128,-0.0617547))
s1 <- crossing(seasonality, seas_name)

# SECOND
seas_name <- 'seasonal'
seasonality <- list(c(0.285505,-0.325352,-0.0109352,0.0779865,-0.132815,0.104675,-0.013919))
s2 <- crossing(seasonality, seas_name)

# THIRD
seas_name <- 'perennial'
seasonality <- list(c(0.2852770,-0.0248801,-0.0529426,-0.0168910,-0.0216681,-0.0242904,-0.0073646))
s3 <- crossing(seasonality, seas_name)

stable <- rbind(s1, s2, s3)

# loop over malariasimulation runs
init_EIR <- c(0.01, 0.1, seq(1,9,1), seq(10, 100, by=5), seq(110,300, by=10),350,400) # set EIR values
ITN <- c('pyr')
ITNuse <- c(0,0.25,0.50,0.75)
treatment <- c(0.30, 0.45, 0.60)
combo <- crossing(stable, init_EIR, ITN, ITNuse, treatment) %>%
  mutate(name = paste0('EIR', "_", seas_name, "_", init_EIR, "_", ITNuse, "_", treatment))

combo <- combo %>%
         filter(!(seas_name == 'seasonal' & init_EIR > 100)) %>%
         filter(!(seas_name == 'perennial' & init_EIR > 100))


# Run tasks -------------------------------------------------------------------
combo <- combo %>% mutate(f = paste0("./03_output/PR_EIR/",combo$name,".rds")) %>%
  mutate(exist=case_when(file.exists(f) ~ 1, !file.exists(f) ~ 0)) %>%
  filter(exist==0) %>%
  select(-f, -exist)

t <- obj$enqueue_bulk(combo, PRmatch)
t$status()


# < run case-study matches -----------------------------------------------------
ITNuse <- c(0.3, 0.5, 0.6)
treatment <- c(0.45)

combo <- crossing(stable, init_EIR, ITN, ITNuse, treatment) %>%
  mutate(name = paste0('EIR', "_", seas_name, "_", init_EIR, "_", ITNuse, "_", treatment))

combo <- combo %>%
  filter(!(seas_name == 'seasonal')) %>%
  filter(!(seas_name == 'perennial' & init_EIR > 100))

combo <- combo %>% mutate(f = paste0("./03_output/PR_EIR/",combo$name,".rds")) %>%
  mutate(exist=case_when(file.exists(f) ~ 1, !file.exists(f) ~ 0)) %>%
  filter(exist==0) %>%
  select(-f, -exist)

t <- obj$enqueue_bulk(combo, PRmatch)
t$status()


# DHS matches ------------------------------------------------------------------
ITNuse <- c(0.22, 0.36, 0.39, 0.44, 0.58, 0.60) # runs for admin1
combo <- crossing(stable, init_EIR, ITN, ITNuse) %>% mutate(name = paste0('EIRadmin1', "_", seas_name, "_", init_EIR, "_", ITNuse))

combo <- combo %>%
  filter(!(seas_name == 'seasonal' & init_EIR > 70)) %>%
  filter(!(seas_name == 'perennial' & init_EIR > 100)) %>%
  filter(!(ITNuse %in% c(0.22, 0.44) & seas_name != 'seasonal')) %>%
  filter(!(ITNuse %in% c(0.36, 0.58) & seas_name != 'highly seasonal')) %>%
  filter(!(ITNuse %in% c(0.39, 0.60) & seas_name != 'perennial'))


combo <- combo %>% mutate(f = paste0("./03_output/PR_EIR/",combo$name,".rds")) %>%
  mutate(exist=case_when(file.exists(f) ~ 1, !file.exists(f) ~ 0)) %>%
  filter(exist==0) %>%
  select(-f, -exist)

t <- obj$enqueue_bulk(combo, PRmatch_u5)
t$status()


# Results ----------------------------------------------------------------------
library(tidyverse)
library(data.table)
library(fuzzyjoin)

# plot EIR / PfPR relationships
# < general settings -----
files <- list.files(path = "./03_output/PR_EIR", pattern = "EIR_", full.names = TRUE)
dat_list <- lapply(files, function (x) readRDS(x))
EIR_prev <-  do.call("rbind", dat_list) %>% as_tibble() %>%
  mutate(init_EIR = as.numeric(init_EIR),
         prev = as.numeric(prev),
         scenario = paste(seas_name, ITN, ITNuse, treatment, sep = '_'))

p <- ggplot(data=EIR_prev) +
  geom_point(aes(x=init_EIR, y=prev), color='black') +
  stat_smooth(aes(x=init_EIR, y=prev), color='red', method = 'gam', n=1000) +
  labs(x='initial EIR', y=expression(paste(italic(Pf),"PR"[2-10]))) +
  facet_wrap('scenario', nrow=10, ncol=4) +
  theme_classic()

p

ggsave('C:/Users/htopazia/OneDrive - Imperial College London/Github/GF-RTSS-CE/03_output/prev_EIR_plot_ITN.pdf', height=7, width=7)


# < admin1s ----
files <- list.files(path = "./03_output/PR_EIR", pattern = "EIRadmin1_", full.names = TRUE)
dat_list <- lapply(files, function (x) readRDS(x))
EIR_prev <-  do.call("rbind", dat_list) %>% as_tibble() %>%
  mutate(init_EIR = as.numeric(init_EIR),
         prev = as.numeric(prev),
         scenario = paste(seas_name, ITN, ITNuse, sep = '_'))

q <- ggplot(data=EIR_prev) +
  geom_point(aes(x=init_EIR, y=prev), color='black') +
  stat_smooth(aes(x=init_EIR, y=prev), color='red', method = 'gam', n=1000) +
  labs(x='initial EIR', y=expression(paste(italic(Pf),"PR"[5-59]))) +
  facet_wrap('scenario', nrow=3, ncol=2) +
  theme_classic()

q

ggsave('C:/Users/htopazia/OneDrive - Imperial College London/Github/GF-RTSS-CE/03_output/prev_EIR_plot_ITN_admin1.pdf', height=7, width=7)


# Match EIR and PfPR -----------------------------------------------------------

# < generic settings ----
# extract points from stat_smooth
p2 <- ggplot_build(p)
p2 <- p2$data[[2]]

pnames <- EIR_prev %>% select(scenario) %>% distinct() %>%
  arrange(scenario) %>% mutate(PANEL=row_number())

p2 <- p2 %>% as_tibble() %>% select(x,y,PANEL) %>% rename(init_EIR=x, pred=y) %>%
  mutate(PANEL = as.numeric(PANEL)) %>%
  left_join(pnames, by='PANEL')

# Pre-intervention baseline PfPR2-10
PfPR <- as_tibble_col(c(.10, .20, .40), column_name="pfpr")

# match via stat_smooth predictions
match <- PfPR %>%
  fuzzyjoin::difference_left_join(p2, by=c("pfpr"="pred"),
                                  max_dist=1, distance_col="dist") %>%
  group_by(pfpr, scenario) %>% slice_min(dist)

match <- match %>% rename(starting_EIR = init_EIR) %>%
  separate(col = scenario, into = c("seas_name", "ITN", "ITNuse", "treatment"), sep="_") %>%
  mutate(ITNuse = as.numeric(ITNuse),
         treatment = as.numeric(treatment)) %>%
  select(pfpr, starting_EIR, seas_name, ITN, ITNuse, treatment)

summary(match)

# save for other HPC runs
saveRDS(match, "./02_code/HPC/EIRestimates.rds")


# < admin1s ----
# extract points from stat_smooth
q2 <- ggplot_build(q)
q2 <- q2$data[[2]]

pnames <- EIR_prev %>% select(scenario) %>% distinct() %>%
  arrange(scenario) %>% mutate(PANEL=row_number())

q2 <- q2 %>% as_tibble() %>% select(x,y,PANEL) %>% rename(init_EIR=x, pred=y) %>%
  mutate(PANEL = as.numeric(PANEL)) %>%
  left_join(pnames, by='PANEL')

# Pre-intervention baseline PfPR
PfPR <- as_tibble_col(c(0.05, 0.07, 0.12, 0.18, 0.39, 0.40), column_name="pfpr")

# match via stat_smooth predictions
match <- PfPR %>%
  fuzzyjoin::difference_left_join(q2, by=c("pfpr"="pred"),
                                  max_dist=1, distance_col="dist") %>%
  group_by(pfpr, scenario) %>% slice_min(dist)

match <- match %>% rename(starting_EIR = init_EIR) %>%
  separate(col = scenario, into = c("seas_name", "ITN", "ITNuse"), sep="_") %>%
  mutate(ITNuse = as.numeric(ITNuse)) %>%
  select(pfpr, starting_EIR, seas_name, ITN, ITNuse)

summary(match)

# save for other HPC runs
saveRDS(match, "./02_code/HPC/EIRestimates_admin1.rds")
