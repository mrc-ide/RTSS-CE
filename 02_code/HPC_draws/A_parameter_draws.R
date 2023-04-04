# This script sets up draws of malariasimulation parameters

# set-up parameter draws for uncertainty runs
require(data.table)
library(tidyverse)

data.dir <- './01_data/parameter_draws/'


# get draws --------------------------------------------------------------------
# pull all parameter draw .txt files and combine
files <- list.files(path = data.dir, pattern = "*.txt", full.names = TRUE)
dat_list <- lapply(files, function (x) read.delim(x, header = F, col.names = c('var', 'value')))

dat <- rbindlist(dat_list, fill = TRUE, idcol = "file")

d <- dat |> group_by(file) |> pivot_wider(names_from = 'var', values_from = 'value')
head(d)
str(d)

saveRDS(d, './02_code/HPC_draws/parameter_draws.rds')
saveRDS(d, 'M:/Hillary/RTSS-CE/02_code/HPC_draws/parameter_draws.rds')
