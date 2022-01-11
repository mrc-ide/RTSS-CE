# packages
library(tidyverse)
library(fuzzyjoin)
library(plotly)
library(kableExtra)
library(malariasimulation)
library(data.table)

library(patchwork)
library(scales)
library(LaCroixColoR)

# devtools::install_github('mrc-ide/malariasimulation@dev', force=TRUE)
# devtools::install_github('johannesbjork/LaCroixColoR')
knitr::opts_knit$set(message=FALSE, warning=FALSE)

# Look here: https://github.com/htopazian/rtss_malariasimulation
# find ideas in rmarkdowns for figures and insert code if helpful
