## Install the ButterflyLifespan package

library(remotes)
remotes::install_github("grosed/ButterflyLifespan/R-package")

## load the library

library(ButterflyLifespan)

## lazy load the data

data(dark_green_fritillary_weekly)
dgf_week <- dark_green_fritillary_weekly


## for demonstration purposes, take a random subset of the data

set.seed(1)
dgf_week <- dgf_week %>% filter(SITENO%in%sample(unique(dgf_week$SITENO),30)) %>% filter(YEAR > 2008)

## analyse the data 

output <- analysis_multiyear(dgf_week,"weekly","slope")

lifespan <- output[[1]][["lifespan"]]






