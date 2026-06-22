## Install the ButterflyLifespan package

library(remotes)
remotes::install_github("grosed/ButterflyLifespan/R-package")

## load the library

library(ButterflyLifespan)

## lazy load the data

data(dark_green_fritillary_weekly)
dgf_week <- dark_green_fritillary_weekly


## analyse the data for 2016

output <- analysis(dgf_week,2016,"weekly")

lifespan <- output[[1]][["lifespan"]]

