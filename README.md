# ButterflyLifespan

## Description

This R package provides methods and data to enable reproduction of the results presented in "Estimation of adult butterfly longevity using long-term citizen-science count data" (Clarke, Dennis and McCrea 2026)

## Installing the package

The package can be installed directly from this github repository from within an active R session by using the **remotes** package. If you do not yet have **remotes** it can be
installed from **CRAN** using


```R
install.packages("remotes")
```

You can now install **ButterflyLifespan** using

```R
library(remotes)
remotes::install_github("grosed/ButterflyLifespan/R-package")
```

## Documentation

The documentation for the methods and data in the package can be accessed in the usual way using the **help** function.

Additionally, the **documentation** directory in this repository contains a pdf version of the package documentation. 

## Examples

An example is provided that demonstrates the process of fitting a GAI multiyear stopover model using the package.
The output is then processed to obtain lifespan estimates with associated uncertainties. To minimise fitting time a random subsample of 30 sites is used which
spans 10 years from 2009-2018. Dark Green Fritillary data are used in the weekly format and lifespan is modelled as a linear trend across time using
the 'phi_type="slope"' option. The data frame produced can be used to create a plot showing how lifespan varies by year.

The example is located in the examples directory of this repository and is available in the form of a jupyter notebook or as a standalone R script.

