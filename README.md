
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PESCAR

<!-- badges: start -->

<!-- badges: end -->

## Installation

# Windows

Install `PESCAR` from GitHub using `pak`:

``` r
install.packages("pak")
pak::pak("Fred-White94/PESCAR")
```

# Linux / Mac

``` r
# Installation
options(repos = c(CRAN = "https://cloud.r-project.org"))
install.packages("remotes")

# Install non-CRAN dependency
remotes::install_gitlab("uvabda/RpESCA", upgrade = "never")

# Install PESCAR
remotes::install_github(
  "Fred-White94/PESCAR",
  dependencies = TRUE,
  upgrade = "never",
  build_vignettes = FALSE
)
```

Then load the package:

``` r
library(PESCAR)
```

## Getting started

A worked simulation tutorial is provided in:

[Rendered simulation tutorial](tutorials/PESCAR_simulation_tutorial.md)

With code: [Simulation
tutorial](vignettes/PESCAR_simulation_tutorial.Rmd)

This tutorial introduces the basic workflow for simulating data, fitting
a PESCAR model, selecting tuning parameters, and inspecting the model
output.

To load a non rendered version of the tutorial in Rstudio:

``` r
file.show(system.file("PESCAR_simulation_tutorial.Rmd", package = "PESCAR"))
```
