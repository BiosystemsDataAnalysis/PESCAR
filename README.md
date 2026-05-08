
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PESCAR

<!-- badges: start -->

<!-- badges: end -->

## Installation

Install `PESCAR` from GitHub using `pak`:

``` r
install.packages("pak")
pak::pak("Fred-White94/PESCAR")
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
