
<!-- README.md is generated from README.Rmd. Please edit that file -->

# posologyr [<img src="man/figures/logo_120.png"  align="right" />](https://github.com/levenc/posologyr/)

<!-- badges: start -->

[![R-CMD-check](https://github.com/levenc/posologyr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/levenc/posologyr/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Overview

Determine individual pharmacokinetic (and
pharmacokinetic-pharmacodynamic) profiles and use them to personalise
drug regimens. You provide the data (observations from therapeutic drug
monitoring, or TDM) and a population pharmacokinetic model, `posologyr`
provides the individual a posteriori estimate and allows you to
determine the optimal dosing.

`posologyr` provides the following functions for dosage optimization:

- `poso_dose_conc()` estimates the optimal dose to reach a target
  concentration at any given time
- `poso_dose_auc()` estimates the optimal dose to reach a target area
  under the concentration/time curve (AUC)
- `poso_time_cmin()` estimates the time needed to reach a target trough
  concentration (Cmin)
- `poso_inter_cmin()` estimates the optimal inter-dose interval to
  reliably achieve a target Cmin between each administration

Individual pharmacokinetic (PK) profiles can be estimated with or
without data from therapeutic drug monitoring (TDM):

- `poso_estim_map()` computes the Maximum A Posteriori (MAP), aka
  Empirical Bayes Estimates (EBE), of individual PK parameters from the
  results of TDM
- `poso_estim_sir()` estimates the full posterior distributions of
  individual PK parameters by Sequential Importance Resampling (SIR)
  from the results of TDM
- `poso_simu_pop()` samples from the the a priori distributions of PK
  parameters

`posologyr` takes advantage of the simulation framework provided by the
[rxode2](https://github.com/nlmixr2/rxode2) package.

## Installation

You can install the development version of `posologyr` from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("levenc/posologyr")
```

## Example of MAP-EBE estimation

`posologyr` allows the adaptation of dosage from two elements: a data
set, and a prior population PK model.

``` r
library("posologyr")
```

Data for input into `posologyr` is the same type of data input for
rxode2, see `vignette("patient_data_input")` for minimal examples.

``` r
patient_data <- data.frame(ID=1,
                           TIME=c(0.0,1.0,14.0),
                           DV=c(NA,25.0,5.5),
                           AMT=c(2000,0,0),
                           DUR=c(0.5,NA,NA),
                           EVID=c(1,0,0),
                           CLCREAT=80,
                           WT=65)
```

A `posologyr` prior ppk model is a named R list. Its structure is
described in `vignette("posologyr_user_defined_models")`.

``` r
mod_run001 <- list(
ppk_model = rxode2::rxode({
  TVCl = THETA_Cl;
  TVVc = THETA_Vc;
  TVKa = THETA_Ka;

  Cl = TVCl*exp(ETA_Cl);
  Vc = TVVc*exp(ETA_Vc);
  Ka = TVKa*exp(ETA_Ka);

  K20 = Cl/Vc;
  Cc = centr/Vc;

  d/dt(depot) = -Ka*depot;
  d/dt(centr) = Ka*depot - K20*centr;
  d/dt(AUC) = Cc;
}),
error_model = function(f,sigma) {
  dv <- cbind(f,1)
  g <- diag(dv%*%sigma%*%t(dv))
  return(sqrt(g))
},
theta = c(THETA_Cl=4.0, THETA_Vc=70.0, THETA_Ka=1.0),
omega = lotri::lotri({ETA_Cl + ETA_Vc + ETA_Ka ~
    c(0.2,
      0, 0.2,
      0, 0, 0.2)}),
sigma = lotri::lotri({prop + add ~ c(0.05,0.0,0.00)}))
```

With these two elements, one can estimate the patient’s MAP-EBE PK
parameters, and the individual concentrations over time.

``` r
poso_estim_map(patient_data,mod_run001)
#> $eta
#>     ETA_Cl     ETA_Vc     ETA_Ka 
#>  0.6019038 -0.4291736  0.1278476 
#> 
#> $model
#> ── Solved rxode2 object ──
#> ── Parameters ($params): ──
#>   THETA_Cl   THETA_Vc   THETA_Ka     ETA_Cl     ETA_Vc     ETA_Ka 
#>  4.0000000 70.0000000  1.0000000  0.6019038 -0.4291736  0.1278476 
#> ── Initial Conditions ($inits): ──
#> depot centr   AUC 
#>     0     0     0 
#> ── First part of data (object): ──
#> # A tibble: 151 × 12
#>    time  TVCl  TVVc  TVKa    Cl    Vc    Ka   K20     Cc depot centr    AUC
#>   <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>  <dbl> <dbl> <dbl>  <dbl>
#> 1   0       4    70     1  7.30  45.6  1.14 0.160  0        0    0   0     
#> 2   0.1     4    70     1  7.30  45.6  1.14 0.160  0.478  378.  21.8 0.0161
#> 3   0.2     4    70     1  7.30  45.6  1.14 0.160  1.83   716.  83.5 0.125 
#> 4   0.3     4    70     1  7.30  45.6  1.14 0.160  3.95  1017. 180.  0.408 
#> 5   0.4     4    70     1  7.30  45.6  1.14 0.160  6.75  1286. 307.  0.938 
#> 6   0.5     4    70     1  7.30  45.6  1.14 0.160 10.1   1526. 461.  1.78  
#> # ℹ 145 more rows
#> 
#> $event
#>      id time  amt evid dur
#>   1:  1  0.0   NA    0  NA
#>   2:  1  0.0 2000    1 0.5
#>   3:  1  0.1   NA    0  NA
#>   4:  1  0.2   NA    0  NA
#>   5:  1  0.3   NA    0  NA
#>  ---                      
#> 148:  1 14.6   NA    0  NA
#> 149:  1 14.7   NA    0  NA
#> 150:  1 14.8   NA    0  NA
#> 151:  1 14.9   NA    0  NA
#> 152:  1 15.0   NA    0  NA
```

The individual profile can be plotted easily

``` r
patient_001 <- poso_estim_map(patient_data,mod_run001)
plot(patient_001$model,Cc)
```

<img src="man/figures/README-map_plot-1.png" width="100%" />

Using `ggplot2` the observed data points can be added to the plot

``` r
# get the observations from the data
indiv_obs           <- patient_data[,c("DV","TIME")]

# set the names to match the names in the rxode2 model
names(indiv_obs)    <- c("value","time")

# call ggplot2
plot(patient_001$model,Cc) +
  ggplot2::geom_point(data=indiv_obs, size= 3, na.rm=TRUE)
```

<img src="man/figures/README-map_plot_dv-1.png" width="100%" />

## Performance of the MAP-EBE algorithm in posologyr

`posologyr` showed comparable performance to NONMEM MAP estimation with
option `MAXEVAL=0`:

- Pharmaceutics **2022**, 14(2), 442;
  [doi:10.3390/pharmaceutics14020442](https://doi.org/10.3390/pharmaceutics14020442)
- Supporting data: <https://github.com/levenc/posologyr-pharmaceutics>
