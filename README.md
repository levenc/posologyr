
<!-- README.md is generated from README.Rmd. Please edit that file -->

# posologyr [<img src="man/figures/logo_120.png"  align="right" />](https://github.com/levenc/posologyr/)

<!-- badges: start -->

[![R-CMD-check](https://github.com/levenc/posologyr/workflows/R-CMD-check/badge.svg)](https://github.com/levenc/posologyr/actions)
<!-- badges: end -->

## Overview

`posologyr` is an R package for the individualization of drug dosage,
taking advantage of population pharmacokinetics, patient
characteristics, and the results of therapeutic drug monitoring.

`posologyr` provides the following functions for dosage optimization:

- `poso_dose_conc()` estimates the optimal dose to reach a target
  concentration at any given time
- `poso_dose_auc()` estimates the optimal dose to reach a target AUC
- `poso_time_cmin()` computes the time needed to reach a target trough
  concentration (Cmin)
- `poso_inter_cmin()` estimates the optimal inter-dose interval to
  reliably achieve a target Cmin between each administration

Individual pharmacokinetic (PK) profiles can be estimated with or
without data from therapeutic drug monitoring (TDM):

- `poso_estim_map()` computes the Maximum A Posteriori (MAP), aka
  Empirical Bayes Estimates (EBE), of individual PK parameters from the
  results of TDM
- `poso_simu_pop()` samples from the the a priori distributions of PK
  parameters
- `poso_estim_sir()` estimates the full posterior distributions of
  individual PK parameters by Sequential Importance Resampling (SIR)
  from the results of TDM

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

A `posologyr` prior ppk model is a list of six objects. Its structure is
described in `vignette("posologyr_user_defined_models")`.

``` r
tobramycin_fictional <- list(
  ppk_model   = rxode2::rxode({
    centr(0) = 0;
    TVke  = THETA_ke*(CLCREAT/67.8)^0.89*(WT/66.4)^(-1.09);
    TVV   = THETA_V*(WT/66.4)^0.80;
    TVk12 = THETA_k12;
    TVk21 = THETA_k21;
    ke    = TVke*exp(ETA_ke);
    V     = TVV*exp(ETA_V);
    k12   = TVk12;
    k21   = TVk21;
    Cc    = centr/V;
    d/dt(centr)  = - ke*centr - k12*centr + k21*periph;
    d/dt(periph) =            + k12*centr - k21*periph;
    d/dt(AUC)    =   Cc;
  }),
  error_model = function(f,sigma){
    g <- sigma[1] + sigma[2]*f
    return(g)
  },
  theta = c(THETA_ke=0.21, THETA_V=19.8,THETA_k12=0.041, THETA_k21=0.12),
  omega = lotri::lotri({ETA_ke + ETA_V + ETA_k12 + ETA_k21 ~
      c(0.08075,
        0      , 0.01203,
        0      , 0      ,  0,
        0      , 0      ,  0, 0)}),
  covariates  = c("CLCREAT","WT"),
  sigma       = c(additive_a = 0, proportional_b = 0.198))
```

With these two elements, one can estimate the patient’s MAP-EBE PK
parameters, and the individual concentrations over time.

``` r
poso_estim_map(patient_data,tobramycin_fictional)

#> $eta
#>      ETA_ke       ETA_V     ETA_k12     ETA_k21 
#> -0.09253437  0.06705227  0.00000000  0.00000000 
#> 
#> $model
#> ── Solved rxode2 object ──
#> ── Parameters ($params): ──
#>    THETA_ke     THETA_V   THETA_k12   THETA_k21      ETA_ke       ETA_V 
#>  0.21000000 19.80000000  0.04100000  0.12000000 -0.09253437  0.06705227 
#> ── Initial Conditions ($inits): ──
#>  centr periph    AUC 
#>      0      0      0 
#> ── First part of data (object): ──
#> # A tibble: 151 × 15
#>    time  TVke   TVV TVk12 TVk21    ke     V   k12   k21    Cc centr periph    AUC CLCREAT    WT
#>   <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>  <dbl>  <dbl>   <dbl> <dbl>
#> 1   0   0.249  19.5 0.041  0.12 0.227  20.8 0.041  0.12   0      0   0      0          80    65
#> 2   0.1 0.249  19.5 0.041  0.12 0.227  20.8 0.041  0.12  19.0  395.  0.809  0.952      80    65
#> 3   0.2 0.249  19.5 0.041  0.12 0.227  20.8 0.041  0.12  37.4  779.  3.20   3.78       80    65
#> 4   0.3 0.249  19.5 0.041  0.12 0.227  20.8 0.041  0.12  55.4 1153.  7.10   8.42       80    65
#> 5   0.4 0.249  19.5 0.041  0.12 0.227  20.8 0.041  0.12  72.9 1517. 12.5   14.8        80    65
#> 6   0.5 0.249  19.5 0.041  0.12 0.227  20.8 0.041  0.12  89.9 1872. 19.2   23.0        80    65
#> # … with 145 more rows
```

Examples of dosage adjustment with posologyr are available in the
following vignettes:

- `vignette("case_study_amikacin")`
- `vignette("case_study_vancomycin")`

## Performance of the MAP-EBE algorithm in posologyr

Posologyr showed comparable performance to NONMEM:

- Pharmaceutics **2022**, 14(2), 442;
  [doi:10.3390/pharmaceutics14020442](https://doi.org/10.3390/pharmaceutics14020442)
- Supporting data: <https://github.com/levenc/posologyr-pharmaceutics>

## Acknowledgments

`posologyr`’s estimation functions were based on [Marc Lavielle’s shiny
app](http://shiny.webpopix.org/mcmc/bayes1/)
