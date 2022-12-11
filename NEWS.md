# posologyr v1.0.0.9999

## Additional features
* `poso_time_cmin()` can now estimate time needed to reach a selected trough 
concentration (Cmin) using the data from TDM directly
* `poso_dose_conc()` can now estimate an optimal dose to reach a target 
concentration following the events from TDM
* `poso_dose_auc()` can now estimate an optimal dose to reach a target auc 
following the events from TDM

# posologyr v1.0.0

## Breaking changes
* `posologyr()` is now an internal function, all exported functions take 
patient data and a prior model as input parameters
* The adaptive MAP forecasting option is removed

## Additional features
* `poso_estim_map()` provides an rxode2 model using MAP-EBE and the input dataset, 
with interpolation of covariates, to make plotting easier

## Internal changes
* RxODE import is updated to rxode2
* All tests are updated to take into account the internalization of the 
`posologyr()` function

## Bug fixes
* `poso_time_cmin()`, `poso_dose_auc()`, `poso_dose_conc()`, and 
`poso_inter_cmin()` no longer fail for models with IOV

# posologyr v0.2.0

* `poso_estim_sir()` estimates the posterior distribution of individual 
parameters by Sequential Importance Resampling (SIR). It is roughly 25 times 
faster than `poso_estim_mcmc()` for 1000 samples.
* `poso_estim_map()` allows the estimation of the individual parameters by 
adaptive MAP forecasting (cf. doi: 10.1007/s11095-020-02908-7) with 
`adapt=TRUE`.
* `poso_simu_pop()`, `poso_estim_map()`, and `poso_estim_sir()` now support 
models with both inter-individual (IIV) and inter-occasion variability (IOV).
* `MASS:mvrnorm` is replaced by `mvtnorm::rmvnorm` for multivariate normal 
distributions. 
* Input validation is added to all exported functions.
* `poso_estim_map()` now uses method="L-BFGS-B" in optim for better convergence 
of the algorithm.
* `poso_inter_cmin()` now uses method="L-BFGS-B" in optim for better convergence 
of the algorithm.
* `poso_dose_conc()` is the new name of `poso_dose_ctime()`.
* Issues #5 and #6 are fixed: `poso_time_cmin()`, `poso_dose_auc()`, 
`poso_dose_conc()`, and `poso_inter_cmin()` now work with prior and posterior 
distributions of ETA, and not only with point estimates (such as the MAP).
* A new `nocb` parameter is added to `posologyr()`. The interpolation method for
time-varying covariates can be either last observation carried forward (locf, 
the RxODE default), or next observation carried backward (nocb, the NONMEM 
default).
* `vignette("uncertainty_estimates")` is removed.
* The built-in models are removed.

# posologyr v0.1.1

* `poso_time_cmin()`, `poso_dose_ctime()`, and `poso_dose_auc()` now work for 
multiple dose regimen.
* `poso_inter_cmin()` allows the optimization of the inter-dose interval for 
multiple dose regimen.
* `vignette("case_study_vancomycin")` illustrates AUC-based optimal dosing, 
multiple dose regimen, and continuous intravenous infusion.

# posologyr v0.1.0

First public release.
