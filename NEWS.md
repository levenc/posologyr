# posologyr v0.1.1.9000

* `poso_estim_map()` allows the estimation of the individual parameters by 
adaptive MAP forecasting (cf. doi: 10.1007/s11095-020-02908-7) with 
`adapt=TRUE`.
* `poso_estim_sir()` allows the estimation of the posterior distribution 
of individual parameters by Sequential Importance Resampling (SIR). It is 
roughly 25 times faster than `poso_estim_mcmc()` for 1000 samples.
* `MASS:mvrnorm` is replaced by `mvtnorm::rmvnorm` for multivariate normal 
distributions. 
* Input validation is added to all exported functions.
* `poso_dose_conc` is the new name of `poso_dose_ctime`.
* `poso_dose_auc()`, and `poso_dose_conc()` now work with distributions of ETA, 
and not only with point estimates (such as the MAP).
* `vignette("uncertainty_estimates")` is removed.

# posologyr v0.1.1

* `poso_time_cmin()`, `poso_dose_ctime()`, and `poso_dose_auc()` now work for 
multiple dose regimen.
* `poso_inter_cmin()` allows the optimization of the inter-dose interval for 
multiple dose regimen.
* `vignette("case_study_vancomycin")` illustrates AUC-based optimal dosing, 
multiple dose regimen, and continuous intravenous infusion.

# posologyr v0.1.0

First public release.
