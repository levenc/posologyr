# posologyr v1.2.4.9000

## Additional feature
* The route of administration (i.e. the compartment in which the drug is to be administered) can now be specified in `poso_time_cmin()`, `poso_dose_conc()`, `poso_dose_auc()` and `poso_inter_cmin()`.

## Documentation
* The README illustrates a simple example of dose adaptation
* `vignette("route_of_administration")` shows how to select a route of administration for optimal dosing

## Bug fix
* Fix a bug where `poso_estim_map()`, `poso_estim_sir()` and `poso_simu_pop()` failed for models featuring a single parameter with IIV.

# posologyr v1.2.4

* Add ability to use rxode2 ui models for the `poso_*` functions.  Once the model has been parsed by `rxode2()` with this package the `model$posologyr` gives the list needed for `poso_*` functions

# posologyr v1.2.3

## Bug fix
* Fix a bug in `poso_dose_conc()`, `poso_dose_auc()` and `poso_inter_cmin()` where the returned estimate of the target value to be optimized against was always equal to zero.

## Documentation
* The documentation for `poso_time_cmin()`, `poso_dose_conc()`, and `poso_dose_auc()` now explicitly states the consequences of setting `tdm` to `TRUE`: which parameters are required, which parameters are ignored, and which parameters behave differently.
* The functions `poso_time_cmin()`, `poso_dose_conc()`, and `poso_dose_auc()` now return a warning if any of the input parameters are ignored.
* Fix incorrect information regarding the duration of the AUC in the documentation of `poso_dose_auc()`

# posologyr v1.2.2
* Relax the requirements of the NONMEM comparison test for time-varying covariates to account for computational differences observed with the alternative BLAS ATLAS on CRAN.

# posologyr v1.2.1
* Add a reference to Kang et al. (2012) <doi:10.4196/kjpp.2012.16.2.97> in the DESCRIPTION (as requested by CRAN)
* Fix messages to the console in the internal function `posologyr()` (as requested by CRAN)
* Fix assignment to parent environment in dose optim functions, using `parent.frame()` (as requested by CRAN)

# posologyr v1.2.0

## Additional features
* `poso_estim_map()`, `poso_estim_sir()` and `poso_estim_mcmc()` can now estimate individual PK profiles for multiple endpoints models (eg. PK-PD, parent-metabolite, blood-CSF...), using a different residual error model for each endpoint.
* `poso_time_cmin()`, `poso_dose_conc()`, `poso_dose_auc()` and `poso_inter_cmin()` now allow you to select the end point of interest for which you want to optimise, provided it is defined in the model.

## Documentation
* `vignette("a_priori_dosing")` illustrates a priori dose selection
* `vignette("a_posteriori_dosing")` illustrates a posteriori dose selection, using TDM data
* `vignette("auc_based_dosing")` shows how to select an optimal dose for a given target AUC using data from TDM
* `vignette("multiple_endpoints")` introduces the new multiple endpoints feature

## Internal changes
* The description of the package is updated

# posologyr v1.1.0

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
