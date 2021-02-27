# Population profile + individual covariates + individual observations

#' Plot typical population model prediction over individual
#' observations
#'
#' Plots the typical population pharmacokinetic model, using the
#' the individual covariates, over the individual observations
#'
#' @param prior_model A posologyr prior population pharmacokinetics model, a
#'    list of seven elements (see 'Details' for the description of the
#'    object)
#' @param dat Dataframe. An individual subject dataset following the
#'     structure of NONMEM/RxODE event records
#' @param covar a named vector of the individual covariates. If omitted,
#'     the covariates from the individual event record `dat` will be
#'     used
#'
#' @details
#' The default values of the arguments `prior_model` and `dat` correspond
#' to the objects created by the convenience function
#' \code{\link{load_ppk_model}}
#'
#' The posologyr prior population pharmacokinetics model is a list of
#' five elements:
#' \describe{
#'  \item{$ppk_model}{A RxODE model implementing the structural
#'      population pharmacokinetics model with no inter-individual
#'      variability, or residual error model}
#'  \item{$error_model}{A function of the residual error model}
#'  \item{$pk_prior}{A list of 3. `name`: a character vector of the names
#'      of the population pharmacokinetc paramters, `reference`: a named
#'      vector of the prior typical value of the population paramaters,
#'      `Omega`: a square variance-covariance matrix of the population
#'      parameters inter-individual variability}
#'  \item{$covariates}{A character vector of the covariates of
#'      the model}
#'  \item{$xi}{The estimates of the parameters of the residual error model}
#' }
#'
#' @return A \code{ggplot2} plot of the typical population time-concentration curve,
#'     using the individual covariates, over the individual observations
#'
#' @examples
#' plot_pop_cov()
#'
#' @export
poso_plot_pop <- function(prior_model=prior_ppk_model,dat=dat_posology,covar=NULL){

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package \"ggplot2\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if (is.null(covar)){
    #get covariates from the prior model and the dataset
    covar <- t(dat[1,prior_model$covariates]) #results in a matrix
    names(covar) <- prior_model$covariates
  }

  indiv_obs_plot <- dat[,c("DV","TIME")]
  names(indiv_obs_plot) <- c("value","time")

  pars_pop <- prior_model$pk_prior$reference

  plot_ppk_model <- RxODE::rxSolve(prior_model$ppk_model,
                                   c(pars_pop,covar),dat)

  plot_ppk_model$time<-seq(0,max(dat$TIME),by=0.2)

  plot(plot_ppk_model,Cc) + ggplot2::ylab("Central concentration") +
    ggplot2::geom_point(data=indiv_obs_plot, na.rm=TRUE)
}

#' Plot the individual Maximum A Posteriori prediction over individual
#' observations
#'
#' Plots the individual Maximum A Posteriori prediction (or Empirical
#' Bayes Estimate), using the individual covariates, over the individual
#' observations
#'
#' @param prior_model A posologyr prior population pharmacokinetics model, a
#'    list of seven elements (see 'Details' for the description of the
#'    object)
#' @param dat Dataframe. An individual subject dataset following the
#'     structure of NONMEM/RxODE event records
#' @param param_psi_map A vector of individual parameters. May be omitted
#'     if a `solved_model` is provided, in which case the
#'     \code{\link{poso_estim_map}}
#'     function will be called
#' @param solved_model An \code{\link[RxODE]{rxSolve}} solve object, created
#'     with the prior RxODE structural population pharmacokinetics model and
#'     the prior typical values of the population parameters from the
#'     `prior_model`, using `dat` as the event record. May be omitted
#'     if a vector of individual parameters `param_psi_map` is provided
#'
#' @details
#' The default values of the arguments `solved_model`, `prior_model` and
#' `dat` correspond to the objects created by the convenience function
#' \code{\link{load_ppk_model}}
#'
#' The posologyr prior population pharmacokinetics model is a list of
#' five elements:
#' \describe{
#'  \item{$ppk_model}{A RxODE model implementing the structural
#'      population pharmacokinetics model with no inter-individual
#'      variability, or residual error model}
#'  \item{$error_model}{A function of the residual error model}
#'  \item{$pk_prior}{A list of 3. `name`: a character vector of the names
#'      of the population pharmacokinetc paramters, `reference`: a named
#'      vector of the prior typical value of the population paramaters,
#'      `Omega`: a square variance-covariance matrix of the population
#'      parameters inter-individual variability}
#'  \item{$covariates}{A character vector of the covariates of
#'      the model}
#'  \item{$xi}{The estimates of the parameters of the residual error model}
#' }
#'
#' @return A \code{ggplot2} plot of the individual Maximum A Posteriori
#'     time-concentration curve, using the individual covariates, over the
#'     individual observations
#'
#' @examples
#' plot_indiv_map()
#'
#' @export
poso_plot_map <- function(prior_model=prior_ppk_model,dat=dat_posology,
                           param_psi_map=NULL,solved_model=solved_ppk_model){

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package \"ggplot2\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  if (is.null(param_psi_map)){#MAP estimates of the individual parameters
    if (!is.null(solved_model)){
      param_psi_map <- poso_estim_map(solved_model,prior_model,dat)
    } else {
      stop("Either param_psi_map or solved_model is needed for this function to work",
           call. = FALSE)
    }
  }

  covar <- t(dat[1,prior_model$covariates]) #results in a matrix
  names(covar) <- prior_model$covariates

  plot_ppk_model <- RxODE::rxSolve(prior_model$ppk_model,
                                   c(param_psi_map,covar),dat)

  plot_ppk_model$time <- seq(0,max(dat$TIME),by=0.2)

  indiv_obs_plot <- dat[,c("DV","TIME")]
  names(indiv_obs_plot) <- c("value","time")

  plot(plot_ppk_model,Cc, ylab="Central compartment") +
    ggplot2::geom_point(data=indiv_obs_plot, na.rm=TRUE)
}

#' Plot the individual posterior distribution over individual
#' observations
#'
#' Plots the posterior distribution, showing the uncertainty around
#' the prediction of the individual profile, using the individual
#' covariates, and over the individual observations
#'
#' @param prior_model A posologyr prior population pharmacokinetics model, a
#'    list of seven elements (see 'Details' for the description of the
#'    object)
#' @param dat Dataframe. An individual subject dataset following the
#'     structure of NONMEM/RxODE event records
#' @param param_psi_mcmc A dataframe of parameters from the posterior
#'     distribution, estimated by Markov Chain Monte Carlo. May be omitted
#'     if a `solved_model` is provided, in which case the
#'     \code{\link{poso_estim_mcmc}}
#'     function will be called
#' @param solved_model An \code{\link[RxODE]{rxSolve}} solve object, created
#'     with the prior RxODE structural population pharmacokinetics model and
#'     the prior typical values of the population parameters from the
#'     `prior_model`, using `dat` as the event record. May be omitted
#'     if a dataframe of parameters from the posterior distribution
#'     `param_psi_mcmc` is provided
#'
#' @details
#' The default values of the arguments `solved_model`, `prior_model` and
#' `dat` correspond to the objects created by the convenience function
#' \code{\link{load_ppk_model}}
#'
#' The posologyr prior population pharmacokinetics model is a list of
#' five elements:
#' \describe{
#'  \item{$ppk_model}{A RxODE model implementing the structural
#'      population pharmacokinetics model with no inter-individual
#'      variability, or residual error model}
#'  \item{$error_model}{A function of the residual error model}
#'  \item{$pk_prior}{A list of 3. `name`: a character vector of the names
#'      of the population pharmacokinetc paramters, `reference`: a named
#'      vector of the prior typical value of the population paramaters,
#'      `Omega`: a square variance-covariance matrix of the population
#'      parameters inter-individual variability}
#'  \item{$covariates}{A character vector of the covariates of
#'      the model}
#'  \item{$xi}{The estimates of the parameters of the residual error model}
#' }
#'
#' @return A \code{ggplot2} plot of the posterior distribution of the
#'     time-concentration curves, using the individual covariates, over the
#'     individual observations
#'
#' @examples
#' plot_indiv_mcmc()
#'
#' @export
poso_plot_mcmc <- function(prior_model=prior_ppk_model,dat=dat_posology,
                            param_psi_mcmc=NULL,solved_model=solved_ppk_model){

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package \"ggplot2\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if (is.null(param_psi_mcmc)){#MCMC estimates of the individual parameters
    if (!is.null(solved_model)){
      param_psi_mcmc <- poso_estim_mcmc(solved_model=solved_model,
                                        prior_model=prior_model,dat=dat_posology)
    } else {
      stop("Either param_psi_mcmc or solved_model is needed for this function to work",
           call. = FALSE)
    }
  }

  covar <- data.frame(dat[1,prior_model$covariates])
  names(covar) <- prior_model$covariates

  plot_ppk_model <- RxODE::rxSolve(prior_model$ppk_model,
                                   cbind(param_psi_mcmc,covar,row.names = NULL),dat)

  plot_ppk_model$time <- seq(0,max(dat$TIME),by=0.2)

  indiv_obs_plot <- dat[,c("DV","TIME")]
  names(indiv_obs_plot) <- c("eff","time")

  plot(confint(plot_ppk_model,"Cc", level=0.95),
       ylab="Central compartment") +
    ggplot2::geom_point(data=indiv_obs_plot, na.rm=TRUE)
}
