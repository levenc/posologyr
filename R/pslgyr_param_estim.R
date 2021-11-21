#-------------------------------------------------------------------------
# posologyr: individual dose optimisation using population PK
# Copyright (C) 2021  Cyril Leven
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Affero General Public License as
#  published by the Free Software Foundation, either version 3 of the
#  License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.
#-------------------------------------------------------------------------

#-------------------------------------------------------------------------
#  Adapted from: http://shiny.webpopix.org/mcmc/bayes1/
#  Copyright (C) 2015 Marc Lavielle, Inria Saclay, CeCILL-B
#
#  Modifications:
#   - interfacing with RxODE
#   - deletion of shiny-specific parts
#   - variable names changed to snake_case
#   - square matrix taken as input, not diagonal
#   - functions return values for both etas and theta
#   - adaptive MAP forecasting
#   - inter-occasion variability (IOV)
#-------------------------------------------------------------------------

#' Estimate the prior distribution of population parameters
#'
#' Estimates the prior distribution of population parameters by Monte Carlo
#' simulations
#'
#' @param object A posologyr list, created by the \code{\link{posologyr}}
#'    function.
#' @param n_simul An integer, the number of simulations to be run. For `n_simul
#'   =0`, all ETAs are set to 0.
#' @param return_model A boolean. Returns a RxODE model using the simulated
#'    ETAs if set to `TRUE`.
#'
#' @return If `return_model` is set to `FALSE`, a list of one element: a
#' dataframe `$eta` of the individual values of ETA.
#' If `return_model` is set to `TRUE`, a list of the dataframe of the
#' individual values of ETA, and a RxODE model using the simulated ETAs.
#'
#' @examples
#' # df_patient01: event table for Patient01, following a 30 minutes intravenous
#' # infusion of tobramycin
#' df_patient01 <- data.frame(ID=1,
#'                         TIME=c(0.0,1.0,14.0),
#'                         DV=c(NA,25.0,5.5),
#'                         AMT=c(2000,0,0),
#'                         EVID=c(1,0,0),
#'                         DUR=c(0.5,NA,NA),
#'                         CLCREAT=80,WT=65)
#' # loading a tobramycin model and Patient01's event record
#' patient01_tobra <- posologyr(prior_model=mod_tobramycin_2cpt_fictional,
#'                              dat=df_patient01)
#' # estimate the prior distribution of population parameters
#' poso_simu_pop(patient01_tobra,n_simul=100)
#'
#' @export
poso_simu_pop <- function(object,n_simul=1000,
                          return_model=TRUE){
  validate_priormod(object)
  validate_dat(object$tdm_data)
  no_covariates <- is.null(object$covariates)

  omega      <- object$omega
  ind_eta    <- which(diag(omega)>0)          # only parameters with IIV
  omega_eta  <- omega[ind_eta,ind_eta]
  eta_mat    <- matrix(0,nrow=1,ncol=ncol(omega))

  if (n_simul > 0) {
    eta_mat <- matrix(0,nrow=n_simul,ncol=ncol(omega))
    eta_sim <- mvtnorm::rmvnorm(n_simul,mean=rep(0,ncol(omega_eta)),
                                sigma=omega_eta)
    eta_mat[,ind_eta] <- eta_sim
  }

  eta_df             <- data.frame(eta_mat)
  names(eta_df)      <- attr(omega,"dimnames")[[1]]

  eta_pop            <- list(eta=eta_df)

  # outputs
  if(return_model){
    model_pop         <- object$solved_ppk_model
    theta             <- rbind(object$theta)

    if(no_covariates){
      params <- cbind(theta,eta_df,row.names=NULL)
    } else {
      covar             <- as.data.frame(object$tdm_data[1,object$covariates])
      names(covar)      <- object$covariates
      params <- cbind(theta,eta_df,covar,row.names=NULL)
    }

    if (!is.null(object$pi_matrix)){
      kappa_mat         <- matrix(0,nrow=1,ncol=ncol(omega))
      kappa_df          <- data.frame(kappa_mat)
      names(kappa_df)   <- attr(object$pi_matrix,"dimnames")[[1]]
      params            <- cbind(params,kappa_df)
    }

    model_pop$params  <- params
    eta_pop$model     <- model_pop
  }

  return(eta_pop)
}

#' Estimate the Maximum A Posteriori individual parameters
#'
#' Estimates the Maximum A Posteriori (MAP) individual parameters,
#' also known as Empirical Bayes Estimates (EBE).
#'
#' @param object A posologyr list, created by the \code{\link{posologyr}}
#'    function.
#' @param adapt A boolean. Should the estimation be performed with the
#'    adaptive MAP method (as opposed to the standard MAP)? A column
#'    `AMS` is required in the patient record to define the segments for
#'    the adaptive MAP approach.
#' @param return_model A boolean. Returns a RxODE model using the estimated
#'    ETAs if set to `TRUE`. If `adapt=TRUE`, the model is solved using the
#'    parameters estimated on the last segment.
#' @param return_ofv A boolean. Returns a the Objective Function Value (OFV)
#'    if set to `TRUE`. Always considered `FALSE` if `adapt=TRUE`.
#' @param return_AMS_models A boolean. Returns a RxODE model using the estimated
#'    ETAs for each Adaptive MAP Segment (AMS) if set to `TRUE`. Ignored if
#'    `adapt=FALSE`.
#'
#' @return A named list consisting of one or more of the following elements
#' depending on the input parameters of the function: `$eta` a named vector
#' of the MAP estimates of the individual values of ETA, `$model` an RxODE
#' model using the estimated ETAs, `AMS_models` a list of RxODE models, one for
#' each Adaptive MAP Segment (AMS).
#'
#' @import data.table
#' @examples
#' # df_patient01: event table for Patient01, following a 30 minutes intravenous
#' # infusion of tobramycin
#' df_patient01 <- data.frame(ID=1,
#'                         TIME=c(0.0,1.0,14.0),
#'                         DV=c(NA,25.0,5.5),
#'                         AMT=c(2000,0,0),
#'                         EVID=c(1,0,0),
#'                         DUR=c(0.5,NA,NA),
#'                         CLCREAT=80,WT=65)
#' # loading a tobramycin model and Patient01's event record
#' patient01_tobra <- posologyr(prior_model=mod_tobramycin_2cpt_fictional,
#'                              dat=df_patient01)
#' # estimate the Maximum A Posteriori individual parameters
#' poso_estim_map(patient01_tobra)
#'
#' @export
poso_estim_map <- function(object,adapt=FALSE,return_model=TRUE,return_ofv=FALSE,
                           return_AMS_models=FALSE){
  validate_priormod(object)
  validate_dat(object$tdm_data)
  estim_with_iov <- check_for_iov(object)
  no_covariates  <- is.null(object$covariates)

  dat           <- object$tdm_data
  solved_model  <- object$solved_ppk_model
  omega         <- object$omega
  theta         <- object$theta
  sigma         <- object$sigma
  error_model   <- object$error_model
  interpolation <- object$interpolation

  ind_eta      <- which(diag(omega)>0)          # only parameters with IIV
  omega_eta    <- omega[ind_eta,ind_eta]        # only variances > 0
  solve_omega  <- try(solve(omega_eta))         # inverse of omega_eta

  eta_map      <- diag(omega)*0

  diag_varcovar_matrix <- diag(omega_eta)

  # initialize the list of outputs
  estim_map    <- list(eta=eta_map)

  if(adapt){ #adaptive MAP estimation  doi: 10.1007/s11095-020-02908-7
    if (is.null(dat$AMS)){
      stop("The AMS column is required in the patient record to define the
      segments for adaptive MAP forecasting")
    }

    start_eta       <- diag(omega_eta)*0    # get a named vector of zeroes

    adaptive_output <- adaptive_map(return_AMS_models=return_AMS_models,
                                    dat=dat,
                                    solved_model=solved_model,
                                    theta=theta,
                                    omega=omega,
                                    start_eta=start_eta,
                                    errpred=errpred,
                                    run_model=run_model,
                                    ind_eta=ind_eta,
                                    sigma=sigma,
                                    solve_omega=solve_omega,
                                    omega_dim=omega_dim,
                                    iov_col=iov_col,
                                    pimat=pimat,
                                    eta_map=eta_map,
                                    error_model=error_model,
                                    estim_with_iov=estim_with_iov,
                                    interpolation=interpolation,
                                    adapt=adapt)
    if(return_AMS_models){
    AMS_models <- adaptive_output$AMS_models
    }
    eta_df     <- adaptive_output$eta_df
    eta_map <- unlist(utils::tail(eta_df,1))

    if(!no_covariates){
      covar            <- t(utils::tail(dat[,object$covariates]))
      names(covar)     <- object$covariates
    }
  }
  else{ #standard MAP estimation

    # avoid empty (NULL) arguments for stats::optim()
    omega_dim <- 0
    iov_col   <- 0
    pimat     <- 0
    eta_df    <- 0

    if (estim_with_iov){
      data_iov     <- dat
      pimat        <- object$pi_matrix

      ind_kappa    <- which(diag(pimat)>0)
      pimat_kappa  <- pimat[ind_kappa,ind_kappa]

      omega_dim    <- ncol(omega_eta)
      pimat_dim    <- ncol(pimat_kappa)

      iov_col      <- init_iov_col(dat=dat,pimat=pimat)
      all_the_mat  <- merge_covar_matrices(omega_eta=omega_eta,
                                           omega_dim=omega_dim,
                                           pimat_dim=pimat_dim,
                                           pimat_kappa=pimat_kappa,
                                           dat=dat)

      solve_omega   <- try(solve(all_the_mat))
      diag_varcovar_matrix <- diag(all_the_mat)
    }

    start_eta        <- init_eta(object,estim_with_iov,omega_iov=all_the_mat)
    model_init       <- 0                             # to appease run_model()
    y_obs            <- dat$DV[dat$EVID == 0]         # only observations

    # initial bounds for the optimization
    bfgs_bounds       <- stats::qnorm(25e-3,0,sqrt(diag_varcovar_matrix),
                                      lower.tail = F)

    optim_attempt     <- 1
    max_attempt       <- 40
    one_more_time     <- TRUE

    # create a table to log the estimates after each attempt
    optim_attempt_log        <- matrix(Inf,nrow=max_attempt,
                                       ncol=(1+length(start_eta)))
    optim_attempt_log        <- data.table(optim_attempt_log)

    data.table::setnames(optim_attempt_log,
                         1:(length(start_eta)+1),
                         c("OFV",names(start_eta)),
                         skip_absent=TRUE)

    while(one_more_time & optim_attempt <= max_attempt){
      r <- try(stats::optim(start_eta,errpred,
                        gr=optim_gradient,
                        run_model=run_model,
                        y_obs=y_obs,
                        theta=theta,
                        ind_eta=ind_eta,
                        sigma=sigma,
                        solve_omega=solve_omega,
                        omega=omega,
                        omega_dim=omega_dim,
                        iov_col=iov_col,
                        pimat=pimat,
                        dat=dat,
                        eta_df=eta_df,
                        model_init=model_init,
                        solved_model=solved_model,
                        error_model=error_model,
                        estim_with_iov=estim_with_iov,
                        interpolation=interpolation,
                        adapt=adapt,
                        method="L-BFGS-B",
                        upper=bfgs_bounds,
                        lower=-bfgs_bounds),
               silent=TRUE)

      if(class(r) != 'try-error'){

        # Objection Function Value: OFV
        OFV_current       <- r$value
        best_attempt_ofv  <- min(optim_attempt_log$OFV)
        second_best_ofv   <- min(sort(optim_attempt_log$OFV)[-1])

        # detection of anomalous estimates calling for a new attempt
        stuck_on_bound    <- TRUE %in% (abs(r$par) >= bfgs_bounds)
        all_eta_are_zero  <- !(FALSE %in% (r$par == 0))
        identical_abs_eta <- isTRUE(length(unique(abs(r$par))) < length(r$par))
        sky_high_ofv      <- OFV_current >= 1e10
        not_the_best      <- isTRUE(OFV_current - best_attempt_ofv > 1e-7)
        no_better_than_2nd_best <- isTRUE(second_best_ofv - best_attempt_ofv > 1e-5)

        need_a_new_start  <- isTRUE(optim_attempt == 1|
                                      all_eta_are_zero|
                                      identical_abs_eta|
                                      sky_high_ofv|
                                      not_the_best|
                                      no_better_than_2nd_best)

        # log the detection of the anomalous estimates, OFV, and ETA estimates
        optim_attempt_log[optim_attempt,] <- data.table("OFV"=OFV_current,
                                                        rbind(r$par))

        if(optim_attempt < max_attempt){

          one_more_time <- FALSE

          # check conditions calling for a new attempt at minimizing the OFV

          if(stuck_on_bound){
            one_more_time <- TRUE
            bfgs_bounds   <- bfgs_bounds + 1

          } else if(need_a_new_start){
            one_more_time <- TRUE
            start_eta     <- init_eta(object,estim_with_iov,
                                      omega_iov=all_the_mat)
          }

        } else if(optim_attempt == max_attempt){

          # if all fails, the "less bad" solution is probably
          # the estimation with the lowest OFV

          OFV   <- NULL    # avoid undefined global variables

          r$par <- unlist(optim_attempt_log[OFV==min(OFV),
                                              2:(length(start_eta)+1)][1,])
          r$value <- min(optim_attempt_log[,OFV])

        }
      } else{ # class(r) == 'try-error'
        one_more_time <- TRUE
        start_eta     <- init_eta(object,estim_with_iov,
                                  omega_iov=all_the_mat)
      }

      optim_attempt     <- optim_attempt + 1
    }


    if (estim_with_iov){
      eta_map[ind_eta] <- r$par[1:omega_dim]
      iov_col <- iov_proposition_as_cols(iov_col=iov_col,dat=dat,pimat=pimat,
                                         omega_dim=omega_dim,
                                         eta_estim=r$par)
      data_iov     <- data.frame(dat,iov_col)
      solved_model <- RxODE::rxSolve(solved_model,c(theta,eta_map),data_iov,
                                     covs_interpolation=interpolation)

      estim_map$data   <- data_iov
    } else {
      eta_map[ind_eta] <- r$par
    }

    if(!no_covariates){
      covar            <- t(dat[1,object$covariates]) #results in a matrix
      names(covar)     <- object$covariates
    }
  }

  # list of all outputs
  estim_map$eta      <- eta_map

  if(return_model){
    model_map        <- solved_model
    if(no_covariates){
      model_map$params <- c(theta,eta_map)
    } else {
      model_map$params <- c(theta,eta_map,covar)
    }
    estim_map$model  <- model_map
  }

  if(return_ofv & !adapt){
    estim_map$ofv <- r$value
  }

  if(adapt & return_AMS_models){
    estim_map$AMS_models <- AMS_models
  }

  return(estim_map)
}

# poso_estim_mcmc: distribution of the individual parameters
#
# Copyright (C) 2017-2021 Emmanuelle Comets <emmanuelle.comets@inserm.fr>,
# Audrey Lavenu, Marc Lavielle (authors of the saemix R package)
# Copyright (C) 2021 Cyril Leven
#
# The saemix package is free software; licensed under the terms of the GNU
# General Public License as published by the Free Software Foundation;
# either version 2 of the License, or (at your option) any later version.
#
# (GPLv2+)

#' Estimate the posterior distribution of individual parameters by MCMC
#'
#' Estimates the posterior distribution of individual parameters by Markov
#' Chain Monte Carlo (using a Metropolis-Hastings algorithm)
#'
#' @param object A posologyr list, created by the \code{\link{posologyr}}
#' function.
#' @param return_model A boolean. Returns a RxODE model using the estimated
#'    ETAs if set to `TRUE`.
#' @param burn_in Number of burn-in iterations for the Metropolis-Hastings
#'    algorithm.
#' @param n_iter Total number of iterations (following the burn-in iterations)
#'  for the Metropolis-Hastings algorithm.
#' @param control A list of parameters controlling the Metropolis-Hastings
#' algorithm.
#'
#' @return If `return_model` is set to `FALSE`, a list of one element: a
#' dataframe `$eta` of ETAs from the posterior distribution, estimated by
#' Markov Chain Monte Carlo.
#' If `return_model` is set to `TRUE`, a list of the dataframe of the posterior
#' distribution of ETA, and a RxODE model using the estimated distributions of ETAs.
#'
#' @author Emmanuelle Comets, Audrey Lavenu, Marc Lavielle, Cyril Leven
#'
#' @references Comets  E, Lavenu A, Lavielle M. Parameter estimation in nonlinear
#' mixed effect models using saemix, an R implementation of the SAEM algorithm.
#' Journal of Statistical Software 80, 3 (2017), 1-41.
#'
#' @examples
#' # df_patient01: event table for Patient01, following a 30 minutes intravenous
#' # infusion of tobramycin
#' df_patient01 <- data.frame(ID=1,
#'                         TIME=c(0.0,1.0,14.0),
#'                         DV=c(NA,25.0,5.5),
#'                         AMT=c(2000,0,0),
#'                         EVID=c(1,0,0),
#'                         DUR=c(0.5,NA,NA),
#'                         CLCREAT=80,WT=65)
#' # loading a tobramycin model and Patient01's event record
#' patient01_tobra <- posologyr(prior_model=mod_tobramycin_2cpt_fictional,
#'                              dat=df_patient01)
#' # estimate the posterior distribution of population parameters
#' \donttest{poso_estim_mcmc(patient01_tobra,n_iter=100)}
#'
#' @export
poso_estim_mcmc <- function(object,return_model=TRUE,burn_in=50,
                            n_iter=1000,control=list(n_kernel=c(2,2,2),
                            stepsize_rw=0.4,proba_mcmc=0.3,nb_max=3)){
  validate_priormod(object)
  validate_dat(object$tdm_data)
  no_covariates <- is.null(object$covariates)

  # Update model predictions with a new set of parameters, for all obs-----
  run_model <- function(x,model=solved_model){
    model$params <- x
    return(model$Cc)
  }

  dat          <- object$tdm_data
  solved_model <- object$solved_ppk_model
  omega        <- object$omega
  sigma        <- object$sigma
  error_model  <- object$error_model

  y_obs        <- dat$DV[dat$EVID == 0]     # only observations
  ind_eta      <- which(diag(omega)>0)      # only parameters with IIV
  nb_etas      <- length(ind_eta)
  omega_eta    <- omega[ind_eta,ind_eta]    # only variances > 0
  solve_omega  <- try(solve(omega_eta))     # inverse of omega_eta
  chol_omega   <- chol(omega_eta)
  rw_init      <- 0.5                       #initial variance parameter for kernels
  d_omega      <- diag(omega_eta)*rw_init
  VK           <- rep(c(1:nb_etas),2)
  n_iter       <- n_iter + burn_in

  # Metropolis-Hastings algorithm------------------------------------------
  theta    <- object$theta
  eta      <- diag(omega_eta)*0
  f        <- do.call(run_model,list(c(theta,eta)))
  g        <- error_model(f,sigma)
  U_y      <- sum(0.5 * ((y_obs - f)/g)^2 + log(g))
  U_eta    <- 0.5 * eta %*% solve_omega %*% eta

  eta_mat     <- matrix(0,nrow=n_iter+1,ncol=ncol(omega))
  eta_mat[1,] <- diag(omega)*0

  for (k_iter in 1:n_iter)
  {
    if (control$n_kernel[1] > 0)
    {
      for (u in 1:control$n_kernel[1])
      {
        etac          <- as.vector(chol_omega%*%stats::rnorm(nb_etas))
        names(etac)   <- attr(omega_eta,"dimnames")[[1]]
        f             <- do.call(run_model,list(c(theta,etac)))
        g             <- error_model(f,sigma)
        Uc_y          <- sum(0.5 * ((y_obs - f)/g)^2 + log(g))
        deltu         <- Uc_y - U_y
        if(deltu < (-1) * log(stats::runif(1)))
        {
          eta        <- etac
          U_y        <- Uc_y
        }
      }
    }
    if (control$n_kernel[2] > 0)
    {
      nb_max        <- min(nb_etas,control$nb_max)
      nbc2          <- nt2     <- replicate(nb_etas,0)
      U_eta         <- 0.5 * eta %*% solve_omega %*% eta
      for (u in 1:control$n_kernel[2])
      {
        for (nrs2 in 1:nb_max)
        {
          for (j in 1:nb_etas)
          {
            jr            <- sample(c(1:nb_etas), nrs2)
            jr            <- jr -jr[1] + j
            vk2           <- jr%%nb_etas + 1
            etac          <- eta
            etac[vk2]     <- eta[vk2] + stats::rnorm(nrs2)*d_omega[vk2]
            f             <- do.call(run_model,list(c(theta,etac)))
            g             <- error_model(f,sigma)
            Uc_y          <- sum(0.5 * ((y_obs - f)/g)^2 + log(g))
            Uc_eta        <- 0.5 * etac %*% solve_omega %*% etac
            deltu         <- Uc_y - U_y + Uc_eta - U_eta
            if(deltu < (-1) * log(stats::runif(1)))
            {
              eta         <- etac
              U_y         <- Uc_y
              U_eta       <- Uc_eta
              nbc2[vk2]   <- nbc2[vk2]+1
            }
            nt2[vk2]      <- nt2[vk2] + 1
          }
        }
      }
      d_omega <- d_omega*(1 + control$stepsize_rw*(nbc2/nt2 - control$proba_mcmc))
    }
    if(control$n_kernel[3]>0) {
      nt2          <- nbc2     <- matrix(data=0,nrow=nb_etas,ncol=1)
      nrs2         <- k_iter%%(nb_etas-1)+2
      for (u in 1:control$n_kernel[3]) {
        if(nrs2<nb_etas) {
          vk        <- c(0,sample(c(1:(nb_etas-1)),nrs2-1))
          nb_iter2  <- nb_etas
        } else {
          vk        <- 0:(nb_etas-1)
          nb_iter2  <- 1
        }
        for(k2 in 1:nb_iter2) {
          vk2             <- VK[k2+vk]
          etac            <- eta
          etac[vk2]       <- eta[vk2]+matrix(stats::rnorm(nrs2),
                                             ncol=nrs2)%*%diag(d_omega[vk2])
          f               <- do.call(run_model,list(c(theta,etac)))
          g               <- error_model(f,sigma)
          Uc_y            <- sum(0.5 * ((y_obs - f)/g)^2 + log(g))
          Uc_eta          <- 0.5*rowSums(etac*(etac%*%solve(omega_eta)))
          deltu           <- Uc_y-U_y+Uc_eta-U_eta
          ind             <- which(deltu<(-log(stats::runif(1))))
          eta[ind]        <- etac[ind]
          U_y[ind]        <- Uc_y[ind]
          U_eta[ind]      <- Uc_eta[ind]
          nbc2[vk2]       <- nbc2[vk2]+length(ind)
          nt2[vk2]        <- nt2[vk2]+1
        }
      }
      d_omega <- d_omega*(1+control$stepsize_rw * (nbc2/nt2-control$proba_mcmc))
    }
    eta_mat[k_iter+1,ind_eta]   <- eta
  }
  eta_df_mcmc            <- data.frame(eta_mat[(burn_in+1):n_iter,])
  names(eta_df_mcmc)     <- attr(omega,"dimnames")[[1]]

  estim_mcmc             <- list(eta=eta_df_mcmc)

  if(return_model){
    model_mcmc        <- solved_model
    theta_return      <- rbind(theta)
    if(no_covariates){
      model_mcmc$params <- cbind(theta_return,eta_df_mcmc,row.names=NULL)
    } else {
      covar             <- as.data.frame(dat[1,object$covariates])
      names(covar)      <- object$covariates
      model_mcmc$params <- cbind(theta_return,eta_df_mcmc,covar,row.names=NULL)
    }
    estim_mcmc$model  <- model_mcmc
  }

  return(estim_mcmc)
}

#' Estimate the posterior distribution of individual parameters by SIR
#'
#' Estimates the posterior distribution of individual parameters by
#' Sequential Importance Resampling (SIR)
#'
#' @param object A posologyr list, created by the \code{\link{posologyr}}
#' function.
#' @param n_sample Number of samples from the S-step
#' @param n_resample Number of samples from the R-step
#' @param return_model A boolean. Returns a RxODE model using the estimated
#'    ETAs if set to `TRUE`.
#'
#' @return If `return_model` is set to `FALSE`, a list of one element: a
#' dataframe `$eta` of ETAs from the posterior distribution, estimated by
#' Sequential Importance Resampling.
#' If `return_model` is set to `TRUE`, a list of the dataframe of the posterior
#' distribution of ETA, and a RxODE model using the estimated distributions of ETAs.
#'
#' @import data.table
#' @examples
#' # df_patient01: event table for Patient01, following a 30 minutes intravenous
#' # infusion of tobramycin
#' df_patient01 <- data.frame(ID=1,
#'                         TIME=c(0.0,1.0,14.0),
#'                         DV=c(NA,25.0,5.5),
#'                         AMT=c(2000,0,0),
#'                         EVID=c(1,0,0),
#'                         DUR=c(0.5,NA,NA),
#'                         CLCREAT=80,WT=65)
#' # loading a tobramycin model and Patient01's event record
#' patient01_tobra <- posologyr(prior_model=mod_tobramycin_2cpt_fictional,
#'                              dat=df_patient01)
#' # estimate the posterior distribution of population parameters
#' poso_estim_sir(patient01_tobra,n_sample=1e4,n_resample=1e3)
#'
#' @export
poso_estim_sir <- function(object,n_sample=1e4,n_resample=1e3,return_model=TRUE){
  validate_priormod(object)
  validate_dat(object$tdm_data)
  estim_with_iov <- check_for_iov(object)
  no_covariates  <- is.null(object$covariates)

  dat           <- object$tdm_data
  solved_model  <- object$solved_ppk_model
  omega         <- object$omega
  sigma         <- object$sigma
  error_model   <- object$error_model
  interpolation <- object$interpolation

  y_obs        <- dat$DV[dat$EVID == 0]     # only observations
  ind_eta      <- which(diag(omega)>0)      # only parameters with IIV
  nb_etas      <- length(ind_eta)
  omega_eta    <- omega[ind_eta,ind_eta]    # only variances > 0
  omega_sim    <- omega_eta
  omega_dim    <- ncol(omega_eta)
  solve_omega  <- try(solve(omega_eta))     # inverse of omega_eta

  theta        <- rbind(object$theta)

  if (estim_with_iov){
    pimat        <- object$pi_matrix

    ind_kappa    <- which(diag(pimat)>0)
    pimat_kappa  <- pimat[ind_kappa,ind_kappa]

    pimat_names  <- attr(pimat_kappa,"dimnames")[[1]]

    pimat_dim    <- ncol(pimat_kappa)
    n_occ        <- length(unique(dat$OCC))

    all_the_mat  <- merge_covar_matrices(omega_eta=omega_eta,
                                         omega_dim=omega_dim,
                                         pimat_dim=pimat_dim,
                                         pimat_kappa=pimat_kappa,
                                         dat=dat)
    omega_sim   <- all_the_mat
    solve_omega <- solve(all_the_mat)
  }

  #SIR algorithm
  # doi: 10.1002/psp4.12492; doi: 10.1007/s10928-016-9487-8

  #S-step
  eta_sim       <- mvtnorm::rmvnorm(n_sample,mean=rep(0,ncol(omega_sim)),
                                    sigma=omega_sim)
  eta_df        <- data.frame(eta_sim)
  names(eta_df) <- attr(omega_eta,"dimnames")[[1]]
  eta_dt        <- data.table::data.table(eta_df)

  param_cols    <- attr(omega_eta,"dimnames")[[1]]
  params        <- cbind(ID=1,eta_dt[,param_cols,with=F],theta)

  if (estim_with_iov){
    eta_dt[,ID:=(1:n_sample)]                           # 1:n_samples ID

    dat_dt <- data.table::data.table(dat)
    dat_dt <- dat_dt[rep(dat_dt[,.I],n_sample)]         # one table per sample
    dat_dt[,ID:=rep(1:n_sample,1,each=nrow(dat))]       # 1:n_samples IDs

    # bind random effects to patient data.table
    ID     <- NULL    # avoid undefined global variables
    dat_dt <- dat_dt[eta_dt,on = list(ID = ID), roll = TRUE]

    # everything but ID, OCC and kappas
    tdm_dt            <- data.table::data.table(dat)
    names_tdm_dt_drop <- names(tdm_dt[,!c("ID","OCC")])
    names_tdm_dt_full <- names(tdm_dt)
    drop_cols         <- c(names_tdm_dt_drop,attr(omega_eta,"dimnames")[[1]])

    # reshape IOV to colums
    iov_col <- t(apply(dat_dt[,!drop_cols,with=F],
                       MARGIN=1,
                       FUN=link_kappa_to_occ,
                       pimat_dim=pimat_dim,
                       pimat_names=pimat_names))
    iov_col_dt <- data.table::data.table(iov_col)

    data_iov   <- cbind(dat_dt[,names_tdm_dt_full,with=F],
                        iov_col_dt[,!c("ID","OCC")])

    param_cols <- c("ID",attr(omega,"dimnames")[[1]])
    params     <- cbind(eta_dt[,param_cols,with=F],theta)
  }

  #I-step
  if(estim_with_iov){
    group_index   <- data.frame(cbind(c(1:10),10,n_sample,nrow(dat)))

    loads_omodels <- apply(group_index,
                          pkmodel=object$solved_ppk_model,
                          params=params,
                          dat=data_iov,
                          interpolation=interpolation,
                          MARGIN=1,FUN=solve_by_groups)
    solved_model  <- do.call(rbind,loads_omodels)

    wide_cc  <- tidyr::pivot_wider(solved_model,
                                   id_cols = "id",
                                   names_from = "time",
                                   values_from = "Cc")
  } else {
    solved_model$params  <- cbind(theta,eta_dt,row.names=NULL)
    wide_cc  <- tidyr::pivot_wider(solved_model,
                                   id_cols = "sim.id",
                                   names_from = "time",
                                   values_from = "Cc")
  }

  LL_func  <- function(simu_obs){ #doi: 10.4196/kjpp.2012.16.2.97
    eta_id   <- simu_obs[1]
    eta      <- eta_sim[eta_id,]
    f        <- simu_obs[-1]
    g        <- error_model(f,sigma)
    minus_LL <- 0.5*objective_function(y_obs=y_obs,f=f,g=g,eta=eta,
                                       solve_omega=solve_omega)
    return(-minus_LL)
  }

  lf        <- apply(wide_cc,MARGIN=1,FUN=LL_func)
  lp        <- mvtnorm::dmvnorm(eta_sim,mean=rep(0,ncol(omega_sim)),
                                sigma=omega_sim,log=TRUE)
  md        <- max(lf - lp)
  wt        <- exp(lf - lp - md)
  probs     <- wt/sum(wt)

  #R-step
  indices   <- sample(1:n_sample, size = n_resample, prob = probs,
                      replace = TRUE)
  if (nb_etas > 1) {
    eta_sim <- eta_sim[indices, ]
  }
  else {
    eta_sim <- eta_sim[indices]
  }

  eta_mat           <- matrix(0,nrow=n_resample,ncol=ncol(omega_eta))
  eta_mat[,ind_eta] <- eta_sim[,1:omega_dim]
  eta_df            <- data.frame(eta_mat)
  names(eta_df)     <- attr(omega_eta,"dimnames")[[1]]

  estim_sir         <- list(eta=eta_df)

  if(return_model){
    if(estim_with_iov){
      params_resample   <- cbind(data.frame(ID=1:n_resample),eta_df,theta)

      # return a list of data.tables of resampled IDs
      loads_otables     <- lapply(indices,
                                  data_iov,
                                  FUN=function(indices,dat)
                                    {dat[ID == indices,]})

      # bind the data.tables nicely
      dat_resample      <- do.call(rbind,loads_otables)

      # overwrite IDs to avoid duplicates, and solve the model once again
      dat_resample[,ID:=rep(1:n_resample,1,each=nrow(dat))]
      estim_sir$model   <- RxODE::rxSolve(object$solved_ppk_model,
                                        params_resample,
                                        dat_resample,
                                        covs_interpolation=interpolation)
    } else {
      params_resample   <- cbind(eta_df,theta)
      model_sir         <- solved_model
      if(no_covariates){
        model_sir$params  <- cbind(params_resample,row.names=NULL)
      } else {
        covar             <- as.data.frame(dat[1,object$covariates])
        names(covar)      <- object$covariates
        model_sir$params  <- cbind(params_resample,covar,row.names=NULL)
      }
      estim_sir$model   <- model_sir
    }
  }
  return(estim_sir)
}
