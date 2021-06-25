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
#   - relative standard error of MAP estimates
#   - adaptive MAP forecasting
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
#'                                 dat=df_patient01)
#' # estimate the prior distribution of population parameters
#' poso_simu_pop(patient01_tobra,n_simul=100)
#'
#' @export
poso_simu_pop <- function(object,n_simul=1000,
                          return_model=TRUE){
  validate_priormod(object)
  validate_dat(object$tdm_data)

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

  if(return_model){
    model_pop         <- object$solved_ppk_model
    theta             <- rbind(object$theta)
    covar             <- object$tdm_data[1,object$covariates]
    names(covar)      <- object$covariates
    model_pop$params  <- cbind(theta,eta_df,covar,row.names=NULL)
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
#'    ETAs if set to `TRUE`.
#' @param return_fim A boolean. Returns the Fisher Information Matrix
#'    (FIM) if set to `TRUE`.
#' @param return_rse A boolean. Returns the relative standard errors
#'    (RSE) of the MAP estimates of ETA if set to `TRUE`.
#'
#' @return A named list consisting of one or more of the following elements
#' depending on the input parameters of the function: `$eta` a named vector
#' of the MAP estimates of the individual values of ETA, `$model` an RxODE
#' model using the estimated ETAs, `$fim` the Fisher information matrix,
#' `$rse` a named vector of RSEs of the MAP estimates of ETAs.
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
#'                                 dat=df_patient01)
#' # estimate the Maximum A Posteriori individual parameters
#' poso_estim_map(patient01_tobra)
#'
#' @export
poso_estim_map <- function(object,adapt=FALSE,return_model=TRUE,
                           return_fim=FALSE,return_rse=FALSE){
  validate_priormod(object)
  validate_dat(object$tdm_data)

  # Update model predictions with a new set of parameters, for all obs-----
  run_model <- function(x,init=model_init,model=solved_model){
    model$params <- x
    if(adapt){
      model$inits  <- init
    }
    return(model$Cc)
    }

  errpred <- function(eta_estim,run_model,y,theta,ind_eta,sigma,solve_omega,
                      adapt=FALSE){
    eta          <- diag(omega)*0
    eta[ind_eta] <- eta_estim

    if(adapt){
      eta        <- eta + eta_df[i,]
      #simulated concentrations with the proposed eta estimates
      f <- do.call(run_model,list(c(theta,eta)))
      g <- error_model(f,sigma)
    }
    else {
      #simulated concentrations with the proposed eta estimates
      f <- do.call(run_model,list(c(theta,eta)))
      g <- error_model(f,sigma)
    }

    #objective function for the Empirical Bayes Estimates
    #doi: 10.4196/kjpp.2012.16.2.97
    U_y   <- sum(0.5 * ((y_obs - f)/g)^2 + log(g))
    #the transpose of a diagonal matrix is itself
    U_eta <- 0.5 * eta_estim %*% solve_omega %*% eta_estim

    optimize_me <- U_y + U_eta
    return(optimize_me)
    }

  dat          <- object$tdm_data
  solved_model <- object$solved_ppk_model
  omega        <- object$omega
  theta        <- object$theta
  sigma        <- object$sigma
  error_model  <- object$error_model

  ind_eta      <- which(diag(omega)>0)          # only parameters with IIV
  omega_eta    <- omega[ind_eta,ind_eta]        # only variances > 0
  solve_omega  <- try(solve(omega_eta))         # inverse of omega_eta
  start_eta    <- diag(omega_eta)*0             # get a named vector of zeroes
  eta_map      <- diag(omega)*0

  if(adapt){ #adaptive MAP estimation  doi: 10.1007/s11095-020-02908-7
    if (is.null(dat$AMS)){
      stop("The AMS column is required in the patient record to define the
      segments for adaptive MAP forecasting")
    }

    segment_id     <- unique(dat$AMS)
    n_segment      <- length(segment_id)
    dat_segment    <- dat

    eta_mat        <- matrix(0,nrow=n_segment+1,ncol=ncol(omega))
    eta_df         <- data.frame(eta_mat)
    names(eta_df ) <- attr(omega,"dimnames")[[1]]

    init_mat       <- matrix(0,nrow=n_segment+1,
                             ncol=length(solved_model$inits))
    init_df        <- data.frame(init_mat)
    init_names     <- names(solved_model$inits)
    names(init_df) <- init_names

    for(i in 1:n_segment){

      if(i>1){ # set TIME == 0 when a segment starts
        dat_segment$TIME  <- dat$TIME -
          utils::tail(which(dat$AMS == segment_id[i-1]),1)$TIME
      }

      dat_segment  <- dat_segment[which(dat_segment$AMS == segment_id[i]),]

      # solved_model for the current segment
      solved_model <- RxODE::rxSolve(object$solved_ppk_model,
                                     c(object$theta,
                                       diag(object$omega)*0),
                                     dat_segment)

      y_obs        <- dat_segment$DV[dat_segment$EVID == 0]

      model_init <- init_df[i,]

      r <- stats::optim(start_eta,errpred,run_model=run_model,y=y_obs,
                        theta=theta,ind_eta=ind_eta,sigma=sigma,
                        solve_omega=solve_omega,hessian=TRUE)

      eta_map[ind_eta]   <- r$par
      eta_df[i+1,]       <- eta_map + eta_df[i,]

      # get the ODE compartment states at the end of the segment
      solved_model$params <- c(object$theta,unlist(eta_df[i+1,]))
      solved_model$inits  <- init_df[i,]
      init_df[i+1,]       <- utils::tail(solved_model[,init_names],1)
    }
    eta_map <- unlist(utils::tail(eta_df,1))

    covar            <- t(dat_segment[1,object$covariates])
    names(covar)     <- object$covariates
  }
  else{ #standard MAP estimation
    y_obs            <- dat$DV[dat$EVID == 0]         # only observations

    r <- stats::optim(start_eta,errpred,run_model=run_model,y=y_obs,theta=theta,
                      ind_eta=ind_eta,sigma=sigma,solve_omega=solve_omega,hessian=TRUE)

    eta_map[ind_eta] <- r$par

    covar            <- t(dat[1,object$covariates]) #results in a matrix
    names(covar)     <- object$covariates
  }

  estim_map          <- list(eta=eta_map)

  if(return_model){
    model_map        <- solved_model
    model_map$params <- c(theta,eta_map,covar)
    estim_map$model  <- model_map
  }
  if(return_fim){
    estim_map$fim    <- r$hessian #the objective function minimized is -LogLikelihood
                                  # hence the hessian is the Fisher Information Matrix
  }
  if(return_rse){
    map_se           <- sqrt(diag(solve(r$hessian))) #the inverse of the fim is the
                                                     # variance-covariance matrix
    map_rse          <- map_se/abs(eta_map[ind_eta])
    estim_map$rse    <- map_rse
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
#'                                 dat=df_patient01)
#' # estimate the posterior distribution of population parameters
#' \donttest{poso_estim_mcmc(patient01_tobra,n_iter=100)}
#'
#' @export
poso_estim_mcmc <- function(object,return_model=TRUE,burn_in=50,
                            n_iter=1000,control=list(n_kernel=c(2,2,2),
                            stepsize_rw=0.4,proba_mcmc=0.3,nb_max=3)){
  validate_priormod(object)
  validate_dat(object$tdm_data)

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
        etac <- as.vector(chol_omega%*%stats::rnorm(nb_etas))
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
      nt2          <- nbc2     <-matrix(data=0,nrow=nb_etas,ncol=1)
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
          etac[vk2]       <- eta[vk2]+matrix(stats::rnorm(nrs2), ncol=nrs2)%*%diag(d_omega[vk2])
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
    covar             <- dat[1,object$covariates]
    names(covar)      <- object$covariates
    model_mcmc$params <- cbind(theta_return,eta_df_mcmc,covar,row.names=NULL)
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
#'                                 dat=df_patient01)
#' # estimate the posterior distribution of population parameters
#' poso_estim_sir(patient01_tobra,n_sample=1e4,n_resample=1e3)
#'
#' @export
poso_estim_sir <- function(object,n_sample=1e5,n_resample=1e3,return_model=TRUE){
  validate_priormod(object)
  validate_dat(object$tdm_data)

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

  theta    <- rbind(object$theta)

  #S-step
  eta_sim  <- mvtnorm::rmvnorm(n_sample,mean=rep(0,ncol(omega_eta)),
                            sigma=omega_eta)
  eta_df        <- data.frame(eta_sim)
  names(eta_df) <- attr(omega_eta,"dimnames")[[1]]

  #I-step
  solved_model$params  <- cbind(theta,eta_df,row.names=NULL)
  wide_cc  <- tidyr::pivot_wider(solved_model,
                                id_cols = "sim.id",
                                names_from = "time",
                                values_from = "Cc")

  LL_func  <- function(simu_obs){
    eta_id   <- simu_obs[1]
    eta      <- eta_sim[eta_id,]
    f        <- simu_obs[-1]
    g        <- error_model(f,sigma)
    U_y      <- sum(0.5 * ((y_obs - f)/g)^2 + log(g))
    U_eta    <- 0.5 * eta %*% solve_omega %*% eta
    minus_LL <- U_y + U_eta
    return(-minus_LL)
  }

  lf        <- apply(wide_cc,MARGIN=1,FUN=LL_func)
  lp        <- mvtnorm::dmvnorm(eta_sim,mean=rep(0,ncol(omega_eta)),
                             sigma=omega_eta,log=TRUE)
  md        <- max(lf - lp)
  wt        <- exp(lf - lp - md)
  probs     <- wt/sum(wt)

  #R-step
  indices   <- sample(1:n_sample, size = n_resample, prob = probs,
                      replace = TRUE)
  if (nb_etas > 1) {
    eta_sim   <- eta_sim[indices, ]
  }
  else {
    eta_sim <- eta_sim[indices]
  }

  eta_mat           <- matrix(0,nrow=n_resample,ncol=ncol(omega))
  eta_mat[,ind_eta] <- eta_sim
  eta_df            <- data.frame(eta_mat)
  names(eta_df)     <- attr(omega,"dimnames")[[1]]

  estim_sir         <- list(eta=eta_df)

  if(return_model){
    model_sir         <- solved_model
    theta_return      <- rbind(theta)
    covar             <- dat[1,object$covariates]
    names(covar)      <- object$covariates
    model_sir$params  <- cbind(theta_return,eta_df,covar,row.names=NULL)
    estim_sir$model   <- model_sir
  }

  return(estim_sir)
}
