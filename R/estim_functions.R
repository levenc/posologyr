#-------------------------------------------------------------------------
#  Adapted from: http://shiny.webpopix.org/mcmc/bayes1/
#  Marc Lavielle, Inria Saclay (June 28th, 2015) CeCILL-B
#
#  Modifications:
#   - interfacing with RxODE
#   - deletion of shiny-specific parts
#   - variable names changed to snake_case
#   - square matrix taken as input, not diagonal
#-------------------------------------------------------------------------

#' Estimate the prior distribution of population parameters
#'
#' Estimates the prior distribution of population parameters by Monte Carlo
#' simulations
#'
#' @param prior_model A posologyr prior population pharmacokinetics model, a
#'    list of seven elements (see 'Details' for the description of the
#'    object)
#' @param n_simul An integer, the number of simulations to be run
#'
#' @details
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
#' @return A list of 2: a dataframe of the simulated population
#'     pharmacokinetic parameters, and a vector of the prior typical
#'     values of the population pharmacokinetic parameters
#'
#' @examples
#' # df_michel: event table for Michel, following a 30 minutes intravenous
#' # infusion of tobramycin
#' df_michel <- data.frame(ID=1,
#'                         TIME=c(0.0,0.5,1.0,14.0),
#'                         DV=c(NA,NA,25.0,5.5),
#'                         AMT=c(1000,-1000,0,0),
#'                         EVID=c(10102,10102,0,0),
#'                         DUR=c(0.5,0.5,NA,NA),
#'                         CLCREAT=80,WT=65)
#' # loading a tobramycin model and Michel's event record
#' load_ppk_model(prior_model=mod_tobramycin_2cpt_fictional,dat=df_michel)
#' # estimate the prior distribution of population parameters
#' poso_simu_pop()
#'
#' @export
poso_simu_pop <- function(solved_model=solved_ppk_model,
                          prior_model=prior_ppk_model,
                          dat=dat_posologyr,n_simul=1000,
                          return_model = TRUE){

  Omega      <- prior_model$pk_prior$Omega
  eta_mat    <- matrix(0,nrow=n_simul,ncol=ncol(Omega))

  for (k in (1:n_simul)){
    eta_sim     <- MASS::mvrnorm(1,mu=rep(0,ncol(Omega)),
                             Sigma=Omega)
    #faster than asking mvrnorm for n_simul samples
    eta_mat[k,] <- eta_sim
  }

  eta_df        <- data.frame(eta_mat)
  names(eta_df) <- attr(Omega,"dimnames")[[1]]

  if(return_model){
    model_pop         <- solved_model
    psi               <- rbind(prior_model$pk_prior$psi)
    covar             <- dat[1,prior_model$covariates]
    names(covar)      <- prior_model$covariates
    model_pop$params  <- cbind(psi,eta_df,covar,row.names = NULL)
    eta_pop           <- list(eta_df,model_pop)
  } else {
    eta_pop           <- eta_df
  }

  return(eta_pop)
}

#' Estimate the Maximum A Posteriori individual parameters
#'
#' Estimates the Maximum A Posteriori (MAP) individual parameters,
#' also known as Empirical Bayes Estimates (EBE).
#'
#' @param solved_model An \code{\link[RxODE]{rxSolve}} solve object, created
#'     with the prior RxODE structural population pharmacokinetics model and the
#'     prior typical values of the population parameters from the
#'     `prior_model`, using `dat` as the event record
#' @param prior_model A posologyr prior population pharmacokinetics model, a
#'    list of seven elements (see 'Details' for the description of the
#'    object)
#' @param dat Dataframe. An individual subject dataset following the
#'     structure of NONMEM/RxODE event records
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
#' @return A named vector of the MAP estimates of the individual parameters
#'
#' @examples
#' # df_michel: event table for Michel, following a 30 minutes intravenous
#' # infusion of tobramycin
#' df_michel <- data.frame(ID=1,
#'                         TIME=c(0.0,0.5,1.0,14.0),
#'                         DV=c(NA,NA,25.0,5.5),
#'                         AMT=c(1000,-1000,0,0),
#'                         EVID=c(10102,10102,0,0),
#'                         DUR=c(0.5,0.5,NA,NA),
#'                         CLCREAT=80,WT=65)
#' # loading a tobramycin model and Michel's event record
#' load_ppk_model(prior_model=mod_tobramycin_2cpt_fictional,dat=df_michel)
#' # estimate the Maximum A Posteriori individual parameters
#' poso_estim_map()
#'
#' @export
poso_estim_map <- function(solved_model=solved_ppk_model,
                              prior_model=prior_ppk_model,
                              dat=dat_posologyr,return_model = TRUE)
{
  # Update model predictions with a new set of parameters, for all obs-----
  run_model <- function(x,model=solved_model){
    model$params <- x
    return(model$Cc)
    }

  errpred <- function(eta_estim,run_model,y,psi,ind_eta,xi,solve_omega){
    eta          <- diag(prior$Omega)*0
    eta[ind_eta] <- eta_estim

    #simulated concentrations with the proposed eta estimates
    f <- do.call(run_model,list(c(psi,eta)))
    g <- error_model(f,xi)

    #http://sia.webpopix.org/nlme.html#estimation-of-the-individual-parameters
    #doi:10.1006/jbin.2001.1033
    U_y <- sum(0.5 * ((y_obs - f)/g)^2 + log(g))
    #the transpose of a diagonal matrix is itself
    U_eta <- 0.5 * eta_estim %*% solve_omega %*% eta_estim

    optimize_me <- U_y + U_eta
    return(optimize_me)
    }

  prior       <- prior_model$pk_prior
  xi          <- prior_model$xi
  error_model <- prior_model$error_model

  y_obs       <- dat$DV[dat$EVID == 0]         # only observations
  ind_eta     <- which(diag(prior$Omega)>0)    # only parameters with IIV
  omega_eta   <- prior$Omega[ind_eta,ind_eta]  # only variances > 0
  solve_omega <- try(solve(omega_eta))         # inverse of omega_eta
  psi         <- prior$psi
  start_eta   <- diag(omega_eta)*0             # get a named vector of zeroes

  r <- optim(start_eta,errpred,run_model=run_model,y=y_obs,psi=psi,
             ind_eta=ind_eta,xi=xi,solve_omega=solve_omega,hessian=TRUE)

  eta_map            <- diag(prior$Omega)*0
  eta_map[ind_eta]   <- r$par

  if(return_model){
    model_map        <- solved_model
    covar            <- t(dat[1,prior_model$covariates]) #results in a matrix
    names(covar)     <- prior_model$covariates
    model_map$params <- c(psi,eta_map,covar)
    estim_map        <- list(eta_map,model_map)
  } else {
    estim_map        <- eta_map
  }

  return(estim_map)
}

#' Estimate the posterior distribution of population parameters
#'
#' Estimates the posterior distribution of population parameters by Markov
#' Chain Monte Carlo (using a Metropolis-Hastings algorithm)
#'
#' @param solved_model An \code{\link[RxODE]{rxSolve}} solve object, created
#'     with the prior RxODE structural population pharmacokinetics model and the
#'     prior typical values of the population parameters from the
#'     `prior_model`, using `dat` as the event record
#' @param prior_model A posologyr prior population pharmacokinetics model, a
#'    list of seven elements (see 'Details' for the description of the
#'    object)
#' @param dat Dataframe. An individual subject dataset following the
#'     structure of NONMEM/RxODE event records
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
#' @return A dataframe of parameters from the posterior distribution,
#'     estimated by Markov Chain Monte Carlo
#'
#' @examples
#' # df_michel: event table for Michel, following a 30 minutes intravenous
#' # infusion of tobramycin
#' df_michel <- data.frame(ID=1,
#'                         TIME=c(0.0,0.5,1.0,14.0),
#'                         DV=c(NA,NA,25.0,5.5),
#'                         AMT=c(1000,-1000,0,0),
#'                         EVID=c(10102,10102,0,0),
#'                         DUR=c(0.5,0.5,NA,NA),
#'                         CLCREAT=80,WT=65)
#' # loading a tobramycin model and Michel's event record
#' load_ppk_model(prior_model=mod_tobramycin_2cpt_fictional,dat=df_michel)
#' # estimate the posterior distribution of population parameters
#' poso_estim_mcmc()
#'
#' @export
poso_estim_mcmc <- function(solved_model=solved_ppk_model,
                            prior_model=prior_ppk_model,dat=dat_posologyr,
                            return_model = TRUE,burn_in=20,n_iter=219,
                            control=list(n_kernel=c(2,2),stepsize_rw=0.4,
                            proba_mcmc=0.4,nb_max=3)){
  # Update model predictions with a new set of parameters, for all obs-----
  run_model <- function(x,model=solved_model){
    model$params <- x
    return(model$Cc)
  }

  prior       <- prior_model$pk_prior
  xi          <- prior_model$xi
  error_model <- prior_model$error_model

  y_obs       <- dat$DV[dat$EVID == 0]        # only observations
  ind_eta     <- which(diag(prior$Omega)>0)   # only parameters with IIV
  nb_etas     <- length(ind_eta)
  omega_eta   <- prior$Omega[ind_eta,ind_eta] # only variances > 0
  solve_omega <- try(solve(omega_eta))        # inverse of omega_eta
  d_omega     <- diag(omega_eta)*0.3

  # Metropolis-Hastings algorithm------------------------------------------
  psi      <- prior$psi
  eta      <- diag(omega_eta)*0
  f        <- do.call(run_model,list(c(psi,eta)))
  g        <- error_model(f,xi)
  U_y      <- sum(0.5 * ((y_obs - f)/g)^2 + log(g))

  eta_mat     <- matrix(0,nrow=n_iter+1,ncol=ncol(prior$Omega))
  eta_mat[1,] <- diag(prior$Omega)*0

  for (k_iter in 1:n_iter)
  {
    if (control$n_kernel[1] > 0)
    {
      for (u in 1:control$n_kernel[1])
      {
        etac          <- MASS::mvrnorm(1,mu=rep(0,nb_etas),
                                                 Sigma=omega_eta)
        names(etac)   <- attr(omega_eta,"dimnames")[[1]]
        f             <- do.call(run_model,list(c(psi,etac)))
        g             <- error_model(f,xi)
        Uc_y          <- sum(0.5 * ((y_obs - f)/g)^2 + log(g))
        deltu         <- Uc_y - U_y
        if(deltu < (-1) * log(runif(1)))
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
            etac[vk2]     <- eta[vk2] + rnorm(nrs2)*d_omega[vk2]
            f             <- do.call(run_model,list(c(psi,etac)))
            g             <- error_model(f,xi)
            Uc_y          <- sum(0.5 * ((y_obs - f)/g)^2 + log(g))
            Uc_eta        <- 0.5 * etac %*% solve_omega %*% etac
            deltu         <- Uc_y - U_y + Uc_eta - U_eta
            if(deltu < (-1) * log(runif(1)))
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
    eta_mat[k_iter+1,ind_eta]   <- eta
  }
  eta_df_mcmc            <- data.frame(eta_mat[burn_in:n_iter,])
  names(eta_df_mcmc)     <- attr(prior$Omega,"dimnames")[[1]]

  if(return_model){
    model_mcmc        <- solved_model
    psi_return        <- rbind(psi)
    covar             <- dat[1,prior_model$covariates]
    names(covar)      <- prior_model$covariates
    model_mcmc$params <- cbind(psi_return,eta_df_mcmc,covar,row.names = NULL)
    estim_mcmc           <- list(eta_df_mcmc,model_mcmc)
  } else {
    estim_mcmc           <- eta_df_mcmc
  }

  return(estim_mcmc)
}
