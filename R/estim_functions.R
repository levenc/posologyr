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
poso_simu_pop <- function(prior_model=prior_ppk_model,n_simul=1000){
  #muphi      <- log(prior_model$pk_prior$reference)
  Omega      <- prior_model$pk_prior$Omega
  eta_mat    <- matrix(0,nrow=n_simul,ncol=ncol(Omega))

  for (k in (1:n_simul)){
    #mvrnorm simule des valeurs de eta mais l'addition au log de psi
    # suppose un modèle individuel exponentiel : trop restrictif
    # simuler uniquement les valeurs de eta, avec mvrnorm par
    # exemple, mais renvoyer les valeurs de eta en sortie de la
    # fonction
    #phi      <- muphi + MASS::mvrnorm(1,mu=rep(0,length(muphi)),
    #                             Sigma=Omega)
    eta_sim  <- MASS::mvrnorm(1,mu=rep(0,ncol(Omega)),
                             Sigma=Omega)
    #faster than asking mvrnorm for n_simul samples
    eta_mat[k,]  <- eta_sim
  }
  eta_df    <- data.frame(eta_mat)
  names(eta_df) <- attr(Omega,"dimnames")[[1]]
  return(eta_df)
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
                              dat=dat_posologyr)
{
  # Update model predictions with a new set of parameters, for all obs-----
  run_model <- function(x,model=solved_model){
    model$params <- x
    return(model$Cc)
    }

  errpred <- function(eta_estim,run_model,y,psi,ind_eta,xi,solve_omega){
    #psi <- exp(muphi) #modifier la fonction pour mettre eta en entrée, et pas phi
    eta          <- diag(prior$Omega)*0
    eta[ind_eta] <- eta_estim
    #psi[ind_eta] <- exp(phi) # même principe mais pour eta

    #eta = log(proposed ind. parameter estimates)-log(prior pop. estimates)
    #eta <- phi-muphi[ind_eta] #inutile si eta en entrée ?

    #simulated concentrations with the proposed ind. parameters estimates
    f <- do.call(run_model,list(c(psi,eta))) # compléter la liste pour ajouter eta + theta (=reference)
    g <- error_model(f,xi)            # pas de changement

    #http://sia.webpopix.org/nlme.html#estimation-of-the-individual-parameters
    #doi:10.1006/jbin.2001.1033
    U_y <- sum(0.5 * ((y_obs - f)/g)^2 + log(g))
    #the transpose of a diagonal matrix is itself
    U_eta <- 0.5 * eta_estim %*% solve_omega %*% eta_estim

    optimize_me <- U_y + U_eta
    return(optimize_me)
    }

  model_map   <- solved_model
  prior       <- prior_model$pk_prior
  xi          <- prior_model$xi
  error_model <- prior_model$error_model

  y_obs       <- dat$DV[dat$EVID == 0]         # only observations
  ind_eta     <- which(diag(prior$Omega)>0)    # only parameters with IIV
  omega_eta   <- prior$Omega[ind_eta,ind_eta]  # only variances > 0
  solve_omega <- try(solve(omega_eta))         # inverse of omega_eta
  #muphi       <- log(prior$psi)                # log prior population estimates
  psi         <- prior$psi
  #l0         <- log(prior$reference[ind_eta]) # starting value for optim # changer pour eta, soit sqrt(diag(Omega))
  start_eta   <- diag(omega_eta)*0         # mais uniquement pour [ind_eta]

  r <- optim(start_eta,errpred,run_model=run_model,y=y_obs,psi=psi,#muphi=muphi,
             ind_eta=ind_eta,xi=xi,solve_omega=solve_omega,hessian=TRUE)

  #ind_est_map <- prior$reference
  #ind_est_map[ind_eta] <- exp(r$par)
  eta_map          <- diag(prior$Omega)*0
  eta_map[ind_eta] <- r$par
  model_map$params <- c(psi,eta_map)
  return(model_map)
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
                            control=list(n_iter=100,n_kernel=c(2,2),
                            stepsize_rw=0.4,proba_mcmc=0.4,nb_max=3)){
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
  mean_phi <- log(prior$reference)
  psi      <- prior$reference
  phi      <- log(psi)
  f        <- do.call(run_model,list(psi))
  g        <- error_model(f,xi)
  U_y      <- sum(0.5 * ((y_obs - f)/g)^2 + log(g))
  eta      <- phi[ind_eta] - mean_phi[ind_eta]
  phic     <- phi

  PSI      <- matrix(0,nrow=control$n_iter+1,ncol=length(psi))
  PSI[1,]  <- psi

  for (k_iter in 1:control$n_iter)
  {
    if (control$n_kernel[1] > 0)
    {
      for (u in 1:control$n_kernel[1])
      {
        etac          <- as.vector(MASS::mvrnorm(1,mu=rep(0,nb_etas),
                                                 Sigma=omega_eta))
        phic[ind_eta] <- mean_phi[ind_eta] + etac
        psic          <- exp(phic)
        f             <- do.call(run_model,list(psic))
        g             <- error_model(f,xi)
        Uc_y          <- sum(0.5 * ((y_obs - f)/g)^2 + log(g))
        deltu         <- Uc_y - U_y
        if(deltu < (-1) * log(runif(1)))
        {
          eta        <- etac
          phi        <- phic
          U_y        <- Uc_y
        }
      }
    }
    if (control$n_kernel[2] > 0)
    {
      nb_max        <- min(nb_etas,control$nb_max)
      nbc2<-nt2     <- replicate(nb_etas,0)
      U_eta         <- 0.5 * eta %*% solve_omega %*% eta
      for (u in 1:control$n_kernel[2])
      {
        for (nrs2 in 1:nb_max)
        {
          for (j in 1:nb_etas)
          {
            jr            <-  sample(c(1:nb_etas), nrs2)
            jr            <- jr -jr[1] + j
            vk2           <- jr%%nb_etas + 1
            etac          <- eta
            etac[vk2]     <- eta[vk2] + rnorm(nrs2)*d_omega[vk2]
            phic[ind_eta] <- mean_phi[ind_eta] + etac
            psic          <- exp(phic)
            f             <- do.call(run_model,list(psic))
            g             <- error_model(f,xi)
            Uc_y          <- sum(0.5 * ((y_obs - f)/g)^2 + log(g))
            Uc_eta        <- 0.5 * etac %*% solve_omega %*% etac
            deltu         <- Uc_y - U_y + Uc_eta - U_eta
            if(deltu < (-1) * log(runif(1)))
            {
              eta         <- etac
              U_y         <- Uc_y
              U_eta       <- Uc_eta
              phi[ind_eta]<-phic[ind_eta]
              nbc2[vk2]   <- nbc2[vk2]+1
            }
            nt2[vk2]      <- nt2[vk2] + 1
          }
        }
      }
      d_omega <- d_omega*(1 + control$stepsize_rw*(nbc2/nt2 - control$proba_mcmc))
    }
    PSI[k_iter+1,]       <- exp(phi)
  }
  ind_est_mcmc           <- data.frame(PSI)
  names(ind_est_mcmc)    <- prior$name
  return(ind_est_mcmc)
}
