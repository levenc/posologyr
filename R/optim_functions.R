#' Predict time to a selected trough concentration
#'
#' Predicts the time needed to reach a selected trough concentration
#' (Cmin) given a population pharmacokinetic model, a set of individual
#' parameters, a dose, and a target Cmin.
#'
#' @param solved_model An \code{\link[RxODE]{rxSolve}} solve object,
#'     created with the prior RxODE structural population pharmacokinetics
#'     model and the prior typical values of the population parameters
#'     from the `prior_model`, using `dat` as the event record. May be omitted
#'     if a vector of individual parameters `param_psi_map` is provided
#' @param prior_model A posologyr prior population pharmacokinetics model,
#'    a list of seven elements (see 'Details' for the description of the
#'    object)
#' @param dat Dataframe. An individual subject dataset following the
#'     structure of NONMEM/RxODE event records
#' @param param_psi_map A vector of individual parameters. May be omitted
#'     if a `solved_model` and an individual event record `dat` are
#'     provided, in which case the \code{\link{poso_estim_map}} function
#'     will be called
#' @param covar a named vector of the individual covariates. May be
#'     omitted if an individual event record `dat` is provided
#' @param from a numeric starting time for the simulation of the
#'     individual time-concentration profile. The default value is
#'     0.2
#' @param last_time a numeric ending time for the simulation of the
#'     individual time-concentration profile. The default value is
#'     72
#' @param dose a numeric dose administered
#' @param duration a numeric duration of infusion, for zero-order
#'     administrations
#' @param target_cmin a numeric target trough concentration
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
#' @return A numeric time to the selected trough concentration, from the
#'     time of administration
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
#' # predict the time needed to reach a concentration of 2.5 mg/l
#' # after the administration of a 2500 mg dose
#' poso_time_cmin(dose=2500,target_cmin=2.5)
#'
#' @export
poso_time_cmin <- function(solved_model=solved_ppk_model,
                           prior_model=prior_ppk_model,dat=dat_posologyr,
                           param_psi_map=NULL,covar=NULL,from=0.2,
                           last_time=72,dose=NULL,duration=NULL,target_cmin=NULL){

  if (is.null(param_psi_map)){#MAP estimates of the individual parameters
    if (!is.null(solved_model)){
      param_psi_map <- poso_estim_map(solved_model,prior_model,dat)
    } else {
      stop("Either param_psi_map or solved_model is needed for this function to work",
           call. = FALSE)
    }
  }

  if (is.null(covar)){
    #get covariates from the prior model and the dataset
    if (!is.null(data)){
      covar <- t(dat[1,prior_model$covariates]) #results in a matrix
      names(covar) <- prior_model$covariates
    } else {
      stop("Either covar or dat is needed for this function to work",
           call. = FALSE)
    }
  }

  #compute the individual time-concentration profile
  event_table_cmin <- RxODE::et(time=0,amt=dose,dur=duration)
  event_table_cmin$add.sampling(seq(from,last_time,by=0.1))

  cmin_ppk_model <- RxODE::rxSolve(object=prior_model$ppk_model,
                                   params=c(param_psi_map,covar),
                                   event_table_cmin)

  time_to_target <-
    min(cmin_ppk_model$time[which(cmin_ppk_model$Cc < target_cmin)])
  return(time_to_target)
}

#' Estimate the optimal dose for a selected target area under the
#' time-concentration curve (AUC)
#'
#' Estimates the optimal dose for a selected target area under the
#' time-concentration curve (AUC) given a population pharmacokinetic
#' model, a set of individual parameters, and a target AUC. Might be
#' useful for non-linear pharmacokinetics.
#'
#' @param solved_model An \code{\link[RxODE]{rxSolve}} solve object, created
#'     with the prior RxODE structural population pharmacokinetics model and
#'     the prior typical values of the population parameters from the
#'     `prior_model`, using `dat` as the event record. May be omitted
#'     if a vector of individual parameters `param_psi_map` is provided
#' @param prior_model A posologyr prior population pharmacokinetics model,
#'    a list of seven elements (see 'Details' for the description of the
#'    object)
#' @param dat Dataframe. An individual subject dataset following the
#'     structure of NONMEM/RxODE event records
#' @param param_psi_map A vector of individual parameters. May be omitted
#'     if a `solved_model` and an individual event record `dat` are
#'     provided, in which case the \code{\link{poso_estim_map}} function will be
#'     called
#' @param covar a named vector of the individual covariates. May be
#'     omitted if an individual event record `dat` is provided
#' @param time_auc a numeric last point in time of the AUC for which the dose
#'     is to be optimized. The AUC is computed from 0 to `time_auc`
#' @param starting_dose numeric starting dose for the optimization
#'     algorithm
#' @param target_auc a numeric target AUC
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
#' @return A numeric optimal dose to reach the target AUC
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
#' # estimate the optimal dose to reach an AUC(0-12h) of 45 h.mg/l
#' poso_dose_auc(time_auc=12,target_auc=45)
#'
#' @export
poso_dose_auc <- function(solved_model=solved_ppk_model,
                          prior_model=prior_ppk_model,dat=dat_posologyr,
                          param_psi_map=NULL,covar=NULL,time_auc=NULL,
                          starting_dose=100,target_auc=NULL){

  if (is.null(param_psi_map)){#MAP estimates of the individual parameters
    if (!is.null(solved_model)){
      param_psi_map <- poso_estim_map(solved_model,prior_model,dat)
    } else {
      stop("Either param_psi_map or solved_model is needed for this function to work",
           call. = FALSE)
    }
  }

  if (is.null(covar)){
    #get covariates from the prior model and the dataset
    if (!is.null(dat)){
      covar <- t(dat[1,prior_model$covariates]) #results in a matrix
      names(covar) <- prior_model$covariates
    } else {
      stop("Either covar or dat is needed for this function to work",
           call. = FALSE)
    }
  }

 err_dose <- function(dose,time_auc,target_auc,prior_model,
                      param_psi_map,covar){

   #compute the individual time-concentration profile
   event_table_auc <- RxODE::et(time=0,amt=dose)
   event_table_auc$add.sampling(time_auc)

   auc_ppk_model <- RxODE::rxSolve(object=prior_model$ppk_model,
                                    params=c(param_psi_map,covar),
                                    event_table_auc)
   #return the difference between the computed AUC and the target
   delta_auc = abs(target_auc - max(auc_ppk_model$auc))
   return(delta_auc)
 }

 optim_dose_auc <- optim(starting_dose,err_dose,time_auc=time_auc,
                         target_auc=target_auc,prior_model=prior_model,
                         param_psi_map=param_psi_map,covar=covar,
                         method="Brent",lower=0, upper=1e5)

 return(optim_dose_auc$par)
}

#' Estimate the optimal dose for a selected target concentration at a
#' selected time point
#'
#' Estimates the optimal dose for a selected target concentration at a
#' selected time point given a population pharmacokinetic model, a set
#' of individual parameters, a selected point in time, and a target
#' concentration.
#'
#' @param solved_model An \code{\link[RxODE]{rxSolve}} solve object, created
#'     with the prior RxODE structural population pharmacokinetics model and
#'     the prior typical values of the population parameters from the
#'     `prior_model`, using `dat` as the event record. May be omitted
#'     if a vector of individual parameters `param_psi_map` is provided
#' @param prior_model A posologyr prior population pharmacokinetics model,
#'    a list of seven elements (see 'Details' for the description of the
#'    object)
#' @param dat Dataframe. An individual subject dataset following the
#'     structure of NONMEM/RxODE event records
#' @param param_psi_map A vector of individual parameters. May be omitted
#'     if a `solved_model` and an individual event record `dat` are
#'     provided, in which case the \code{\link{poso_estim_map}} function will be
#'     called
#' @param covar a named vector of the individual covariates. May be
#'     omitted if an individual event record `dat` is provided
#' @param time_c a numeric point in time for which the dose is to be
#'     optimized.
#' @param starting_dose numeric starting dose for the optimization
#'     algorithm
#' @param duration a numeric duration of infusion, for zero-order
#'     administrations
#' @param target_ctime a numeric target concentration
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
#' @return A numeric optimal dose to reach the target concentration
#'     at the selected point in time
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
#' # estimate the optimal dose to reach a concentration of 80 mg/l
#' # one hour after starting the infusion
#' poso_dose_ctime(time_c=1,target_ctime=80)
#'
#' @export
poso_dose_ctime <- function(solved_model=solved_ppk_model,
                          prior_model=prior_ppk_model,dat=dat_posologyr,
                          param_psi_map=NULL,covar=NULL,time_c=NULL,
                          starting_dose=100,duration=NULL,target_ctime=NULL){

  if (is.null(param_psi_map)){#MAP estimates of the individual parameters
    if (!is.null(solved_model)){
      param_psi_map <- poso_estim_map(solved_model,prior_model,dat)
    } else {
      stop("Either param_psi_map or solved_model is needed for this function to work",
           call. = FALSE)
    }
  }
  if (is.null(covar)){
    #get covariates from the prior model and the dataset
    if (!is.null(data)){
      covar <- t(dat[1,prior_model$covariates]) #results in a matrix
      names(covar) <- prior_model$covariates
    } else {
      stop("Either covar or dat is needed for this function to work",
           call. = FALSE)
    }
  }

  err_dose <- function(dose,time_c,target_ctime,prior_model,
                       duration=duration,param_psi_map,covar){

    #compute the individual time-concentration profile
    event_table_ctime <- RxODE::et(time=0,amt=dose,dur=duration)
    event_table_ctime$add.sampling(time_c)

    ctime_ppk_model <- RxODE::rxSolve(object=prior_model$ppk_model,
                                    params=c(param_psi_map,covar),
                                    event_table_ctime)
    #return the difference between the computed ctime and the target
    delta_ctime = abs(target_ctime - ctime_ppk_model$Cc)
    return(delta_ctime)
  }

  optim_dose_ctime <- optim(starting_dose,err_dose,time_c=time_c,
                          target_ctime=target_ctime,prior_model=prior_model,
                          duration=duration,param_psi_map=param_psi_map,
                          covar=covar,method="Brent",lower=0, upper=1e5)

  return(optim_dose_ctime$par)
}
