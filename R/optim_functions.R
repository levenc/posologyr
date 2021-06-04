#-------------------------------------------------------------------------
# posologyr: individual dose optimisation using population PK
# Copyright (C) 2021  Cyril Leven
#
#    This program is free software: you can redistribute it and/or modify
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

#' Predict time to a selected trough concentration
#'
#' Predicts the time needed to reach a selected trough concentration
#' (Cmin) given a population pharmacokinetic model, a set of individual
#' parameters, a dose, and a target Cmin.
#'
#' @param object A posologyr list, created by the \code{\link{posologyr}}
#'     function.
#' @param param_map A vector of individual parameters. May be omitted,
#'     in which case the \code{\link{poso_estim_map}} function
#'     will be called.
#' @param from Numeric. Starting time for the simulation of the
#'     individual time-concentration profile. The default value is
#'     0.2
#' @param last_time Numeric. Ending time for the simulation of the
#'     individual time-concentration profile. The default value is
#'     72.
#' @param dose Numeric. Dose administered.
#' @param add_dose Numeric. Additional doses administered at
#'     inter-dose interval after the first dose. Optional.
#' @param interdose_interval Numeric. Time for the inter-dose interval
#'     for multiple dose regimen. Must be provided when add_dose is used.
#' @param duration Numeric. Duration of infusion, for zero-order
#'     administrations.
#' @param target_cmin Numeric. Target trough concentration (Cmin).
#'
#' @return A numeric time to the selected trough concentration, from the
#'     time of administration.
#'
#' @examples
#' # df_patient01: event table for Patient01, following a 30 minutes intravenous
#' # infusion of tobramycin
#' df_patient01 <- data.frame(ID=1,
#'                         TIME=c(0.0,0.5,1.0,14.0),
#'                         DV=c(NA,NA,25.0,5.5),
#'                         AMT=c(1000,-1000,0,0),
#'                         EVID=c(10102,10102,0,0),
#'                         DUR=c(0.5,0.5,NA,NA),
#'                         CLCREAT=80,WT=65)
#' # loading a tobramycin model and Patient01's event record
#' patient01_tobra <- posologyr(prior_model=mod_tobramycin_2cpt_fictional,
#'                                 dat=df_patient01)
#' # predict the time needed to reach a concentration of 2.5 mg/l
#' # after the administration of a 2500 mg dose over a 30 minutes
#' # infusion
#' poso_time_cmin(patient01_tobra,dose=2500,duration=0.5,target_cmin=2.5)
#'
#' @export
poso_time_cmin <- function(object=NULL,param_map=NULL,from=0.2,
                           last_time=72,dose=NULL,add_dose=NULL,
                           interdose_interval=NULL,duration=NULL,
                           target_cmin=NULL){

  if (is.null(param_map)){ #theta_pop + MAP estimates of eta + covariates
      model_map <- poso_estim_map(object,return_model=TRUE)
      param_map <- model_map[[2]]$params
  }
  if (!is.null(add_dose)){
    if (is.null(interdose_interval)){
      stop("interdose_interval is mandatory when add_dose is used.",
           call.=FALSE)
    }
  }

  #compute the individual time-concentration profile
  if (!is.null(add_dose)){
    event_table_cmin <- RxODE::et(amt=dose,dur=duration,
                                   ii=interdose_interval,
                                   addl=add_dose)
    time_last_dose   <- add_dose*interdose_interval
    event_table_cmin$add.sampling(seq(time_last_dose+from,
                                      time_last_dose+last_time,
                                      by=0.1))
  }
  else {
    event_table_cmin <- RxODE::et(amt=dose,dur=duration)
    event_table_cmin$add.sampling(seq(from,last_time,by=0.1))
    time_last_dose   <- 0
  }


  cmin_ppk_model <- RxODE::rxSolve(object=object$ppk_model,
                                   params=param_map,
                                   event_table_cmin)

  time_to_target <-
    min(cmin_ppk_model$time[which(cmin_ppk_model$Cc < target_cmin)]) -
    time_last_dose
  return(time_to_target)
}

#' Estimate the optimal dose for a selected target area under the
#' time-concentration curve (AUC)
#'
#' Estimates the optimal dose for a selected target area under the
#' time-concentration curve (AUC) given a population pharmacokinetic
#' model, a set of individual parameters, and a target AUC.
#'
#' @param object A posologyr list, created by the \code{\link{posologyr}}
#'     function.
#' @param param_map A vector of individual parameters. May be omitted,
#'     in which case the \code{\link{poso_estim_map}} function
#'     will be called.
#' @param time_auc Numeric. Last point in time of the AUC for which the dose
#'     is to be optimized. The AUC is computed from 0 to `time_auc`.
#' @param starting_time Numeric. First point in time of the AUC, for multiple
#' dose regimen. The default is zero.
#' @param starting_dose Numeric. Starting dose for the optimization
#'     algorithm.
#' @param duration Numeric. Duration of infusion, for zero-order
#'     administrations.
#' @param add_dose Numeric. Additional doses administered at
#'     inter-dose interval after the first dose. Optional.
#' @param interdose_interval Numeric. Time for the interdose interval
#'     for multiple dose regimen. Must be provided when add_dose is used.
#' @param target_auc Numeric. Target AUC
#'
#' @return A numeric optimal dose to reach the target AUC.
#'
#' @examples
#' # df_patient01: event table for Patient01, following a 30 minutes intravenous
#' # infusion of tobramycin
#' df_patient01 <- data.frame(ID=1,
#'                         TIME=c(0.0,0.5,1.0,14.0),
#'                         DV=c(NA,NA,25.0,5.5),
#'                         AMT=c(1000,-1000,0,0),
#'                         EVID=c(10102,10102,0,0),
#'                         DUR=c(0.5,0.5,NA,NA),
#'                         CLCREAT=80,WT=65)
#' patient01_tobra <- posologyr(prior_model=mod_tobramycin_2cpt_fictional,
#'                                 dat=df_patient01)
#' # estimate the optimal dose to reach an AUC(0-12h) of 45 h.mg/l
#' poso_dose_auc(patient01_tobra,time_auc=12,target_auc=45)
#'
#' @export
poso_dose_auc <- function(object=NULL,param_map=NULL,time_auc=NULL,
                          starting_time=0,starting_dose=100,
                          interdose_interval=NULL,add_dose=NULL,
                          duration=NULL,target_auc=NULL){

  if (is.null(param_map)){ #theta_pop + MAP estimates of eta + covariates
    model_map <- poso_estim_map(object,return_model=TRUE)
    param_map <- model_map[[2]]$params
  }
  if (!is.null(add_dose)){
    if (is.null(interdose_interval)){
      stop("interdose_interval is mandatory when add_dose is used.",
           call.=FALSE)
    }
    if (starting_time+time_auc>interdose_interval*add_dose){
      stop("The auc time window is outside of the dosing time range:
           starting_time+time_auc>interdose_interval*add_dose.",
           call.=FALSE)
    }
  }


  err_dose <- function(dose,time_auc,starting_time,target_auc,
                       interdose_interval,add_dose,prior_model,
                       duration=duration,param_map){

   #compute the individual time-concentration profile
  if (!is.null(add_dose)){
    event_table_auc <- RxODE::et(amt=dose,dur=duration,
                                ii=interdose_interval,
                                addl=add_dose)
  }
  else {
    event_table_auc <- RxODE::et(amt=dose,dur=duration)
  }
   event_table_auc$add.sampling(starting_time)
   event_table_auc$add.sampling(starting_time+time_auc)

   auc_ppk_model <- RxODE::rxSolve(object=prior_model$ppk_model,
                                   params=param_map,
                                   event_table_auc)

   auc_proposal  <- max(auc_ppk_model$AUC)-min(auc_ppk_model$AUC)

   #return the difference between the computed AUC and the target
   delta_auc<-(target_auc - auc_proposal)^2
   return(delta_auc)
 }

 optim_dose_auc <- optim(starting_dose,err_dose,time_auc=time_auc,
                         starting_time=starting_time,add_dose=add_dose,
                         interdose_interval=interdose_interval,
                         target_auc=target_auc,prior_model=object,
                         duration=duration,param_map=param_map,
                         method="Brent",lower=0,upper=1e5)

 return(optim_dose_auc$par)
}

#' Estimate the optimal dose for a selected target concentration at a
#' selected point in time
#'
#' Estimates the optimal dose for a selected target concentration at a
#' selected point in time given a population pharmacokinetic model, a set
#' of individual parameters, a selected point in time, and a target
#' concentration.
#'
#' @param object A posologyr list, created by the \code{\link{posologyr}}
#'     function.
#' @param param_map A vector of individual parameters. May be omitted,
#'     in which case the \code{\link{poso_estim_map}} function
#'     will be called.
#' @param time_c Numeric. Point in time for which the dose is to be
#'     optimized.
#' @param starting_dose Numeric. Starting dose for the optimization
#'     algorithm.
#' @param add_dose Numeric. Additional doses administered at
#'     inter-dose interval after the first dose. Optional.
#' @param interdose_interval Numeric. Time for the interdose interval
#'     for multiple dose regimen. Must be provided when add_dose is used.
#' @param duration Numeric. Duration of infusion, for zero-order
#'     administrations.
#' @param target_conc Numeric. Target concentration.
#'
#' @return A numeric optimal dose to reach the target concentration
#'     at the selected point in time.
#'
#' @examples
#' # df_patient01: event table for Patient01, following a 30 minutes intravenous
#' # infusion of tobramycin
#' df_patient01 <- data.frame(ID=1,
#'                         TIME=c(0.0,0.5,1.0,14.0),
#'                         DV=c(NA,NA,25.0,5.5),
#'                         AMT=c(1000,-1000,0,0),
#'                         EVID=c(10102,10102,0,0),
#'                         DUR=c(0.5,0.5,NA,NA),
#'                         CLCREAT=80,WT=65)
#' # loading a tobramycin model and Patient01's event record
#' patient01_tobra <- posologyr(prior_model=mod_tobramycin_2cpt_fictional,
#'                                 dat=df_patient01)
#' # estimate the optimal dose to reach a concentration of 80 mg/l
#' # one hour after starting the 30-minutes infusion
#' poso_dose_ctime(patient01_tobra,time_c=1,duration=0.5,target_conc=80)
#'
#' @export
poso_dose_ctime <- function(object=NULL,param_map=NULL,time_c=NULL,
                            starting_dose=100,interdose_interval=NULL,
                            add_dose=NULL,duration=NULL,
                            target_conc=NULL){

  if (is.null(param_map)){ #theta_pop + MAP estimates of eta + covariates
    model_map <- poso_estim_map(object,return_model=TRUE)
    param_map <- model_map[[2]]$params
  }
  if (!is.null(add_dose)){
    if (is.null(interdose_interval)){
      stop("interdose_interval is mandatory when add_dose is used.",
           call.=FALSE)
    }
    if (time_c>(add_dose*interdose_interval)){
      stop("The target time is outside of the dosing time range:
           time_c>(add_dose*interdose_interval).",
           call.=FALSE)
    }
  }

  err_dose <- function(dose,time_c,target_conc,prior_model,
                       add_dose,interdose_interval,
                       duration=duration,param_map){

    #compute the individual time-concentration profile
    if (!is.null(add_dose)){
      event_table_ctime <- RxODE::et(amt=dose,dur=duration,
                                   ii=interdose_interval,
                                   addl=add_dose)
    }
    else {
      event_table_ctime <- RxODE::et(amt=dose,dur=duration)
    }
    event_table_ctime$add.sampling(time_c)

    ctime_ppk_model <- RxODE::rxSolve(object=prior_model$ppk_model,
                                    params=param_map,
                                    event_table_ctime)
    #return the difference between the computed ctime and the target
    delta_ctime <- (target_conc - ctime_ppk_model$Cc)^2
    return(delta_ctime)
  }

  optim_dose_ctime <- optim(starting_dose,err_dose,time_c=time_c,
                          target_conc=target_conc,prior_model=object,
                          add_dose=add_dose,
                          interdose_interval=interdose_interval,
                          duration=duration,param_map=param_map,
                          method="Brent",lower=0, upper=1e5)

  return(optim_dose_ctime$par)
}

#' Estimate the optimal inter-dose interval for a given dose and a
#' selected target trough concentration
#'
#' Estimates the optimal inter-dose interval for a selected target
#' trough concentration (Cmin), given a dose, a population
#' pharmacokinetic model, a set of individual parameters, and a
#' target concentration.
#'
#' @param object A posologyr list, created by the \code{\link{posologyr}}
#'     function.
#' @param param_map A vector of individual parameters. May be omitted,
#'     in which case the \code{\link{poso_estim_map}} function
#'     will be called.
#' @param dose Numeric. The dose given.
#' @param starting_interval Numeric. Starting inter-dose interval for
#'     the optimization algorithm.
#' @param add_dose Numeric. Additional doses administered at
#'     inter-dose interval after the first dose.
#' @param duration Numeric. Duration of infusion, for zero-order
#'     administrations.
#' @param target_cmin Numeric. Target trough concentration (Cmin).
#'
#' @return A numeric inter-dose interval to reach the target trough
#'     concentration before each dosing of a multiple dose regimen.
#'
#' @examples
#' # df_patient01: event table for Patient01, following a 30 minutes intravenous
#' # infusion of tobramycin
#' df_patient01 <- data.frame(ID=1,
#'                         TIME=c(0.0,0.5,1.0,14.0),
#'                         DV=c(NA,NA,25.0,5.5),
#'                         AMT=c(1000,-1000,0,0),
#'                         EVID=c(10102,10102,0,0),
#'                         DUR=c(0.5,0.5,NA,NA),
#'                         CLCREAT=80,WT=65)
#' # loading a tobramycin model and Patient01's event record
#' patient01_tobra <- posologyr(prior_model=mod_tobramycin_2cpt_fictional,
#'                              dat=df_patient01)
#' # estimate the optimal interval to reach a cmin of of 2.5 mg/l
#' # before each administration
#' poso_inter_cmin(patient01_tobra,dose=1500,duration=0.5,target_cmin=2.5)
#'
#' @export
poso_inter_cmin <- function(object=NULL,param_map=NULL,dose=NULL,
                            starting_interval=12,add_dose=10,
                            duration=NULL,target_cmin=NULL){

  if (is.null(param_map)){ #theta_pop + MAP estimates of eta + covariates
    model_map <- poso_estim_map(object,return_model=TRUE)
    param_map <- model_map[[2]]$params
  }

  err_inter <- function(interdose_interval,dose,target_cmin,
                        prior_model,add_dose,duration=duration,
                        param_map){
    #compute the individual time-concentration profile
    event_table_cmin <- RxODE::et(amt=dose,dur=duration,
                                  ii=interdose_interval,
                                  addl=add_dose)
    event_table_cmin$add.sampling(interdose_interval*add_dose-0.1)

    cmin_ppk_model <- RxODE::rxSolve(object=prior_model$ppk_model,
                                     params=param_map,
                                     event_table_cmin)
    #return the difference between the computed cmin and the target
    delta_cmin <- (target_cmin - cmin_ppk_model$Cc)^2
    return(delta_cmin)
  }

  optim_dose_cmin <- optim(starting_interval,err_inter,dose=dose,
                           target_cmin=target_cmin,prior_model=object,
                           add_dose=add_dose,duration=duration,
                           param_map=param_map,method="Brent",
                           lower=0,upper=1e5)

  return(optim_dose_cmin$par)
}
