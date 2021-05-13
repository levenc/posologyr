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
#' @param from a numeric starting time for the simulation of the
#'     individual time-concentration profile. The default value is
#'     0.2
#' @param last_time a numeric ending time for the simulation of the
#'     individual time-concentration profile. The default value is
#'     72.
#' @param dose a numeric dose administered.
#' @param duration a numeric duration of infusion, for zero-order
#'     administrations.
#' @param target_cmin a numeric target trough concentration.
#'
#' @return A numeric time to the selected trough concentration, from the
#'     time of administration
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
#' # after the administration of a 2500 mg dose
#' poso_time_cmin(patient01_tobra,dose=2500,target_cmin=2.5)
#'
#' @export
poso_time_cmin <- function(object=NULL,param_map=NULL,from=0.2,
                           last_time=72,dose=NULL,duration=NULL,
                           target_cmin=NULL){

  if (is.null(param_map)){ #theta_pop + MAP estimates of eta + covariates
      model_map <- poso_estim_map(object,return_model = TRUE)
      param_map <- model_map[[2]]$params
  }

  #compute the individual time-concentration profile
  event_table_cmin <- RxODE::et(time=0,amt=dose,dur=duration)
  event_table_cmin$add.sampling(seq(from,last_time,by=0.1))

  cmin_ppk_model <- RxODE::rxSolve(object=object$ppk_model,
                                   params=param_map,
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
#' model, a set of individual parameters, and a target AUC.
#'
#' @param object A posologyr list, created by the \code{\link{posologyr}}
#'     function.
#' @param param_map A vector of individual parameters. May be omitted,
#'     in which case the \code{\link{poso_estim_map}} function
#'     will be called
#' @param time_auc a numeric last point in time of the AUC for which the dose
#'     is to be optimized. The AUC is computed from 0 to `time_auc`
#' @param starting_dose numeric starting dose for the optimization
#'     algorithm
#' @param duration a numeric duration of infusion, for zero-order
#'     administrations
#' @param target_auc a numeric target AUC
#'
#' @return A numeric optimal dose to reach the target AUC
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
                          starting_dose=100,duration=NULL,
                          target_auc=NULL){

  if (is.null(param_map)){ #theta_pop + MAP estimates of eta + covariates
    model_map <- poso_estim_map(object,return_model = TRUE)
    param_map <- model_map[[2]]$params
  }

  err_dose <- function(dose,time_auc,target_auc,prior_model,
                       duration=duration,param_map){
   #compute the individual time-concentration profile
   event_table_auc <- RxODE::et(time=0,amt=dose,dur=duration)
   event_table_auc$add.sampling(time_auc)

   auc_ppk_model <- RxODE::rxSolve(object=prior_model$ppk_model,
                                    params=param_map,
                                    event_table_auc)
   #return the difference between the computed AUC and the target
   delta_auc = (target_auc - max(auc_ppk_model$AUC))^2
   return(delta_auc)
 }

 optim_dose_auc <- optim(starting_dose,err_dose,time_auc=time_auc,
                         target_auc=target_auc,prior_model=object,
                         duration=duration,param_map=param_map,
                         method="Brent",lower=0,upper=1e5)

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
#' @param object A posologyr list, created by the \code{\link{posologyr}}
#'     function.
#' @param param_map A vector of individual parameters. May be omitted,
#'     in which case the \code{\link{poso_estim_map}} function
#'     will be called
#' @param time_c a numeric point in time for which the dose is to be
#'     optimized
#' @param starting_dose numeric starting dose for the optimization
#'     algorithm
#' @param duration a numeric duration of infusion, for zero-order
#'     administrations
#' @param target_conc a numeric target concentration
#'
#' @return A numeric optimal dose to reach the target concentration
#'     at the selected point in time
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
#' # one hour after starting the infusion
#' poso_dose_ctime(patient01_tobra,time_c=1,target_conc=80)
#'
#' @export
poso_dose_ctime <- function(object=NULL,param_map=NULL,time_c=NULL,
                            starting_dose=100,duration=NULL,
                            target_conc=NULL){

  if (is.null(param_map)){ #theta_pop + MAP estimates of eta + covariates
    model_map <- poso_estim_map(object,return_model = TRUE)
    param_map <- model_map[[2]]$params
  }

  err_dose <- function(dose,time_c,target_conc,prior_model,
                       duration=duration,param_map){

    #compute the individual time-concentration profile
    event_table_ctime <- RxODE::et(time=0,amt=dose,dur=duration)
    event_table_ctime$add.sampling(time_c)

    ctime_ppk_model <- RxODE::rxSolve(object=prior_model$ppk_model,
                                    params=param_map,
                                    event_table_ctime)
    #return the difference between the computed ctime and the target
    delta_ctime = (target_conc - ctime_ppk_model$Cc)^2
    return(delta_ctime)
  }

  optim_dose_ctime <- optim(starting_dose,err_dose,time_c=time_c,
                          target_conc=target_conc,prior_model=object,
                          duration=duration,param_map=param_map,
                          method="Brent",lower=0, upper=1e5)

  return(optim_dose_ctime$par)
}
