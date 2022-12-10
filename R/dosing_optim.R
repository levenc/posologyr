#-------------------------------------------------------------------------
# posologyr: individual dose optimisation using population PK
# Copyright (C) 2022  Cyril Leven
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
#' @param dat Dataframe. An individual subject dataset following the
#'     structure of NONMEM/rxode2 event records.
#' @param prior_model A \code{posologyr} prior population pharmacokinetics
#'    model, a list of six objects.
#' @param tdm A boolean. If `TRUE`: computes the predicted time to reach the
#'    target trough concentration (Cmin) following the last event from `dat`,
#'    and using Maximum A Posteriori estimation. If `FALSE` : performs the
#'    estimation for a simulated scenario defined by the remaining parameters.
#' @param target_cmin Numeric. Target trough concentration (Cmin).
#' @param dose Numeric. Dose administered.
#' @param estim_method A character string. An estimation method to be used for
#'    the individual parameters. The default method "map" is the Maximum A
#'    Posteriori estimation, the method "prior" simulates from the prior
#'    population model, and "sir" uses the Sequential Importance Resampling
#'    algorithm to estimate the a posteriori distribution of the individual
#'    parameters. This argument is ignored if `indiv_param` is provided.
#' @param nocb A boolean. For time-varying covariates: the next observation
#'     carried backward (nocb) interpolation style, similar to NONMEM.  If
#'     `FALSE`, the last observation carried forward (locf) style will be used.
#'     Defaults to `FALSE`.
#' @param p Numeric. The proportion of the distribution of cmin to consider for
#'    the estimation. Mandatory for `estim_method=sir`.
#' @param greater_than A boolean. If `TRUE`: targets a time leading to a
#'    proportion `p` of the cmins to be greater than `target_cmin`.
#'    Respectively, lower if `FALSE`.
#' @param from Numeric. Starting time for the simulation of the
#'     individual time-concentration profile. The default value is
#'     0.2
#' @param last_time Numeric. Ending time for the simulation of the
#'     individual time-concentration profile. The default value is
#'     72.
#' @param add_dose Numeric. Additional doses administered at
#'     inter-dose interval after the first dose. Optional.
#' @param interdose_interval Numeric. Time for the inter-dose interval
#'     for multiple dose regimen. Must be provided when add_dose is used.
#' @param duration Numeric. Duration of infusion, for zero-order
#'     administrations.
#' @param indiv_param Optional. A set of individual parameters : THETA,
#'     estimates of ETA, and covariates.
#'
#' @return A numeric time to the selected trough concentration, from the
#'     time of administration.
#'
#' @examples
#' # model
#' mod_run001 <- list(
#' ppk_model = rxode2::rxode({
#'   centr(0) = 0;
#'   depot(0) = 0;
#'
#'   TVCl = THETA_Cl;
#'   TVVc = THETA_Vc;
#'   TVKa = THETA_Ka;
#'
#'   Cl = TVCl*exp(ETA_Cl);
#'   Vc = TVVc*exp(ETA_Vc);
#'   Ka = TVKa*exp(ETA_Ka);
#'
#'   K20 = Cl/Vc;
#'   Cc = centr/Vc;
#'
#'   d/dt(depot) = -Ka*depot;
#'   d/dt(centr) = Ka*depot - K20*centr;
#'   d/dt(AUC) = Cc;
#' }),
#' error_model = function(f,sigma) {
#'   dv <- cbind(f,1)
#'   g <- diag(dv%*%sigma%*%t(dv))
#'   return(sqrt(g))
#' },
#' theta = c(THETA_Cl=4.0, THETA_Vc=70.0, THETA_Ka=1.0),
#' omega = lotri::lotri({ETA_Cl + ETA_Vc + ETA_Ka ~
#'     c(0.2,
#'       0, 0.2,
#'       0, 0, 0.2)}),
#' sigma = lotri::lotri({prop + add ~ c(0.05,0.0,0.00)}))
#' # df_patient01: event table for Patient01, following a 30 minutes intravenous
#' # infusion
#' df_patient01 <- data.frame(ID=1,
#'                         TIME=c(0.0,1.0,14.0),
#'                         DV=c(NA,25.0,5.5),
#'                         AMT=c(2000,0,0),
#'                         EVID=c(1,0,0),
#'                         DUR=c(0.5,NA,NA))
#' # predict the time needed to reach a concentration of 2.5 mg/l
#' # after the administration of a 2500 mg dose over a 30 minutes
#' # infusion
#' poso_time_cmin(dat=df_patient01,prior_model=mod_run001,
#' dose=2500,duration=0.5,from=0.5,target_cmin=2.5)
#'
#' @export
poso_time_cmin <- function(dat=NULL,prior_model=NULL,tdm=FALSE,
                           target_cmin,dose=NULL,estim_method="map",nocb=FALSE,
                           p=NULL,greater_than=TRUE,from=0.2,last_time=72,
                           add_dose=NULL,interdose_interval=NULL,
                           duration=NULL,indiv_param=NULL){
  object <- posologyr(prior_model,dat,nocb)

  #Get or simulate the individual scenario: events and observations------------
  if(tdm){ #using TDM data
    #clear all unused input
    estim_method <- NULL
    dose <- NULL
    p <- NULL
    greater_than <- NULL
    add_dose <- NULL
    interdose_interval <- NULL
    duration <- NULL
    indiv_param <- NULL
    select_proposal_from_distribution <- FALSE
    #MAP estimation of the individual PK profile from TDM data
    cmin_map <- poso_estim_map(dat=dat,prior_model=prior_model,nocb=nocb)
    #individual parameters as ETA + covariates
    covar <- utils::tail(dat[,prior_model$covariates],1)
    indiv_param <- cbind(cmin_map$model$params,covar)
    #time of the last dose
    time_last_dose <- utils::tail(cmin_map$event[cmin_map$event$evid == 1,],1)$time
    #starting time following the last dose
    from <- time_last_dose + from
    #last observation of the poso_estim_map output
    lobs_map <- utils::tail(cmin_map$event,1)$time
    #last observation desired
    lobs <- lobs_map + last_time
    #bind the inital eventTable with enough repetitions of the last row
    # allows to keep EVID==0 and the last known values of the covariates for every obs
    # where .N is a shortcut for "last row"
    extended_et <- rbind(cmin_map$event,
                         cmin_map$event[rep(.N,(lobs-lobs_map)/0.1)])
    #fill the time column with the desired observation times: seq from the last
    # observation of the MAP output to the last observation needed
    time <- NULL    # avoid undefined global variables
    extended_et[time>=lobs_map,time:=seq(lobs_map,lobs,
                                         length.out=1+(lobs-lobs_map)/0.1)]
    #Solve the model with the extended et
    #rxSolve needs:
    # an rxode2 model (eg. tobramycin_fictional$ppk_model)
    # an eventTable
    # parameters THETA from the prior model, ETA from the MAP estimation
    # optionnal: an interpolation method for time-varying covariates
    cmin_ppk_model <- rxode2::rxSolve(prior_model$ppk_model,extended_et,
                                      c(prior_model$theta,cmin_map$eta),
                                      covsInterpolation =
                                        ifelse(nocb,"nocb","locf"))
  } else { #tdm == FALSE: simulation of the individual scenario
    #Read and process the parameters required for the simulation
    read_input  <- read_optim_distribution_input(dat=dat,
                                                 prior_model=prior_model,
                                                 nocb=nocb,object=object,p=p,
                                                 estim_method=estim_method,
                                                 indiv_param=indiv_param)
    indiv_param <- read_input[[1]]
    select_proposal_from_distribution <- read_input[[2]]
    #Input validation
    if (!is.null(add_dose)){
      if (is.null(interdose_interval)){
        stop("interdose_interval is mandatory when add_dose is used.",
             call.=FALSE)
      }
    }

    #simulation of the individual scenario
    # create an event table with the required number of administrations
    if (!is.null(add_dose)){ #more than one dose is needed
      event_table_cmin <- rxode2::et(amt=dose,dur=duration,
                                     ii=interdose_interval,
                                     addl=add_dose)
      time_last_dose   <- add_dose*interdose_interval
      #add observations from the time of the last dose onward only
      event_table_cmin$add.sampling(seq(time_last_dose+from,
                                        time_last_dose+last_time,
                                        by=0.1))
    }
    else { #only one dose is needed: simulation of a single administration
      event_table_cmin <- rxode2::et(amt=dose,dur=duration)
      event_table_cmin$add.sampling(seq(from,last_time,by=0.1))
      time_last_dose   <- 0
    }

    #solve the model with custom event table to simulate the individual scenario
    cmin_ppk_model <- rxode2::rxSolve(object=object$ppk_model,
                                      params=indiv_param,
                                      event_table_cmin,
                                      nDisplayProgress=1e5)
  } #end of the individual scenario estimation

  #Estimate the time to Cmin----------------------------------------------------
  if (select_proposal_from_distribution){ #simulated scenario distribution of ETA
    wide_cmin      <- tidyr::pivot_wider(cmin_ppk_model,
                                        id_cols = "time",
                                        names_from = "sim.id",
                                        values_from = "Cc")

    get_p_cmin <- function(all_cmin){
      all_cmin     <- all_cmin[-1]
      sorted_cmin  <- sort(all_cmin)
      n_cmin       <- length(sorted_cmin)
      cmin_index   <- ceiling(p * n_cmin)

      if (greater_than){
        cmin_proposal <- sorted_cmin[n_cmin - cmin_index]
      } else {
        cmin_proposal <- sorted_cmin[cmin_index]
      }

      return(cmin_proposal)
    }

    cmin_p            <- apply(wide_cmin,MARGIN=1,FUN=get_p_cmin)
    wide_cmin$cmin_p  <- cmin_p

    #compute time_to_target taking into account the multiple doses if needed
    # time_to_target: time after the last dose for which cmin_p is lower than
    # target_cmin
    time_to_target <-
      min(wide_cmin$time[wide_cmin$cmin_p < target_cmin]) - time_last_dose

    #unlist cmin_distribution to get a vector of numeric values
    cmin_distribution <- wide_cmin[wide_cmin$time == (time_to_target +
                                     time_last_dose),]
    cmin_distribution$time <- cmin_distribution$cmin_p <- NULL
    cmin_distribution      <- unlist(cmin_distribution,use.names=FALSE)
  } else { #select_proposal_from_distribution == FALSE, TDM or not
    #compute time_to_target taking into account the multiple doses if needed
    # time_to_target: time after the last dose for which Cc is lower than
    # target_cmin
    time_to_target <-
      min(cmin_ppk_model[cmin_ppk_model$time > from &
                         cmin_ppk_model$Cc < target_cmin,]$time) - time_last_dose
    #here cmin_distribution is a point estimate of Cc
    cmin_distribution <<-
      cmin_ppk_model$Cc[cmin_ppk_model$time == (time_to_target +
                                                  time_last_dose)]
  }

  ifelse(select_proposal_from_distribution,
         type_of_estimate <-"distribution",
         type_of_estimate <- "point estimate")

  time_cmin <- list(time=time_to_target,
                   type_of_estimate=type_of_estimate,
                   cmin_estimate=cmin_distribution,
                   indiv_param=indiv_param)
  return(time_cmin)
}

#' Estimate the optimal dose for a selected target area under the
#' time-concentration curve (AUC)
#'
#' Estimates the optimal dose for a selected target area under the
#' time-concentration curve (AUC) given a population pharmacokinetic
#' model, a set of individual parameters, and a target AUC.
#'
#' @param dat Dataframe. An individual subject dataset following the
#'     structure of NONMEM/rxode2 event records.
#' @param prior_model A \code{posologyr} prior population pharmacokinetics
#'    model, a list of six objects.
#' @param time_auc Numeric. The target AUC is computed from 0 to `time_auc`.
#' @param target_auc Numeric. The target AUC.
#' @param estim_method A character string. An estimation method to be used for
#'    the individual parameters. The default method "map" is the Maximum A
#'    Posteriori estimation, the method "prior" simulates from the prior
#'    population model, and "sir" uses the Sequential Importance Resampling
#'    algorithm to estimate the a posteriori distribution of the individual
#'    parameters. This argument is ignored if `indiv_param` is provided.
#' @param nocb A boolean. for time-varying covariates: the next observation
#'     carried backward (nocb) interpolation style, similar to NONMEM.  If
#'     `FALSE`, the last observation carried forward (locf) style will be used.
#'     Defaults to `FALSE`.
#' @param p Numeric. The proportion of the distribution of AUC to consider for
#'    the optimization. Mandatory for `estim_method=sir`.
#' @param greater_than A boolean. If `TRUE`: targets a dose leading to a
#'    proportion `p` of the AUCs to be greater than `target_auc`. Respectively,
#'    lower if `FALSE`.
#' @param starting_time Numeric. First point in time of the AUC, for multiple
#'     dose regimen. The default is zero.
#' @param interdose_interval Numeric. Time for the interdose interval for
#'     multiple dose regimen. Must be provided when add_dose is used.
#' @param add_dose Numeric. Additional doses administered at
#'     inter-dose interval after the first dose. Optional.
#' @param duration Numeric. Duration of infusion, for zero-order
#'     administrations.
#' @param starting_dose Numeric. Starting dose for the optimization
#'     algorithm.
#' @param indiv_param Optional. A set of individual parameters : THETA,
#'     estimates of ETA, and covariates.
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{dose}{Numeric. An optimal dose for the selected target AUC.}
#'   \item{type_of_estimate}{Character string. The type of estimate of the
#'   individual parameters. Either a point estimate, or a distribution.}
#'   \item{auc_estimate}{A vector of numeric estimates of the AUC. Either a
#'   single value (for a point estimate of ETA), or a distribution.}
#'   \item{indiv_param}{A `data.frame`. The set of individual parameters used
#'   for the determination of the optimal dose : THETA, estimates of ETA, and
#'   covariates}
#' }
#'
#' @examples
#' # model
#' mod_run001 <- list(
#' ppk_model = rxode2::rxode({
#'   centr(0) = 0;
#'   depot(0) = 0;
#'
#'   TVCl = THETA_Cl;
#'   TVVc = THETA_Vc;
#'   TVKa = THETA_Ka;
#'
#'   Cl = TVCl*exp(ETA_Cl);
#'   Vc = TVVc*exp(ETA_Vc);
#'   Ka = TVKa*exp(ETA_Ka);
#'
#'   K20 = Cl/Vc;
#'   Cc = centr/Vc;
#'
#'   d/dt(depot) = -Ka*depot;
#'   d/dt(centr) = Ka*depot - K20*centr;
#'   d/dt(AUC) = Cc;
#' }),
#' error_model = function(f,sigma) {
#'   dv <- cbind(f,1)
#'   g <- diag(dv%*%sigma%*%t(dv))
#'   return(sqrt(g))
#' },
#' theta = c(THETA_Cl=4.0, THETA_Vc=70.0, THETA_Ka=1.0),
#' omega = lotri::lotri({ETA_Cl + ETA_Vc + ETA_Ka ~
#'     c(0.2,
#'       0, 0.2,
#'       0, 0, 0.2)}),
#' sigma = lotri::lotri({prop + add ~ c(0.05,0.0,0.00)}))
#' # df_patient01: event table for Patient01, following a 30 minutes intravenous
#' # infusion
#' df_patient01 <- data.frame(ID=1,
#'                         TIME=c(0.0,1.0,14.0),
#'                         DV=c(NA,25.0,5.5),
#'                         AMT=c(2000,0,0),
#'                         EVID=c(1,0,0),
#'                         DUR=c(0.5,NA,NA))
#' # estimate the optimal dose to reach an AUC(0-12h) of 45 h.mg/l
#' poso_dose_auc(dat=df_patient01,prior_model=mod_run001,
#' time_auc=12,target_auc=45)
#'
#' @export
poso_dose_auc <- function(dat=NULL,prior_model=NULL,time_auc,target_auc,
                          estim_method="map",nocb=FALSE,p=NULL,
                          greater_than=TRUE,starting_time=0,
                          interdose_interval=NULL,add_dose=NULL,
                          duration=NULL,starting_dose=100,indiv_param=NULL){

  # Input validation -----------------------------------------------------------
  object <- posologyr(prior_model,dat,nocb)

  read_input  <- read_optim_distribution_input(dat=dat,
                                               prior_model=prior_model,
                                               nocb=nocb,object=object,p=p,
                                               estim_method=estim_method,
                                               indiv_param=indiv_param)
  indiv_param <- read_input[[1]]
  select_proposal_from_distribution <- read_input[[2]]

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
  # Optimization ---------------------------------------------------------------
  err_dose <- function(dose,time_auc,starting_time,target_auc,
                       interdose_interval,add_dose,prior_model,
                       duration,indiv_param){

    #compute the individual time-concentration profile
    if (!is.null(add_dose)){
      event_table_auc <- rxode2::et(amt=dose,dur=duration,
                                   ii=interdose_interval,
                                   addl=add_dose)
    } else {
      event_table_auc <- rxode2::et(amt=dose,dur=duration)
    }

    event_table_auc$add.sampling(starting_time)
    event_table_auc$add.sampling(starting_time+time_auc)

    auc_ppk_model <- rxode2::rxSolve(object=prior_model$ppk_model,
                                    params=indiv_param,
                                    event_table_auc,
                                    nDisplayProgress=1e5)

    if (select_proposal_from_distribution){
      wide_auc      <- tidyr::pivot_wider(auc_ppk_model,
                                          id_cols = "sim.id",
                                          names_from = "time",
                                          values_from = "AUC")

      get_auc_difference <- function(auc_pairs){
        auc_diff <- auc_pairs[3] - auc_pairs[2]
        return(auc_diff)
      }

      auc_difference <- apply(wide_auc,MARGIN=1,FUN=get_auc_difference)
      sorted_auc     <- sort(auc_difference)
      n_auc          <- length(sorted_auc)
      auc_index      <- ceiling(p * n_auc)

      # assign the distribution of auc to the parent environment
      auc_distribution <<- sorted_auc

      if (greater_than){
        auc_proposal <- sorted_auc[n_auc - auc_index]
      } else {
        auc_proposal <- sorted_auc[auc_index]
      }
    } else {
      auc_proposal  <- max(auc_ppk_model$AUC)-min(auc_ppk_model$AUC)
      auc_distribution <<- auc_proposal
    }

    #return the difference between the computed AUC and the target
    delta_auc<-(target_auc - auc_proposal)^2
    return(delta_auc)
  }

  #initialization of auc_distribution to avoid a global variable
  auc_distribution <- 0

  optim_dose_auc <- stats::optim(starting_dose,err_dose,time_auc=time_auc,
                                 starting_time=starting_time,add_dose=add_dose,
                                 interdose_interval=interdose_interval,
                                 target_auc=target_auc,prior_model=object,
                                 duration=duration,indiv_param=indiv_param,
                                 method="Brent",lower=0,upper=1e5)

  ifelse(select_proposal_from_distribution,
         type_of_estimate <-"distribution",
         type_of_estimate <- "point estimate")

  dose_auc <- list(dose=optim_dose_auc$par,
                   type_of_estimate=type_of_estimate,
                   auc_estimate=auc_distribution,
                   indiv_param=indiv_param)
  return(dose_auc)
}

#' Estimate the optimal dose for a selected target concentration
#'
#' Estimates the optimal dose for a selected target concentration at a
#' selected point in time given a population pharmacokinetic model, a set
#' of individual parameters, a selected point in time, and a target
#' concentration.
#'
#' @param dat Dataframe. An individual subject dataset following the
#'     structure of NONMEM/rxode2 event records.
#' @param prior_model A \code{posologyr} prior population pharmacokinetics
#'    model, a list of six objects.
#' @param tdm A boolean. If `TRUE`: estimates the optimal dose for a selected
#'    target concentration at a selected point in time following the events from
#'     `dat`, and using Maximum A Posteriori estimation. If `FALSE` : performs
#'    the estimation in a simulated scenario defined by the remaining parameters.
#' @param time_c Numeric. Point in time for which the dose is to be
#'     optimized.
#' @param time_dose Numeric. Time when the dose is to be given.
#' @param target_conc Numeric. Target concentration.
#' @param estim_method A character string. An estimation method to be used for
#'    the individual parameters. The default method "map" is the Maximum A
#'    Posteriori estimation, the method "prior" simulates from the prior
#'    population model, and "sir" uses the Sequential Importance Resampling
#'    algorithm to estimate the a posteriori distribution of the individual
#'    parameters. This argument is ignored if `indiv_param` is provided.
#' @param nocb A boolean. for time-varying covariates: the next observation
#'     carried backward (nocb) interpolation style, similar to NONMEM.  If
#'     `FALSE`, the last observation carried forward (locf) style will be used.
#'     Defaults to `FALSE`.
#' @param p Numeric. The proportion of the distribution of concentrations to
#'    consider for the optimization. Mandatory for `estim_method=sir`.
#' @param greater_than A boolean. If `TRUE`: targets a dose leading to a
#'    proportion `p` of the concentrations to be greater than `target_conc`.
#'    Respectively, lower if `FALSE`.
#' @param starting_dose Numeric. Starting dose for the optimization
#'     algorithm.
#' @param add_dose Numeric. Additional doses administered at
#'     inter-dose interval after the first dose. Optional.
#' @param interdose_interval Numeric. Time for the interdose interval
#'     for multiple dose regimen. Must be provided when add_dose is used.
#' @param duration Numeric. Duration of infusion, for zero-order
#'     administrations.
#' @param indiv_param Optional. A set of individual parameters : THETA,
#'     estimates of ETA, and covariates.
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{dose}{Numeric. An optimal dose for the selected target concentration.}
#'   \item{type_of_estimate}{Character string. The type of estimate of the
#'   individual parameters. Either a point estimate, or a distribution.}
#'   \item{conc_estimate}{A vector of numeric estimates of the conc. Either a
#'   single value (for a point estimate of ETA), or a distribution.}
#'   \item{indiv_param}{A `data.frame`. The set of individual parameters used
#'   for the determination of the optimal dose : THETA, estimates of ETA, and
#'   covariates}
#' }
#'
#' @examples
#' # model
#' mod_run001 <- list(
#' ppk_model = rxode2::rxode({
#'   centr(0) = 0;
#'   depot(0) = 0;
#'
#'   TVCl = THETA_Cl;
#'   TVVc = THETA_Vc;
#'   TVKa = THETA_Ka;
#'
#'   Cl = TVCl*exp(ETA_Cl);
#'   Vc = TVVc*exp(ETA_Vc);
#'   Ka = TVKa*exp(ETA_Ka);
#'
#'   K20 = Cl/Vc;
#'   Cc = centr/Vc;
#'
#'   d/dt(depot) = -Ka*depot;
#'   d/dt(centr) = Ka*depot - K20*centr;
#'   d/dt(AUC) = Cc;
#' }),
#' error_model = function(f,sigma) {
#'   dv <- cbind(f,1)
#'   g <- diag(dv%*%sigma%*%t(dv))
#'   return(sqrt(g))
#' },
#' theta = c(THETA_Cl=4.0, THETA_Vc=70.0, THETA_Ka=1.0),
#' omega = lotri::lotri({ETA_Cl + ETA_Vc + ETA_Ka ~
#'     c(0.2,
#'       0, 0.2,
#'       0, 0, 0.2)}),
#' sigma = lotri::lotri({prop + add ~ c(0.05,0.0,0.00)}))
#' # df_patient01: event table for Patient01, following a 30 minutes intravenous
#' # infusion
#' df_patient01 <- data.frame(ID=1,
#'                         TIME=c(0.0,1.0,14.0),
#'                         DV=c(NA,25.0,5.5),
#'                         AMT=c(2000,0,0),
#'                         EVID=c(1,0,0),
#'                         DUR=c(0.5,NA,NA))
#' # estimate the optimal dose to reach a concentration of 80 mg/l
#' # one hour after starting the 30-minutes infusion
#' poso_dose_conc(dat=df_patient01,prior_model=mod_run001,
#' time_c=1,duration=0.5,target_conc=80)
#'
#' @export
poso_dose_conc <- function(dat=NULL,prior_model=NULL,tdm=FALSE,
                           time_c,time_dose=NULL,target_conc,estim_method="map",
                           nocb=FALSE,p=NULL,greater_than=TRUE,
                           starting_dose=100,interdose_interval=NULL,
                           add_dose=NULL,duration=NULL,indiv_param=NULL){
  object <- posologyr(prior_model,dat,nocb)

  if(tdm){ #using TDM data
    #input validation
    if (time_c<=time_dose){
      stop("time_c cannot be before time_dose.",
           call.=FALSE)
    }
    #clear all unused input
    estim_method <- NULL
    p <- NULL
    greater_than <- NULL
    add_dose <- NULL
    interdose_interval <- NULL
    indiv_param <- NULL
    select_proposal_from_distribution <- FALSE
    #MAP estimation of the individual PK profile from TDM data
    ctime_map <- poso_estim_map(dat=dat,prior_model=prior_model,nocb=nocb)
    #individual parameters as ETA + covariates
    covar <- utils::tail(dat[,prior_model$covariates],1)
    indiv_param <- cbind(ctime_map$model$params,covar)
    #remove all observations from the eventTable, keep only dosing
    extended_et <- ctime_map$event[ctime_map$event$evid == 1,]
    #last event of the poso_estim_map output
    last_event <- utils::tail(extended_et,1)$time
    #more input validation
    if (time_dose<=last_event){
      stop("time_dose must occur after the last recorded dosing.",
           call.=FALSE)
    }
    #bind the inital eventTable with enough repetitions of the last row
    # one row for dosing
    # one row for the observation
    extended_et <- rbind(extended_et,
                         extended_et[rep(.N,2)])
    #fill the time column with the desired observation times:
    # the simulated administration, and subsequent concentration
    time <- NULL    # avoid undefined global variables
    extended_et[time>=last_event,time:=c(last_event,time_dose,time_c)]

    err_dose_tdm <- function(dose,time_dose,target_conc,time_c,prior_model,
                             extended_et,nocb){
      #New dose at time_dose in the extended_et
      # This shorter syntax with "list" allows for setting several variable at once
      extended_et[time==time_dose,c("evid","amt","dur"):=list(1,dose,duration)]
      #New observation at time_c
      extended_et[time==time_c,c("evid","amt","dur"):=list(0,NA,NA)]
      #Solve the model with the extended et
      ctime_ppk_model <- rxode2::rxSolve(prior_model$ppk_model,extended_et,
                                         c(prior_model$theta,ctime_map$eta),
                                         covsInterpolation =
                                           ifelse(nocb,"nocb","locf"))

      conc_proposal     <- ctime_ppk_model$Cc
      conc_distribution <<- conc_proposal #to parent environment
      #return the difference between the computed ctime and the target
      delta_conc <- (target_conc - conc_proposal)^2
      return(delta_conc)
    }
    #initialization of conc_distribution to avoid a global variable
    conc_distribution <- 0

    optim_dose_conc <- stats::optim(starting_dose,err_dose_tdm,
                                    time_dose=time_dose,target_conc=target_conc,
                                    time_c=time_c,prior_model=prior_model,
                                    extended_et=extended_et,nocb=nocb,
                                    method="Brent",lower=0, upper=1e5)

  } else { #tdm == FALSE: simulation of the individual scenario
    #time_dose is only used when tdm=TRUE
    time_dose <- NULL
    #Read and process the parameters required for the simulation
    read_input  <- read_optim_distribution_input(dat=dat,
                                                 prior_model=prior_model,
                                                 nocb=nocb,object=object,p=p,
                                                 estim_method=estim_method,
                                                 indiv_param=indiv_param)
    indiv_param <- read_input[[1]]
    select_proposal_from_distribution <- read_input[[2]]
    #Input validation
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
                         duration=duration,indiv_param){

      #compute the individual time-concentration profile
      if (!is.null(add_dose)){
        event_table_ctime <- rxode2::et(amt=dose,dur=duration,
                                        ii=interdose_interval,
                                        addl=add_dose)
      }
      else {
        event_table_ctime <- rxode2::et(amt=dose,dur=duration)
      }
      event_table_ctime$add.sampling(time_c)

      ctime_ppk_model <- rxode2::rxSolve(object=prior_model$ppk_model,
                                         params=indiv_param,
                                         event_table_ctime,
                                         nDisplayProgress=1e5)

      if (select_proposal_from_distribution){

        sorted_conc     <- sort(ctime_ppk_model$Cc)
        n_conc          <- length(sorted_conc)
        conc_index      <- ceiling(p * n_conc)

        # assign the distribution of concentrations to the parent environment
        conc_distribution <<- sorted_conc

        if (greater_than){
          conc_proposal <- sorted_conc[n_conc - conc_index]
        } else {
          conc_proposal <- sorted_conc[conc_index]
        }
      } else {
        conc_proposal     <- ctime_ppk_model$Cc
        conc_distribution <<- conc_proposal
      }

      #return the difference between the computed ctime and the target
      delta_conc <- (target_conc - conc_proposal)^2
      return(delta_conc)
    }
    #initialization of conc_distribution to avoid a global variable
    conc_distribution <- 0

    optim_dose_conc <- stats::optim(starting_dose,err_dose,time_c=time_c,
                                    target_conc=target_conc,prior_model=object,
                                    add_dose=add_dose,
                                    interdose_interval=interdose_interval,
                                    duration=duration,indiv_param=indiv_param,
                                    method="Brent",lower=0, upper=1e5)
  }

  ifelse(select_proposal_from_distribution,
         type_of_estimate <-"distribution",
         type_of_estimate <- "point estimate")

  dose_conc <- list(dose=optim_dose_conc$par,
                   type_of_estimate=type_of_estimate,
                   conc_estimate=conc_distribution,
                   indiv_param=indiv_param)
  return(dose_conc)
}

#' Estimate the optimal inter-dose interval for a given dose and a
#' selected target trough concentration
#'
#' Estimates the optimal inter-dose interval for a selected target
#' trough concentration (Cmin), given a dose, a population
#' pharmacokinetic model, a set of individual parameters, and a
#' target concentration.
#'
#' @param dat Dataframe. An individual subject dataset following the
#'     structure of NONMEM/rxode2 event records.
#' @param prior_model A \code{posologyr} prior population pharmacokinetics
#'    model, a list of six objects.
#' @param target_cmin Numeric. Target trough concentration (Cmin).
#' @param dose Numeric. The dose given.
#' @param estim_method A character string. An estimation method to be used for
#'    the individual parameters. The default method "map" is the Maximum A
#'    Posteriori estimation, the method "prior" simulates from the prior
#'    population model, and "sir" uses the Sequential Importance Resampling
#'    algorithm to estimate the a posteriori distribution of the individual
#'    parameters. This argument is ignored if `indiv_param` is provided.
#' @param nocb A boolean. for time-varying covariates: the next observation
#'     carried backward (nocb) interpolation style, similar to NONMEM.  If
#'     `FALSE`, the last observation carried forward (locf) style will be used.
#'     Defaults to `FALSE`.
#' @param p Numeric. The proportion of the distribution of concentrations to
#'    consider for the optimization. Mandatory for `estim_method=sir`.
#' @param greater_than A boolean. If `TRUE`: targets a dose leading to a
#'    proportion `p` of the concentrations to be greater than `target_conc`.
#'    Respectively, lower if `FALSE`.
#' @param starting_interval Numeric. Starting inter-dose interval for
#'     the optimization algorithm.
#' @param add_dose Numeric. Additional doses administered at
#'     inter-dose interval after the first dose.
#' @param duration Numeric. Duration of infusion, for zero-order
#'     administrations.
#' @param indiv_param Optional. A set of individual parameters : THETA,
#'     estimates of ETA, and covariates.
#'
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{interval}{Numeric. An inter-dose interval to reach the target trough
#'     concentration before each dosing of a multiple dose regimen.}
#'   \item{type_of_estimate}{Character string. The type of estimate of the
#'   individual parameters. Either a point estimate, or a distribution.}
#'   \item{conc_estimate}{A vector of numeric estimates of the conc. Either a
#'   single value (for a point estimate of ETA), or a distribution.}
#'   \item{indiv_param}{A `data.frame`. The set of individual parameters used
#'   for the determination of the optimal dose : THETA, estimates of ETA, and
#'   covariates}
#' }
#'
#' @examples
#' # model
#' mod_run001 <- list(
#' ppk_model = rxode2::rxode({
#'   centr(0) = 0;
#'   depot(0) = 0;
#'
#'   TVCl = THETA_Cl;
#'   TVVc = THETA_Vc;
#'   TVKa = THETA_Ka;
#'
#'   Cl = TVCl*exp(ETA_Cl);
#'   Vc = TVVc*exp(ETA_Vc);
#'   Ka = TVKa*exp(ETA_Ka);
#'
#'   K20 = Cl/Vc;
#'   Cc = centr/Vc;
#'
#'   d/dt(depot) = -Ka*depot;
#'   d/dt(centr) = Ka*depot - K20*centr;
#'   d/dt(AUC) = Cc;
#' }),
#' error_model = function(f,sigma) {
#'   dv <- cbind(f,1)
#'   g <- diag(dv%*%sigma%*%t(dv))
#'   return(sqrt(g))
#' },
#' theta = c(THETA_Cl=4.0, THETA_Vc=70.0, THETA_Ka=1.0),
#' omega = lotri::lotri({ETA_Cl + ETA_Vc + ETA_Ka ~
#'     c(0.2,
#'       0, 0.2,
#'       0, 0, 0.2)}),
#' sigma = lotri::lotri({prop + add ~ c(0.05,0.0,0.00)}))
#' # df_patient01: event table for Patient01, following a 30 minutes intravenous
#' # infusion
#' df_patient01 <- data.frame(ID=1,
#'                         TIME=c(0.0,1.0,14.0),
#'                         DV=c(NA,25.0,5.5),
#'                         AMT=c(2000,0,0),
#'                         EVID=c(1,0,0),
#'                         DUR=c(0.5,NA,NA))
#' # estimate the optimal interval to reach a cmin of of 2.5 mg/l
#' # before each administration
#' poso_inter_cmin(dat=df_patient01,prior_model=mod_run001,
#' dose=1500,duration=0.5,target_cmin=2.5)
#'
#' @export
poso_inter_cmin <- function(dat=NULL,prior_model=NULL,dose,target_cmin,
                            estim_method="map",nocb=FALSE,p=NULL,
                            greater_than=TRUE,starting_interval=12,add_dose=10,
                            duration=NULL,indiv_param=NULL){
  object <- posologyr(prior_model,dat,nocb)

  read_input  <- read_optim_distribution_input(dat=dat,
                                               prior_model=prior_model,
                                               nocb=nocb,object=object,p=p,
                                               estim_method=estim_method,
                                               indiv_param=indiv_param)
  indiv_param <- read_input[[1]]
  select_proposal_from_distribution <- read_input[[2]]

  err_inter <- function(interdose_interval,dose,target_cmin,
                        prior_model,add_dose,duration=duration,
                        indiv_param){
    #compute the individual time-concentration profile
    event_table_cmin <- rxode2::et(amt=dose,dur=duration,
                                  ii=interdose_interval,
                                  addl=add_dose)
    event_table_cmin$add.sampling(interdose_interval*add_dose-0.1)

    cmin_ppk_model <- rxode2::rxSolve(object=prior_model$ppk_model,
                                     params=indiv_param,
                                     event_table_cmin,
                                     nDisplayProgress=1e5)

    if (select_proposal_from_distribution){

      sorted_cmin     <- sort(cmin_ppk_model$Cc)
      n_cmin          <- length(sorted_cmin)
      cmin_index      <- ceiling(p * n_cmin)

      # assign the distribution of cmin to the parent environment
      cmin_distribution <<- sorted_cmin

      if (greater_than){
        cmin_proposal <- sorted_cmin[n_cmin - cmin_index]
      } else {
        cmin_proposal <- sorted_cmin[cmin_index]
      }
    } else {
      cmin_proposal     <-  cmin_ppk_model$Cc
      cmin_distribution <<- cmin_proposal
    }

    #return the difference between the computed cmin and the target,
    # normalized by the computed cmin to avoid divergence of the algorithm
    delta_cmin <- ((target_cmin - cmin_proposal)/cmin_proposal)^2
    return(delta_cmin)
  }

  #initialization of cmin_distribution to avoid a global variable
  cmin_distribution <- 0

  #cf. optim documentation: optim will work with one-dimensional pars, but the
  # default method does not work well (and will warn). Method "Brent" uses
  # optimize and needs bounds to be available; "BFGS" often works well enough if
  # not.
  optim_inter_cmin <- stats::optim(starting_interval,err_inter,dose=dose,
                                  target_cmin=target_cmin,prior_model=object,
                                  add_dose=add_dose,duration=duration,
                                  indiv_param=indiv_param,method="L-BFGS-B")

  ifelse(select_proposal_from_distribution,
         type_of_estimate <-"distribution",
         type_of_estimate <- "point estimate")

  inter_cmin <- list(interval=optim_inter_cmin$par,
                    type_of_estimate=type_of_estimate,
                    conc_estimate=cmin_distribution,
                    indiv_param=indiv_param)

  return(inter_cmin)
}

#-------------------------------------------------------------------------
#
# Internal functions for optimal dosing
#
#-------------------------------------------------------------------------

# get the parameters for distribution-based optimal dosing
read_optim_distribution_input <- function(dat,prior_model,
                                          nocb,object,p,
                                          estim_method,
                                          indiv_param){
  if (is.null(indiv_param)){ #theta_pop + estimates of eta + covariates
    if (estim_method=="map"){
      model_map   <- poso_estim_map(dat,prior_model,nocb=nocb,return_model=TRUE)
      if(is.null(object$covariates)){
        indiv_param <- model_map[[2]]$params
      } else {
        covar <- as.data.frame(object$tdm_data[length(object$tdm_data[,1]),
                                                     object$covariates])
        names(covar) <- object$covariates
        indiv_param <- cbind(model_map[[2]]$params,covar,row.names=NULL)
      }
      if (!is.null(object$pi_matrix)){
        kappa_mat         <- matrix(0,nrow=1,ncol=ncol(object$pi_matrix))
        kappa_df          <- data.frame(kappa_mat)
        names(kappa_df)   <- attr(object$pi_matrix,"dimnames")[[1]]
        indiv_param       <- cbind(indiv_param,kappa_df)
      }
      select_proposal_from_distribution <- FALSE
      if (!is.null(p)){
        warning('p is not needed with estim_method="map", p is ignored')
      }
    } else  if (estim_method=="prior"){
      if (!is.null(p)){
        if (p < 0 || p >= 1){
          stop('p must be between 0 and 1')
        }
        model_pop   <- poso_simu_pop(dat,prior_model,n_simul=1e5,return_model=TRUE)
        select_proposal_from_distribution <- TRUE
      } else {
        model_pop   <- poso_simu_pop(dat,prior_model,n_simul=0,return_model=TRUE)
        select_proposal_from_distribution <- FALSE
      }
      if(is.null(object$covariates)){
        indiv_param <- model_pop[[2]]$params
      } else {
        covar <- as.data.frame(object$tdm_data[length(object$tdm_data[,1]),
                                               object$covariates])
        names(covar) <- object$covariates
        indiv_param <- cbind(model_pop[[2]]$params,covar,row.names=NULL)
      }
    } else if (estim_method=="sir"){
      if (p < 0 || p >= 1){
        stop('p must be between 0 and 1')
      }
      model_sir   <- poso_estim_sir(dat,prior_model,nocb=nocb,
                                    n_sample=1e5,n_resample=1e4,
                                    return_model=TRUE)
      if(is.null(object$covariates)){
        indiv_param <- model_sir[[2]]$params
      } else {
        covar <- as.data.frame(object$tdm_data[length(object$tdm_data[,1]),
                                               object$covariates])
        names(covar) <- object$covariates
        indiv_param <- cbind(model_sir[[2]]$params,covar,row.names=NULL)
      }
      if (!is.null(object$pi_matrix)){
        kappa_mat         <- matrix(0,nrow=1,ncol=ncol(object$pi_matrix))
        kappa_df          <- data.frame(kappa_mat)
        names(kappa_df)   <- attr(object$pi_matrix,"dimnames")[[1]]
        indiv_param       <- cbind(indiv_param,kappa_df)
      }
      select_proposal_from_distribution <- TRUE
    } else {
      print(estim_method)
      stop("'estim_method' not recognized")
    }
  } else {
    if (FALSE %in% (c(names(object$solved_ppk_model$params))
                    %in% names(indiv_param))){
      stop("The names of indiv_param do not match the parameters of the object")
    }
    if (!is.null(p) && (length(rbind(indiv_param[,1])) < 1000)){
      warn_1000 <-
        sprintf("In order to perform the optimization using a parameter distribution, you
need at least 1000 parameter samples. Only the first set of parameters will
be used.")
      warning(warn_1000)
      indiv_param <- rbind(indiv_param)[1,]
      select_proposal_from_distribution <- FALSE
    } else if (!is.null(p) && (length(rbind(indiv_param)[,1]) >= 1000)){
      select_proposal_from_distribution <- TRUE
    } else { # p==NULL, using one set of parameters
      indiv_param <- rbind(indiv_param)[1,]
      select_proposal_from_distribution <- FALSE
    }
  }
  return(list(indiv_param,select_proposal_from_distribution))
}
