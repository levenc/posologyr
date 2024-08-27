#-------------------------------------------------------------------------
# posologyr: individual dose optimization using population PK
# Copyright (C) Cyril Leven
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

#' Estimate the time required to reach a target trough concentration (Cmin)
#'
#' Estimates the time required to reach a target trough concentration (Cmin)
#' given a population pharmacokinetic model, a set of individual
#' parameters, a dose, and a target Cmin.
#'
#' @param dat Dataframe. An individual subject dataset following the
#'    structure of NONMEM/rxode2 event records.
#' @param prior_model A \code{posologyr} prior population pharmacokinetics
#'    model, a list of six objects.
#' @param tdm A boolean. If `TRUE`: computes the predicted time to reach the
#'    target trough concentration (Cmin) following the last event from `dat`,
#'    and using Maximum A Posteriori estimation. Setting `tdm` to `TRUE` causes
#'    the following to occur:
#'
#'    * the simulation  starts at the time of the last recorded dose (from the
#'    TDM data) plus `from`;
#'    * the simulation stops at the time of the last recorded dose (from the TDM
#'     data) plus `last_time`;
#'    * the arguments `dose`, `duration`, `estim_method`, `p`, `greater_than`,
#'    `interdose_interval`, `add_dose`, `indiv_param` and `starting_time` are
#'    ignored.
#'
#' @param target_cmin Numeric. Target trough concentration (Cmin).
#' @param dose Numeric. Dose administered. This argument is ignored if `tdm` is
#'    set to `TRUE`.
#' @param cmt_dose Character or numeric. The compartment in which the dose is
#'    to be administered. Must match one of the compartments in the prior model.
#'    Defaults to 1.
#' @param endpoint Character. The endpoint of the prior model to be optimised
#'    for. The default is "Cc", which is the central concentration.
#' @param estim_method A character string. An estimation method to be used for
#'    the individual parameters. The default method "map" is the Maximum A
#'    Posteriori estimation, the method "prior" simulates from the prior
#'    population model, and "sir" uses the Sequential Importance Resampling
#'    algorithm to estimate the a posteriori distribution of the individual
#'    parameters. This argument is ignored if `indiv_param` is provided, or if
#'    `tdm` is set to `TRUE`.
#' @param nocb A boolean. For time-varying covariates: the next observation
#'    carried backward (nocb) interpolation style, similar to NONMEM. If
#'    `FALSE`, the last observation carried forward (locf) style will be used.
#'    Defaults to `FALSE`.
#' @param p Numeric. The proportion of the distribution of Cmin to consider for
#'    the estimation. Mandatory for `estim_method=sir`. This argument is ignored
#'    if `tdm` is set to `TRUE`.
#' @param greater_than A boolean. If `TRUE`: targets a time leading to a
#'    proportion `p` of the cmins to be greater than `target_cmin`.
#'    Respectively, lower if `FALSE`. This argument is ignored if `tdm` is set
#'    to `TRUE`.
#' @param from Numeric. Starting time for the simulation of the individual
#'    time-concentration profile. The default value is 0.2. When `tdm` is set
#'    to `TRUE` the simulation starts at the time of the last recorded dose plus
#'    `from`.
#' @param last_time Numeric. Ending time for the simulation of the individual
#'    time-concentration profile. The default value is 72. When `tdm` is set to
#'    `TRUE` the simulation stops at the time of the last recorded dose plus
#'    `last_time`.
#' @param add_dose Numeric. Additional doses administered at inter-dose interval
#'    after the first dose. Optional. This argument is ignored if `tdm` is
#'    set to `TRUE`.
#' @param interdose_interval Numeric. Time for the inter-dose interval for
#'    multiple dose regimen. Must be provided when add_dose is used. This
#'    argument is ignored if `tdm` is set to `TRUE`.
#' @param duration Numeric. Duration of infusion, for zero-order
#'    administrations. This argument is ignored if `tdm` is set to `TRUE`.
#' @param indiv_param Optional. A set of individual parameters : THETA,
#'    estimates of ETA, and covariates.
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{time}{Numeric. Time needed to reach the selected Cmin.}
#'   \item{type_of_estimate}{Character string. The type of estimate of the
#'   individual parameters. Either a point estimate, or a distribution.}
#'   \item{cmin_estimate}{A vector of numeric estimates of the Cmin. Either a
#'   single value (for a point estimate of ETA), or a distribution.}
#'   \item{indiv_param}{A `data.frame`. The set of individual parameters used
#'   for the determination of the time needed to reach a selected Cmin: THETA,
#'   estimates of ETA, and covariates}
#'   }
#'
#' @examples
#' rxode2::setRxThreads(2L) # limit the number of threads
#'
#' # model
#' mod_run001 <- function() {
#'   ini({
#'     THETA_Cl <- 4.0
#'     THETA_Vc <- 70.0
#'     THETA_Ka <- 1.0
#'     ETA_Cl ~ 0.2
#'     ETA_Vc ~ 0.2
#'     ETA_Ka ~ 0.2
#'     prop.sd <- sqrt(0.05)
#'   })
#'   model({
#'     TVCl <- THETA_Cl
#'     TVVc <- THETA_Vc
#'     TVKa <- THETA_Ka
#'
#'     Cl <- TVCl*exp(ETA_Cl)
#'     Vc <- TVVc*exp(ETA_Vc)
#'     Ka <- TVKa*exp(ETA_Ka)
#'
#'     K20 <- Cl/Vc
#'     Cc <- centr/Vc
#'
#'     d/dt(depot) = -Ka*depot
#'     d/dt(centr) = Ka*depot - K20*centr
#'     Cc ~ prop(prop.sd)
#'   })
#' }
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
                           target_cmin,dose=NULL,cmt_dose=1,endpoint="Cc",
                           estim_method="map",nocb=FALSE,p=NULL,
                           greater_than=TRUE,from=0.2,last_time=72,
                           add_dose=NULL,interdose_interval=NULL,
                           duration=0,indiv_param=NULL){

  prior_model <- get_prior_model(prior_model)
  object <- posologyr(prior_model,dat,nocb)

  #Get or simulate the individual scenario: events and observations-------------
  if(tdm){ #using TDM data
    #input validation
    if (estim_method != "map"){
      warning("estim_method is ignored when tdm=TRUE")
    }
    if (!is.null(dose)){
      warning("dose is ignored when tdm=TRUE")
    }
    if (duration != 0){
      warning("duration is ignored when tdm=TRUE")
    }
    if (!is.null(interdose_interval)){
      warning("interdose_interval is ignored when tdm=TRUE")
    }
    if (!is.null(add_dose)){
      warning("add_dose is ignored when tdm=TRUE")
    }
    if (!is.null(indiv_param)){
      warning("indiv_param is ignored when tdm=TRUE")
    }
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
    #MAP estimation of the individual PK profile from the TDM data
    cmin_map <- poso_estim_map(dat=dat,prior_model=prior_model,nocb=nocb)
    #individual parameters as ETA + covariates
    covar <- utils::tail(dat[,prior_model$covariates],1)
    indiv_param <- cbind(cmin_map$model$params,covar)
    #time of the last dose
    time_last_dose <- utils::tail(cmin_map$event[cmin_map$event$evid == 1,],
                                  1)$time
    #starting time following the last dose
    from <- time_last_dose + from
    #last observation of the poso_estim_map output
    lobs_map <- utils::tail(cmin_map$event,1)$time
    #last observation desired
    lobs <- lobs_map + last_time
    #bind the inital eventTable with enough repetitions of the last row
    # allows to keep EVID==0 and the last known values of the covariates
    # for every obs
    # where .N is a shortcut for "last row"
    extended_et <- rbind(cmin_map$event,
                         cmin_map$event[rep(.N,round((lobs-lobs_map)/0.1))])
    #fill the time column with the desired observation times: seq from the last
    # observation of the MAP output to the last observation needed
    time <- NULL    # avoid undefined global variables
    extended_et[time==lobs_map,time:=seq(from=lobs_map,to=lobs,by=0.1)]
    #Solve the model with the extended et
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
        stop("interdose_interval is mandatory when add_dose is used.")
      }
    }

    #simulation of the individual scenario
    # create an event table with the required number of administrations
    if (!is.null(add_dose)){ #more than one dose is needed
      event_table_cmin <- rxode2::et(amt=dose,dur=duration,cmt=cmt_dose,
                                     ii=interdose_interval,
                                     addl=add_dose)
      time_last_dose   <- add_dose*interdose_interval
      #add observations from the time of the last dose onward only
      event_table_cmin$add.sampling(seq(time_last_dose+from,
                                        time_last_dose+last_time,
                                        by=0.1))
    }
    else { #only one dose is needed: simulation of a single administration
      event_table_cmin <- rxode2::et(amt=dose,dur=duration,cmt=cmt_dose)
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
  if (select_proposal_from_distribution){ #simu. scenario distribution of ETA
    wide_cmin      <- tidyr::pivot_wider(cmin_ppk_model,
                                        id_cols = "time",
                                        names_from = "sim.id",
                                        values_from = endpoint)

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
    # time_to_target: time after the last dose for which tne endpoint is lower
    # than target_cmin
    time_to_target <-
      min(cmin_ppk_model[cmin_ppk_model$time > from &
                         cmin_ppk_model[,endpoint] < target_cmin,
                         ]$time) - time_last_dose
    #here cmin_distribution is a point estimate of the endpoint
    cmin_distribution <-
      cmin_ppk_model[cmin_ppk_model$time == (time_to_target +
                                                  time_last_dose),endpoint]
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

#' Estimate the dose needed to reach a target area under the concentration-time
#' curve (AUC)
#'
#' estimates the dose needed to reach a target area under the concentration-time
#' curve (AUC) given a population pharmacokinetic model, a set of individual
#' parameters, and a target AUC.
#'
#' @param dat Dataframe. An individual subject dataset following the
#'    structure of NONMEM/rxode2 event records.
#' @param prior_model A \code{posologyr} prior population pharmacokinetics
#'    model, a list of six objects.
#' @param tdm A boolean. If `TRUE`: estimates the optimal dose for a selected
#'    target auc over a selected duration following the events from `dat`, and
#'    using Maximum A Posteriori estimation. Setting `tdm` to `TRUE` causes the
#'    following to occur:
#'
#'    * the `time_dose` argument is required and is used as the starting point
#'    for the AUC calculation instead of `starting_time`;
#'    * the arguments `estim_method`, `p`, `greater_than`, `interdose_interval`,
#'    `add_dose`, `indiv_param` and `starting_time` are ignored.
#'
#' @param time_auc Numeric. A duration. The target AUC is computed from
#'    `starting_time` to `starting_time` + `time_auc`.
#'    When `tdm` is set to `TRUE` the target AUC is computed from `time_dose` to
#'    `time_dose` + `time_auc` instead.
#' @param time_dose Numeric. Time when the dose is to be given. Only used and
#'    mandatory, when `tdm` is set to `TRUE`.
#' @param cmt_dose Character or numeric. The compartment in which the dose is
#'    to be administered. Must match one of the compartments in the prior model.
#'    Defaults to 1.
#' @param target_auc Numeric. The target AUC.
#' @param estim_method A character string. An estimation method to be used for
#'    the individual parameters. The default method "map" is the Maximum A
#'    Posteriori estimation, the method "prior" simulates from the prior
#'    population model, and "sir" uses the Sequential Importance Resampling
#'    algorithm to estimate the a posteriori distribution of the individual
#'    parameters. This argument is ignored if `indiv_param` is provided, or if
#'    `tdm` is set to `TRUE`.
#' @param nocb A boolean. for time-varying covariates: the next observation
#'    carried backward (nocb) interpolation style, similar to NONMEM.  If
#'    `FALSE`, the last observation carried forward (locf) style will be used.
#'    Defaults to `FALSE`.
#' @param p Numeric. The proportion of the distribution of AUC to consider for
#'    the optimization. Mandatory for `estim_method=sir`. This argument is
#'    ignored if `tdm` is set to `TRUE`.
#' @param greater_than A boolean. If `TRUE`: targets a dose leading to a
#'    proportion `p` of the AUCs to be greater than `target_auc`. Respectively,
#'    lower if `FALSE`. This argument is ignored if `tdm` is set to `TRUE`.
#' @param starting_time Numeric. First point in time of the AUC, for multiple
#'    dose regimen. The default is zero. This argument is ignored if `tdm` is
#'    set to `TRUE`, and `time_dose` is used as a starting point instead.
#' @param interdose_interval Numeric. Time for the interdose interval for
#'    multiple dose regimen. Must be provided when add_dose is used. This
#'    argument is ignored if `tdm` is set to `TRUE`.
#' @param add_dose Numeric. Additional doses administered at inter-dose interval
#'    after the first dose. Optional. This argument is ignored if `tdm` is set
#'    to `TRUE`.
#' @param duration Numeric. Duration of infusion, for zero-order
#'    administrations.
#' @param starting_dose Numeric. Starting dose for the optimization
#'    algorithm.
#' @param indiv_param Optional. A set of individual parameters : THETA,
#'    estimates of ETA, and covariates. This argument is ignored if `tdm` is
#'    set to `TRUE`.
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
#' rxode2::setRxThreads(2L) # limit the number of threads
#'
#' # model
#' mod_run001 <- function() {
#'   ini({
#'     THETA_Cl <- 4.0
#'     THETA_Vc <- 70.0
#'     THETA_Ka <- 1.0
#'     ETA_Cl ~ 0.2
#'     ETA_Vc ~ 0.2
#'     ETA_Ka ~ 0.2
#'     prop.sd <- sqrt(0.05)
#'   })
#'   model({
#'     TVCl <- THETA_Cl
#'     TVVc <- THETA_Vc
#'     TVKa <- THETA_Ka
#'
#'     Cl <- TVCl*exp(ETA_Cl)
#'     Vc <- TVVc*exp(ETA_Vc)
#'     Ka <- TVKa*exp(ETA_Ka)
#'
#'     K20 <- Cl/Vc
#'     Cc <- centr/Vc
#'
#'     d/dt(depot) = -Ka*depot
#'     d/dt(centr) = Ka*depot - K20*centr
#'     Cc ~ prop(prop.sd)
#'   })
#' }
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
poso_dose_auc <- function(dat=NULL,prior_model=NULL,tdm=FALSE,
                          time_auc,time_dose=NULL,cmt_dose=1,target_auc,
                          estim_method="map",nocb=FALSE,
                          p=NULL,greater_than=TRUE,starting_time=0,
                          interdose_interval=NULL,add_dose=NULL,
                          duration=0,starting_dose=100,indiv_param=NULL){
  prior_model <- get_prior_model(prior_model)
  object <- posologyr(prior_model,dat,nocb)

  #dedicated environment to retrieve variables created inside functions
  hand_made_env <- new.env()

  if(tdm){ #using TDM data
    #input validation
    if (is.null(time_dose)){
      stop("time_dose is mandatory when tdm=TRUE")
    }
    if (estim_method != "map"){
      warning("estim_method is ignored when tdm=TRUE")
    }
    if (starting_time != 0){
      warning("starting_time is ignored when tdm=TRUE")
    }
    if (!is.null(interdose_interval)){
      warning("interdose_interval is ignored when tdm=TRUE")
    }
    if (!is.null(add_dose)){
      warning("add_dose is ignored when tdm=TRUE")
    }
    if (!is.null(indiv_param)){
      warning("indiv_param is ignored when tdm=TRUE")
    }
    #clear all unused input
    estim_method <- NULL
    p <- NULL
    greater_than <- NULL
    add_dose <- NULL
    starting_time <- NULL
    interdose_interval <- NULL
    indiv_param <- NULL
    select_proposal_from_distribution <- FALSE
    #auc starts immediately after the optimized dose
    # starting_time <- time_dose
    ending_time <- time_dose + time_auc
    #MAP estimation of the individual PK profile from the TDM data
    auc_map <- poso_estim_map(dat=dat,prior_model=prior_model,nocb=nocb)
    #individual parameters as ETA + covariates
    covar <- utils::tail(dat[,prior_model$covariates],1)
    indiv_param <- cbind(auc_map$model$params,covar)
    #remove all observations from the eventTable, keep only dosing
    extended_et <- auc_map$event[auc_map$event$evid == 1,]
    #last event of the poso_estim_map output
    last_event <- utils::tail(extended_et,1)$time
    #more input validation
    if (time_dose<=last_event){
      stop("time_dose must occur after the last recorded dosing.")
    }
    #bind the inital eventTable with enough repetitions of the last row
    # one row for dosing
    # one row for the observation at time_dose
    # one row for the observation at ending_time
    extended_et <- rbind(extended_et,
                         extended_et[rep(.N,3)])
    #fill the time column with the desired observation times:
    # the simulated administration, and subsequent concentration
    time <- NULL    # avoid undefined global variables
    extended_et[time>=last_event,time:=c(last_event,time_dose,time_dose,
                                         ending_time)]
    err_dose_tdm <- function(dose,time_dose,target_auc,ending_time,prior_model,
                             extended_et,nocb){
      #time_dose and starting_time are two duplicated rows in the eventTable
      # get their row numbers to assign the correct information to each
      dose_and_start <- extended_et[,which(time==time_dose)]
      #New dose at time_dose in the extended_et
      # Shorter syntax with "list" allows for setting several variable at once
      extended_et[dose_and_start[1],
                  c("evid","amt","cmt","dur"):=list(1,dose,cmt_dose,duration)]
      #New observation at time_dose and ending_time
      extended_et[dose_and_start[2],c("evid","amt","cmt","dur"):=list(0,NA,NA,NA)]
      extended_et[time==ending_time,c("evid","amt","cmt","dur"):=list(0,NA,NA,NA)]
      #Solve the model with the extended et
      auc_ppk_model <- rxode2::rxSolve(prior_model$ppk_model,extended_et,
                                         c(prior_model$theta,auc_map$eta),
                                         covsInterpolation =
                                           ifelse(nocb,"nocb","locf"))

      auc_proposal     <- max(auc_ppk_model$AUC)-min(auc_ppk_model$AUC)
      #assign the proposed AUC to a dedicated environment
      assign("auc_estimate",auc_proposal,
             envir = hand_made_env,
             inherits = FALSE)
      #return the difference between the computed auc and the target
      delta_auc <- (target_auc - auc_proposal)^2
      return(delta_auc)
    }

    optim_dose_auc <- stats::optim(starting_dose,err_dose_tdm,
                                    time_dose=time_dose,target_auc=target_auc,
                                    ending_time=ending_time,
                                    prior_model=prior_model,
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

    if (!is.null(add_dose)){
      if (is.null(interdose_interval)){
        stop("interdose_interval is mandatory when add_dose is used.")
      }
      if (starting_time+time_auc>interdose_interval*add_dose){
        stop("The auc time window is outside of the dosing time range:
           starting_time+time_auc>interdose_interval*add_dose.")
      }
    }
    # Optimization -------------------------------------------------------------
    err_dose <- function(dose,time_auc,starting_time,target_auc,
                         interdose_interval,add_dose,prior_model,
                         duration,indiv_param){

      #compute the individual time-concentration profile
      if (!is.null(add_dose)){
        event_table_auc <- rxode2::et(amt=dose,dur=duration,cmt=cmt_dose,
                                      ii=interdose_interval,
                                      addl=add_dose)
      } else {
        event_table_auc <- rxode2::et(amt=dose,dur=duration,cmt=cmt_dose)
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

        #assign the distribution of auc to a dedicated environment
        assign("auc_estimate",sorted_auc,
               envir = hand_made_env,
               inherits = FALSE)

        if (greater_than){
          auc_proposal <- sorted_auc[n_auc - auc_index]
        } else {
          auc_proposal <- sorted_auc[auc_index]
        }
      } else {
        auc_proposal  <- max(auc_ppk_model$AUC)-min(auc_ppk_model$AUC)
        #assign the proposed AUC to a dedicated environment
        assign("auc_estimate",auc_proposal,
               envir = hand_made_env,
               inherits = FALSE)
      }

      #return the difference between the computed AUC and the target
      delta_auc<-(target_auc - auc_proposal)^2
      return(delta_auc)
    }

    optim_dose_auc <- stats::optim(starting_dose,err_dose,time_auc=time_auc,
                                   starting_time=starting_time,
                                   add_dose=add_dose,
                                   interdose_interval=interdose_interval,
                                   target_auc=target_auc,prior_model=object,
                                   duration=duration,indiv_param=indiv_param,
                                   method="Brent",lower=0,upper=1e5)
  }

  ifelse(select_proposal_from_distribution,
         type_of_estimate <-"distribution",
         type_of_estimate <- "point estimate")

  dose_auc <- list(dose=optim_dose_auc$par,
                   type_of_estimate=type_of_estimate,
                   auc_estimate=hand_made_env$auc_estimate,
                   indiv_param=indiv_param)
  return(dose_auc)
}

#' Estimate the optimal dose to achieve a target concentration at any given time
#'
#' Estimates the optimal dose to achieve a target concentration at any given
#' time given a population pharmacokinetic model, a set of individual
#' parameters, a selected point in time, and a target concentration.
#'
#' @param dat Dataframe. An individual subject dataset following the
#'     structure of NONMEM/rxode2 event records.
#' @param prior_model A \code{posologyr} prior population pharmacokinetics
#'    model, a list of six objects.
#' @param tdm A boolean. If `TRUE`: estimates the optimal dose for a selected
#'    target concentration at a selected point in time following the events from
#'     `dat`, and using Maximum A Posteriori estimation. Setting `tdm` to `TRUE` causes the
#'    following to occur:
#'
#'    * the arguments `estim_method`, `p`, `greater_than`, `interdose_interval`,
#'    `add_dose`, `indiv_param` and `starting_time` are ignored.
#'
#' @param time_c Numeric. Point in time for which the dose is to be
#'     optimized.
#' @param time_dose Numeric. Time when the dose is to be given.
#' @param target_conc Numeric. Target concentration.
#' @param cmt_dose Character or numeric. The compartment in which the dose is
#'    to be administered. Must match one of the compartments in the prior model.
#'    Defaults to 1.
#' @param endpoint Character. The endpoint of the prior model to be optimised
#'    for. The default is "Cc", which is the central concentration.
#' @param estim_method A character string. An estimation method to be used for
#'    the individual parameters. The default method "map" is the Maximum A
#'    Posteriori estimation, the method "prior" simulates from the prior
#'    population model, and "sir" uses the Sequential Importance Resampling
#'    algorithm to estimate the a posteriori distribution of the individual
#'    parameters. This argument is ignored if `indiv_param` is provided or if
#'    `tdm` is set to `TRUE`.
#' @param nocb A boolean. for time-varying covariates: the next observation
#'     carried backward (nocb) interpolation style, similar to NONMEM.  If
#'     `FALSE`, the last observation carried forward (locf) style will be used.
#'     Defaults to `FALSE`.
#' @param p Numeric. The proportion of the distribution of concentrations to
#'    consider for the optimization. Mandatory for `estim_method=sir`. This
#'    argument is ignored if `tdm` is set to `TRUE`.
#' @param greater_than A boolean. If `TRUE`: targets a dose leading to a
#'    proportion `p` of the concentrations to be greater than `target_conc`.
#'    Respectively, lower if `FALSE`. This argument is ignored if `tdm` is
#'    set to `TRUE`.
#' @param starting_dose Numeric. Starting dose for the optimization
#'     algorithm.
#' @param add_dose Numeric. Additional doses administered at inter-dose interval
#'     after the first dose. Optional. This argument is ignored if `tdm` is set
#'     to `TRUE`.
#' @param interdose_interval Numeric. Time for the interdose interval
#'     for multiple dose regimen. Must be provided when add_dose is used. This
#'     argument is ignored if `tdm` is set to `TRUE`.
#' @param duration Numeric. Duration of infusion, for zero-order
#'     administrations.
#' @param indiv_param Optional. A set of individual parameters : THETA,
#'     estimates of ETA, and covariates. This argument is ignored if `tdm` is
#'     set to `TRUE`.
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{dose}{Numeric. An optimal dose for the selected target
#'   concentration.}
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
#' rxode2::setRxThreads(2L) # limit the number of threads
#'
#' # model
#' mod_run001 <- function() {
#'   ini({
#'     THETA_Cl <- 4.0
#'     THETA_Vc <- 70.0
#'     THETA_Ka <- 1.0
#'     ETA_Cl ~ 0.2
#'     ETA_Vc ~ 0.2
#'     ETA_Ka ~ 0.2
#'     prop.sd <- sqrt(0.05)
#'   })
#'   model({
#'     TVCl <- THETA_Cl
#'     TVVc <- THETA_Vc
#'     TVKa <- THETA_Ka
#'
#'     Cl <- TVCl*exp(ETA_Cl)
#'     Vc <- TVVc*exp(ETA_Vc)
#'     Ka <- TVKa*exp(ETA_Ka)
#'
#'     K20 <- Cl/Vc
#'     Cc <- centr/Vc
#'
#'     d/dt(depot) = -Ka*depot
#'     d/dt(centr) = Ka*depot - K20*centr
#'     Cc ~ prop(prop.sd)
#'   })
#' }
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
                           time_c,time_dose=NULL,target_conc,cmt_dose=1,
                           endpoint="Cc",estim_method="map",nocb=FALSE,p=NULL,
                           greater_than=TRUE,starting_dose=100,
                           interdose_interval=NULL,add_dose=NULL,duration=0,
                           indiv_param=NULL){
  prior_model <- get_prior_model(prior_model)
  object <- posologyr(prior_model,dat,nocb)

  #dedicated environment to retrieve variables created inside functions
  hand_made_env <- new.env()

  if(tdm){ #using TDM data
    #input validation
    if (time_c<=time_dose){
      stop("time_c cannot be before time_dose.")
    }
    if (is.null(time_dose)){
      stop("time_dose is mandatory when tdm=TRUE")
    }
    if (estim_method != "map"){
      warning("estim_method is ignored when tdm=TRUE")
    }
    if (!is.null(interdose_interval)){
      warning("interdose_interval is ignored when tdm=TRUE")
    }
    if (!is.null(add_dose)){
      warning("add_dose is ignored when tdm=TRUE")
    }
    if (!is.null(indiv_param)){
      warning("indiv_param is ignored when tdm=TRUE")
    }
    #clear all unused input
    estim_method <- NULL
    p <- NULL
    greater_than <- NULL
    add_dose <- NULL
    interdose_interval <- NULL
    indiv_param <- NULL
    select_proposal_from_distribution <- FALSE
    #MAP estimation of the individual PK profile from the TDM data
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
      stop("time_dose must occur after the last recorded dosing.")
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
      # Shorter syntax with "list" allows for setting several variable at once
      extended_et[time==time_dose,
                  c("evid","amt","cmt","dur"):=list(1,dose,cmt_dose,duration)]
      #New observation at time_c
      extended_et[time==time_c,c("evid","amt","cmt","dur"):=list(0,NA,NA,NA)]
      #Solve the model with the extended et
      ctime_ppk_model <- rxode2::rxSolve(prior_model$ppk_model,extended_et,
                                         c(prior_model$theta,ctime_map$eta),
                                         covsInterpolation =
                                           ifelse(nocb,"nocb","locf"))

      conc_proposal     <- ctime_ppk_model[,endpoint]
      #assign the proposed concentration to a dedicated environment
      assign("conc_distribution",conc_proposal,
             envir = hand_made_env,
             inherits = FALSE)
      #return the difference between the computed ctime and the target
      delta_conc <- (target_conc - conc_proposal)^2
      return(delta_conc)
    }

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
        stop("interdose_interval is mandatory when add_dose is used.")
      }
      if (time_c>(add_dose*interdose_interval)){
        stop("The target time is outside of the dosing time range:
           time_c>(add_dose*interdose_interval).")
      }
    }
    err_dose <- function(dose,time_c,target_conc,prior_model,
                         add_dose,interdose_interval,
                         duration=duration,indiv_param){

      #compute the individual time-concentration profile
      if (!is.null(add_dose)){
        event_table_ctime <- rxode2::et(amt=dose,dur=duration,cmt=cmt_dose,
                                        ii=interdose_interval,
                                        addl=add_dose)
      }
      else {
        event_table_ctime <- rxode2::et(amt=dose,dur=duration,cmt=cmt_dose)
      }
      event_table_ctime$add.sampling(time_c)

      ctime_ppk_model <- rxode2::rxSolve(object=prior_model$ppk_model,
                                         params=indiv_param,
                                         event_table_ctime,
                                         nDisplayProgress=1e5)

      if (select_proposal_from_distribution){

        sorted_conc     <- sort(ctime_ppk_model[,endpoint])
        n_conc          <- length(sorted_conc)
        conc_index      <- ceiling(p * n_conc)

        #assign the distribution of concentrations to a dedicated environment
        assign("conc_distribution",sorted_conc,
               envir = hand_made_env,
               inherits = FALSE)

        if (greater_than){
          conc_proposal <- sorted_conc[n_conc - conc_index]
        } else {
          conc_proposal <- sorted_conc[conc_index]
        }
      } else {
        conc_proposal     <- ctime_ppk_model[,endpoint]
        #assign the proposed concentration to a dedicated environment
        assign("conc_distribution",conc_proposal,
               envir = hand_made_env,
               inherits = FALSE)
      }

      #return the difference between the computed ctime and the target
      delta_conc <- (target_conc - conc_proposal)^2
      return(delta_conc)
    }

    optim_dose_conc <- stats::optim(starting_dose,err_dose,time_c=time_c,
                                    target_conc=target_conc,prior_model=object,
                                    add_dose=add_dose,
                                    interdose_interval=interdose_interval,
                                    duration=duration,indiv_param=indiv_param,
                                    method="Brent",lower=0, upper=1e5)
  }

  conc_distribution <- hand_made_env$conc_distribution

  ifelse(select_proposal_from_distribution,
         type_of_estimate <-"distribution",
         type_of_estimate <- "point estimate")

  dose_conc <- list(dose=optim_dose_conc$par,
                   type_of_estimate=type_of_estimate,
                   conc_estimate=conc_distribution,
                   indiv_param=indiv_param)
  return(dose_conc)
}

#' Estimate the optimal dosing interval to consistently achieve a target trough
#' concentration (Cmin)
#'
#' Estimates the optimal dosing interval to consistently achieve a target Cmin,
#' given a dose, a population pharmacokinetic model, a set of individual
#' parameters, and a target concentration.
#'
#' @param dat Dataframe. An individual subject dataset following the
#'     structure of NONMEM/rxode2 event records.
#' @param prior_model A \code{posologyr} prior population pharmacokinetics
#'    model, a list of six objects.
#' @param target_cmin Numeric. Target trough concentration (Cmin).
#' @param dose Numeric. The dose given.
#' @param cmt_dose Character or numeric. The compartment in which the dose is
#'    to be administered. Must match one of the compartments in the prior model.
#'    Defaults to 1.
#' @param endpoint Character. The endpoint of the prior model to be optimised
#'    for. The default is "Cc", which is the central concentration.
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
#' rxode2::setRxThreads(2L) # limit the number of threads
#'
#' # model
#' mod_run001 <- function() {
#'   ini({
#'     THETA_Cl <- 4.0
#'     THETA_Vc <- 70.0
#'     THETA_Ka <- 1.0
#'     ETA_Cl ~ 0.2
#'     ETA_Vc ~ 0.2
#'     ETA_Ka ~ 0.2
#'     prop.sd <- sqrt(0.05)
#'   })
#'   model({
#'     TVCl <- THETA_Cl
#'     TVVc <- THETA_Vc
#'     TVKa <- THETA_Ka
#'
#'     Cl <- TVCl*exp(ETA_Cl)
#'     Vc <- TVVc*exp(ETA_Vc)
#'     Ka <- TVKa*exp(ETA_Ka)
#'
#'     K20 <- Cl/Vc
#'     Cc <- centr/Vc
#'
#'     d/dt(depot) = -Ka*depot
#'     d/dt(centr) = Ka*depot - K20*centr
#'     Cc ~ prop(prop.sd)
#'   })
#' }
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
                            cmt_dose=1,endpoint="Cc",estim_method="map",
                            nocb=FALSE,p=NULL,greater_than=TRUE,
                            starting_interval=12,add_dose=10,duration=0,
                            indiv_param=NULL){

  prior_model <- get_prior_model(prior_model)

  object <- posologyr(prior_model,dat,nocb)

  #dedicated environment to retrieve variables created inside functions
  hand_made_env <- new.env()

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
    event_table_cmin <- rxode2::et(amt=dose,dur=duration,cmt=cmt_dose,
                                  ii=interdose_interval,
                                  addl=add_dose)
    event_table_cmin$add.sampling(interdose_interval*add_dose-0.1)

    cmin_ppk_model <- rxode2::rxSolve(object=prior_model$ppk_model,
                                     params=indiv_param,
                                     event_table_cmin,
                                     nDisplayProgress=1e5)

    if (select_proposal_from_distribution){

      sorted_cmin     <- sort(cmin_ppk_model[,endpoint])
      n_cmin          <- length(sorted_cmin)
      cmin_index      <- ceiling(p * n_cmin)

      #assign the estimated cmin to a dedicated environment
      assign("cmin_estimate",sorted_cmin,
             envir = hand_made_env,
             inherits = FALSE)

      if (greater_than){
        cmin_proposal <- sorted_cmin[n_cmin - cmin_index]
      } else {
        cmin_proposal <- sorted_cmin[cmin_index]
      }
    } else {
      cmin_proposal     <-  cmin_ppk_model[,endpoint]
      #assign the estimated cmin to a dedicated environment
      assign("cmin_estimate",cmin_proposal,
             envir = hand_made_env,
             inherits = FALSE)
    }

    #return the difference between the computed cmin and the target,
    # normalized by the computed cmin to avoid divergence of the algorithm
    delta_cmin <- ((target_cmin - cmin_proposal)/cmin_proposal)^2
    return(delta_cmin)
  }

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
                    conc_estimate=hand_made_env$cmin_estimate,
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
  prior_model <- get_prior_model(prior_model)
  if (is.null(indiv_param)){ #theta_pop + estimates of eta + covariates
    if (estim_method=="map"){
      model_map <- poso_estim_map(dat,prior_model,nocb=nocb,return_model=TRUE)
      if(is.null(object$covariates)){
        indiv_param <- model_map$model$params
      } else {
        covar <- as.data.frame(object$tdm_data[length(object$tdm_data[,1]),
                                                     object$covariates])
        names(covar) <- object$covariates
        indiv_param <- cbind(model_map$model$params,covar,row.names=NULL)
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
        model_pop   <- poso_simu_pop(dat,prior_model,
                                     n_simul=1e5,return_model=TRUE)
        select_proposal_from_distribution <- TRUE
      } else {
        model_pop   <- poso_simu_pop(dat,prior_model,
                                     n_simul=0,return_model=TRUE)
        select_proposal_from_distribution <- FALSE
      }
      if(is.null(object$covariates)){
        indiv_param <- model_pop$model$params
      } else {
        covar <- as.data.frame(object$tdm_data[length(object$tdm_data[,1]),
                                               object$covariates])
        names(covar) <- object$covariates
        indiv_param <- cbind(model_pop$model$params,covar,row.names=NULL)
      }
    } else if (estim_method=="sir"){
      if (p < 0 || p >= 1){
        stop('p must be between 0 and 1')
      }
      model_sir   <- poso_estim_sir(dat,prior_model,nocb=nocb,
                                    n_sample=1e5,n_resample=1e4,
                                    return_model=TRUE)
      if(is.null(object$covariates)){
        indiv_param <- model_sir$model$params
      } else {
        covar <- as.data.frame(object$tdm_data[length(object$tdm_data[,1]),
                                               object$covariates])
        names(covar) <- object$covariates
        indiv_param <- cbind(model_sir$model$params,covar,row.names=NULL)
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
        sprintf("In order to perform the optimization using a parameter
        distribution, you need at least 1000 parameter samples. Only the first
        set of parameters will be used.")
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
