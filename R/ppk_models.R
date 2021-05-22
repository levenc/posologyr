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

#' Creates a posologyr list from a prior model and an individual event
#' record
#'
#' Creates a list for a \code{posologyr} prior model, an individual event
#' record, and an \code{\link[RxODE]{rxSolve}} solve object,
#' created from the prior ppk model and the individual event record.
#'
#' @param prior_model A \code{posologyr} prior population pharmacokinetics model, a
#'    list of six objects.
#' @param dat Dataframe. An individual subject dataset following the
#'     structure of NONMEM/RxODE event records.
#'
#' \code{posologyr} will check the validity of the compiled RxODE
#' model. If \code{prior_model$ppk_model$isValid()} returns \code{FALSE},
#' \code{posologyr} will call \code{\link[RxODE]{RxODE}}
#' to recompile the model before solving it using \code{prior_model}
#' and \code{dat}.
#'
#' @return A list of eight objects: the posologyr prior population
#' pharmacokinetics model given as `prior_model` parameter (6 objects),
#' the individual event record (\code{tdm_data}) given as `dat` parameter,
#' and the solved model (\code{solved_ppk_model}).
#' \describe{
#'  \item{ppk_model}{A RxODE model implementing the structural
#'      population pharmacokinetics model with the individual model
#'      (i.e. the model of inter-individual variability) and the
#'      covariates}
#'  \item{error_model}{A function of the residual error model}
#'  \item{theta}{A named vector of the population estimates of the
#'      fixed effects parameters (called THETAs, following NONMEM
#'      terminology)}
#'  \item{omega}{A named square variance-covariance matrix of the
#'      population parameters inter-individual variability}
#'  \item{covariates}{A character vector of the covariates of
#'      the model}
#'  \item{sigma}{The estimates of the parameters of the residual error model}
#'  \item{tdm_data}{A dataframe. The individual subject dataset
#'   given as `dat` parameter}
#'  \item{solved_ppk_model}{An \code{\link[RxODE]{rxSolve}} solve object,
#'   created with `prior_ppk_model` and using `dat` as the event record.}
#' }
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
#' posologyr(prior_model=mod_tobramycin_2cpt_fictional,
#'                                 dat=df_patient01)
#'
#' @export
posologyr <- function(prior_model=NULL,dat=NULL){

  # check the validity of the compiled model and call RxODE::RxODE
  # on invalid models
    if (!prior_model$ppk_model$isValid()){
    cat("Invalid RxODE model, trying to recompile...")
    prior_model$ppk_model <-
      try(RxODE::RxODE(prior_model$ppk_model), silent=TRUE)
    if (prior_model$ppk_model$isValid()){
      cat("Success","\n")
    } else {
      stop("Failed. The RxODE model is still invalid. Aborting",
           call. = FALSE)
    }
  }
  solved_ppk_model <- RxODE::rxSolve(prior_model$ppk_model,
                                     c(prior_model$theta,
                                       diag(prior_model$omega)*0),
                                     dat)

  # assign the objects to a single list
  prior_model$tdm_data    <- dat
  prior_model$solved_ppk_model <- solved_ppk_model

  return(prior_model)
}

#' Residual error model combined 1
#'
#' Residual error model combined 1. Constant error model
#' if no proportional coefficient is provided. Proportional
#' error model if no constant (or additive) error
#' coefficient is provided.
#'
#' @param f Numeric vector, output of a pharmacokinetic model
#' @param sigma Numeric vector of the coefficients for the
#' residual error model
#'
#' @details Implements the following function:
#' \code{g <- sigma[1] + sigma[2]*f}
#'
#' @return Numeric vector, residual error
#' @export
error_model_comb1 <- function(f,sigma){
  g <- sigma[1] + sigma[2]*f
  return(g)
}

#' Residual error model combined 2
#'
#' Residual error model combined 2.
#'
#' @param f Numeric vector, output of a pharmacokinetic model
#' @param sigma Numeric vector of the coefficients for the
#' residual error model
#'
#' @details Implements the following function:
#' \code{g <- sqrt(sigma[1]^2 + sigma[2]^2*f^2)}
#'
#' @return Numeric vector, residual error
#' @export
error_model_comb2 <- function(f,sigma){
  g <- sqrt(sigma[1]^2 + sigma[2]^2*f^2)
  return(g)
}
