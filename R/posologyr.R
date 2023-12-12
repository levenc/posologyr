#-------------------------------------------------------------------------
# posologyr: individual dose optimisation using population PK
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

#' Get prior model making sure to account for rxode2 ui
#'
#' @param prior_model prior model list or ui
#' @return  prior model lis
#' @noRd
#' @author Matthew L. Fidler
get_prior_model <- function(prior_model, auc=FALSE) {
  if (!is.list(prior_model)) {
    rxui <- try(rxode2::as.rxUi(prior_model), silent=FALSE)
    if (inherits(rxui, "rxUi")) {
      if (auc) {
        prior_model <- rxui$posologyr_auc
      } else {
        prior_model <- rxui$posologyr
      }
    }
  }
  prior_model
}

#' Creates a posologyr list from a prior model and an individual event
#' record
#'
#' Creates a list for a \code{posologyr} prior model, an individual event
#' record, and an \code{\link[rxode2]{rxSolve}} solve object,
#' created from the prior ppk model and the individual event record.
#'
#' @param prior_model A \code{posologyr} prior population pharmacokinetics
#'    model, a list of six objects.
#' @param dat Dataframe. An individual subject dataset following the
#'     structure of NONMEM/rxode2 event records.
#' @param nocb A boolean. for time-varying covariates: the next observation
#'     carried backward (nocb) interpolation style, similar to NONMEM.  If
#'     `FALSE`, the last observation carried forward (locf) style will be used.
#'     Defaults to `FALSE`.
#'
#' \code{posologyr} will check the validity of the compiled rxode2
#' model. If \code{prior_model$ppk_model$isValid()} returns \code{FALSE},
#' \code{posologyr} will call \code{\link[rxode2]{rxode2}}
#' to recompile the model before solving it using \code{prior_model}
#' and \code{dat}.
#'
#' @return A list of eight objects: the posologyr prior population
#' pharmacokinetics model given as `prior_model` parameter (6 objects),
#' the individual event record (\code{tdm_data}) given as `dat` parameter,
#' and the solved model (\code{solved_ppk_model}).
#' \describe{
#'  \item{ppk_model}{A rxode2 model implementing the structural
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
#'  \item{solved_ppk_model}{An \code{\link[rxode2]{rxSolve}} solve object,
#'   created with `prior_ppk_model` and using `dat` as the event record.}
#' }
#'
#' @keywords internal
#' @noRd
posologyr <- function(prior_model=NULL,dat=NULL,nocb=FALSE){
  prior_model <- get_prior_model(prior_model)
  validate_priormod(prior_model)
  validate_dat(dat)

  # check the validity of the compiled model and call rxode2::rxode
  # on invalid models
    if (!prior_model$ppk_model$isValid()){
      warning("Invalid rxode2 model, trying to recompile...")
      prior_model$ppk_model <- try(rxode2::rxode(prior_model$ppk_model),
                                   silent=TRUE)
      if (prior_model$ppk_model$isValid()){
        message("Success","\n")
       } else {
         stop("Failed. The rxode2 model is still invalid. Aborting",
              call. = FALSE)
       }
    }

  # interpolation method for time-varying covariates, nocb or locf
  interpolation <- ifelse(nocb,"nocb","locf")

  solved_ppk_model <- rxode2::rxSolve(prior_model$ppk_model,
                                     c(prior_model$theta,
                                       diag(prior_model$omega)*0,
                                       diag(prior_model$pi_matrix)*0),
                                     dat,covsInterpolation=interpolation)

  # assign the objects to a single list
  prior_model$tdm_data         <- as.data.frame(dat)
  prior_model$solved_ppk_model <- solved_ppk_model
  prior_model$interpolation    <- interpolation

  return(prior_model)
}
