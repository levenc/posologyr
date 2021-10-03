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

validate_dat <- function(tdm_data){
  if (is.null(tdm_data$TIME)) stop("TIME is missing from the patient record")
  if (is.null(tdm_data$EVID)) stop("EVID is missing from the patient record")
  if (is.null(tdm_data$AMT)) stop("AMT is missing from the patient record")
  if (is.null(tdm_data$DV)) stop("DV is missing from the patient record")
}

validate_priormod <- function(priormod){
  if (is.null(priormod$ppk_model)) stop("ppk_model is missing from the prior
                                        posologyr model")
  if (is.null(priormod$error_model)) stop("error_model is missing from the prior
                                          posologyr model")
  if (is.null(priormod$theta)) stop("theta is missing from the prior
                                    posologyr model")
  if (is.null(priormod$omega)) stop("omega is missing from the prior
                                    posologyr model")
  if (is.null(priormod$covariates)) stop("covariates is missing from the prior
                                         posologyr model")
  if (is.null(priormod$sigma)) stop("sigma is missing from the prior
                                    posologyr model")

  # ETA == 0 are not needed in ppk_model$params, only check for ETA > 0
  ind_eta      <- which(diag(priormod$omega)>0)
  omega_eta    <- priormod$omega[ind_eta,ind_eta]

  if (FALSE %in% (attr(omega_eta,"dimnames")[[1]] %in% priormod$ppk_model$params)){
    stop("The names of the omega matrix do not match the parameters of ppk_model")
  }
  if (FALSE %in% (attr(priormod$theta,"names") %in% priormod$ppk_model$params)){
    stop("The names of the theta vector do not match the parameters of ppk_model")
  }
}

check_for_iov <- function(object){
  if (is.null(object$tdm_data$OCC) & !is.null(object$pi_matrix)){
    stop("OCC is missing from the patient record")
  } else if (is.null(object$tdm_data$OCC) | is.null(object$pi_matrix)){
    estim_with_iov <- FALSE
  } else {
    estim_with_iov <- TRUE
  }
  return(estim_with_iov)
}
