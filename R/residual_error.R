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

#' Residual error model mixed (idem NONMEM)
#'
#' Mixed residual error model, similar to NONMEM implementation.
#'
#' @param f Numeric vector, output of a pharmacokinetic model
#' @param sigma Matrix of the coefficients for the
#' residual error model
#'
#' @return Numeric vector, residual error
#' @export
error_model_mixednm <- function(f,sigma){
  dv <- cbind(f,1)
  g  <- diag(dv%*%sigma%*%t(dv))
  return(sqrt(g))
}
