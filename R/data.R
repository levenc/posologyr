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

#' Fictionnal IV tobramycin model for test purpose
#'
#' A fictional bicompartmental population pharmacokinetic
#' model of intravenous tobramycin.
#'
#' @format A  list of six objects:
#' \describe{
#'  \item{$ppk_model}{A RxODE model implementing the structural
#'      population pharmacokinetics model with the individual model
#'      (i.e. the model of inter-individual variability) and the
#'      covariates}
#'  \item{$error_model}{A function of the residual error model}
#'  \item{$psi}{A named vector of the population estimates of the
#'      fixed effects parameters (called THETAs, following NONMEM
#'      terminology)}
#'  \item{$omega}{A named square variance-covariance matrix of the
#'      population parameters inter-individual variability}
#'  \item{$covariates}{A character vector of the covariates of
#'      the model}
#'  \item{$xi}{The estimates of the parameters of the residual error model}
#' }
#' @source \url{https://www.page-meeting.org/pdf_assets/1954-2017_05_01_poster_Tobramycin.pdf}
"mod_tobramycin_2cpt_fictional"

#' Fictionnal oral amoxicillin model for test purpose
#'
#' A fictional one-compartment population pharmacokinetic
#' model of oral amoxicillin, with non-saturable oral
#' absorption, and non-zero covariance between the
#' inter-individual variability of V and Cl.
#'
#' @format A  list of six objects:
#' \describe{
#'  \item{$ppk_model}{A RxODE model implementing the structural
#'      population pharmacokinetics model with the individual model
#'      (i.e. the model of inter-individual variability) and the
#'      covariates}
#'  \item{$error_model}{A function of the residual error model}
#'  \item{$psi}{A named vector of the population estimates of the
#'      fixed effects parameters (called THETAs, following NONMEM
#'      terminology)}
#'  \item{$omega}{A named square variance-covariance matrix of the
#'      population parameters inter-individual variability}
#'  \item{$covariates}{A character vector of the covariates of
#'      the model}
#'  \item{$xi}{The estimates of the parameters of the residual error model}
#' }
#' @source In-house model created specifically for posologyr
"mod_amoxicillin_oral_1cpt_fictional"
