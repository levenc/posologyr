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

#' FICTIONAL IV tobramycin model for test purpose
#'
#' A FICTIONAL two-compartment population pharmacokinetic
#' model of intravenous tobramycin.
#'
#' @format A  list of six objects:
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
#' }
#' @source \url{https://www.page-meeting.org/pdf_assets/1954-2017_05_01_poster_Tobramycin.pdf}
"mod_tobramycin_2cpt_fictional"

#' Two-compartment model of IV amikacin (Burdet 2015)
#'
#' Population pharmacokinetics model from a study of single-dose amikacin
#' in critically ill patients with suspected ventilator-associated pneumonia
#'
#' @format A  list of six objects:
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
#' }
#' @details
#' Covariates included in the model:
#' \describe{
#'  \item{CLCREAT4H}{4-h creatinine clearance in ml/min}
#'  \item{TBW}{Total body weight in kg}
#'  \item{PoverF}{PaO2/FIO2 ratio in mmHg}
#' }
#' @source \doi{10.1007/s00228-014-1766-y}
"mod_amikacin_2cpt_Burdet2015"

#' Two-compartment model of IV vancomycin (Goti 2018)
#'
#' Population pharmacokinetics model from a study of hospitalized
#' patients with and without hemodialysis
#'
#' @format A  list of six objects:
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
#' }
#' @details
#' Covariates included in the model:
#' \describe{
#'  \item{CLCREAT}{Creatinine clearance in ml/min calculated by the
#'  Cockroft-Gault formula}
#'  \item{WT}{Body weight in kg}
#'  \item{DIAL}{Hemodialysis status: 0 in the case of non-dialysis
#'  subjects, 1 in the case of dialysis subjects}
#' }
#' @source \doi{10.1097/FTD.0000000000000490}
"mod_vancomycin_2cpt_Goti2018"
