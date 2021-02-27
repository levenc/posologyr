#' Fictionnal IV tobramycin model for test purposes
#'
#' A fictional bicompartmental population pharmacokinetic
#' model of intravenous tobramycin.
#'
#' @format A list of five elements:
#' \describe{
#'  \item{$ppk_model}{A RxODE model implementing the structural
#'      population pharmacokinetics model with no inter-individual
#'      variability, or residual error model}
#'  \item{$error_model}{A function of the residual error model}
#'  \item{$pk_prior}{A list of 3. `name`: a character vector of the names
#'      of the population pharmacokinetc paramters, `reference`: a named
#'      vector of the prior typical value of the population paramaters,
#'      `Omega`: a square variance-covariance matrix of the population
#'      parameters inter-individual variability}
#'  \item{$covariates}{A character vector of the covariates of
#'      the model}
#'  \item{$xi}{The estimates of the parameters of the residual error model}
#' }
#' @source \url{https://www.page-meeting.org/pdf_assets/1954-2017_05_01_poster_Tobramycin.pdf}
"mod_tobramycin_2cpt_fictional"

#' Fictionnal oral amoxicillin model for test purposes
#'
#' A fictional one-compartment population pharmacokinetic
#' model of oral amoxicillin, with non-saturable oral
#' absorption, and non-zero covariance between the
#' inter-individual variability of V and Cl.
#'
#' @format A list of five elements:
#' \describe{
#'  \item{$ppk_model}{A RxODE model implementing the structural
#'      population pharmacokinetics model with no inter-individual
#'      variability, or residual error model}
#'  \item{$error_model}{A function of the residual error model}
#'  \item{$pk_prior}{A list of 3. `name`: a character vector of the names
#'      of the population pharmacokinetc paramters, `reference`: a named
#'      vector of the prior typical value of the population paramaters,
#'      `Omega`: a square variance-covariance matrix of the population
#'      parameters inter-individual variability}
#'  \item{$covariates}{A character vector of the covariates of
#'      the model}
#'  \item{$xi}{The estimates of the parameters of the residual error model}
#' }
#' @source In-house model created specifically for posologyr
"mod_amox_oral_1cpt_fictional"
