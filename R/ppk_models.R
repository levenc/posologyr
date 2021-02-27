#' Easy loading of model and event record for posologyr
#'
#' Creates, or renames, and loads in the global environment the
#' objects needed for posologyr functions: prior model, solved
#' model, and individual event record. The objects are assigned
#' with names corresponding to the defaults of posologyr functions.
#'
#' @param prior_model A posologyr prior population pharmacokinetics model, a
#'    list of seven elements (see 'Details' for the description of the
#'    object)
#' @param dat Dataframe. An individual subject dataset following the
#'     structure of NONMEM/RxODE event records
#'
#' @details
#'
#' The posologyr prior population pharmacokinetics model is a list of
#' seven elements:
#' \describe{
#'  \item{$description}{A brief description of the population
#'      pharmacokinetics model}
#'  \item{$reference}{A named character vector. Bibliographic reference
#'      of the model (DOI)}
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
#'
#' @return Assigns 3 objects to \code{.GlobalEnv}, the global environment:
#' \describe{
#'  \item{prior_ppk_model}{The posologyr prior population pharmacokinetics
#'   model given as `prior_model` parameter}
#'  \item{dat_posology}{A dataframe. The individual subject dataset
#'   given as `dat` parameter}
#'  \item{solved_ppk_model}{An \code{\link[RxODE]{rxSolve}} solve object,
#'   created with `prior_ppk_model` and using `dat` as the event record.}
#' }
#'
#' @examples
#' load_ppk_model(prior_model=mod_amikacin_2cpt_delattre,dat=df_michel)
#'
#' @export
load_ppk_model <- function(prior_model=NULL,dat=NULL){
  solved_ppk_model <- RxODE::rxSolve(prior_model$ppk_model,
                                     prior_model$pk_prior$reference,
                                     dat)

  assign("prior_ppk_model", prior_model, envir = .GlobalEnv)
  assign("solved_ppk_model", solved_ppk_model, envir = .GlobalEnv)
  assign("dat_posology", dat, envir = .GlobalEnv)

  if (requireNamespace("crayon", quietly = TRUE)) {
    bold_green <- crayon::combine_styles("bold","green")
    cat(" Full model + prior information loaded as",
        bold_green("prior_ppk_model"),"\n","Solved model created as",
        bold_green("solved_ppk_model"),"\n","Dataset loaded as",
        bold_green("dat_posology"),"\n")
  } else {
    cat(" Full model + prior information loaded as",
        "prior_ppk_model","\n","Solved model created as",
        "solved_ppk_model","\n","Dataset loaded as",
        "dat_posology","\n")
  }
}

#' Residuel error model combined 1
#'
#' Residual error model combined 1. Constant error model
#' if no proportional coefficient is provided. Proportional
#' error model if no constant (or additive) error
#' coefficient is provided.
#'
#' @param f Numeric vector, output of a pharmacokinetic model
#' @param xi Numeric vector of the coefficients for the
#' residual error model
#'
#' @details Implements the following function:
#' \code{g <- xi[1] + xi[2]*f}
#'
#' @return Numeric vector, residual error
#' @export
error_model_comb1 <- function(f,xi){
  g <- xi[1] + xi[2]*f
  return(g)
}

#' Residuel error model combined 2
#'
#' Residual error model combined 2.
#'
#' @param f Numeric vector, output of a pharmacokinetic model
#' @param xi Numeric vector of the coefficients for the
#' residual error model
#'
#' @details Implements the following function:
#' \code{g <- sqrt(xi[1]^2 + xi[2]^2*f^2)}
#'
#' @return Numeric vector, residual error
#' @export
error_model_comb2 <- function(f,xi){
  g <- sqrt(xi[1]^2 + xi[2]^2*f^2)
  return(g)
}
