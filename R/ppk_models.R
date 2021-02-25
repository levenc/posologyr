#' Easy loading of model and event record for posologyr
#'
#' Creates and loads in the global environment the suitable objects for the
#' use of the posologyr functions.
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
#' @return A list of 2: a dataframe of the simulated population
#'     pharmacokinetic parameters, and a vector of the prior typical
#'     values of the population pharmacokinetic parameters
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

# error models -----------------------------------------------------------

# combined1 (also prop and const error models)
error_model_comb1 <- function(f,xi){
  g <- xi[1] + xi[2]*f
  return(g)
}

# combined2
error_model_comb2 <- function(f,xi){
  g <- sqrt(xi[1]^2 + xi[2]^2*f^2)
  return(g)
}

# popPK structural models-------------------------------------------------
# Fictionnal tobramycin 2cpt: example
mod_tobramycin_2cpt_fictionnal <- list(
  description = c("Fictionnal tobramycin model for test purposes,
                  based on
  https://www.page-meeting.org/pdf_assets/1954-2017_05_01_poster_Tobramycin.pdf"),
  reference = NULL,
  ppk_model   = RxODE::RxODE({
    centr(0) = 0
    ke = TVke*(CLCREAT/67.8)^0.89*(WT/66.4)^-1.09
    V  = TVV*(WT/66.4)^0.80
    Cc  = centr/V;
    d/dt(centr)  = - ke*centr - k12*centr + k21*periph;
    d/dt(periph) =            + k12*centr - k21*periph;
    d/dt(AUC)    =   Cc;
  }),
  error_model = error_model_comb1,
  pk_prior    = list( name = c('TVke','TVV','k12','k21'),
                      reference = c(TVke=0.21, TVV=19.8, k12=0.041, k21=0.12),
                      Omega = matrix(c(0.08075, 0      ,  0, 0,
                                       0      , 0.01203,  0, 0,
                                       0      , 0      ,  0, 0,
                                       0      , 0      ,  0, 0),
                                     ncol=4,byrow=TRUE)),
  covariates  = c("CLCREAT","WT"),
  xi          = c(additive_a = 0, proportional_b = 0.198))

# Amikacin 2cpt critically ill patients with sepsis
mod_amikacin_2cpt_delattre <- list(
  description = c("Amikacin: two-compartment model with first order elimination,
                  in critically ill patients with sepsis,
                  with (CLCREAT) according to Cockroft-Gault equation"),
  reference = c(doi = "10.1097/FTD.0b013e3181f675c2"),
  ppk_model   = RxODE::RxODE({
    centr(0) = 0
    Cl  = TCl+(1.42*CLCREAT);
    k12 = Q/V1;
    k21 = Q/V2;
    ke  = Cl/V1;
    Cc  = centr/V1;
    d/dt(centr)  = - ke*centr - k12*centr + k21*periph;
    d/dt(periph) =            + k12*centr - k21*periph;
    d/dt(AUC)    =   Cc;
  }),
  error_model = error_model_comb1,
  pk_prior    = list( name = c('V1','V2','Q','TCl'),
                      reference = c(V1=19.2, V2=9.38, Q=4.38, TCl=0.77),
                      Omega = matrix(c(0.1423, 0     ,      0,      0,
                                       0     , 0.1740,      0,      0,
                                       0     ,      0, 0.0275,      0,
                                       0     ,      0,      0, 0.2467),
                                     ncol=4,byrow=TRUE)),
  covariates  = c("CLCREAT"),
  xi          = c(additive_a = 1.03, proportional_b = 0.268))

# Imatinib 1cpt CML or GIST patients
mod_imatinib_oral_1cpt_Widmer <- list(
  description = c("Imatinib: one compartment model with first order elimination,
                  in patients with CML (PATHO = 0) or GIST (PATHO = 1),
                  SEX = 1 for male patients."),
  reference = c(doi = "10.1111/j.1365-2125.2006.02719.x"),
  ppk_model   = RxODE::RxODE({
    depot(0) = 0
    centr(0) = 0
    Cl = TCl+54.2*(WT-70)/70+1.49*SEX-1.49*(1-SEX)-5.81*(AGE-50)/50-0.806*PATHO+0.806*(1-PATHO);
    V  = TV+46.2*SEX-46.2*(1-SEX);
    ke = Cl/V;
    Cc = centr/V;
    d/dt(depot) = - ka*depot;
    d/dt(centr) = + ka*depot - ke*centr;
    d/dt(AUC)   =   Cc;
  }),
  error_model = error_model_comb1,
  pk_prior    = list( name = c('TCl','TV','ka'),
                      reference = c(TCl=14.3, TV=347, ka=0.61),
                      Omega = matrix(c(0.1219, 0.1790,      0,
                                       0.1790, 0.3343,      0,
                                       0     ,      0,      0),
                                     ncol=3,byrow=TRUE)),
  covariates  = c("WT", "SEX", "AGE", "PATHO"),
  xi          = c(additive_a = 0, proportional_b = 0.310))

# Vancomycin 2cpt adult patients
mod_vancomycin_2cpt_Thomson <- list(
  description = c("Vancomycin: two-compartment model with first order elimination,
                  in adult patients total, body weight as covariate (WT),
                  (CLCREAT) calculated using TBW and Cockroft equation."),
  reference = c(doi = "10.1093/jac/dkp085"),
  ppk_model   = RxODE::RxODE({
    centr(0) = 0
    Cl  = TCl + (CLCREAT/66) * 0.0154;
    V1  = TV1 * WT;
    V2  = TV2 * WT;
    k12 = Q/V1;
    k21 = Q/V2;
    ke  = Cl/V1;
    Cc  = centr/V1;
    d/dt(centr)  = - ke*centr - k12*centr + k21*periph;
    d/dt(periph) =            + k12*centr - k21*periph;
    d/dt(AUC)    =   Cc;
  }),
  error_model = error_model_comb1,
  pk_prior    = list( name = c('TCl','TV1','TV2','Q'),
                      reference = c(TCl=2.99, TV1=0.675, TV2=0.732, Q=2.28),
                      Omega = matrix(c(0.0704, 0     ,      0,      0,
                                       0     , 0.0223,      0,      0,
                                       0     ,      0, 0.9895,      0,
                                       0     ,      0,      0, 0.2152),
                                     ncol=4,byrow=TRUE)),
  covariates  = c("CLCREAT", "WT"),
  xi          = c(additive_a = 1.6, proportional_b = 0.15))

# Linezolid 2cpt critically ill patients
mod_linezolid_2cpt_Soraluce <- list(
  description = c("Linezolid: two-compartment model with first order elimination,
                  in critically ill patients, with measured creatinin clearance,
                  U*V/P as covariate (CLCREAT)"),
  reference = c(doi = "10.3390/pharmaceutics12010054"),
  ppk_model   = RxODE::RxODE({
    centr(0) = 0
    Cl  = TCl + (CLCREAT/44) * 4.35;
    k12 = Q/V1;
    k21 = Q/V2;
    ke  = Cl/V1;
    Cc  = centr/V1;
    d/dt(centr)  = - ke*centr - k12*centr + k21*periph;
    d/dt(periph) =            + k12*centr - k21*periph;
    d/dt(AUC)    =   Cc;
  }),
  error_model = error_model_comb1,
  pk_prior    = list( name = c('TCl','V1','V2','Q'),
                      reference = c(TCl=2.62, V1=16.2, V2=29.0, Q=71.7),
                      Omega = matrix(c(0.3208, 0     ,      0,      0,
                                       0     , 0.3607,      0,      0,
                                       0     ,      0,      0,      0,
                                       0     ,      0,      0,      0),
                                     ncol=4,byrow=TRUE)),
  covariates  = c("CLCREAT"),
  xi          = c(additive_a = 0.266, proportional_b = 0.159))
