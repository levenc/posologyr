validate_dat <- function(tdm_data){
  if (is.null(tdm_data$TIME)) stop("TIME is missing from the patient record")
  if (is.null(tdm_data$EVID)) stop("EVID is missing from the patient record")
  if (is.null(tdm_data$AMT)) stop("TIME is missing from the patient record")
  if (is.null(tdm_data$DV)) stop("TIME is missing from the patient record")
}

validate_priormod <- function(priormod){
  if (is.null(priormod$ppk_model)) stop("ppk_model is missing from the prior
                                        posologyr model")
  if (is.null(priormod$error_model)) stop("error_model is missing from the prior
                                        posologyr model")
  if (is.null(priormod$theta)) stop("theta is missing from the prior
                                        posologyr model")
  if (is.null(priormod$omega)) stop("error_model is missing from the prior
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
