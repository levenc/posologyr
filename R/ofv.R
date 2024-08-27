#-------------------------------------------------------------------------
# posologyr: individual dose optimization using population PK
# Copyright (C) Cyril Leven
#
#  This program is free software: you can redistribute it and/or modify
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

# Update model predictions with a new set of parameters, for all obs
run_model <- function(x,solved_model=NULL,estim_with_iov=NULL,endpoints=NULL){
  if (!estim_with_iov){ #rxode2 already updated in errpred() if estim_with_iov
    solved_model$params <- x
  }
  return(solved_model[endpoints])
}

# Objective function for the Empirical Bayes Estimates
# doi: 10.4196/kjpp.2012.16.2.97
objective_function <- function(y_obs=NULL,f=NULL,g=NULL,
                               eta=NULL,solve_omega=NULL){

  # 1) When the prediction f is zero, g can be zero (depending on the residual
  # error model).
  # 2) log(0) is NaN, the limit of log(x) when x approaches zero is -Inf,
  # 3) log(1) is zero, and 1^2 is 1
  g[which(g == 0)] <- 1

  U_y   <-  sum(((y_obs - f)/g)^2 + log(g^2))

  # the transpose of a diagonal matrix is itself
  U_eta <- eta %*% solve_omega %*% eta

  if (TRUE %in% is.na(f)){
    # if rxode2 fails to solve the model, the proposed ETA is not optimal,
    # assign a large value to OFV to divert the algorithm from this area
    OFV <- 10^10
  } else {
    OFV <- U_y + U_eta
  }

  return(OFV)
}

# Prediction error to optimize for MAP-EBE
errpred <- function(eta_estim=NULL,
                    run_model=NULL,
                    y_obs=NULL,
                    endpoints=NULL,
                    theta=NULL,
                    ind_eta=NULL,
                    sigma=NULL,
                    solve_omega=NULL,
                    omega=NULL,
                    omega_dim=NULL,
                    iov_col=NULL,
                    pimat=NULL,
                    dat=NULL,
                    solved_model=NULL,
                    error_model=NULL,
                    estim_with_iov=NULL,
                    interpolation=NULL){

  eta          <- diag(omega)*0

  if (estim_with_iov){
    eta[ind_eta] <- eta_estim[1:omega_dim]
    iov_col <- iov_proposition_as_cols(iov_col=iov_col,dat=dat,pimat=pimat,
                                       omega_dim=omega_dim,
                                       eta_estim=eta_estim)
    dat <- data.frame(dat,iov_col)
    solved_model <- rxode2::rxSolve(solved_model,c(theta,eta),dat,
                                   covsInterpolation=interpolation)
  } else {
    eta[ind_eta] <- eta_estim
  }
  #simulated concentrations with the proposed eta estimates
  f_all_endpoints <- do.call(run_model,list(c(theta,eta),
                                            solved_model=solved_model,
                                            estim_with_iov=estim_with_iov,
                                            endpoints=endpoints))

  obs_res <- residual_error_all_endpoints(f_all_endpoints=f_all_endpoints,
                                          y_obs=y_obs,
                                          error_model=error_model,
                                          sigma=sigma,
                                          endpoints=endpoints)

  optimize_me <- objective_function(y_obs=y_obs[,"DV"],
                                    f=obs_res$f_all_endpoints$f,
                                    g=obs_res$g_all_endpoints$g,
                                    eta=eta_estim,
                                    solve_omega=solve_omega)
  return(optimize_me)
}
