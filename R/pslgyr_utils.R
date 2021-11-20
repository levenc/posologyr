#-------------------------------------------------------------------------
# posologyr: individual dose optimisation using population PK
# Copyright (C) 2021  Cyril Leven
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

# solve a large data set group by group, workaround for
# https://github.com/nlmixrdevelopment/RxODE/issues/459
solve_by_groups <- function(index,pkmodel,params,dat,interpolation){
  number_of_observ   <- index[4]
  number_of_subjects <- index[3]
  number_of_groups   <- index[2]
  group_number       <- index[1]

  group_size         <- number_of_subjects/number_of_groups

  start_eta          <- group_size*(group_number-1)+1
  stop_eta           <- start_eta+group_size-1

  start_dat          <- group_size*number_of_observ*(group_number-1)+1
  stop_dat           <- start_dat+(group_size)*number_of_observ-1

  group_model <- RxODE::rxSolve(pkmodel,params[start_eta:stop_eta,],
                                dat[start_dat:stop_dat,],
                                covs_interpolation=interpolation,
                                returnType="data.table")
  return(group_model)
}

# pracma gradient of errpred for L-BFGS-B in optim
optim_gradient <- function(x,
                           run_model=run_model,
                           y_obs=y_obs,
                           theta=theta,
                           ind_eta=ind_eta,
                           sigma=sigma,
                           solve_omega=solve_omega,
                           omega=omega,
                           omega_dim=omega_dim,
                           iov_col=iov_col,
                           pimat=pimat,
                           dat=dat,
                           eta_df=eta_df,
                           model_init=model_init,
                           solved_model=solved_model,
                           error_model=error_model,
                           estim_with_iov=estim_with_iov,
                           interpolation=interpolation,
                           adapt=adapt){
  pracma::grad(errpred,x,
               run_model=run_model,
               y_obs=y_obs,
               theta=theta,
               ind_eta=ind_eta,
               sigma=sigma,
               solve_omega=solve_omega,
               omega=omega,
               omega_dim=omega_dim,
               iov_col=iov_col,
               pimat=pimat,
               dat=dat,
               eta_df=eta_df,
               model_init=model_init,
               solved_model=solved_model,
               error_model=error_model,
               estim_with_iov=estim_with_iov,
               interpolation=interpolation,
               adapt=adapt)
}
