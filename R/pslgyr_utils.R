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


# select a suitable vector of ETA to start the optimization
init_eta <- function(object,estim_with_iov,omega_iov=NULL){

  dat           <- object$tdm_data
  solved_model  <- object$solved_ppk_model

  theta         <- rbind(object$theta)

  omega         <- object$omega
  ind_eta       <- which(diag(omega)>0)
  ifelse(estim_with_iov,
         omega_eta<-omega_iov,
         omega_eta<-omega[ind_eta,ind_eta])
  solve_omega   <- try(solve(omega_eta))

  pimat        <- object$pi_matrix
  ind_kappa    <- which(diag(pimat)>0)
  pimat_kappa  <- pimat[ind_kappa,ind_kappa]
  pimat_names  <- attr(pimat_kappa,"dimnames")[[1]]
  pimat_dim    <- ncol(pimat_kappa)

  sigma         <- object$sigma
  interpolation <- object$interpolation

  y_obs         <- dat$DV[dat$EVID == 0]     # only observations
  error_model   <- object$error_model

  n_sample      <- 1e2

  #simulate n ETAs
  eta_sim       <- mvtnorm::rmvnorm(n_sample,mean=rep(0,ncol(omega_eta)),
                                    sigma=omega_eta)

  eta_df        <- data.frame(eta_sim)
  names(eta_df) <- attr(omega[ind_eta,ind_eta],"dimnames")[[1]]
  eta_dt        <- data.table::data.table(eta_df)

  param_cols    <- attr(omega[ind_eta,ind_eta],"dimnames")[[1]]
  params        <- cbind(ID=1,eta_dt[,param_cols,with=F],theta)

  if (estim_with_iov){
    eta_dt[,ID:=(1:n_sample)]                           # 1:n_samples ID

    dat_dt <- data.table::data.table(dat)
    dat_dt <- dat_dt[rep(dat_dt[,.I],n_sample)]         # one table per sample
    dat_dt[,ID:=rep(1:n_sample,1,each=nrow(dat))]       # 1:n_samples IDs

    # bind random effects to patient data.table
    ID     <- NULL    # avoid undefined global variables
    dat_dt <- dat_dt[eta_dt,on = list(ID = ID), roll = TRUE]

    # everything but ID, OCC and kappas
    tdm_dt            <- data.table::data.table(dat)
    names_tdm_dt_drop <- names(tdm_dt[,!c("ID","OCC")])
    names_tdm_dt_full <- names(tdm_dt)
    drop_cols         <- c(names_tdm_dt_drop,
                           attr(omega[ind_eta,ind_eta],"dimnames")[[1]])

    # reshape IOV to colums
    iov_col <- t(apply(dat_dt[,!drop_cols,with=F],
                       MARGIN=1,
                       FUN=link_kappa_to_occ,
                       pimat_dim=pimat_dim,
                       pimat_names=pimat_names))
    iov_col_dt <- data.table::data.table(iov_col)

    data_iov   <- cbind(dat_dt[,names_tdm_dt_full,with=F],
                        iov_col_dt[,!c("ID","OCC")])

    param_cols <- c("ID",attr(omega,"dimnames")[[1]])
    params     <- cbind(eta_dt[,param_cols,with=F],theta)

    #group_index   <- data.frame(cbind(c(1:10),10,n_sample,nrow(dat)))

    long_cc <- RxODE::rxSolve(solved_model,
                              cbind(theta,eta_dt,row.names=NULL),
                              data_iov,covs_interpolation=interpolation,
                              returnType="data.table")

    if(length(unique(dat$TIME))!=(length(dat$TIME))){
      wide_cc <- dcast(long_cc, formula = id ~ time, value.var = "Cc",
                       fun.aggregate = list(parent=function(x){x[1]},
                                            metabolite=function(x){x[2]}))
    } else {
      wide_cc <- dcast(long_cc, formula = id ~ time, value.var = "Cc")
    }

    wide_cc <- wide_cc[complete.cases(wide_cc[,2])]
  } else {
    long_cc <- RxODE::rxSolve(solved_model,
                              cbind(theta,eta_dt,row.names=NULL),
                              dat,covs_interpolation=interpolation,
                              returnType="data.table")

    if(length(unique(dat$TIME))!=(length(dat$TIME))){
      wide_cc <- dcast(long_cc, formula = sim.id ~ time, value.var = "Cc",
                       fun.aggregate = list(parent=function(x){x[1]},
                                            metabolite=function(x){x[2]}))
    } else {
      wide_cc <- dcast(long_cc, formula = sim.id ~ time, value.var = "Cc")
    }

    wide_cc <- wide_cc[complete.cases(wide_cc[,2])]
  }

  LL_func  <- function(simu_obs){ #doi: 10.4196/kjpp.2012.16.2.97
    eta_id   <- simu_obs[1]
    eta      <- eta_sim[eta_id,]
    f        <- simu_obs[-1]
    g        <- error_model(f,sigma)
    minus_LL <- 0.5*objective_function(y_obs=y_obs,f=f,g=g,eta=eta,
                                       solve_omega=solve_omega)
    return(-minus_LL)
  }

  log_likelihood  <- unlist(apply(wide_cc,MARGIN=1,FUN=LL_func))
  start_eta       <- eta_sim[which(log_likelihood == max(log_likelihood)),
                             1:ncol(omega_eta)]

  if(estim_with_iov){
    names(start_eta) <- c(colnames(omega[ind_eta,ind_eta]),
                          seq(1,length(start_eta)-ncol(omega[ind_eta,ind_eta])))
  } else{
    names(start_eta) <- colnames(omega_eta)
  }

  return(start_eta)
}
