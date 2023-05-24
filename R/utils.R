#-------------------------------------------------------------------------
# posologyr: individual dose optimisation using population PK
# Copyright (C) 2022  Cyril Leven
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

  group_model <- rxode2::rxSolve(pkmodel,params[start_eta:stop_eta,],
                                dat[start_dat:stop_dat,],
                                covsInterpolation=interpolation,
                                returnType="data.table")
  return(group_model)
}

# select a suitable vector of ETA to start the optimization
init_eta <- function(object,estim_with_iov,omega_iov=NULL,endpoints=NULL){

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

  if(length(ind_kappa)==1){
    pimat_kappa  <- pimat
    pimat_names  <- attr(pimat_kappa,"dimnames")[[1]]
    pimat_dim    <- 1
  } else{
    pimat_kappa  <- pimat[ind_kappa,ind_kappa]
    pimat_names  <- attr(pimat_kappa,"dimnames")[[1]]
    pimat_dim    <- ncol(pimat_kappa)
  }

  sigma         <- object$sigma
  interpolation <- object$interpolation

  ifelse(setequal(endpoints,"Cc"),
         y_obs <- data.frame(DV=dat[dat$EVID==0,"DV"],DVID="Cc"),
         y_obs <- dat[dat$EVID==0,c("DV","DVID")])

  error_model   <- object$error_model

  n_sample      <- 10

  #simulate n ETAs
  eta_sim       <- mvtnorm::rmvnorm(n_sample,mean=rep(0,ncol(omega_eta)),
                                    sigma=omega_eta)

  if (estim_with_iov){
    # matrix large enough for omega + pi
    eta_mat           <- matrix(0,nrow=n_sample,
                                ncol=ncol(omega_iov)-
                                ncol(omega[ind_eta,ind_eta])+
                                ncol(omega))

    # IIV
    eta_mat[,ind_eta] <-
      eta_sim[,1:length(ind_eta)]

    # IOV
    eta_mat[,(ncol(omega)+1):ncol(eta_mat)] <-
      eta_sim[,(length(ind_eta)+1):ncol(eta_sim)]

  } else{
    eta_mat           <- matrix(0,nrow=n_sample,ncol=ncol(omega))
    eta_mat[,ind_eta] <- eta_sim

  }

  eta_df            <- data.frame(eta_mat)

  names(eta_df) <- attr(omega,"dimnames")[[1]]
  eta_dt        <- data.table::data.table(eta_df)

  param_cols    <- attr(omega,"dimnames")[[1]]
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
                           attr(omega,"dimnames")[[1]])

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

    f_all_endpoints <- rxode2::rxSolve(solved_model,
                                       cbind(theta,eta_dt,row.names=NULL),
                                       data_iov,covsInterpolation=interpolation,
                                       returnType="data.table")

    data.table::setnames(f_all_endpoints,"id","sim.id")
  } else {
    f_all_endpoints <- rxode2::rxSolve(solved_model,
                                       cbind(theta,eta_dt,row.names=NULL),
                                       dat,covsInterpolation=interpolation,
                                       returnType="data.table")
  }

  # binding the simulated observations with DVID, the DVID column is recycled
  f_all_endpoints <- data.table::data.table(f_all_endpoints,DVID=y_obs$DVID)
  g_all_endpoints <- f_all_endpoints

  if (setequal(endpoints,"Cc")){    # for retro-compatibility purposes
    g_all_endpoints$Cc <- error_model(f_all_endpoints$Cc,sigma)
  } else {
    for (edp in endpoints){
      g_all_endpoints[,edp] <- error_model[[edp]](as.matrix(f_all_endpoints[,get(edp)]),
                                                 sigma[[edp]])
    }
  }

  f <- g <- DVID <- NULL # avoid undefined global variables

  f_all_endpoints[, f := get(as.character(DVID)),by = seq_len(nrow(f_all_endpoints))]
  g_all_endpoints[, g := get(as.character(DVID)),by = seq_len(nrow(g_all_endpoints))]

  f_all_sim <- dcast(f_all_endpoints, formula = sim.id ~ rowid(sim.id), value.var = "f")
  g_all_sim <- dcast(g_all_endpoints, formula = sim.id ~ rowid(sim.id), value.var = "g")

  LL_func  <- function(simu_obs){ #doi: 10.4196/kjpp.2012.16.2.97
    eta_id   <- simu_obs[1]
    eta      <- eta_sim[eta_id,]
    f        <- simu_obs[-1]
    g        <- g_all_sim[eta_id,-1]
    minus_LL <- 0.5*objective_function(y_obs=y_obs$DV,f=f,g=g,eta=eta,
                                       solve_omega=solve_omega)
    return(-minus_LL)
  }

  log_likelihood  <- unlist(apply(f_all_sim,MARGIN=1,FUN=LL_func))

  start_eta       <- eta_sim[which(log_likelihood == max(log_likelihood)),
                             1:ncol(omega_eta)]

  if(!is.null(dim(start_eta))){ # if all proposals are equally bad, start_eta
                                # is not a vector, but an array, and the output
                                # of dim() exists, in that case return the first
                                # row
    start_eta <- start_eta[1,]
  }

  if(estim_with_iov){
    names(start_eta) <- c(colnames(omega[ind_eta,ind_eta]),
                          seq(1,length(start_eta)-ncol(omega[ind_eta,ind_eta])))
  } else{
    names(start_eta) <- colnames(omega_eta)
  }

  return(start_eta)
}

#extrapolate covariates to allow more sampling times in rxode2 solved models
extrapol_cov <- function(x=NULL,dat=NULL,covar=NULL,interpol_approx=NULL,
                         f=NULL,event_table=NULL){
  x <- which(covar == x) #input char, return position in the covar vector
  approx_cov <- stats::approxfun(x = dat$TIME,
                          y = dat[[covar[x]]],
                          yleft = dat[[covar[x]]][x],
                          yright = dat[[covar[x]]][length(dat[[covar[x]]])],
                          method = interpol_approx,
                          f = f)
  approx_cov(event_table$time)
}

#extrapolate kappas to allow more sampling times in rxode2 solved models
extrapol_iov <- function(x=NULL,dat=NULL,iov_kappa=NULL,event_table=NULL){
  x <- which(iov_kappa == x) #input char, return position in the kappa vector
  approx_iov <- stats::approxfun(x = dat$TIME,
                          y = dat[[iov_kappa[x]]],
                          yleft = dat[[iov_kappa[x]]][x],
                          yright = dat[[iov_kappa[x]]][length(dat[[iov_kappa[x]]])],
                          method = "constant")
  approx_iov(event_table$time)
}

