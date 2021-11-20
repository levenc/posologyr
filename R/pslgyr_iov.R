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

# Get propositions for values of kappa and put them in colums to be added
#  to the dataset for RxODE
iov_proposition_as_cols <- function(iov_col=NULL,
                                    dat=NULL,
                                    pimat=NULL,
                                    omega_dim=NULL,
                                    eta_estim=NULL){
  n_iov        <- ncol(iov_col)-1   #minus one because of $OCC

  for (i in seq_along(unique(dat$OCC))){
    start_estim_iov <- omega_dim+1+n_iov*(i-1)
    end_estim_iov   <- start_estim_iov+n_iov-1
    iov_vector_i    <- eta_estim[start_estim_iov:end_estim_iov]

    occ_size   <- length(iov_col[iov_col$OCC == i,1])
    iov_mat_i  <- matrix(iov_vector_i,
                         nrow=occ_size,
                         ncol=n_iov,
                         byrow=TRUE)

    iov_col[iov_col$OCC == i,attr(pimat,"dimnames")[[1]]] <- iov_mat_i
  }
  return(iov_col)
}

# make a single matrix of omega and pi_matrix
merge_covar_matrices <- function(omega_eta=NULL,
                                 omega_dim=NULL,
                                 pimat_dim=NULL,
                                 pimat_kappa=NULL,
                                 dat=NULL){
  matrix_dim   <- omega_dim+pimat_dim*(length(unique(dat$OCC)))
  all_the_mat  <- matrix(0,nrow=matrix_dim,ncol=matrix_dim)
  all_the_mat[1:omega_dim,1:omega_dim] <- omega_eta
  for (i in unique(dat$OCC)){
    if (TRUE){
      start_pi_mat <- omega_dim+pimat_dim*(i-1)+1
      end_pi_mat   <- omega_dim+pimat_dim*(i)
      all_the_mat[start_pi_mat:end_pi_mat,
                  start_pi_mat:end_pi_mat] <- pimat_kappa
    }
  }
  return(all_the_mat)
}

# create colums to store the estimations of KAPPA
init_iov_col <- function(dat=NULL,
                         pimat=NULL){
  iov_col        <- matrix(0,nrow=nrow(dat),ncol=nrow(pimat))
  iov_col        <- data.frame(iov_col,dat$OCC)
  names(iov_col) <- c(attr(pimat,"dimnames")[[1]],"OCC")
  return(iov_col)
}

# reshape kappa from wide to long
link_kappa_to_occ <- function(input,pimat_dim,pimat_names){
  current_id    <- input[1]
  current_occ   <- input[2]
  start_kappa   <- (current_occ*pimat_dim)+1 #+2: after ID and OCC, then -1
  end_kappa     <- start_kappa+pimat_dim-1
  range_kappa   <- start_kappa:end_kappa
  kappas        <- input[range_kappa]
  names(kappas) <- pimat_names
  return(c(current_id,current_occ,kappas))
}
