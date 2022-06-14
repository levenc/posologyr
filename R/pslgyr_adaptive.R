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

# iterative adaptive MAP estimation
adaptive_map <- function(return_AMS_models=NULL,
                         dat=NULL,
                         solved_model=NULL,
                         theta=NULL,
                         omega=NULL,
                         start_eta=NULL,
                         errpred=NULL,
                         run_model=NULL,
                         ind_eta=NULL,
                         sigma=NULL,
                         solve_omega=NULL,
                         omega_dim=NULL,
                         iov_col=NULL,
                         pimat=NULL,
                         eta_map=NULL,
                         error_model=NULL,
                         estim_with_iov=NULL,
                         interpolation=NULL,
                         adapt=NULL){

  segment_id     <- unique(dat$AMS)
  n_segment      <- length(segment_id)
  dat_segment    <- dat

  eta_mat        <- matrix(0,nrow=n_segment+1,ncol=ncol(omega))
  eta_df         <- data.frame(eta_mat)
  names(eta_df)  <- attr(omega,"dimnames")[[1]]

  init_mat       <- matrix(0,nrow=n_segment+1,
                           ncol=length(solved_model$inits))
  init_df        <- data.frame(init_mat)
  init_names     <- names(solved_model$inits)
  names(init_df) <- init_names

  if(return_AMS_models){
    AMS_models <- list()
  }

  # For each segment, split the data set, run the estimation, and log the
  # the ODE compartment states at the end, to be used as init for the
  # next segment
  for(index_segment in 1:n_segment){

    dat_segment  <- dat[which(dat$AMS == segment_id[index_segment]),]

    if(index_segment>1){ # set TIME == 0 when a segment starts
      dat_segment$TIME  <- dat_segment$TIME -
        utils::tail(dat[dat$AMS == segment_id[index_segment-1],]$TIME,1)
    }

    # solved_model for the current segment
    solved_model <- rxode2::rxSolve(solved_model,c(theta,start_eta),
                                   dat_segment,
                                   covsInterpolation=interpolation)

    y_obs        <- dat_segment$DV[dat_segment$EVID == 0]

    model_init   <- init_df[index_segment,]

    r <- stats::optim(start_eta,errpred,
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
                      adapt=adapt,
                      index_segment=index_segment,
                      method="L-BFGS-B")

    eta_map[ind_eta]              <- r$par
    eta_df[index_segment+1,]      <- eta_map + eta_df[index_segment,]

    # get the ODE compartment states at the end of the segment
    solved_model$params <- c(theta,unlist(eta_df[index_segment+1,]))
    solved_model$inits  <- init_df[index_segment,]
    init_df[index_segment+1,]       <- utils::tail(solved_model[,init_names],1)

    if(return_AMS_models){
      # get the rxode2 model for the current segment
      AMS_models[[index_segment]]   <- solved_model
    }
  }
  adaptive_output <- list(eta_df=eta_df)

  if(return_AMS_models){
    adaptive_output$AMS_models <- AMS_models
  }
  return(adaptive_output)
}
