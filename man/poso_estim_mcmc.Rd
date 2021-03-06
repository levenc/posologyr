% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estim_functions.R
\name{poso_estim_mcmc}
\alias{poso_estim_mcmc}
\title{Estimate the posterior distribution of individual parameters by MCMC}
\usage{
poso_estim_mcmc(
  object,
  return_model = TRUE,
  burn_in = 50,
  n_iter = 1000,
  control = list(n_kernel = c(2, 2, 2), stepsize_rw = 0.4, proba_mcmc = 0.3, nb_max =
    3)
)
}
\arguments{
\item{object}{A posologyr list, created by the \code{\link{posologyr}}
function.}

\item{return_model}{A boolean. Returns a RxODE model using the estimated
ETAs if set to \code{TRUE}.}

\item{burn_in}{Number of burn-in iterations for the Metropolis-Hastings
algorithm.}

\item{n_iter}{Total number of iterations (following the burn-in iterations)
for the Metropolis-Hastings algorithm.}

\item{control}{A list of parameters controlling the Metropolis-Hastings
algorithm.}
}
\value{
If \code{return_model} is set to \code{FALSE}, a list of one element: a
dataframe \verb{$eta} of ETAs from the posterior distribution, estimated by
Markov Chain Monte Carlo.
If \code{return_model} is set to \code{TRUE}, a list of the dataframe of the posterior
distribution of ETA, and a RxODE model using the estimated distributions of ETAs.
}
\description{
Estimates the posterior distribution of individual parameters by Markov
Chain Monte Carlo (using a Metropolis-Hastings algorithm)
}
\examples{
# df_patient01: event table for Patient01, following a 30 minutes intravenous
# infusion of tobramycin
df_patient01 <- data.frame(ID=1,
                        TIME=c(0.0,1.0,14.0),
                        DV=c(NA,25.0,5.5),
                        AMT=c(2000,0,0),
                        EVID=c(1,0,0),
                        DUR=c(0.5,NA,NA),
                        CLCREAT=80,WT=65)
# loading a tobramycin model and Patient01's event record
patient01_tobra <- posologyr(prior_model=mod_tobramycin_2cpt_fictional,
                                dat=df_patient01)
# estimate the posterior distribution of population parameters
\donttest{poso_estim_mcmc(patient01_tobra,n_iter=100)}

}
\references{
Comets  E, Lavenu A, Lavielle M. Parameter estimation in nonlinear
mixed effect models using saemix, an R implementation of the SAEM algorithm.
Journal of Statistical Software 80, 3 (2017), 1-41.
}
\author{
Emmanuelle Comets, Audrey Lavenu, Marc Lavielle, Cyril Leven
}
