#-------------------------------------------------------------------------
# posologyr: individual dose optimisation using population PK
# Copyright (C) 2021  Cyril Leven
#
#    This program is free software: you can redistribute it and/or modify
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

mod_tobramycin_2cpt_fictional <- list(
  ppk_model   = RxODE::RxODE({
    centr(0) = 0;
    TVke  = THETA_ke*(CLCREAT/67.8)^0.89*(WT/66.4)^-1.09;
    TVV   = THETA_V*(WT/66.4)^0.80;
    TVk12 = THETA_k12;
    TVk21 = THETA_k21;
    ke    = TVke*exp(eta_ke);
    V     = TVV*exp(eta_V);
    k12   = TVk12;
    k21   = TVk21;
    Cc    = centr/V;
    d/dt(centr)  = - ke*centr - k12*centr + k21*periph;
    d/dt(periph) =            + k12*centr - k21*periph;
    d/dt(AUC)    =   Cc;
  }),
  error_model = function(f,xi){
    g <- xi[1] + xi[2]*f
    return(g)
  },
  pk_prior    = list(psi = c(THETA_ke=0.21, THETA_V=19.8,
                                    THETA_k12=0.041, THETA_k21=0.12),
                     Omega = lotri::lotri({eta_ke + eta_V + eta_k12 + eta_k21 ~
                          c(0.08075,
                            0      , 0.01203,
                            0      , 0      ,  0,
                            0      , 0      ,  0, 0)})),
  covariates  = c("CLCREAT","WT"),
  xi          = c(additive_a = 0, proportional_b = 0.198))

mod_amoxicillin_oral_1cpt_fictional <- list(
  ppk_model   = RxODE::RxODE({
    depot(0) = 0;
    centr(0) = 0;
    TVka = THETA_ka;
    ka   = TVka*exp(eta_ka);
    TVV  = THETA_V;
    V    = TVV*exp(eta_V);
    TVCl = THETA_Cl*(CLCREAT/90)^0.62;
    Cl   = TVCl*exp(eta_Cl);
    ke   = Cl/V;
    Cc   = centr/V;
    d/dt(depot) = -ka*depot;
    d/dt(centr) =  ka*depot - ke*centr;
    d/dt(AUC)   =  Cc;
  }),
  error_model = function(f,xi){
    g <- xi[1] + xi[2]*f
    return(g)
  },
  pk_prior    = list(psi = c(THETA_ka=0.52, THETA_V=0.56,
                             THETA_Cl=7.33),
                      Omega = lotri::lotri({eta_ka + eta_V + eta_Cl ~
                          c(0.0289,
                           0     , 0.0256 ,
                           0     ,-0.01792,  0.0400)})),
  covariates  = c("CLCREAT"),
  xi          = c(additive_a = 0.24, proportional_b = 0.27))

## save the models in a .rda data file
save(mod_tobramycin_2cpt_fictional,
     mod_amoxicillin_oral_1cpt_fictional,
     file="data/lib_ppk_model.rda")
