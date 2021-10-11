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
    tTVke  = log(THETA_ke)+log(CLCREAT/67.8)*0.89+log(WT/66.4)*(-1.09);
    tTVV   = log(THETA_V)+log(WT/66.4)*0.80;
    tTVk12 = log(THETA_k12);
    tTVk21 = log(THETA_k21);
    ke     = exp(tTVke+ETA_ke);
    V      = exp(tTVV+ETA_V);
    k12    = exp(tTVk12);
    k21    = exp(tTVk21);
    Cc     = centr/V;
    d/dt(centr)  = - ke*centr - k12*centr + k21*periph;
    d/dt(periph) =            + k12*centr - k21*periph;
    d/dt(AUC)    =   Cc;
  }),
  error_model = function(f,sigma){
    g <- sigma[1] + sigma[2]*f
    return(g)
  },
  theta = c(THETA_ke=0.21, THETA_V=19.8,THETA_k12=0.041, THETA_k21=0.12),
  omega = lotri::lotri({ETA_ke + ETA_V + ETA_k12 + ETA_k21 ~
                          c(0.08075,
                            0      , 0.01203,
                            0      , 0      ,  0,
                            0      , 0      ,  0, 0)}),
  covariates  = c("CLCREAT","WT"),
  sigma       = c(additive_a = 0, proportional_b = 0.198))

mod_vancomycin_2cpt_Goti2018 <- list(
  ppk_model   = RxODE::RxODE({
    centr(0) = 0;
    TVCl  = THETA_Cl*(CLCREAT/120)^0.8*(0.7^DIAL);
    TVVc  = THETA_Vc*(WT/70)          *(0.5^DIAL);
    TVVp  = THETA_Vp;
    TVQ   = THETA_Q;
    Cl    = TVCl*exp(ETA_Cl);
    Vc    = TVVc*exp(ETA_Vc);
    Vp    = TVVp*exp(ETA_Vp);
    Q     = TVQ;
    ke    = Cl/Vc;
    k12   = Q/Vc;
    k21   = Q/Vp;
    Cc    = centr/Vc;
    d/dt(centr)  = - ke*centr - k12*centr + k21*periph;
    d/dt(periph) =            + k12*centr - k21*periph;
    d/dt(AUC)    =   Cc;
  }),
  error_model = function(f,sigma){
    g <- sigma[1] + sigma[2]*f
    return(g)
  },
  theta = c(THETA_Cl=4.5, THETA_Vc=58.4, THETA_Vp=38.4,THETA_Q=6.5),
  omega = lotri::lotri({ETA_Cl + ETA_Vc + ETA_Vp + ETA_Q ~
      c(0.147,
        0     ,   0.510,
        0     ,       0,   0.282,
        0     ,       0,       0,    0)}),
  covariates  = c("CLCREAT","WT","DIAL"),
  sigma       = c(additive_a = 3.4, proportional_b = 0.227))

mod_amikacin_2cpt_Burdet2015 <- list(
  ppk_model   = RxODE::RxODE({
    centr(0) = 0;
    TVCl  = THETA_Cl*(CLCREAT4H/82)^0.7;
    TVVc  = THETA_Vc*(TBW/78)^0.9*(PoverF/169)^0.4;
    TVVp  = THETA_Vp;
    TVQ   = THETA_Q;
    Cl    = TVCl*exp(ETA_Cl);
    Vc    = TVVc*exp(ETA_Vc);
    Vp    = TVVp*exp(ETA_Vp);
    Q     = TVQ *exp(ETA_Q);
    ke    = Cl/Vc;
    k12   = Q/Vc;
    k21   = Q/Vp;
    Cc    = centr/Vc;
    d/dt(centr)  = - ke*centr - k12*centr + k21*periph;
    d/dt(periph) =            + k12*centr - k21*periph;
    d/dt(AUC)    =   Cc;
  }),
  error_model = function(f,sigma){
    g <- sigma[1] + sigma[2]*f
    return(g)
  },
  theta = c(THETA_Cl=4.3, THETA_Vc=15.9, THETA_Vp=21.4,THETA_Q=12.1),
  omega = lotri::lotri({ETA_Cl + ETA_Vc + ETA_Vp + ETA_Q ~
      c(0.1,
        0.01     ,   0.05 ,
        0.01     ,   0.02 ,   0.2  ,
        -0.06    ,   0.004,   0.003,    0.08)}),
  covariates  = c("CLCREAT4H","TBW","PoverF"),
  sigma       = c(additive_a = 0.2, proportional_b = 0.1))

mod_piperacillin_2cpt_Roberts2010 <- list(
  ppk_model   = RxODE::RxODE({
    centr(0)   = 0;
    # Time lag from dose infuser to patient
    TVLAGTIME  = THETA_LAGTIME
    LAGTIME    = TVLAGTIME*exp(ETA_LAGTIME)
    # -------------------------------------
    TVCl  = THETA_Cl*(TBW/70);
    TVV1  = THETA_V1;
    TVV2  = THETA_V2;
    TVQ   = THETA_Q;
    Cl    = TVCl*exp(ETA_Cl + KAPPA_Cl);
    Cli   = TVCl*exp(ETA_Cl);
    V1    = TVV1*exp(ETA_V1 + KAPPA_V1);
    V1i   = TVV1*exp(ETA_V1);
    V2    = TVV2*exp(ETA_V2);
    Q     = TVQ *exp(ETA_Q);
    ke    = Cl/V1;
    k12   = Q/V1;
    k21   = Q/V2;
    Cc    = centr/V1;
    d/dt(centr)  = - ke*centr - k12*centr + k21*periph;
    d/dt(periph) =            + k12*centr - k21*periph;
    lag(centr)   =   LAGTIME;
    d/dt(AUC)    =   Cc;
  }),
  error_model = function(f,sigma){
    g <- sigma[1] + sigma[2]*f
    return(g)
  },
  theta = c(THETA_LAGTIME=0.07,THETA_Cl=17.1,THETA_V1=7.2,THETA_V2=17.8,
            THETA_Q=52.0),
  omega = lotri::lotri({ETA_LAGTIME + ETA_Cl + ETA_V1 + ETA_V2 + ETA_Q ~
      c(0.1747673,
        0.00     ,  0.08507985,
        0.00     ,  0.00      ,  0.0673745,
        0.00     ,  0.00      ,  0.00     ,  0.429067,
        0.00     ,  0.00      ,  0.00     ,  0.00    ,  0.2247455)}),
  pi_matrix = lotri::lotri({KAPPA_Cl + KAPPA_V1 ~
      c(0.1934626,
        0.00     ,  0.05783106)}),
  covariates  = c("TBW"),
  sigma       = c(additive_a = 3.2, proportional_b = 0.253))

## save the models in a .rda data file
save(mod_tobramycin_2cpt_fictional,
     mod_vancomycin_2cpt_Goti2018,
     mod_amikacin_2cpt_Burdet2015,
     mod_piperacillin_2cpt_Roberts2010,
     file="data/lib_ppk_model.rda")
