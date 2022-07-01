mod_piperacillin_2cpt_Roberts2010 <- list(
  ppk_model   = rxode2::rxode({
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

df_patient32 <- data.frame(ID=6,
                         TIME=c(0.0,8,16,24,32,40,47.9,48,56,64,71.9),
                         DV=c(NA,NA,NA,NA,NA,NA,63.0,NA,NA,NA,19.7),
                         AMT=c(4000,4000,4000,4000,4000,4000,NA,4000,4000,4000,
                               NA),
                         DUR=c(8,8,8,8,8,8,NA,8,8,8,NA),
                         EVID=c(1,1,1,1,1,1,0,1,1,1,0),
                         TBW=92,OCC=c(1,1,1,1,1,1,1,2,2,2,2))

test_that("Dosing optim functions can use models with IOV", {
  expect_equal(poso_time_cmin(dat=df_patient32,
                              prior_model=mod_piperacillin_2cpt_Roberts2010,
                              from=4,
                              dose=1500,
                              duration=4,
                              target_cmin=10)$time,4.8)
  expect_equal(poso_dose_auc(dat=df_patient32,
                             prior_model=mod_piperacillin_2cpt_Roberts2010,
                             time_auc=24,
                             duration=4,
                             target_auc=200)$dose,3793,tolerance=1e-1)
  expect_equal(poso_dose_conc(dat=df_patient32,
                              prior_model=mod_piperacillin_2cpt_Roberts2010,
                              time_c=6,
                              duration=4,
                              target_conc=10)$dose,3388,tolerance=1e-1)
  expect_equal(poso_inter_cmin(dat=df_patient32,
                               prior_model=mod_piperacillin_2cpt_Roberts2010,
                               dose=5000,
                               add_dose=10,
                               duration=4,
                               target_cmin=20)$interval,5.66,tolerance=1e-1)
  })
