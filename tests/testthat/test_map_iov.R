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

df_patient06_pipera <- data.frame(ID=6,
                                  TIME=c(0.0,8,16,24,32,40,47.9,48,
                                         56,64,71.9,72,80,88,95.9,
                                         96,104,112,120),
                                  DV=c(NA,NA,NA,NA,NA,NA,63.0,NA,NA,
                                       NA,19.7,NA,NA,NA,31,NA,NA,NA,16),
                                  AMT=c(4000,4000,4000,4000,4000,4000,NA,
                                       4000,4000,4000,NA,4000,4000,4000,
                                       NA,4000,4000,4000,NA),
                                  DUR=c(8,8,8,8,8,8,NA,8,8,8,NA,8,8,8,NA,
                                        8,8,8,NA),
                                  EVID=c(1,1,1,1,1,1,0,1,1,1,0,1,1,1,0,
                                        1,1,1,0),
                                  TBW=92,
                                  OCC=c(1,1,1,1,1,1,1,2,2,2,2,3,3,3,3,
                                       4,4,4,4))

patient06_pipera <- posologyr(prior_model=mod_piperacillin_2cpt_Roberts2010,
                             dat=df_patient06_pipera)

patient06_pipera_map_iov <- poso_estim_map(patient06_pipera,
                                      return_model=TRUE)

test_that("MAP estimates match Monolix MAP estimates", {
  expect_equal(patient06_pipera_map_iov$model$LAGTIME[1], 0.07, tolerance=1e-3)
  expect_equal(unique(patient06_pipera_map_iov$model$Cl)[1], 10.58, tolerance=1e-3)
  expect_equal(unique(patient06_pipera_map_iov$model$Cl)[2], 24.0, tolerance=1e-3)
  expect_equal(unique(patient06_pipera_map_iov$model$Cl)[3], 18.4, tolerance=1e-3)
  expect_equal(unique(patient06_pipera_map_iov$model$Cl)[4], 26.3, tolerance=1e-3)
  expect_equal(patient06_pipera_map_iov$model$V1[1], 7.2, tolerance=1e-3)
  expect_equal(patient06_pipera_map_iov$model$Q[1], 52, tolerance=1e-3)
  expect_equal(patient06_pipera_map_iov$model$V2[1], 17.8, tolerance=1e-3)
})
