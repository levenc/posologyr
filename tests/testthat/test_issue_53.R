test_that("poso_simu_pop() outputs a model with interpolated covariates", {
  patient_data <- data.frame(ID=1,
                             TIME=c(0.0,1.0,14.0),
                             DV=c(NA,25.0,5.5),
                             AMT=c(2000,0,0),
                             DUR=c(0.5,NA,NA),
                             EVID=c(1,0,0),
                             CLCREAT=c(60,60,20),
                             WT=65)

  mod_run001 <- function() {
    ini({
      THETA_Cl <- 4.0
      THETA_Vc <- 70.0
      THETA_Ka <- 1.0
      ETA_Cl ~ 0.2
      prop.sd <- sqrt(0.05)
    })
    model({
      TVCl <- THETA_Cl*(CLCREAT/90)
      TVVc <- THETA_Vc
      TVKa <- THETA_Ka

      Cl <- TVCl*exp(ETA_Cl)
      Vc <- TVVc
      Ka <- TVKa

      K20 <- Cl/Vc
      Cc <- centr/Vc

      d/dt(depot) = -Ka*depot
      d/dt(centr) = Ka*depot - K20*centr
      Cc ~ prop(prop.sd)
    })
  }

  f <- poso_simu_pop(patient_data, mod_run001, n_simul=1000,return_model=TRUE)

  #expect length n_simul*length(seq(0,15,by=.1))
  expect_equal(length(f$model$time), 1000*151)

  #unique values of time
  expect_equal(length(unique(f$model$time)), 151)

  #unique values of Cl: depends on unique values of CLCREAT and n_simul
  expect_equal(length(unique(f$model$Cl)), 1000*2)

})

test_that("poso_simu_pop() outputs a model with interpolated covariates (with IOV)", {
mod_piperacillin <- list(
  ppk_model   = rxode2::rxode({
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

df_patient <- data.frame(ID=6,
                         TIME=c(0.0,8,16,24,32,40,47.9,48,56,64,71.9),
                         DV=c(NA,NA,NA,NA,NA,NA,63.0,NA,NA,NA,19.7),
                         AMT=c(4000,4000,4000,4000,4000,4000,NA,
                               4000,4000,4000,NA),
                         DUR=c(8,8,8,8,8,8,NA,8,8,8,NA),
                         EVID=c(1,1,1,1,1,1,0,1,1,1,0),
                         TBW=c(92,92,92,60,60,60,60,60,60,60,60),
                         OCC=c(1,1,1,1,1,1,1,2,2,2,2))

f <- poso_simu_pop(df_patient,mod_piperacillin,n_simul=1000,return_model=TRUE)

#expect length n_simul*length(seq(0,15,by=.1))
expect_equal(length(f$model$time), 1000*730)

#unique values of time
expect_equal(length(unique(f$model$time)), 730)

#unique values of Cl: depends on unique values of TBW and n_simul
expect_equal(length(unique(f$model$Cl)), 1000*2)

#unique values of KAPPA_Cl and KAPPA_V1: zero
expect_equal(unique(f$model$KAPPA_Cl), 0)
expect_equal(unique(f$model$KAPPA_V1), 0)
})
