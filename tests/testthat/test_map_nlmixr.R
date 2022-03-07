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

df_patient03_amik <- data.frame(ID=1,TIME=c(0,1,6),
                                DV=c(NA,75,32),
                                EVID=c(1,0,0),
                                AMT=c(1600,0,0),
                                DUR=c(0.5,NA,NA),
                                CLCREAT4H=50,TBW=62,PoverF=169,AMS=1)

patient03_amik  <- posologyr(mod_amikacin_2cpt_Burdet2015,
                             dat=df_patient03_amik)

patient03_amik_map  <- poso_estim_map(patient03_amik,
                                      return_model=TRUE)

patient03_amik_map_ad <- poso_estim_map(patient03_amik,
                                        adapt=TRUE,
                                        return_model=TRUE)

test_that("MAP estimates match nlmixr posthoc estimates", {
  expect_equal(patient03_amik_map$model$Cl[1], 2.62, tolerance=1e-2)
  expect_equal(patient03_amik_map$model$Vc[1], 11.01, tolerance=1e-2)
  expect_equal(patient03_amik_map$model$Vp[1], 14.57, tolerance=1e-2)
  expect_equal(patient03_amik_map$model$Q[1], 11.96, tolerance=1e-2)
})

test_that("adaptive MAP estimates match standard MAP estimates for
          a single segment", {
  expect_equal(patient03_amik_map$model$Cl[1],
               patient03_amik_map_ad$model$Cl[1],
               tolerance=1e-3)
  expect_equal(patient03_amik_map$model$Vc[1],
               patient03_amik_map_ad$model$Vc[1],
               tolerance=1e-3)
  expect_equal(patient03_amik_map$model$Vp[1],
               patient03_amik_map_ad$model$Vp[1],
               tolerance=1e-3)
  expect_equal(patient03_amik_map$model$Q[1],
               patient03_amik_map_ad$model$Q[1],
               tolerance=1e-3)
})
