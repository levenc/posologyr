mod_tobramycin_2cpt_fictional <- list(
  ppk_model   = rxode2::rxode({
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
  ppk_model   = rxode2::rxode({
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

df_patient01_tobra <- data.frame(ID=1,TIME=c(0.0,1.0,14.0),
                                 DV=c(NA,25.0,5.5),
                                 AMT=c(500,0,0),
                                 DUR=c(0.5,NA,NA),
                                 EVID=c(1,0,0),
                                 CLCREAT=80,WT=65)

df_patient02_vanco <- data.frame(ID=1,TIME=c(0.0,12.0,22.2,37.5),
                                 DV=c(NA,14.8,NA,22.5),
                                 AMT=c(1900,0,1750,0),
                                 DUR=c(1,NA,1,NA),
                                 EVID=c(1,0,1,0),
                                 CLCREAT=34,WT=62,DIAL=0)

patient01_tobra <- posologyr(prior_model=mod_tobramycin_2cpt_fictional,
                             dat=df_patient01_tobra)

patient02_vanco <- posologyr(mod_vancomycin_2cpt_Goti2018,
                             dat=df_patient02_vanco)

patient01_tobra_map <- poso_estim_map(patient01_tobra,
                                      return_model=TRUE)

patient02_vanco_map <- poso_estim_map(patient02_vanco,
                                      return_model=TRUE)

test_that("MAP estimates match Monolix MAP estimates", {
  expect_equal(patient01_tobra_map$model$ke[1], 0.1258, tolerance=1e-3)
  expect_equal(patient01_tobra_map$model$V[1], 18.21, tolerance=1e-2)
  expect_equal(patient02_vanco_map$model$Cl[1], 1.72, tolerance=1e-2)
  expect_equal(patient02_vanco_map$model$Vc[1], 59.8, tolerance=1e-2)
  expect_equal(patient02_vanco_map$model$Vp[1], 41.6, tolerance=1e-2)
})
