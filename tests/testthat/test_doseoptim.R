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

df_patient01_tobra <- data.frame(ID=1,TIME=c(0.0,1.0,14.0),
                                 DV=c(NA,25.0,5.5),
                                 AMT=c(500,0,0),
                                 DUR=c(0.5,NA,NA),
                                 EVID=c(1,0,0),
                                 CLCREAT=80,WT=65)

df_patient04_vanco <- data.frame(ID=1,TIME=c(0.0,13.0,24.2,48),
                                 DV=c(NA,12,NA,10.5),
                                 AMT=c(2000,0,1400,0),
                                 DUR=c(2,NA,2,NA),
                                 EVID=c(1,0,1,0),
                                 CLCREAT=65,WT=70,DIAL=0)

patient01_tobra <- posologyr(prior_model=mod_tobramycin_2cpt_fictional,
                             dat=df_patient01_tobra)

patient04_tdm_vanco <- posologyr(prior_model=mod_vancomycin_2cpt_Goti2018,
                            dat=df_patient04_vanco)

params_patient01_tobra_map <- c(THETA_ke=0.21,
                             CLCREAT=80,
                             WT=65,
                             THETA_V=19.8,
                             THETA_k12=0.041,
                             THETA_k21=0.12,
                             ETA_ke=-0.6828811,
                             ETA_V=-0.0663349)

params_patient04_vanco_map <- c(THETA_Cl=4.5,
                                CLCREAT=65,
                                DIAL=0,
                                THETA_Vc=58.4,
                                WT=70,
                                THETA_Vp=38.4,
                                THETA_Q=6.5,
                                ETA_Cl=0.10796221,
                                ETA_Vc=0.04643863,
                                ETA_Vp=0.06859768)

test_that("Same optimal dose with or without providing MAP estimates", {
  expect_equal(poso_dose_conc(patient01_tobra,
                               time_c=1,
                               duration=0.5,
                               target_conc=30,
                               indiv_param=params_patient01_tobra_map)$dose,
               poso_dose_conc(patient01_tobra,
                               time_c=1,
                               duration=0.5,
                               target_conc=30)$dose,
               tolerance=1e-3)
  expect_equal(poso_dose_auc(patient01_tobra,
                               time_auc=12,
                               duration=0.5,
                               target_auc=200,
                               indiv_param=params_patient01_tobra_map)$dose,
               poso_dose_auc(patient01_tobra,
                               time_auc=12,
                               duration=0.5,
                               target_auc=200,
                               indiv_param=NULL)$dose,
               tolerance=1e-3)
  expect_equal(poso_time_cmin(patient01_tobra,
                             dose=620,
                             duration=0.5,
                             target_cmin=0.5,
                             indiv_param=params_patient01_tobra_map,)$time,
               poso_time_cmin(patient01_tobra,
                             dose=620,
                             duration=0.5,
                             target_cmin=0.5)$time,
               tolerance=1e-3)
})

test_that("Optimization results do not deviate from known values
          for single dose administration", {
  expect_equal(poso_time_cmin(patient04_tdm_vanco,
                              from=2,
                              dose=1500,
                              duration=2,
                              target_cmin=10,
                              indiv_param=params_patient04_vanco_map)$time,
               11.2)
  expect_equal(poso_dose_auc(patient04_tdm_vanco,
                             time_auc=24,
                             duration=2,
                             target_auc=400,
                             indiv_param=params_patient04_vanco_map)$dose,
               2377.758,
               tolerance=1e-3)
  expect_equal(poso_dose_conc(patient04_tdm_vanco,
                              time_c=24,
                              duration=2,
                              target_conc=11.04931,
                              indiv_param=params_patient04_vanco_map)$dose,
               2530.699,
               tolerance=1e-3)
})

test_that("Optimization results do not deviate from known values
          for multiple dose regimen", {
  expect_equal(poso_time_cmin(patient04_tdm_vanco,
                              from=2,
                              dose=1500,
                              interdose_interval=24,
                              add_dose=10,
                              duration=2,
                              target_cmin=10,
                              indiv_param=params_patient04_vanco_map)$time,
               34.8)
  expect_equal(poso_dose_auc(patient04_tdm_vanco,
                             time_auc=24,
                             starting_time=24*9,
                             interdose_interval=24,
                             add_dose=10,
                             duration=2,
                             target_auc=400,
                             indiv_param=params_patient04_vanco_map)$dose,
               1229.366,
               tolerance=1e-3)
  expect_equal(poso_dose_conc(patient04_tdm_vanco,
                              time_c=24*9.9,
                              interdose_interval=24,
                              add_dose=10,
                              duration=24,
                              target_conc=400/24,
                              indiv_param=params_patient04_vanco_map)$dose,
               1229.426,
               tolerance=1e-3)
  expect_equal(poso_inter_cmin(patient04_tdm_vanco,
                               indiv_param=params_patient04_vanco_map,
                               dose=2000,
                               add_dose=10,
                               duration=2,
                               target_cmin=10)$interval,
               35.9,
               tolerance=1e-3)
})
