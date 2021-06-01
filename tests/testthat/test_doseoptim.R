df_patient01 <- data.frame(ID=1,TIME=c(0.0,0.5,1.0,14.0),
                        DV=c(NA,NA,25.0,5.5),
                        AMT=c(1000,-1000,0,0),
                        EVID=c(10102,10102,0,0),
                        CLCREAT=80,WT=65)

df_patient04_vanco <- data.frame(ID=1,TIME=c(0.0,2,13.0,24.2,26.2,48),
                                 DV=c(NA,NA,12,NA,NA,10.5),
                                 AMT=c(1000,-1000,0,700,-700,0),
                                 EVID=c(10101,10101,0,10101,10101,0),
                                 CLCREAT=65,WT=70,DIAL=0)

patient01_tobra <- posologyr(prior_model=mod_tobramycin_2cpt_fictional,
                             dat=df_patient01)

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
  expect_equal(poso_dose_ctime(patient01_tobra,
                               param_map=params_patient01_tobra_map,
                               time_c=1,
                               duration=0.5,
                               target_conc=30),
               poso_dose_ctime(patient01_tobra,
                               param_map=NULL,
                               time_c=1,
                               duration=0.5,
                               target_conc=30),
               tolerance=1e-3)
  expect_equal(poso_dose_auc(patient01_tobra,
                               param_map=params_patient01_tobra_map,
                               time_auc=12,
                               duration=0.5,
                               target_auc=200),
               poso_dose_auc(patient01_tobra,
                               param_map=NULL,
                               time_auc=12,
                               duration=0.5,
                               target_auc=200),
               tolerance=1e-3)
  expect_equal(poso_time_cmin(patient01_tobra,
                             param_map=params_patient01_tobra_map,
                             dose=620,
                             duration=0.5,
                             target_cmin=0.5),
               poso_time_cmin(patient01_tobra,
                             param_map=NULL,
                             dose=620,
                             duration=0.5,
                             target_cmin=0.5),
               tolerance=1e-3)
})

test_that("Optimization results do not deviate from known values
          for single dose administration", {
  expect_equal(poso_time_cmin(patient04_tdm_vanco,
                              param_map=params_patient04_vanco_map,
                              from=2,
                              dose=1500,
                              duration=2,
                              target_cmin=10),
               11.2)
  expect_equal(poso_dose_auc(patient04_tdm_vanco,
                             param_map=params_patient04_vanco_map,
                             time_auc=24,
                             duration=2,
                             target_auc=400),
               2377.758,
               tolerance=1e-3)
  expect_equal(poso_dose_ctime(patient04_tdm_vanco,
                               param_map=params_patient04_vanco_map,
                               time_c=24,
                               duration=2,
                               target_conc=11.04931),
               2530.699,
               tolerance=1e-3)
})

test_that("Optimization results do not deviate from known values
          for multiple dose regimen", {
  expect_equal(poso_time_cmin(patient04_tdm_vanco,
                              param_map=params_patient04_vanco_map,
                              from=2,
                              dose=1500,
                              interdose_interval=24,
                              add_dose=10,
                              duration=2,
                              target_cmin=10),
               34.8)
  expect_equal(poso_dose_auc(patient04_tdm_vanco,
                             param_map=params_patient04_vanco_map,
                             time_auc=24,
                             starting_time=24*9,
                             interdose_interval=24,
                             add_dose=10,
                             duration=2,
                             target_auc=400),
               1229.366,
               tolerance=1e-3)
  expect_equal(poso_dose_ctime(patient04_tdm_vanco,
                               param_map=params_patient04_vanco_map,
                               time_c=24*9.9,
                               interdose_interval=24,
                               add_dose=10,
                               duration=24,
                               target_conc=400/24),
               1229.426,
               tolerance=1e-3)
})
