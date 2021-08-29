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
