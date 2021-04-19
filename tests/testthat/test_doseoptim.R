df_michel <- data.frame(ID=1,TIME=c(0.0,0.5,1.0,14.0),
                        DV=c(NA,NA,25.0,5.5),
                        AMT=c(1000,-1000,0,0),
                        EVID=c(10102,10102,0,0),
                        CLCREAT=80,WT=65)

mod_tobra <- mod_tobramycin_2cpt_fictional

solved_model_tobra <- RxODE::rxSolve(mod_tobra$ppk_model,
                                   c(mod_tobra$psi,
                                     diag(mod_tobra$omega)*0),
                                   df_michel)

michel_tobra_map <- poso_estim_map(solved_model = solved_model_tobra,
                                   prior_model = mod_tobra,
                                   dat = df_michel,
                                   return_model = TRUE)

params_michel_tobra_map <- c(THETA_ke=0.21,
                             CLCREAT=80,
                             WT=65,
                             THETA_V=19.8,
                             THETA_k12=0.041,
                             THETA_k21=0.12,
                             ETA_ke=-0.6828811,
                             ETA_V=-0.0663349)

test_that("Same optimal dose with or without providing MAP estimates", {
  expect_equal(poso_dose_ctime(solved_model = NULL,
                               prior_model = mod_tobra,
                               dat = NULL,
                               param_map = params_michel_tobra_map,
                               time_c = 1,
                               duration = 0.5,
                               target_conc = 30),
               poso_dose_ctime(solved_model = solved_model_tobra,
                               prior_model = mod_tobra,
                               dat = df_michel,
                               param_map = NULL,
                               time_c = 1,
                               duration = 0.5,
                               target_conc = 30),
               tolerance = 1e-3)
  expect_equal(poso_dose_auc(solved_model = NULL,
                               prior_model = mod_tobra,
                               dat = NULL,
                               param_map = params_michel_tobra_map,
                               time_auc = 12,
                               duration = 0.5,
                               target_auc = 200),
               poso_dose_auc(solved_model = solved_model_tobra,
                               prior_model = mod_tobra,
                               dat = df_michel,
                               param_map = NULL,
                               time_auc = 12,
                               duration = 0.5,
                               target_auc = 200),
               tolerance = 1e-3)
  expect_equal(poso_time_cmin(solved_model = NULL,
                             prior_model = mod_tobra,
                             dat = NULL,
                             param_map = params_michel_tobra_map,
                             dose = 620,
                             duration = 0.5,
                             target_cmin = 0.5),
               poso_time_cmin(solved_model = solved_model_tobra,
                             prior_model = mod_tobra,
                             dat = df_michel,
                             param_map = NULL,
                             dose = 620,
                             duration = 0.5,
                             target_cmin = 0.5),
               tolerance = 1e-3)
})
