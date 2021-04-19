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

test_that("MAP estimates match Monolix estimates", {
  expect_equal(michel_tobra_map[[2]]$ke[1], 0.1258, tolerance = 1e-3)
  expect_equal(michel_tobra_map[[2]]$V[1], 18.21, tolerance = 1e-2)
}) # 3 significant digits
