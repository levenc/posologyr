df_patient01 <- data.frame(ID=1,TIME=c(0.0,0.5,1.0,14.0),
                        DV=c(NA,NA,25.0,5.5),
                        AMT=c(1000,-1000,0,0),
                        EVID=c(10102,10102,0,0),
                        CLCREAT=80,WT=65)

patient01_tobra <- posologyr(prior_model=mod_tobramycin_2cpt_fictional,
                          dat=df_patient01)

patient01_tobra_map <- poso_estim_map(patient01_tobra,
                                   return_model = TRUE)

test_that("MAP estimates match Monolix estimates", {
  expect_equal(patient01_tobra_map[[2]]$ke[1], 0.1258, tolerance = 1e-3)
  expect_equal(patient01_tobra_map[[2]]$V[1], 18.21, tolerance = 1e-2)
}) # 3 significant digits
