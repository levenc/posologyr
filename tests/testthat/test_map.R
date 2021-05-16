df_patient01_tobra <- data.frame(ID=1,TIME=c(0.0,0.5,1.0,14.0),
                                 DV=c(NA,NA,25.0,5.5),
                                 AMT=c(1000,-1000,0,0),
                                 EVID=c(10101,10101,0,0),
                                 CLCREAT=80,WT=65)

df_patient02_vanco <- data.frame(ID=1,TIME=c(0.0,1,12.0,22.2,23.2,37.5),
                                 DV=c(NA,NA,14.8,NA,NA,22.5),
                                 AMT=c(1900,-1900,0,1750,-1750,0),
                                 EVID=c(10101,10101,0,10101,10101,0),
                                 CLCREAT=34,WT=62,DIAL=0)

patient01_tobra <- posologyr(prior_model=mod_tobramycin_2cpt_fictional,
                             dat=df_patient01_tobra)

patient02_vanco <- posologyr(mod_vancomycin_2cpt_Goti2018,
                             dat=df_patient02_vanco)

patient01_tobra_map <- poso_estim_map(patient01_tobra,
                                      return_model = TRUE)

patient02_vanco_map <- poso_estim_map(patient02_vanco,
                                      return_model = TRUE)


test_that("MAP estimates match Monolix estimates", {
  expect_equal(patient01_tobra_map[[2]]$ke[1], 0.1258, tolerance = 1e-3)
  expect_equal(patient01_tobra_map[[2]]$V[1], 18.21, tolerance = 1e-2)
  expect_equal(patient02_vanco_map[[2]]$Cl[1], 1.72, tolerance = 1e-2)
  expect_equal(patient02_vanco_map[[2]]$Vc[1], 59.8, tolerance = 1e-2)
  expect_equal(patient02_vanco_map[[2]]$Vp[1], 41.6, tolerance = 1e-2)
}) # 3 significant digits
