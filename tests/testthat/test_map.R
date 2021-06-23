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

df_patient03_amik <- data.frame(ID=1,TIME=c(0,1,6),
                                DV=c(NA,75,32),
                                EVID=c(1,0,0),
                                AMT=c(1600,0,0),
                                DUR=c(0.5,NA,NA),
                                CLCREAT4H=50,TBW=62,PoverF=169,AMS=1)

patient01_tobra <- posologyr(prior_model=mod_tobramycin_2cpt_fictional,
                             dat=df_patient01_tobra)

patient02_vanco <- posologyr(mod_vancomycin_2cpt_Goti2018,
                             dat=df_patient02_vanco)

patient03_amik  <- posologyr(mod_amikacin_2cpt_Burdet2015,
                             dat=df_patient03_amik)

patient01_tobra_map <- poso_estim_map(patient01_tobra,
                                      return_model=TRUE)

patient02_vanco_map <- poso_estim_map(patient02_vanco,
                                      return_model=TRUE)

patient03_amik_map  <- poso_estim_map(patient03_amik,
                                      return_model=TRUE)

patient03_amik_map_ad <- poso_estim_map(patient03_amik,
                                        adapt=TRUE,
                                        return_model=TRUE)

test_that("MAP estimates match Monolix MAP estimates", {
  expect_equal(patient01_tobra_map$model$ke[1], 0.1258, tolerance=1e-3)
  expect_equal(patient01_tobra_map$model$V[1], 18.21, tolerance=1e-2)
  expect_equal(patient02_vanco_map$model$Cl[1], 1.72, tolerance=1e-2)
  expect_equal(patient02_vanco_map$model$Vc[1], 59.8, tolerance=1e-2)
  expect_equal(patient02_vanco_map$model$Vp[1], 41.6, tolerance=1e-2)
})

test_that("MAP estimates match nlmixr posthoc estimates", {
  expect_equal(patient03_amik_map$model$Cl[1], 2.62, tolerance=1e-2)
  expect_equal(patient03_amik_map$model$Vc[1], 11.01, tolerance=1e-2)
  expect_equal(patient03_amik_map$model$Vp[1], 14.57, tolerance=1e-2)
  expect_equal(patient03_amik_map$model$Q[1], 11.96, tolerance=1e-2)
})

test_that("adaptive MAP estimates match standard MAP estimates for
          a single segment", {
  expect_equal(patient03_amik_map$model$Cl[1],
               patient03_amik_map_ad$model$Cl[1])
  expect_equal(patient03_amik_map$model$Vc[1],
               patient03_amik_map_ad$model$Vc[1])
  expect_equal(patient03_amik_map$model$Vp[1],
               patient03_amik_map_ad$model$Vp[1])
  expect_equal(patient03_amik_map$model$Q[1],
               patient03_amik_map_ad$model$Q[1])
})
