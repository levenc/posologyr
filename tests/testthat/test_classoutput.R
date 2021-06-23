df_patient01_tobra <- data.frame(ID=1,TIME=c(0.0,1.0,14.0),
                                 DV=c(NA,25.0,5.5),
                                 AMT=c(2000,0,0),
                                 DUR=c(0.5,NA,NA),
                                 EVID=c(1,0,0),
                                 CLCREAT=80,WT=65)

patient01_tobra      <- posologyr(prior_model=mod_tobramycin_2cpt_fictional,
                          dat=df_patient01_tobra)

test_that("poso_simu_pop returns the expected objects", {
  p01_pop_mod    <- poso_simu_pop(patient01_tobra)
  p01_pop_nomod  <- poso_simu_pop(patient01_tobra,return_model = FALSE)

  expect_equal("data.frame" %in% attributes(p01_pop_mod$eta)$class,
               TRUE)
  expect_equal("rxSolve" %in% attributes(p01_pop_mod$model)$class,
               TRUE)
  expect_equal("data.frame" %in% attributes(p01_pop_nomod$eta)$class,
               TRUE)
  expect_equal("rxSolve" %in% attributes(p01_pop_nomod$model)$class,
               FALSE)
})

test_that("poso_estim_map returns the expected objects", {
  p01_map_mod    <- poso_estim_map(patient01_tobra)
  p01_map_nomod  <- poso_estim_map(patient01_tobra,return_model = FALSE)

  expect_equal("data.frame" %in% attributes(p01_map_mod$eta)$class,
               FALSE)
  expect_equal(attributes(p01_map_mod$eta)$class,
               NULL)
  expect_equal("rxSolve" %in% attributes(p01_map_mod$model)$class,
               TRUE)
  expect_equal("data.frame" %in% attributes(p01_map_nomod$eta)$class,
               FALSE)
  expect_equal(attributes(p01_map_nomod$eta)$class,
               NULL)
  expect_equal("rxSolve" %in% attributes(p01_map_nomod$model)$class,
               FALSE)
})

test_that("poso_estim_mcmc returns the expected objects", {
  p01_mcmc_mod   <- poso_estim_mcmc(patient01_tobra,burn_in=0,n_iter = 50)
  p01_mcmc_nomod <- poso_estim_mcmc(patient01_tobra,burn_in=0,n_iter = 50,
                                    return_model = FALSE)

  expect_equal("data.frame" %in% attributes(p01_mcmc_mod$eta)$class,
               TRUE)
  expect_equal("rxSolve" %in% attributes(p01_mcmc_mod$model)$class,
               TRUE)
  expect_equal("data.frame" %in% attributes(p01_mcmc_nomod$eta)$class,
               TRUE)
  expect_equal("rxSolve" %in% attributes(p01_mcmc_nomod$model)$class,
               FALSE)
})

test_that("poso_estim_sir returns the expected objects", {
  p01_sir_mod   <- poso_estim_sir(patient01_tobra)
  p01_sir_nomod <- poso_estim_sir(patient01_tobra,return_model = FALSE)

  expect_equal("data.frame" %in% attributes(p01_sir_mod$eta)$class,
               TRUE)
  expect_equal("rxSolve" %in% attributes(p01_sir_mod$model)$class,
               TRUE)
  expect_equal("data.frame" %in% attributes(p01_sir_nomod$eta)$class,
               TRUE)
  expect_equal("rxSolve" %in% attributes(p01_sir_nomod$model)$class,
               FALSE)
})
