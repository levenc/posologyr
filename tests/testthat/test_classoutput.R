df_patient01 <- data.frame(ID=1,TIME=c(0.0,0.5,1.0,14.0),
                        DV=c(NA,NA,25.0,5.5),
                        AMT=c(1000,-1000,0,0),
                        EVID=c(10102,10102,0,0),
                        CLCREAT=80,WT=65)

patient01_tobra      <- posologyr(prior_model=mod_tobramycin_2cpt_fictional,
                          dat=df_patient01)

test_that("poso_simu_pop returns the expected objects", {
  p01_pop_mod    <- poso_simu_pop(patient01_tobra)
  p01_pop_nomod  <- poso_simu_pop(patient01_tobra,return_model = FALSE)

  expect_equal("data.frame" %in% attributes(p01_pop_mod[[1]])$class,
               TRUE)
  expect_equal("rxSolve" %in% attributes(p01_pop_mod[[2]])$class,
               TRUE)
  expect_equal("data.frame" %in% attributes(p01_pop_nomod[[1]])$class,
               FALSE)
  expect_equal(attributes(p01_pop_nomod[[1]])$class,
               NULL)
  expect_equal("rxSolve" %in% attributes(p01_pop_nomod[[2]])$class,
               FALSE)
  expect_equal(attributes(p01_pop_nomod[[2]])$class,
               NULL)
  expect_equal(class(p01_pop_nomod[[1]]),"numeric")
})

test_that("poso_estim_map returns the expected objects", {
  p01_map_mod    <- poso_estim_map(patient01_tobra)
  p01_map_nomod  <- poso_estim_map(patient01_tobra,return_model = FALSE)

  expect_equal(attributes(p01_map_mod[[1]])$class,
               NULL)
  expect_equal(class(p01_map_nomod[[1]]),"numeric")
  expect_equal("rxSolve" %in% attributes(p01_map_mod[[2]])$class,
               TRUE)
  expect_equal(attributes(p01_map_nomod[[1]])$class,
               NULL)
  expect_equal("rxSolve" %in% attributes(p01_map_nomod[[2]])$class,
               FALSE)
  expect_equal(attributes(p01_map_nomod[[2]])$class,
               NULL)
  expect_equal(class(p01_map_nomod[[1]]),"numeric")
})

test_that("poso_estim_mcmc returns the expected objects", {
  p01_mcmc_mod   <- poso_estim_mcmc(patient01_tobra,burn_in=0,n_iter = 50)
  p01_mcmc_nomod <- poso_estim_mcmc(patient01_tobra,burn_in=0,n_iter = 50,
                                    return_model = FALSE)

  expect_equal("data.frame" %in% attributes(p01_mcmc_mod[[1]])$class,
               TRUE)
  expect_equal("rxSolve" %in% attributes(p01_mcmc_mod[[2]])$class,
               TRUE)
  expect_equal("data.frame" %in% attributes(p01_mcmc_nomod[[1]])$class,
               FALSE)
  expect_equal(attributes(p01_mcmc_nomod[[1]])$class,
               NULL)
  expect_equal("rxSolve" %in% attributes(p01_mcmc_nomod[[2]])$class,
               FALSE)
  expect_equal(attributes(p01_mcmc_nomod[[2]])$class,
               NULL)
  expect_equal(class(p01_mcmc_nomod[[1]]),"numeric")
})
