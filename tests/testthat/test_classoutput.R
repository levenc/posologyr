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

df_patient06_pipera <- data.frame(ID=6,TIME=c(0.0,8,16,24,32,40,47.9,48,
                                         56,64,71.9,72,80,88,95.9,
                                         96,104,112,120),
                                  DV=c(NA,NA,NA,NA,NA,NA,63.0,NA,NA,
                                       NA,19.7,NA,NA,NA,31,NA,NA,NA,16),
                                  AMT=c(4000,4000,4000,4000,4000,4000,NA,
                                        4000,4000,4000,NA,4000,4000,4000,
                                        NA,4000,4000,4000,NA),
                                  DUR=c(8,8,8,8,8,8,NA,8,8,8,NA,8,8,8,NA,
                                        8,8,8,NA),
                                  EVID=c(1,1,1,1,1,1,0,1,1,1,0,1,1,1,0,
                                         1,1,1,0),TBW=92,
                                  OCC=c(1,1,1,1,1,1,1,2,2,2,2,3,3,3,3,
                                        4,4,4,4))

patient06_pipera <- posologyr(prior_model=mod_piperacillin_2cpt_Roberts2010,
                              dat=df_patient06_pipera)

test_that("poso_simu_pop returns the expected objects (with IOV)", {
  p06_pop_mod    <- poso_simu_pop(patient06_pipera)
  p06_pop_nomod  <- poso_simu_pop(patient06_pipera,return_model = FALSE)

  expect_equal("data.frame" %in% attributes(p06_pop_mod$eta)$class,
               TRUE)
  expect_equal("rxSolve" %in% attributes(p06_pop_mod$model)$class,
               TRUE)
  expect_equal("data.frame" %in% attributes(p06_pop_nomod$eta)$class,
               TRUE)
  expect_equal("rxSolve" %in% attributes(p06_pop_nomod$model)$class,
               FALSE)
})

test_that("poso_estim_map returns the expected objects (with IOV)", {
  p06_map_mod    <- poso_estim_map(patient06_pipera)
  p06_map_nomod  <- poso_estim_map(patient06_pipera,return_model = FALSE)

  expect_equal("data.frame" %in% attributes(p06_map_mod$eta)$class,
               FALSE)
  expect_equal(attributes(p06_map_mod$eta)$class,
               NULL)
  expect_equal("rxSolve" %in% attributes(p06_map_mod$model)$class,
               TRUE)
  expect_equal("data.frame" %in% attributes(p06_map_nomod$eta)$class,
               FALSE)
  expect_equal(attributes(p06_map_nomod$eta)$class,
               NULL)
  expect_equal("rxSolve" %in% attributes(p06_map_nomod$model)$class,
               FALSE)
})

test_that("poso_estim_sir returns the expected objects (with IOV)", {
  p06_sir_mod   <- poso_estim_sir(patient06_pipera,
                                  n_sample=5e2,
                                  n_resample=1e2)
  p06_sir_nomod <- poso_estim_sir(patient06_pipera,
                                  n_sample=5e2,
                                  n_resample=1e2,
                                  return_model = FALSE)

  expect_equal("data.frame" %in% attributes(p06_sir_mod$eta)$class,
               TRUE)
  expect_equal("rxSolve" %in% attributes(p06_sir_mod$model)$class,
               TRUE)
  expect_equal("data.frame" %in% attributes(p06_sir_nomod$eta)$class,
               TRUE)
  expect_equal("rxSolve" %in% attributes(p06_sir_nomod$model)$class,
               FALSE)
})
