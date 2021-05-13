df_patient01 <- data.frame(ID=1,TIME=c(0.0,0.5,1.0,14.0),
                        DV=c(NA,NA,25.0,5.5),
                        AMT=c(1000,-1000,0,0),
                        EVID=c(10102,10102,0,0),
                        CLCREAT=80,WT=65)

p01_tobra <- posologyr(prior_model=mod_tobramycin_2cpt_fictional,
                          dat=df_patient01)

test_that("poso_mcmc_pop accepts burn_in == n_iter", {
  expect_equal(class(poso_estim_mcmc(p01_tobra,burn_in = 20,
                                     n_iter = 20,
                                     return_model = FALSE)[[1]]),"numeric")
})

test_that("poso_mcmc_pop accepts burn_in < n_iter", {
  expect_equal(class(poso_estim_mcmc(p01_tobra,burn_in = 22,
                                     n_iter = 20,
                                     return_model = FALSE)[[1]]),"numeric")
})
