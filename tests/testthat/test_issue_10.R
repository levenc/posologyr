mod_tobramycin_2cpt_fictional <- list(
  ppk_model   = RxODE::RxODE({
    centr(0) = 0;
    tTVke  = log(THETA_ke)+log(CLCREAT/67.8)*0.89+log(WT/66.4)*(-1.09);
    tTVV   = log(THETA_V)+log(WT/66.4)*0.80;
    tTVk12 = log(THETA_k12);
    tTVk21 = log(THETA_k21);
    ke     = exp(tTVke+ETA_ke);
    V      = exp(tTVV+ETA_V);
    k12    = exp(tTVk12);
    k21    = exp(tTVk21);
    Cc     = centr/V;
    d/dt(centr)  = - ke*centr - k12*centr + k21*periph;
    d/dt(periph) =            + k12*centr - k21*periph;
    d/dt(AUC)    =   Cc;
  }),
  error_model = function(f,sigma){
    g <- sigma[1] + sigma[2]*f
    return(g)
  },
  theta = c(THETA_ke=0.21, THETA_V=19.8,THETA_k12=0.041, THETA_k21=0.12),
  omega = lotri::lotri({ETA_ke + ETA_V + ETA_k12 + ETA_k21 ~
      c(0.08075,
        0      , 0.01203,
        0      , 0      ,  0,
        0      , 0      ,  0, 0)}),
  covariates  = c("CLCREAT","WT"),
  sigma       = c(additive_a = 0, proportional_b = 0.198))

df_patient01_tobra <- data.frame(ID=1,TIME=c(0.0,1.0,14.0),
                                 DV=c(NA,25.0,5.5),
                                 AMT=c(2000,0,0),
                                 DUR=c(0.5,NA,NA),
                                 EVID=c(1,0,0),
                                 CLCREAT=80,WT=65)

p01_tobra <- posologyr(prior_model=mod_tobramycin_2cpt_fictional,
                          dat=df_patient01_tobra)

test_that("poso_mcmc_pop accepts burn_in == n_iter", {
  expect_equal(class(poso_estim_mcmc(p01_tobra,burn_in = 10,
                                     n_iter = 10,
                                     return_model = FALSE)$eta),"data.frame")
})

test_that("poso_mcmc_pop accepts burn_in < n_iter", {
  expect_equal(class(poso_estim_mcmc(p01_tobra,burn_in = 12,
                                     n_iter = 10,
                                     return_model = FALSE)$eta),"data.frame")
})
