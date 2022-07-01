mod_tobramycin_2cpt_fictional <- list(
  ppk_model   = rxode2::rxode({
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

mod_piperacillin_2cpt_Roberts2010 <- list(
  ppk_model   = rxode2::rxode({
    centr(0)   = 0;
    # Time lag from dose infuser to patient
    TVLAGTIME  = THETA_LAGTIME
    LAGTIME    = TVLAGTIME*exp(ETA_LAGTIME)
    # -------------------------------------
    TVCl  = THETA_Cl*(TBW/70);
    TVV1  = THETA_V1;
    TVV2  = THETA_V2;
    TVQ   = THETA_Q;
    Cl    = TVCl*exp(ETA_Cl + KAPPA_Cl);
    Cli   = TVCl*exp(ETA_Cl);
    V1    = TVV1*exp(ETA_V1 + KAPPA_V1);
    V1i   = TVV1*exp(ETA_V1);
    V2    = TVV2*exp(ETA_V2);
    Q     = TVQ *exp(ETA_Q);
    ke    = Cl/V1;
    k12   = Q/V1;
    k21   = Q/V2;
    Cc    = centr/V1;
    d/dt(centr)  = - ke*centr - k12*centr + k21*periph;
    d/dt(periph) =            + k12*centr - k21*periph;
    lag(centr)   =   LAGTIME;
    d/dt(AUC)    =   Cc;
  }),
  error_model = function(f,sigma){
    g <- sigma[1] + sigma[2]*f
    return(g)
  },
  theta = c(THETA_LAGTIME=0.07,THETA_Cl=17.1,THETA_V1=7.2,THETA_V2=17.8,
            THETA_Q=52.0),
  omega = lotri::lotri({ETA_LAGTIME + ETA_Cl + ETA_V1 + ETA_V2 + ETA_Q ~
      c(0.1747673,
        0.00     ,  0.08507985,
        0.00     ,  0.00      ,  0.0673745,
        0.00     ,  0.00      ,  0.00     ,  0.429067,
        0.00     ,  0.00      ,  0.00     ,  0.00    ,  0.2247455)}),
  pi_matrix = lotri::lotri({KAPPA_Cl + KAPPA_V1 ~
      c(0.1934626,
        0.00     ,  0.05783106)}),
  covariates  = c("TBW"),
  sigma       = c(additive_a = 3.2, proportional_b = 0.253))

df_patient01_tobra <- data.frame(ID=1,TIME=c(0.0,1.0,14.0),
                                 DV=c(NA,25.0,5.5),
                                 AMT=c(2000,0,0),
                                 DUR=c(0.5,NA,NA),
                                 EVID=c(1,0,0),
                                 CLCREAT=80,WT=65)

test_that("poso_simu_pop returns the expected objects", {
  p01_pop_mod    <- poso_simu_pop(dat=df_patient01_tobra,
                                  prior_model=mod_tobramycin_2cpt_fictional)
  p01_pop_nomod  <- poso_simu_pop(dat=df_patient01_tobra,
                                  prior_model=mod_tobramycin_2cpt_fictional,
                                  return_model = FALSE)

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
  p01_map_mod    <- poso_estim_map(dat=df_patient01_tobra,
                                   prior_model=mod_tobramycin_2cpt_fictional)
  p01_map_nomod  <- poso_estim_map(dat=df_patient01_tobra,
                                   prior_model=mod_tobramycin_2cpt_fictional
                                   ,return_model = FALSE)

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
  p01_mcmc_mod   <- poso_estim_mcmc(dat=df_patient01_tobra,
                                    prior_model=mod_tobramycin_2cpt_fictional,
                                    burn_in=0,n_iter = 10)
  p01_mcmc_nomod <- poso_estim_mcmc(dat=df_patient01_tobra,
                                    prior_model=mod_tobramycin_2cpt_fictional,
                                    burn_in=0,n_iter = 10,
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
  p01_sir_mod   <- poso_estim_sir(dat=df_patient01_tobra,
                                  prior_model=mod_tobramycin_2cpt_fictional,
                                  n_sample=5e2,
                                  n_resample=1e2)
  p01_sir_nomod <- poso_estim_sir(dat=df_patient01_tobra,
                                  prior_model=mod_tobramycin_2cpt_fictional,
                                  n_sample=5e2,
                                  n_resample=1e2,
                                  return_model = FALSE)

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

test_that("poso_simu_pop returns the expected objects (with IOV)", {
  p06_pop_mod    <- poso_simu_pop(dat=df_patient06_pipera,
                                  prior_model=mod_piperacillin_2cpt_Roberts2010)
  p06_pop_nomod  <- poso_simu_pop(dat=df_patient06_pipera,
                                  prior_model=mod_piperacillin_2cpt_Roberts2010,
                                  return_model = FALSE)

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
  p06_map_mod    <- poso_estim_map(dat=df_patient06_pipera,
                                   prior_model=mod_piperacillin_2cpt_Roberts2010)
  p06_map_nomod  <- poso_estim_map(dat=df_patient06_pipera,
                                   prior_model=mod_piperacillin_2cpt_Roberts2010,
                                   return_model = FALSE)

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
  p06_sir_mod   <- poso_estim_sir(dat=df_patient06_pipera,
                                  prior_model=mod_piperacillin_2cpt_Roberts2010,
                                  n_sample=5e2,
                                  n_resample=1e2)
  p06_sir_nomod <- poso_estim_sir(dat=df_patient06_pipera,
                                  prior_model=mod_piperacillin_2cpt_Roberts2010,
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
