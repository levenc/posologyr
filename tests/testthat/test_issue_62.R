mod_meropenem_Guohua1   <- function() {
  ini({
    #Fixed effects: POP
    THETA_1 = 1.59
    THETA_2 = 3.69
    THETA_3 = 2.39
    THETA_4 = 35.1

    #Random effects: IIV
    ETA_Cl  ~ 0.206
    ETA_V   ~ 0.298

    # residual variability
    add_sd  <- 0.209
    prop_sd <- 0.316
  })
  model({
    #Individual and covariates
    TVCl  = THETA_1 +THETA_2 *(CLCR / 87)*(1-CRRT)+THETA_3*CRRT
    TVV   = THETA_4*(BW/85)
    Cl    = TVCl * exp(ETA_Cl)
    V     = TVV * exp(ETA_V)
    ke    = Cl / V
    Cp    = centr / V

    d / dt(centr)  = -ke * centr
    Cp ~ add(add_sd) + prop(prop_sd)
  })
}

df_patient <- data.frame(
  ID   = 1,
  TIME = c(0, 8, 16, 16, 24),
  AMT  = c(2000, 2000, 2000,0, 0),
  DUR  = c(8, 8, 8, NA, NA),
  DV   = c(NA, NA, NA, 12, NA),
  EVID = c(1, 1, 1, 0, 0),
  CLCR = c(30, 30, 30, 30,30),
  CRRT = c(1,1,1,1,1),
  BW   = c(60, 60, 60,60, 60),
  QFD  = c(2, 2, 2, 2, 2)
)

test_that("objective_function() throws an error when using a data.frame with
          DV==NA and EVID==0 on the same row", {
  expect_error(posologyr:::objective_function(y_obs=NA,f=4,g=1,eta=1,
                                              solve_omega=1))
  expect_error(poso_estim_map(df_patient,mod_meropenem_Guohua1))
          })
