test_that("poso_replace_et() respects nocb and midpoint covariate interpolation", {
  rxode2::setRxThreads(1L)

  mod_cov <- list(
    ppk_model = rxode2::rxode({
      centr(0) = 0
      depot(0) = 0
      Cl = THETA_Cl * (WT / 70) * exp(ETA_Cl)
      Vc = THETA_Vc * exp(ETA_Vc)
      Ka = THETA_Ka * exp(ETA_Ka)
      K20 = Cl / Vc
      Cc = centr / Vc
      d/dt(depot) = -Ka * depot
      d/dt(centr) = Ka * depot - K20 * centr
    }),
    error_model = function(f, sigma) {
      dv <- cbind(f, 1)
      sqrt(diag(dv %*% sigma %*% t(dv)))
    },
    theta = c(THETA_Cl = 4.0, THETA_Vc = 70.0, THETA_Ka = 1.0),
    omega = lotri::lotri({ETA_Cl + ETA_Vc + ETA_Ka ~ c(0.2, 0, 0.2, 0, 0, 0.2)}),
    sigma = lotri::lotri({prop + add ~ c(0.05, 0.0, 0.00)}),
    covariates = c("WT")
  )

  dat_cov <- data.frame(
    ID = 1,
    TIME = c(0, 1, 12),
    DV = c(NA, 20, 5),
    AMT = c(2000, 0, 0),
    EVID = c(1, 0, 0),
    WT = c(70, 70, 140)
  )

  pop <- poso_simu_pop(
    dat = dat_cov,
    prior_model = mod_cov,
    n_simul = 0,
    return_model = TRUE
  )

  new_et <- rxode2::et(amt = 2000)
  new_et$add.sampling(c(0, 1, 6, 12))

  res_locf <- poso_replace_et(pop$model, mod_cov, new_et, interpolation = "locf")
  res_mid  <- poso_replace_et(pop$model, mod_cov, new_et, interpolation = "midpoint")
  res_nocb <- poso_replace_et(pop$model, mod_cov, new_et, interpolation = "nocb")

  expect_equal(res_locf$WT, c(70, 70, 70, 140))
  expect_equal(res_mid$WT,  c(70, 70, 105, 140))
  expect_equal(res_nocb$WT, c(70, 70, 140, 140))

  # Higher WT increases Cl, so Cc at t = 6 and t = 12 must strictly decrease:
  # locf (WT=70) > midpoint (WT=105) > nocb (WT=140)
  expect_gt(res_locf$Cc[3], res_mid$Cc[3])
  expect_gt(res_mid$Cc[3], res_nocb$Cc[3])
  expect_gt(res_locf$Cc[4], res_mid$Cc[4])
  expect_gt(res_mid$Cc[4], res_nocb$Cc[4])

  # Direct rxSolve with the interpolated covariates and matching covsInterpolation
  et_nocb <- cbind(new_et, WT = c(70, 70, 70, 140, 140))
  expected_nocb <- rxode2::rxSolve(
    pop$model,
    et_nocb,
    pop$model$params,
    covsInterpolation = "nocb"
  )
  expect_equal(res_nocb$Cc, expected_nocb$Cc)
})
