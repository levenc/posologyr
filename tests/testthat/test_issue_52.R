test_that("single IIV", {
  patient_data <- data.frame(ID=1,
                             TIME=c(0.0,1.0,14.0),
                             DV=c(NA,25.0,5.5),
                             AMT=c(2000,0,0),
                             DUR=c(0.5,NA,NA),
                             EVID=c(1,0,0),
                             CLCREAT=80,
                             WT=65)

  mod_run001 <- function() {
    ini({
      THETA_Cl <- 4.0
      THETA_Vc <- 70.0
      THETA_Ka <- 1.0
      ETA_Cl ~ 0.2
      prop.sd <- sqrt(0.05)
    })
    model({
      TVCl <- THETA_Cl
      TVVc <- THETA_Vc
      TVKa <- THETA_Ka

      Cl <- TVCl*exp(ETA_Cl)
      Vc <- TVVc
      Ka <- TVKa

      K20 <- Cl/Vc
      Cc <- centr/Vc

      d/dt(depot) = -Ka*depot
      d/dt(centr) = Ka*depot - K20*centr
      Cc ~ prop(prop.sd)
    })
  }

  f <- poso_estim_map(patient_data, mod_run001)

  expect_equal(f$eta, c(ETA_Cl = 0.7299818),tolerance=1e-3)

})
