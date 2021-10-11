df_patient06_pipera <- data.frame(ID=6,
                                  TIME=c(0.0,8,16,24,32,40,47.9,48,
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
                                        1,1,1,0),
                                  TBW=92,
                                  OCC=c(1,1,1,1,1,1,1,2,2,2,2,3,3,3,3,
                                       4,4,4,4))

patient06_pipera <- posologyr(prior_model=mod_piperacillin_2cpt_Roberts2010,
                             dat=df_patient06_pipera)

patient06_pipera_map_iov <- poso_estim_map(patient06_pipera,
                                      return_model=TRUE)

test_that("MAP estimates match Monolix MAP estimates", {
  expect_equal(patient06_pipera_map_iov$model$LAGTIME[1], 0.07, tolerance=1e-3)
  expect_equal(patient06_pipera_map_iov$model$Cl[1], 10.58, tolerance=1e-3)
  expect_equal(patient06_pipera_map_iov$model$Cl[2], 24.0, tolerance=1e-3)
  expect_equal(patient06_pipera_map_iov$model$Cl[3], 18.4, tolerance=1e-3)
  expect_equal(patient06_pipera_map_iov$model$Cl[4], 26.3, tolerance=1e-3)
  expect_equal(patient06_pipera_map_iov$model$V1[1], 7.2, tolerance=1e-3)
  expect_equal(patient06_pipera_map_iov$model$V1[2], 7.2, tolerance=1e-3)
  expect_equal(patient06_pipera_map_iov$model$V1[3], 7.2, tolerance=1e-3)
  expect_equal(patient06_pipera_map_iov$model$V1[4], 7.2, tolerance=1e-3)
  expect_equal(patient06_pipera_map_iov$model$Q[1], 52, tolerance=1e-3)
  expect_equal(patient06_pipera_map_iov$model$V2[1], 17.8, tolerance=1e-3)
})
