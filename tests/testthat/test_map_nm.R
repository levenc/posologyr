mod_run302 <- list(
  ppk_model   = rxode2::rxode({

    TVCl  = THETA_Cl*((BW/75)^1.5)*(0.75^SEX);
    TVVc  = THETA_Vc;
    TVKa  = THETA_Ka;

    Cl    = TVCl*exp(ETA_Cl);
    Vc    = TVVc*exp(ETA_Vc);
    Ka    = TVKa*exp(ETA_Ka);

    K20   = Cl/Vc;
    Cc    = centr/Vc;

    d/dt(depot)  = - Ka*depot;
    d/dt(centr)  =   Ka*depot - K20*centr;
    d/dt(AUC)    =   Cc;
  }),
  error_model = function(f,sigma){
    dv <- cbind(f,1)
    g  <- diag(dv%*%sigma%*%t(dv))
    return(sqrt(g))
  },
  covariates = c("BW","SEX"),
  theta = c(THETA_Cl=4.0, THETA_Vc=70.0, THETA_Ka=1.0),
  omega = lotri::lotri({ETA_Cl + ETA_Vc + ETA_Ka ~
      c(0.2   ,
        0     ,     0.2,
        0     ,       0,   0.2)}),
  sigma       = lotri::lotri({prop + add ~ c(0.05,0.0,0.00)}))

# IV to centr
pat_302_1 <- data.frame(ID=1,
                        TIME=c(0.00000,1.50000,4.40000,7.10000,24.60000,72.00000,
                               96.00000,120.00000,144.00000,168.00000,192.00000,
                               216.00000),
                        EVID=c(1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1),
                        AMT=c(10000,0,0,0,0,10000,10000,10000,10000,10000,10000,10000),
                        CMT=2,
                        II=0,
                        ADDL=0,
                        RATE=c(10000,0,0,0,0,10000,10000,10000,10000,10000,10000,10000),
                        #MDV=c(1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1),
                        DV=c(NA,107.26700,100.99500,54.05460,5.97081,NA,NA,NA,NA,NA,NA,NA),
                        BW=c(77,113,92,132,126,128,41,56,68,45,93,126),
                        SEX=c(0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1))
set.seed(0)
pat302_nm  <- posologyr(mod_run302,dat=pat_302_1,nocb=TRUE)
pat302_nm_map  <- poso_estim_map(pat302_nm,return_model=TRUE)

test_that("MAP estimates match NONMEM posthoc estimates (time varying covariates)", {
  expect_equal(pat302_nm_map$eta[1], c(ETA_Cl=0.398845969), tolerance=1e-3)
  expect_equal(pat302_nm_map$eta[2], c(ETA_Vc=0.0762201167), tolerance=1e-3)
  expect_equal(pat302_nm_map$eta[3], c(ETA_Ka=0), tolerance=1e-3)
})
