# prior model
mod_run003 <- list(
  ppk_model = RxODE::RxODE({
    centr(0) = 0;
    depot(0) = 0;

    TVCl   = THETA_Cl;
    TVVc   = THETA_Vc;
    TVKa   = THETA_Ka;
    TVALAG = THETA_ALAG;

    Cl   = TVCl*exp(ETA_Cl);
    Vc   = TVVc*exp(ETA_Vc);
    Ka   = TVKa*exp(ETA_Ka);
    ALAG = TVALAG*exp(ETA_ALAG);

    K20 = Cl/Vc;
    Cc = centr/Vc;

    d/dt(depot) = -Ka*depot;
    d/dt(centr) = Ka*depot - K20*centr;
    d/dt(AUC)   = Cc;

    alag(depot) = ALAG;

  }),
  error_model = function(f,sigma) {
    dv <- cbind(f,1)
    g <- diag(dv%*%sigma%*%t(dv))
    return(sqrt(g))
  },
  theta = c(THETA_Cl=4.0, THETA_Vc=70.0, THETA_Ka=1.0, THETA_ALAG=1.0),
  omega = lotri::lotri({ETA_Cl + ETA_Vc + ETA_Ka + ETA_ALAG ~
      c(0.2,
        0, 0.2,
        0, 0, 0.2,
        0, 0, 0, 0.2)}),
  sigma = lotri::lotri({prop + add ~ c(0.05,0.0,0.00)})
)

#patient data
lagpat <- data.frame(ID=1,
                     TIME=c(0,1.5,4.4,7.1,24.6,72),
                     EVID=c(1,0,0,0,0,1),
                     AMT=c(10000,0,0,0,0,10000),
                     CMT=c(1,2,2,2,2,1),
                     II=c(0,0,0,0,0,24),
                     ADDL=c(0,0,0,0,0,6),
                     MDV=c(1,0,0,0,0,1),
                     DV=c(NA,67.4275,115.4490,83.7768,21.8323,NA))

set.seed(42)

test_that("poso_estim_map provides estimates even when predicted concentrations
          are zero", {
  expect_equal(poso_estim_map(posologyr(mod_run003,lagpat,nocb=TRUE),
                              return_ofv=TRUE,
                              return_model=F)$eta,
            c(ETA_Cl=0.396185952,ETA_Vc=0.0370173441,
              ETA_Ka=0.0825083501,ETA_ALAG=-0.0967599813),tolerance=1e-3)
})
