# prior model
mod_run202 <- list(
  ppk_model   = rxode2::rxode({
    depot(0) = 0;
    centr(0) = 0;

    TVVmax  = THETA_Vmax;
    TVVc  = THETA_Vc;
    TVKa  = THETA_Ka;
    TVKm  = THETA_Km;

    Vmax  = TVVmax*exp(ETA_Vmax);
    Vc    = TVVc*exp(ETA_Vc);
    Ka    = TVKa*exp(ETA_Ka);
    Km    = TVKm*exp(ETA_Km);
    Cc    = centr/Vc;

    d/dt(depot)  = - Ka*depot;
    d/dt(centr)  =   Ka*depot - Vmax*(centr/Vc)/(Km+(centr/Vc));
    d/dt(AUC)    =   Cc;
  }),
  error_model = function(f,sigma){
    dv <- cbind(f,1)
    g  <- diag(dv%*%sigma%*%t(dv))
    return(sqrt(g))
  },
  theta = c(THETA_Vmax=10000, THETA_Vc=70.0, THETA_Ka=1.0, THETA_Km=2500),
  omega = lotri::lotri({ETA_Vmax + ETA_Vc + ETA_Ka + ETA_Km ~
      c(0.2   ,
        0     ,     0.2,
        0     ,       0,     0,
        0     ,       0,     0,   0.2)}),
  sigma       = lotri::lotri({prop + add ~ c(0.05,0.0,0.00)}))

# patient data
kmpat <- data.frame(ID=1,
                    TIME=c(0,1.5,4.4,7.1,24.6,72),
                    EVID=c(1,0,0,0,0,1),
                    AMT=c(10000,0,0,0,0,10000),
                    CMT=2,II=c(0,0,0,0,0,24),
                    ADDL=c(0,0,0,0,0,6),
                    RATE=c(10000,0,0,0,0,10000),
                    MDV=c(1,0,0,0,0,1),
                    DV=c(NA,111.57400,97.86790,64.29320,9.91143,NA))

set.seed(42)

test_that("poso_estim_map provides estimates even when the IIV is zero", {
  expect_equal(poso_estim_map(kmpat,mod_run202,nocb=TRUE,
                              return_ofv=TRUE,
                              return_model=F)$eta,
            c(ETA_Vmax=0.3799632360,ETA_Vc=8.405155e-02,
              ETA_Ka=0.00000000,ETA_Km=-0.3684082160),tolerance=1e-3)
})
