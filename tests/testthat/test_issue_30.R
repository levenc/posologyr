# model with a single IOV
mod_MTXHD_joerger2011 <- list(
  ppk_model   = rxode2::rxode({
    centr(0) = 0;
    TVCl  = THETA_Cl*(CLCREAT/95)^0.28*(BSA/1.75)^0.15;
    TVVc  = THETA_Vc;
    TVVp  = THETA_Vp;
    TVQ   = THETA_Q;
    Cl    = TVCl*exp(ETA_Cl+KAPPA_Cl);
    Vc    = TVVc;
    Vp    = TVVp*exp(ETA_Vp);
    Q     = TVQ;
    ke    = Cl/Vc;
    k12   = Q/Vc;
    k21   = Q/Vp;
    Cc    = centr/Vc;
    d/dt(centr)  = - ke*centr - k12*centr + k21*periph;
    d/dt(periph) =            + k12*centr - k21*periph;
    d/dt(AUC)    =   Cc;
  }),
  error_model = function(f,sigma){
    dv <- cbind(f,1)
    g  <- diag(dv%*%sigma%*%t(dv))
    return(sqrt(g))
  },
  theta = c(THETA_Cl=10.8, THETA_Vc=34, THETA_Vp=6.3,THETA_Q=0.35),
  omega = lotri::lotri({ETA_Cl + ETA_Vc + ETA_Vp + ETA_Q ~
      c(0.10878,
        0    ,  0 ,
        0    ,  0     ,   0.116808,
        0    ,  0     ,   0    ,    0)}),
  pi_matrix = lotri::lotri({KAPPA_Cl ~
      c(0.017)}),
  covariates  = c("CLCREAT","BSA"),
  sigma       = lotri::lotri({prop + add ~ c(0.318,0.0,0.00)}))

# patient record
df_patientA_tdm <-data.frame(ID=1,
                             TIME=c(0,24,28,34,36,42,48),
                             DV=c(NA,25.0,18.5,11.4,10.32,6.59,4.58),
                             AMT=c(11723,0,0,0,0,0,0),
                             DUR=c(3,NA,NA,NA,NA,NA,NA),
                             EVID=c(1,0,0,0,0,0,0),
                             OCC=c(1,1,1,1,1,1,1),
                             CLCREAT=c(59,120,130,156,175,175,175),
                             BSA=2.049,
                             OCC=1)

# estimation
patAmap <- poso_estim_map(posologyr(mod_MTXHD_joerger2011,df_patientA_tdm))


test_that("poso_estim_map provides estimates for models with a single IOV", {
  expect_equal(patAmap$eta,
            c(ETA_Cl=-0.81668461,ETA_Vc=0.0,
              ETA_Vp=0.08900547,ETA_Q=0.0),tolerance=1e-3)
})
