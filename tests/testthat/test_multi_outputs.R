# Difference from the usual warfarin PKPD model: some IIV are set to zero to
# reduce the dimension of omega and reduce the execution time of the test

if(Sys.getenv("POSOLOGYR_DEV_MACHINE")=="TRUE"){
mod_warfarin_nlmixr <- list(
  ppk_model   = rxode2::rxode({
    ktr <- exp(THETA_ktr + ETA_ktr)
    ka <- exp(THETA_ka + ETA_ka)
    cl <- exp(THETA_cl + ETA_cl)
    v <- exp(THETA_v + ETA_v)
    emax = expit(THETA_emax + ETA_emax)
    ec50 =  exp(THETA_ec50 + ETA_ec50)
    kout = exp(THETA_kout + ETA_kout)
    e0 = exp(THETA_e0 + ETA_e0)
    ##
    DCP = center/v
    PD=1-emax*DCP/(ec50+DCP)
    ##
    effect(0) = e0
    kin = e0*kout
    ##
    d/dt(depot) = -ktr * depot
    d/dt(gut) =  ktr * depot -ka * gut
    d/dt(center) =  ka * gut - cl / v * center
    d/dt(effect) = kin*PD -kout*effect
    ##
    cp = center / v
    pca = effect
  }),
  error_model = list(
    cp = function(f,sigma){
      g <- sigma[1]^2 + (sigma[2]^2)*(f^2)
      return(sqrt(g))
    },
    pca = function(f,sigma){
      g <- sigma[1]^2 + (sigma[2]^2)*(f^2)
      return(sqrt(g))
    }
  ),
  theta = c(THETA_ktr=0.106,
            THETA_ka=-0.087,
            THETA_cl=-2.03,
            THETA_v=2.07,
            THETA_emax=3.4,
            THETA_ec50=0.00724,
            THETA_kout=-2.9,
            THETA_e0=4.57),
  omega = lotri::lotri({ETA_ktr + ETA_ka + ETA_cl + ETA_v + ETA_emax + ETA_ec50 +
      ETA_kout + ETA_e0 ~
      c(1.024695,
        0.00     ,   0.9518403 ,
        0.00     ,   0.00 ,   0.5300943  ,
        0.00    ,    0.00,   0.00,    0,#0.4785394,
        0.00    ,    0.00,   0.00,   0.00,   0.7134424,
        0.00    ,    0.00,   0.00,   0.00,   0.00,   0.7204165,
        0.00    ,    0.00,   0.00,   0.00,   0.00,   0.00,   0,#0.3563706,
        0.00    ,    0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0)}),#0.2660827)}),
  sigma       = list(
    cp=c(additive_a = 0.144, proportional_b = 0.15),
    pca=c(additive_a = 3.91, proportional_b = 0.0)
    )
  )

warf_01 <- data.frame(ID=1,
                      TIME=c(0.0,0.5,1.0,2.0,3.0,6.0,9.0,12.0,24.0,24.0,36.0,
                             36.0,48.0,48.0,72.0,72.0,96.0,120.0,144.0),
                      DV=c(0.0,0.0,1.9,3.3,6.6,9.1,10.8,8.6,5.6,44.0,4.0,27.0,
                           2.7,28.0,0.8,31.0,60.0,65.0,71.0),
                      DVID=c("cp","cp","cp","cp","cp","cp","cp","cp","cp","pca",
                             "cp","pca","cp","pca","cp","pca","pca","pca","pca"),
                      EVID=c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                      AMT=c(100,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
set.seed(1)
map_warf_01  <- poso_estim_map(dat=warf_01,
                               prior_model=mod_warfarin_nlmixr,
                               return_model=TRUE)

test_that("MAP estimates match nlmixr posthoc estimates for multi endpoints models", {
  expect_equal(map_warf_01$model$ktr[1], 0.6392321, tolerance=1e-3)
  expect_equal(map_warf_01$model$ka[1], 0.5294417, tolerance=1e-3)
  expect_equal(map_warf_01$model$cl[1], 0.2852736, tolerance=1e-3)
  expect_equal(map_warf_01$model$v[1], 7.924823, tolerance=1e-3)
  expect_equal(map_warf_01$model$emax[1], 0.9304377, tolerance=1e-3)
  expect_equal(map_warf_01$model$ec50[1], 0.7909184, tolerance=1e-3)
  expect_equal(map_warf_01$model$kout[1], 0.05502322, tolerance=1e-3)
  expect_equal(map_warf_01$model$e0[1], 96.54411, tolerance=1e-3)
})
}
## estimates from nlmixr
#pk.turnover.emax3 <- function() {
#  ini({
#   tktr <- 0.106
#   tka <- -0.087
#   tcl <- -2.03
#   tv <- 2.07
#   ##
#   eta.ktr ~ sqrt(1.05)
#   eta.ka ~ sqrt(0.906)
#   eta.cl ~ 0#sqrt(0.281)
#   eta.v ~ 0#sqrt(0.229)
#   prop.err <- 0.15
#   pkadd.err <- 0.144
#   ##
#   temax <- 3.4
#   tec50 <- 0.00724
#   tkout <- -2.9
#   te0 <- 4.57
#   ##
#   eta.emax ~ sqrt(.509)
#   eta.ec50  ~ sqrt(.519)
#   eta.kout ~ 0#sqrt(.127)
#   eta.e0 ~ 0#sqrt(.0708)
#   ##
#   pdadd.err <- 3.91
# })
# model({
#   ktr <- exp(tktr + eta.ktr)
#   ka <- exp(tka + eta.ka)
#   cl <- exp(tcl + eta.cl)
#   v <- exp(tv + eta.v)
#   emax = expit(temax+eta.emax)
#   ec50 =  exp(tec50 + eta.ec50)
#   kout = exp(tkout + eta.kout)
#   e0 = exp(te0 + eta.e0)
#   ##
#   DCP = center/v
#   PD=1-emax*DCP/(ec50+DCP)
#   ##
#   effect(0) = e0
#   kin = e0*kout
#   ##
#   d/dt(depot) = -ktr * depot
#   d/dt(gut) =  ktr * depot -ka * gut
#  d/dt(center) =  ka * gut - cl / v * center
#   d/dt(effect) = kin*PD -kout*effect
#   ##
#   cp = center / v
#    cp ~ prop(prop.err) + add(pkadd.err)
#    effect ~ add(pdadd.err) | pca
#  })
#}
# nlmixr(pk.turnover.emax3,warfarin,"focei",
#       control=foceiControl(maxOuterIterations=0))
