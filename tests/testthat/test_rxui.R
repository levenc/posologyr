pk_turnover_emax3 <- function() {
  ini({
    tktr <- log(1)
    tka <- log(1)
    tcl <- log(0.1)
    tv <- log(10)
    ##
    eta_ktr ~ 1
    eta_ka ~ 1
    eta_cl ~ 2
    eta_v ~ 1
    prop_err <- 0.1
    pkadd_err <- 0.1
    ##
    temax <- logit(0.8)
    tec50 <- log(0.5)
    tkout <- log(0.05)
    te0 <- log(100)
    ##
    eta_emax ~ .5
    eta_ec50  ~ .5
    eta_kout ~ .5
    eta_e0 ~ .5
    ##
    pdadd_err <- 10
  })
  model({
    ktr <- exp(tktr + eta_ktr)
    ka <- exp(tka + eta_ka)
    cl <- exp(tcl + eta_cl)
    v <- exp(tv + eta_v)
    emax = expit(temax+eta_emax)
    ec50 =  exp(tec50 + eta_ec50)
    kout = exp(tkout + eta_kout)
    e0 = exp(te0 + eta_e0)
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
    cp ~ prop(prop_err) + add(pkadd_err)
    effect ~ add(pdadd_err) | pca
  })
}

pk_turnover_emax3 <- rxode2::rxode2(pk_turnover_emax3)
