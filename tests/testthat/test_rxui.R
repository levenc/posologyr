test_that("readme model equivalent", {
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
      ETA_Vc ~ 0.2
      ETA_Ka ~ 0.2
      prop.sd <- sqrt(0.05)
    })
    model({
      TVCl <- THETA_Cl
      TVVc <- THETA_Vc
      TVKa <- THETA_Ka

      Cl <- TVCl*exp(ETA_Cl)
      Vc <- TVVc*exp(ETA_Vc)
      Ka <- TVKa*exp(ETA_Ka)

      K20 <- Cl/Vc
      Cc <- centr/Vc

      d/dt(depot) = -Ka*depot
      d/dt(centr) = Ka*depot - K20*centr
      Cc ~ prop(prop.sd)
    })
  }

  f <- poso_estim_map(patient_data, mod_run001)

  expect_equal(f$eta, c(ETA_Cl = 0.601903430923289, ETA_Vc = -0.429173629074746, ETA_Ka = 0.127847428673983),
               tolerance=1e-4)

})

test_that("mod_amikacin_Burdet2015", {

  skip_on_cran()

  mod_amikacin_Burdet2015 <- function() {
    ini({
      THETA_Cl=4.3
      THETA_Vc=15.9
      THETA_Vp=21.4
      THETA_Q=12.1
      ETA_Cl + ETA_Vc + ETA_Vp + ETA_Q ~
        c(0.1,
          0.01     ,   0.05 ,
          0.01     ,   0.02 ,   0.2  ,
          -0.06    ,   0.004,   0.003,    0.08)
      add_sd <- 0.2
      prop_sd <- 0.1
    })
    model({
      centr(0) = 0
      TVCl  = THETA_Cl*(CLCREAT4H/82)^0.7
      TVVc  = THETA_Vc*(TBW/78)^0.9*(PoverF/169)^0.4
      TVVp  = THETA_Vp
      TVQ   = THETA_Q
      Cl    = TVCl*exp(ETA_Cl)
      Vc    = TVVc*exp(ETA_Vc)
      Vp    = TVVp*exp(ETA_Vp)
      Q     = TVQ *exp(ETA_Q)
      ke    = Cl/Vc
      k12   = Q/Vc
      k21   = Q/Vp
      Cp    = centr/Vc
      d/dt(centr)  = - ke*centr - k12*centr + k21*periph
      d/dt(periph) =            + k12*centr - k21*periph
      #
      Cp ~ add(add_sd) + prop(prop_sd) + combined1()
    })
  }

  df_patientA <- data.frame(ID=1,TIME=c(0,1,6),
                            DV=c(NA,58,14),
                            EVID=c(1,0,0),
                            AMT=c(2000,0,0),
                            DUR=c(0.5,NA,NA),
                            CLCREAT4H=50,TBW=62,PoverF=169)

  patA_map <- poso_estim_map(dat=df_patientA,
                             prior_model=mod_amikacin_Burdet2015)

  expect_equal(patA_map$eta,
               c(ETA_Cl = 0.449952188793169, ETA_Vc = 0.273058665000001, ETA_Vp = 0.70615858706018, ETA_Q = -0.138843327030614),
               tolerance=1e-4)

  cmin <- poso_time_cmin(dat=df_patientA,
                 prior_model=mod_amikacin_Burdet2015,
                 tdm = TRUE,
                 target_cmin = 2.5)

  expect_equal(cmin$time, 33.9, tolerance=0.1)

  map_dose <- poso_dose_conc(dat=df_patientA,
                             prior_model=mod_amikacin_Burdet2015,
                             tdm=TRUE,
                             time_c = 35,               #target concentration at t = 35 h
                             time_dose = 34,            #dosing at t = 34 h
                             duration = 0.5,
                             target_conc = 80)

  expect_equal(map_dose$dose, 2447.917, tolerance = 1)

  map_interval <- poso_inter_cmin(dat=df_patientA,
                                  prior_model=mod_amikacin_Burdet2015,
                                  dose = map_dose$dose,
                                  duration = 0.5,
                                  target_cmin = 2.5)

  expect_equal(map_interval$interval, 38.57782, tolerance=0.1)

})

test_that("mod_vancomycin_Goti2018", {

  skip_on_cran()

  mod_vancomycin_Goti2018 <- function() {
    ini({
      THETA_Cl <- 4.5
      THETA_Vc <- 58.4
      THETA_Vp <- 38.4
      THETA_Q <- 6.5
      ETA_Cl ~ 0.147
      ETA_Vc ~ 0.510
      ETA_Vp ~ 0.282
      add.sd <- 3.4
      prop.sd <- 0.227
    })
    model({
      centr(0) = 0;
      TVCl  = THETA_Cl*(CLCREAT/120)^0.8*(0.7^DIAL);
      TVVc  = THETA_Vc*(WT/70)          *(0.5^DIAL);
      TVVp  = THETA_Vp;
      TVQ   = THETA_Q;
      Cl    = TVCl*exp(ETA_Cl);
      Vc    = TVVc*exp(ETA_Vc);
      Vp    = TVVp*exp(ETA_Vp);
      Q     = TVQ;
      ke    = Cl/Vc;
      k12   = Q/Vc;
      k21   = Q/Vp;
      Cc    = centr/Vc;
      d/dt(centr)  = - ke*centr - k12*centr + k21*periph;
      d/dt(periph) =            + k12*centr - k21*periph;
      d/dt(AUC) <- Cc
      Cc ~ add(add.sd) + prop(prop.sd) + combined1()
    })
  }

  df_patientB <- data.frame(ID=1,TIME=c(0.0,13.0,24.2,48),
                            DV=c(NA,12,NA,9.5),
                            AMT=c(2000,0,1000,0),
                            DUR=c(2,NA,2,NA),
                            EVID=c(1,0,1,0),
                            CLCREAT=65,WT=70,DIAL=0)

  patB_map <- poso_estim_map(dat=df_patientB,
                             prior_model=mod_vancomycin_Goti2018)

  AUC_map_first_dose <- patB_map$model$AUC[which(patB_map$model$time == 24)]

  expect_equal(AUC_map_first_dose, 337.8963, tolerance=1)

  AUC_map_second_dose <- patB_map$model$AUC[which(patB_map$model$time == 48)] - AUC_map_first_dose
  AUC_map_second_dose

  expect_equal(325.3071, AUC_map_second_dose, tolerance=1)

  a <- poso_dose_auc(dat=df_patientB,
                prior_model=mod_vancomycin_Goti2018,
                tdm=TRUE,
                time_auc=24,            #AUC24
                time_dose = 48,         #48 h: immediately following the last observation
                duration=2,             #infused over 2 h
                target_auc=400)

  expect_equal(a$dose, 1411.593, tolerance=1)


  # could possibly use steady state flags instead; regardless:
  a2 <- poso_dose_auc(dat=df_patientB,
                prior_model=mod_vancomycin_Goti2018,
                time_auc=24,
                starting_time=24*9,
                interdose_interval=24,
                add_dose=10,
                duration=2,
                target_auc=400)

  expect_equal(a$dose, 1198.267, tolerance=1)

  dur <- poso_dose_auc(dat=df_patientB,
                       prior_model=mod_vancomycin_Goti2018,
                       time_auc=24,
                       starting_time=24*9,
                       interdose_interval=24,
                       add_dose=10,
                       duration=24,
                       target_auc=400)

  expect_equal(1198.923, dur$dose, tolerance=1)

})

if(Sys.getenv("POSOLOGYR_DEV_MACHINE")=="TRUE"){ #skip on CRAN or github actions
test_that("warfarin example", {

  mod_warfarin_nlmixr <- function() {
    ini({
      THETA_ktr=0.106
      THETA_ka=-0.087
      THETA_cl=-2.03
      THETA_v=2.07
      THETA_emax=3.4
      THETA_ec50=0.00724
      THETA_kout=-2.9
      THETA_e0=4.57
      ETA_ktr ~ 1.024695
      ETA_ka ~ 0.9518403
      ETA_cl ~ 0.5300943
      ETA_v ~ 0.4785394
      ETA_emax ~ 0.7134424
      ETA_ec50 ~ 0.7204165
      ETA_kout ~ 0.3563706
      ETA_e0 ~ 0.2660827
      cp.sd <- 0.144
      cp.prop.sd <- 0.15
      pca.sd <- 3.91
    })
    model({
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
      cp ~ add(cp.sd) + prop(cp.prop.sd)
      pca ~ add(pca.sd)
    })
  }

  warf_01 <- data.frame(ID=1,
                        TIME=c(0.0,0.5,1.0,2.0,3.0,6.0,9.0,12.0,24.0,24.0,36.0,
                               36.0,48.0,48.0,72.0,72.0,96.0,120.0,144.0),
                        DV=c(0.0,0.0,1.9,3.3,6.6,9.1,10.8,8.6,5.6,44.0,4.0,27.0,
                             2.7,28.0,0.8,31.0,60.0,65.0,71.0),
                        DVID=c("cp","cp","cp","cp","cp","cp","cp","cp","cp","pca",
                               "cp","pca","cp","pca","cp","pca","pca","pca","pca"),
                        EVID=c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                        AMT=c(100,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))

  map_warf_01  <- poso_estim_map(warf_01,mod_warfarin_nlmixr)

  expect_equal(c(ETA_ktr = -0.434704890568383, ETA_ka = -0.694024512420889, ETA_cl = 0.823372753994899, ETA_v = -0.0230522657927063, ETA_emax = 0.0684548537529647, ETA_ec50 = 0.141904594137433, ETA_kout = -0.198702717691606, ETA_e0 = -0.123821088809355),
               map_warf_01$eta, tolerance =1e-4)

})
test_that("+var() tests", {

  library(rxode2)

  mod_warfarin_nlmixr <- function() {
    ini({
      THETA_ktr=0.106
      THETA_ka=-0.087
      THETA_cl=-2.03
      THETA_v=2.07
      THETA_emax=3.4
      THETA_ec50=0.00724
      THETA_kout=-2.9
      THETA_e0=4.57
      ETA_ktr ~ 1.024695
      ETA_ka ~ 0.9518403
      ETA_cl ~ 0.5300943
      ETA_v ~ 0.4785394
      ETA_emax ~ 0.7134424
      ETA_ec50 ~ 0.7204165
      ETA_kout ~ 0.3563706
      ETA_e0 ~ 0.2660827
      cp.sd <- 0.144^2
      cp.prop.sd <- 0.15^2
      pca.sd <- 3.91
    })
    model({
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
      cp ~ add(cp.sd) + prop(cp.prop.sd) + var()
      pca ~ add(pca.sd)
    })
  }

  mod_warfarin_nlmixr <- try(mod_warfarin_nlmixr(), silent=TRUE)

  if (inherits(mod_warfarin_nlmixr, "rxUi")) {

    expect_equal(mod_warfarin_nlmixr$posologyr_sigma,
                 list(cp = c(cp.sd = 0.144, cp.prop.sd = 0.15), pca = c(pca.sd = 3.91)))

    mod_warfarin_nlmixr <- function() {
      ini({
        THETA_ktr=0.106
        THETA_ka=-0.087
        THETA_cl=-2.03
        THETA_v=2.07
        THETA_emax=3.4
        THETA_ec50=0.00724
        THETA_kout=-2.9
        THETA_e0=4.57
        ETA_ktr ~ 1.024695
        ETA_ka ~ 0.9518403
        ETA_cl ~ 0.5300943
        ETA_v ~ 0.4785394
        ETA_emax ~ 0.7134424
        ETA_ec50 ~ 0.7204165
        ETA_kout ~ 0.3563706
        ETA_e0 ~ 0.2660827
        cp.prop.sd <- 0.15^2
        pca.sd <- 3.91^2
      })
      model({
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
        cp ~  prop(cp.prop.sd)  + var()
        pca ~ add(pca.sd) + var()
      })
    }

    mod_warfarin_nlmixr <- mod_warfarin_nlmixr()

    expect_equal(mod_warfarin_nlmixr$posologyr_sigma,
                 list(cp = c(cp.prop.sd = 0.15), pca = c(pca.sd = 3.91)))

  }

})
}
