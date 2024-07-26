mod_ganciclovir_Caldes_AAC2009 <- function() {
  ini({
    THETA_cl  <- 7.49
    THETA_v1  <- 31.90
    THETA_cld <- 10.20
    THETA_v2  <- 32.0
    THETA_ka  <- 0.895
    THETA_baf <- 0.825
    ETA_cl ~ 0.107
    ETA_v1 ~ 0.227
    ETA_ka ~ 0.464
    ETA_baf ~ 0.049
    add.sd <- 0.465
    prop.sd <- 0.143
  })
  model({
    TVcl  = THETA_cl*(ClCr/57);
    TVv1  = THETA_v1;
    TVcld = THETA_cld;
    TVv2  = THETA_v2;
    TVka  = THETA_ka;
    TVbaf = THETA_baf;

    cl  = TVcl*exp(ETA_cl);
    v1  = TVv1*exp(ETA_v1);
    cld = TVcld;
    v2  = TVv2;
    ka  = TVka*exp(ETA_ka);
    baf = TVbaf*exp(ETA_baf);

    k10 = cl/v1;
    k12 = cld / v1;
    k21 = cld / v2;
    Cc = centr/v1;

    d/dt(depot)  = -ka*depot
    d/dt(centr)  =  ka*depot - k10*centr - k12*centr + k21*periph;
    d/dt(periph) =                         k12*centr - k21*periph;
    d/dt(AUC)    = Cc;

    f(depot)=baf;
    alag(depot)=0.382;

    Cc ~ add(add.sd) + prop(prop.sd) + combined1()
  })
}

df_patient <- data.frame(ID=1,TIME=c(0,1,2,6),
                              DV=c(NA,3.4,2.7,1.4),
                              AMT=c(450,0,0,0),
                              EVID=c(1,0,0,0),
                              ClCr=80)

test_that("Optimal dose can be estimated for dosing in any compartment", {
  # poso_dose_conc--------------------------------------------------------------
  # cmt_dose = "depot"
  expect_equal(poso_dose_conc(df_patient,mod_ganciclovir_Caldes_AAC2009,
                              tdm=TRUE,time_c = 25,time_dose = 24.1,
                              target_conc = 6,cmt_dose = "depot")$dose,
               1014.197,tolerance=1e-3)
  expect_equal(poso_dose_conc(dat=df_patient,
                              prior_model=mod_ganciclovir_Caldes_AAC2009,
                              time_c=1,target_conc=6,cmt_dose = "depot")$dose,
               928.35,tolerance=1e-3)
  # cmt_dose="centr"
  expect_equal(poso_dose_conc(df_patient,mod_ganciclovir_Caldes_AAC2009,
                              tdm=TRUE,time_c = 25,time_dose = 24.1,
                              target_conc = 6,cmt_dose = "centr")$dose,
               367.7199,tolerance=1e-3)
  expect_equal(poso_dose_conc(dat=df_patient,
                              prior_model=mod_ganciclovir_Caldes_AAC2009,
                              time_c=1,
                              target_conc=6,cmt_dose = "centr")$dose,
               393.84,tolerance=1e-3)
  # poso_time_cmin -------------------------------------------------------------
  # cmt_dose = 1
  expect_equal(poso_time_cmin(df_patient,mod_ganciclovir_Caldes_AAC2009,
                              tdm=TRUE,target_cmin = 0.5,from=1)$time,
               13,tolerance=1e-3)
  expect_equal(poso_time_cmin(df_patient,mod_ganciclovir_Caldes_AAC2009,
                              target_cmin = 0.5,from=1,dose=500)$time,
               13.7,tolerance=1e-3)
  # cmt_dose="centr"
  expect_equal(poso_time_cmin(df_patient,mod_ganciclovir_Caldes_AAC2009,
                               tdm=TRUE,from=1,target_cmin = 0.5,
                               cmt_dose = "centr")$time,
               13,tolerance=1e-3)
  expect_equal(poso_time_cmin(df_patient,mod_ganciclovir_Caldes_AAC2009,
                              target_cmin = 0.5,from=1,dose=500,
                              cmt_dose = "centr")$time,
               14.6,tolerance=1e-3)
  # poso_dose_auc -------------------------------------------------------------
  # cmt_dose = 1
  expect_equal(poso_dose_auc(df_patient,mod_ganciclovir_Caldes_AAC2009,
                             tdm=TRUE,time_auc = 24,time_dose = 24.1,
                             target_auc = 50)$dose,
               924.458,tolerance=1e-3)
  expect_equal(poso_dose_auc(df_patient,mod_ganciclovir_Caldes_AAC2009,
                             time_auc = 24,target_auc = 50)$dose,
               938.59,tolerance=1e-3)
  # cmt_dose="centr"
  expect_equal(poso_dose_auc(df_patient,mod_ganciclovir_Caldes_AAC2009,
                             tdm=TRUE,time_auc = 24,time_dose = 24.1,
                             target_auc = 50,cmt_dose = "centr")$dose,
               645.857,tolerance=1e-3)
  expect_equal(poso_dose_auc(df_patient,mod_ganciclovir_Caldes_AAC2009,
                             time_auc = 24,target_auc = 50,
                             cmt_dose = "centr")$dose,
               655.729,tolerance=1e-3)
  # poso_dose_auc -------------------------------------------------------------
  # cmt_dose = 1
  expect_equal(poso_inter_cmin(df_patient,mod_ganciclovir_Caldes_AAC2009,
                               dose=500,target_cmin = 1)$interval,
               10.69046,tolerance=1e-3)
  # cmt_dose="centr"
  expect_equal(poso_inter_cmin(df_patient,mod_ganciclovir_Caldes_AAC2009,
                               dose=500,target_cmin = 1,
                               cmt_dose = "centr")$interval,
               11.34967,tolerance=1e-3)
})
