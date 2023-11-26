mod_daptomycin_Dvorchik_AAC2004 <- list(
  ppk_model = rxode2::RxODE({
    centr(0) = 0;

    TVCl = ((THETA_Cl+0.00514*(ClCr-91.2))+0.14*(TEMP-37.2))*0.8^(SEX);
    TVVc = THETA_Vc;
    TVVp = (THETA_Vp+0.0458*(WT-75.1))*1.93;
    TVQ = THETA_Q+0.0593*(WT-75.1);

    Cl = TVCl*exp(ETA_Cl);
    Vc = TVVc*exp(ETA_Vc);
    Vp = TVVp*exp(ETA_Vp);
    Q = TVQ*exp(ETA_Q);

    K10 = Cl/Vc;
    K12 = Q / Vc;
    K21 = Q / Vp;
    Cc = centr/Vc;

    d/dt(centr) = -K10*centr - K12*centr + K21*periph;
    d/dt(periph) = K12*centr - K21*periph;
    d/dt(AUC) = Cc;
  }),
  error_model = function(f,sigma) {
    dv <- cbind(f,1)
    g <- diag(dv%*%sigma%*%t(dv))
    return(sqrt(g))
  },
  theta = c(THETA_Cl=0.807, THETA_Vc=4.80, THETA_Vp=3.13,THETA_Q=3.46),
  omega = lotri::lotri({ETA_Cl + ETA_Vc + ETA_Vp + ETA_Q  ~
      c(sqrt(log(0.306^2+1)),
        0, sqrt(log(0.567^2+1)),
        0, 0, sqrt(log(0.191^2+1)),
        0, 0, 0, sqrt(log(0.652^2+1)))}),
  covariates  = c("ClCr","TEMP","SEX","WT"),
  sigma = lotri::lotri({prop + add ~ c(0.00,0.00,4.72)})
)

df_patient_dap <- data.frame(ID=1,
                           TIME=c(0.0,23.5,24,48.75,49,97.25),
                           DV=c(NA,26.9,NA,48.7,NA,27.7),
                           AMT=c(1000,NA,1000,NA,1000,NA),
                           DUR=c(1,NA,1,NA,1,NA),
                           EVID=c(1,0,1,0,1,0),
                           SEX=1,WT=100,ClCr=53,TEMP=37.2)

test_that("Optimization results do not deviate from known values
          following TDM events when the last time is a float number", {
  expect_equal(poso_time_cmin(df_patient_dap,
                              mod_daptomycin_Dvorchik_AAC2004,
                              tdm=TRUE,target_cmin = 24)$time,
               55.1,
               tolerance=1e-2)
  expect_equal(poso_dose_conc(df_patient_dap,
                              mod_daptomycin_Dvorchik_AAC2004,
                              tdm=TRUE,time_c = 105,time_dose = 104,
                              target_conc = 60,duration = 1)$dose,
               323.538,
               tolerance=1e-3)
  expect_equal(poso_dose_auc(df_patient_dap,
                             mod_daptomycin_Dvorchik_AAC2004,
                             tdm=TRUE,time_auc = 24,time_dose = 104,
                             target_auc = 666,duration = 1)$dose,
               222.902,
               tolerance=1e-3)
})

df_patient_dap_int <- data.frame(ID=1,
                             TIME=c(0.0,23.5,24,48.75,49,98),
                             DV=c(NA,26.9,NA,48.7,NA,27.7),
                             AMT=c(1000,NA,1000,NA,1000,NA),
                             DUR=c(1,NA,1,NA,1,NA),
                             EVID=c(1,0,1,0,1,0),
                             SEX=1,WT=100,ClCr=53,TEMP=37.2)

test_that("Optimization results do not deviate from known values
          following TDM events when the last time is an integer", {
            expect_equal(poso_time_cmin(df_patient_dap_int,
                                        mod_daptomycin_Dvorchik_AAC2004,
                                        tdm=TRUE,target_cmin = 24)$time,
                         55.9,
                         tolerance=1e-2)
            expect_equal(poso_dose_conc(df_patient_dap_int,
                                        mod_daptomycin_Dvorchik_AAC2004,
                                        tdm=TRUE,time_c = 105,time_dose = 104,
                                        target_conc = 60,duration = 1)$dose,
                         322.9,
                         tolerance=1e-3)
            expect_equal(poso_dose_auc(df_patient_dap_int,
                                       mod_daptomycin_Dvorchik_AAC2004,
                                       tdm=TRUE,time_auc = 24,time_dose = 104,
                                       target_auc = 666,duration = 1)$dose,
                         214,
                         tolerance=1e-3)
})

df_patient_dap_tdm <- data.frame(ID=1,
                                 TIME=c(0.0,23.5,24),
                                 DV=c(NA,26.9,NA),
                                 AMT=c(1000,NA,1000),
                                 DUR=c(1,NA,1),
                                 EVID=c(1,0,1),
                                 SEX=1,WT=100,ClCr=53,TEMP=37.2)

test_that("The functions issue a warning when parameters are ignored because
          TDM=true", {
            #test issue #45
            #estim_method is ignored
            expect_warning(poso_time_cmin(df_patient_dap_tdm,
                                          mod_daptomycin_Dvorchik_AAC2004,
                                          tdm=TRUE,target_cmin = 24,
                                          estim_method="sir"))
            #dose is ignored
            expect_warning(poso_time_cmin(df_patient_dap_tdm,
                                          mod_daptomycin_Dvorchik_AAC2004,
                                          tdm=TRUE,target_cmin = 24,
                                          dose=100))
            #duration is ignored
            expect_warning(poso_time_cmin(df_patient_dap_tdm,
                                          mod_daptomycin_Dvorchik_AAC2004,
                                          tdm=TRUE,target_cmin = 24,
                                          duration=0.5))
            #interdose_interval is ignored
            expect_warning(poso_time_cmin(df_patient_dap_tdm,
                                          mod_daptomycin_Dvorchik_AAC2004,
                                          tdm=TRUE,target_cmin = 24,
                                          interdose_interval=12))
            #add_dose is ignored
            expect_warning(poso_time_cmin(df_patient_dap_tdm,
                                          mod_daptomycin_Dvorchik_AAC2004,
                                          tdm=TRUE,target_cmin = 24,
                                          add_dose=10))
            #indiv_param is ignored
            expect_warning(poso_time_cmin(df_patient_dap_tdm,
                                          mod_daptomycin_Dvorchik_AAC2004,
                                          tdm=TRUE,target_cmin = 24,
                                          indiv_param="anything"))
            #estim_method is ignored
            expect_warning(poso_dose_auc(df_patient_dap_tdm,
                                         mod_daptomycin_Dvorchik_AAC2004,
                                         tdm=TRUE,time_dose=122,time_auc=12,
                                         target_auc = 400,
                                         estim_method="sir"))
            #time_dose must be provided
            expect_error(poso_dose_auc(df_patient_dap_tdm,
                                         mod_daptomycin_Dvorchik_AAC2004,
                                         tdm=TRUE,time_auc=12,
                                         target_auc = 400,
                                         time_dose=NULL))
            #starting_time is ignored
            expect_warning(poso_dose_auc(df_patient_dap_tdm,
                                         mod_daptomycin_Dvorchik_AAC2004,
                                         tdm=TRUE,time_dose=122,time_auc=12,
                                         target_auc = 400,
                                         starting_time=24))
            #interdose_interval is ignored
            expect_warning(poso_dose_auc(df_patient_dap_tdm,
                                         mod_daptomycin_Dvorchik_AAC2004,
                                         tdm=TRUE,time_dose=122,time_auc=12,
                                         target_auc = 400,
                                         interdose_interval=12))
            #add_dose is ignored
            expect_warning(poso_dose_auc(df_patient_dap_tdm,
                                         mod_daptomycin_Dvorchik_AAC2004,
                                         tdm=TRUE,time_dose=122,time_auc=12,
                                         target_auc = 400,
                                         add_dose=10))
            #indiv_param is ignored
            expect_warning(poso_dose_auc(df_patient_dap_tdm,
                                         mod_daptomycin_Dvorchik_AAC2004,
                                         tdm=TRUE,time_dose=122,time_auc=12,
                                         target_auc = 400,
                                         indiv_param="anything"))
            #estim_method is ignored
            expect_warning(poso_dose_conc(df_patient_dap_tdm,
                                         mod_daptomycin_Dvorchik_AAC2004,
                                         tdm=TRUE,time_c=123,time_dose=122,
                                         target_conc=50,
                                         estim_method="sir"))
            #time_dose must be provided
            expect_error(poso_dose_conc(df_patient_dap_tdm,
                                         mod_daptomycin_Dvorchik_AAC2004,
                                         tdm=TRUE,time_c=123,target_conc=50,
                                         time_dose=NULL))
            #time_dose must happen before time_c
            expect_error(poso_dose_conc(df_patient_dap_tdm,
                                         mod_daptomycin_Dvorchik_AAC2004,
                                         tdm=TRUE,target_conc=50,
                                         time_c=121,time_dose=122))
            #interdose_interval is ignored
            expect_warning(poso_dose_conc(df_patient_dap_tdm,
                                         mod_daptomycin_Dvorchik_AAC2004,
                                         tdm=TRUE,time_c=123,time_dose=122,
                                         target_conc=50,
                                         interdose_interval=12))
            #add_dose is ignored
            expect_warning(poso_dose_conc(df_patient_dap_tdm,
                                         mod_daptomycin_Dvorchik_AAC2004,
                                         tdm=TRUE,time_c=123,time_dose=122,
                                         target_conc=50,
                                         add_dose=10))
            #indiv_param is ignored
            expect_warning(poso_dose_conc(df_patient_dap_tdm,
                                         mod_daptomycin_Dvorchik_AAC2004,
                                         tdm=TRUE,time_c=123,time_dose=122,
                                         target_conc=50,
                                         indiv_param="anything"))
})
