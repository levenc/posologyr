mod_tobramycin_2cpt_fictional <- list(
  ppk_model   = RxODE::RxODE({
    centr(0) = 0;
    ke = TVke*(CLCREAT/67.8)^0.89*(WT/66.4)^-1.09;
    V  = TVV*(WT/66.4)^0.80;
    Cc  = centr/V;
    d/dt(centr)  = - ke*centr - k12*centr + k21*periph;
    d/dt(periph) =            + k12*centr - k21*periph;
    d/dt(AUC)    =   Cc;
  }),
  error_model = function(f,xi){
    g <- xi[1] + xi[2]*f
    return(g)
  },
  pk_prior    = list( name = c('TVke','TVV','k12','k21'),
                      reference = c(TVke=0.21, TVV=19.8, k12=0.041, k21=0.12),
                      Omega = matrix(c(0.08075, 0      ,  0, 0,
                                       0      , 0.01203,  0, 0,
                                       0      , 0      ,  0, 0,
                                       0      , 0      ,  0, 0),
                                     ncol=4,byrow=TRUE)),
  covariates  = c("CLCREAT","WT"),
  xi          = c(additive_a = 0, proportional_b = 0.198))

mod_amoxicillin_oral_1cpt_fictional <- list(
  ppk_model   = RxODE::RxODE({
    depot(0) = 0
    centr(0) = 0
    Cl = TVCl*(CLCREAT/90)^0.62
    ke = Cl/V;
    Cc = centr/V;
    d/dt(depot) =-ka*depot;
    d/dt(centr) = ka*depot - ke*centr;
    d/dt(AUC)   = Cc;
  }),
  error_model = function(f,xi){
    g <- xi[1] + xi[2]*f
    return(g)
  },
  pk_prior    = list( name = c('ka','V','TVCl'),
                      reference = c(ka=0.52, V=0.56, TVCl=7.33),
                      Omega = matrix(c(0.0289, 0      ,  0,
                                       0     , 0.0256 , -0.01792,
                                       0     ,-0.01792,  0.0400),ncol=3,byrow=TRUE)),
  covariates  = c("CLCREAT"),
  xi          = c(additive_a = 0.24, proportional_b = 0.27))

## creations of the .Rdata files needed to embed the models
save(mod_tobramycin_2cpt_fictional,
     mod_amoxicillin_oral_1cpt_fictional,
     file="data/lib_ppk_model.rda")
