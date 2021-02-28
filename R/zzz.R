.onLoad <- function(libname, pkgname){
  compiled_models <- ls('package:posologyr', pattern='mod_')
  for (i in 1:length(compiled_models)){
    #get and assign allow
    if (!get(compiled_models[i])$ppk_model$isValid()){
      assign(compiled_models[i],
             `[[<-`(get(compiled_models[i]), 'ppk_model',
             value = RxODE::RxODE(get(compiled_models[i]$ppk_model))),
             envir = as.environment('package:posologyr'))
    }
  }
  invisible()
}
