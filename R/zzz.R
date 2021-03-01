#.onLoad <- function(libname, pkgname){
#  compiled_models <- data(package = pkgname)
#  compiled_models <- compiled_models$results[,3]
#  for (i in 1:length(compiled_models)){
#   #if(TRUE){
#   #  print(get(compiled_models[1], envir = parent.env(environment())))
#   #}
#   #if (!get(compiled_models[1],
#   #        envir = parent.env(environment()))$ppk_model$isValid()){
#   #  print(get(compiled_models[1], envir = parent.env(environment())))
#   #  assign(compiled_models[i],
#   #         `[[<-`(get(compiled_models[1],
#   #                    envir = parent.env(environment())), 'ppk_model',
#   #         value = RxODE::RxODE(get(compiled_models[1],
#   #                                  envir = parent.env(environment()))$ppk_model)),
#   #         envir = parent.env(environment()))
#   #}
# }
# invisible()
#}
