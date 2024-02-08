.onLoad <- function(libname, pkgname) {
  rxode2::.s3register("rxode2::rxUiGet", "posologyr")
  rxode2::.s3register("rxode2::rxUiGet", "posologyr_ppk_model")
  rxode2::.s3register("rxode2::rxUiGet", "posologyr_error_model")
  rxode2::.s3register("rxode2::rxUiGet", "posologyr_sigma")
}
