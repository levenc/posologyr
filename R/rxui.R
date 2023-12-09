#' This creates the posologyr error lines from a rxui model
#'
#' @param line line to parse
#'
#' @return error lines for posology
#'
#' @export
#'
#' @keywords internal
#'
#' @author Matthew L. Fidler
#'
posologyr_error_lines <- function(line) {
  UseMethod("posologyr_error_lines")
}

#' @rdname posologyr_error_lines
#' @export
posologyr_error_lines.norm <- function(line) {
  env <- line[[1]]
  pred1 <- line[[2]]
  ret <- vector("list", 6)
  .yj <- as.double(pred1$transform) - 1
  ret[[1]] <- bquote(rx_yj_ ~ .(.yj))
  ret[[2]] <- bquote(rx_lambda_~.(rxode2::.rxGetLambdaFromPred1AndIni(env, pred1)))
  ret[[3]] <- bquote(rx_low_ ~ .(rxode2::.rxGetLowBoundaryPred1AndIni(env, pred1)))
  ret[[4]] <- bquote(rx_hi_ ~ .(rxode2::.rxGetHiBoundaryPred1AndIni(env, pred1)))
  ret[[5]] <- bquote(rx_pred_f_ ~ .(rxode2::.rxGetPredictionF(env, pred1)))
  if (length(env$predDf$cond) == 1L) {
    ret[[6]] <- bquote(Cc <- .(rxode2::.rxGetPredictionFTransform(env, pred1, .yj)))
  } else {
    ret[[6]] <- bquote(.(str2lang(paste0("rxEndpoint", pred1$dvid))) <- .(rxode2::.rxGetPredictionFTransform(env, pred1, .yj)))
  }
  ret
}

#' @rdname posologyr_error_lines
#' @export
posologyr_error_lines.t <- function(line) {
  stop("t isn't supported yet", call.=FALSE)
}

#' @rdname posologyr_error_lines
#' @export
posologyr_error_lines.default  <- function(line) {
  stop("distribution not supported", call.=FALSE)
}

# This handles the errors for focei
create_posologyr_line_object <- function(x, line) {
  pred_df <- x$predDf
  if (line > nrow(pred_df)) {
    return(NULL)
  }
  pred_line <- pred_df[line, ]
  ret <- list(x, pred_line, line)
  class(ret) <- c(paste(pred_line$distribution), "posologyr_error_lines")
  ret
}

#' @rdname posologyr_error_lines
#' @export
posologyr_error_lines.rxUi <- function(line) {
  pred_df <- line$predDf
  lapply(seq_along(pred_df$cond), function(c) {
    mod <- create_posologyr_line_object(line, c)
    posologyr_error_lines(mod)
  })
}


#' @export
rxUiGet.posologyr_ppk_model <- function(x, ...) {
  ui <- x[[1]]
  if (is.null(ui$predDf)) {
    stop("need endpoint defined for now")
  }
  if (length(ui$predDf$cond) == 1L) {
    # Here Cc  needs to be the endpoint
    mv <- rxode2::rxModelVars(ui)
    if (any(mv$lhs == "Cc")) {
      ui <- rxode2::rxRename(ui, rxCc=Cc)
    }
  }
  rxode2::rxCombineErrorLines(ui, errLines=posologyr_error_lines(ui),
                              cmtLines=FALSE, dvidLine=FALSE)
}
attr(rxUiGet.posologyr_ppk_model, "desc") <- "posologyr ppk_model element"

#' @export
rxUiGet.posologyr_error_model <- function(x, ...) {

}
attr(rxUiGet.posologyr_error_model, "desc") <- "posologyr error_model element"


#' @export
rxUiGet.posologyr_sigma <- function(x, ...) {

}

#' @export
rxUiGet.posologyr <- function(x, ...) {
  ui <- x[[1]]
  list(ppk_model=rxUiGet.posologyr_ppk_model(x, ...),
       error_model=rxUiGet.posologyr_error_model(x, ...),
       theta=ui$theta,
       omega=ui$omega,
       sigma=rxUiGet.posologyr_sigma(x, ...))
}
attr(rxUiGet.posologyr, "desc") <- "posologyr model from ui"
