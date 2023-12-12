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
    ret[[6]] <- bquote(.(str2lang(pred1$cond)) <- .(rxode2::.rxGetPredictionFTransform(env, pred1, .yj)))
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
create_posologyr_line_object <- function(x, line, type="posologyr_error_lines") {
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

posologyr_ppk_model_rxui <- function(ui, auc=FALSE) {
}

#' @export
rxUiGet.posologyr_ppk_model <- function(x, ...) {
  ui <- x[[1]]
  if (is.null(ui$predDf)) {
    stop("need endpoint defined for now")
  }
  add_auc <- FALSE
  if (length(ui$predDf$cond) == 1L) {
    # Here Cc  needs to be the endpoint
    mv <- rxode2::rxModelVars(ui)
    if (any(mv$lhs == "Cc")) {
      ui <- rxode2::rxRename(ui, rxCc=Cc)
    }
    if (!any(mv$state == "AUC")) {
      # add AUC for single endpoint models
      states <- mv$state
      add_auc <- TRUE
    }
  }
  mod <- rxode2::rxCombineErrorLines(ui, errLines=posologyr_error_lines(ui),
                                     cmtLines=FALSE, dvidLine=FALSE, useIf = FALSE)
  mod[[1]] <- str2lang("rxode2::rxode2")
  if (add_auc) {
    message("Added AUC to model")
    m2 <- mod[[2]]
    mod[[2]] <- as.call(lapply(seq_len(length(m2) + 1), function(i) {
      if (i > length(m2)) return(str2lang("d/dt(AUC)=Cc"))
      m2[[i]]
    }))
  }
  mod
}
attr(rxUiGet.posologyr_ppk_model, "desc") <- "posologyr ppk_model element"

#'  Get the additive error model function or estimate from ui
#'
#' @param ui environment of rxode2 ui
#' @param pred1 prediction line for endpint
#' @param fun boolean to return function
#' @return function or estimate of error line
#' @noRd
#' @author Matthew L. Fidler
posologyr_get_error_model_add <- function(ui, pred1, fun=TRUE) {
  if (!is.na(pred1$a)) {
    stop("residual variability parameters must be in ini block (not more complex)",
         call.=FALSE)
  } else {
    .cnd <- pred1$cond
    .w <- which(ui$iniDf$err %in% c("add", "lnorm", "logitNorm", "probitNorm") & ui$iniDf$condition == .cnd)
    if (length(.w) == 1L) {
      .p1 <- setNames(ui$iniDf$est[.w], ui$iniDf$name[.w])
    } else {
      stop("cannot find additive standard deviation for '", .cnd, "'",
           ifelse(length(ui$predDf$condition) == 1L, "", "; this parameter could be estimated by another endpoint, to fix move outside of error expression."), call.=FALSE)
    }
  }
  if (!fun) {
    return(.p1)
  }
  # on standard deviation scale
  f <- function(f, sigma) {
    sigma[1]
  }
  f
}

#' Get the proportional error model function or estimate from ui
#'
#' @param ui environment of rxode2 ui
#' @param pred1 prediction line for endpint
#' @param fun boolean to return function
#' @return function or estimate of error line
#' @noRd
#' @author Matthew L. Fidler
posologyr_get_error_model_prop <- function(ui, pred1, fun=TRUE) {
  type <- as.character(pred1$errTypeF)
  if (!(type %in% c("untransformed", "none"))) {
    stop("f can only be untransformed for poslogyr", call.=FALSE)
  }
  if (!fun) {
    if (!is.na(pred1$b)) {
      .p1 <- str2lang(pred1$b)
    } else {
      .cnd <- pred1$cond
      .w <- which(ui$iniDf$err %in% c("prop", "propF", "propT") & ui$iniDf$condition == .cnd)
      if (length(.w) == 1L) {
        .p1 <- setNames(ui$iniDf$est[.w], ui$iniDf$name[.w])
      } else {
        stop("cannot find proportional standard deviation", call.=FALSE)
      }
    }
  }
  if (!fun) {
    return(.p1)
  }
  f <- function(f, sigma) {
    sigma[1] * f
  }
  f
}

#' Get the power error model function or estimate from ui
#'
#' @param ui environment of rxode2 ui
#' @param pred1 prediction line for endpint
#' @param fun boolean to return function
#' @return function or estimate of error line
#' @noRd
#' @author Matthew L. Fidler
posologyr_get_error_model_pow <- function(env, pred1, fun=TRUE) {
  stop("pow not supported with posologyr", call.=FALSE)
}

#' Get the add+prop error model function or estimate from ui
#'
#' @param ui environment of rxode2 ui
#' @param pred1 prediction line for endpint
#' @param fun boolean to return function
#' @return function or estimate of error line
#' @noRd
#' @author Matthew L. Fidler
posologyr_get_error_model_add_prop <- function(ui, pred1, fun=TRUE) {
  if (!is.na(pred1$a)) {
    stop("residual variability parameters must be in ini block (not more complex)",
         call.=FALSE)
  } else {
    .cnd <- pred1$cond
    .w <- which(ui$iniDf$err %in% c("add", "lnorm", "probitNorm", "logitNorm") & ui$iniDf$condition == .cnd)
    if (length(.w) == 1L) {
      .p1 <- setNames(ui$iniDf$est[.w], ui$iniDf$name[.w])
    } else {
      stop("cannot find additive standard deviation", call.=FALSE)
    }
  }
  if (!is.na(pred1$b)) {
    stop("residual variability parameters must be in ini block (not more complex)",
         call.=FALSE)
  } else {
    .cnd <- pred1$cond
    .w <- which(ui$iniDf$err %in% c("prop", "propT", "propF") & ui$iniDf$condition == .cnd)
    if (length(.w) == 1L) {
      .p2 <- setNames(ui$iniDf$est[.w], ui$iniDf$name[.w])
    } else {
      stop("cannot find proportional standard deviation", call.=FALSE)
    }
  }
  if (pred1$addProp == "default") {
    .addProp <- rxode2::rxGetControl(ui, "addProp", getOption("rxode2.addProp", "combined2"))
  } else {
    .addProp <- pred1$addProp
  }
  if (!fun) return(c(.p1, .p2))
  if (.addProp == "combined2") {
    f <- function(f, sigma) {
      sqrt(sigma[1]^2 + f^2 * sigma[2]^2)
    }
  } else {
    # combined1, standard deviations add
    f <- function(f, sigma) {
      sigma[1] + f * sigma[2]
    }
  }
}

#' Get the add+prop error model function or estimate from ui
#'
#' @param ui environment of rxode2 ui
#' @param pred1 prediction line for endpint
#' @param fun boolean to return function
#' @return function or estimate of error line
#' @noRd
#' @author Matthew L. Fidler
posologyr_get_error_model_add_pow <- function(ui, pred1, fun=TRUE) {
  stop("add+pow not supported with posologyr", call.=FALSE)
}
#'  Get error model based on pred line
#'
#' @param ui rxode2 ui environment
#' @param pred1 pred1 line
#' @param fun boolean to return function or ini estimate
#' @return function for error type
#' @noRd
#' @author Matthew L. Fidler
#' @keywords internal
posologyr_get_error_model <- function(ui, pred1, fun=TRUE) {
  switch(as.character(pred1$errType),
         "add"=posologyr_get_error_model_add(ui, pred1, fun), # 1
         "prop"=posologyr_get_error_model_prop(ui, pred1, fun), # 2
         "pow"=posologyr_get_error_model_pow(ui, pred1, fun), # 3
         "add + prop"=posologyr_get_error_model_add_prop(ui, pred1, fun),# 4
         "add + pow"=posologyr_get_error_model_add_pow(ui, pred1, fun) # 5
         )
}

#' @export
rxUiGet.posologyr_error_model <- function(x, ...) {
  ui <- x[[1]]
  pred_df <- ui$predDf
  ret <- lapply(seq_along(pred_df$cond), function(c) {
    posologyr_get_error_model(ui, pred_df[c, ], fun=TRUE)
  })
  if (length(ret) == 1) return(ret[[1]])
  names(ret) <- pred_df$cond
  ret
}
attr(rxUiGet.posologyr_error_model, "desc") <- "posologyr error_model element"


#' @export
rxUiGet.posologyr_sigma <- function(x, ...) {
  ui <- x[[1]]
  pred_df <- ui$predDf
  ret <- lapply(seq_along(pred_df$cond), function(c) {
    posologyr_get_error_model(ui, pred_df[c, ], fun=FALSE)
  })
  if (length(ret) == 1) return(ret[[1]])
  names(ret) <- pred_df$cond
  ret
}

#' @export
rxUiGet.posologyr <- function(x, ...) {
  ui <- x[[1]]
  mod <- eval(rxUiGet.posologyr_ppk_model(x, ...))
  ret <- list(ppk_model=mod,
       error_model=rxUiGet.posologyr_error_model(x, ...),
       theta=ui$theta,
       omega=ui$omega,
       sigma=rxUiGet.posologyr_sigma(x, ...),
       covariates=ui$allCovs)
  if (length(ret$covariates) == 0L) ret$covariates <- NULL
  ret
}
attr(rxUiGet.posologyr, "desc") <- "posologyr model from ui"
