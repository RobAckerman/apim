#' Deviance Difference Test for APIM Models
#'
#' @description
#' Computes and prints the chi-square statistic, degrees of freedom, and
#' p-value for a deviance difference test comparing a reduced and full model.
#' Accepts either fitted \code{gls} or \code{glmmTMB} model objects directly,
#' or raw deviance values and parameter counts for use with output from other
#' software.
#'
#' When comparing models that differ in their fixed effects, both models must
#' be estimated using Maximum Likelihood (ML). When comparing models that
#' differ only in their random effects, REML estimation is acceptable.
#'
#' @param model_r A fitted \code{gls} or \code{glmmTMB} reduced model object.
#'   Required if not supplying raw deviance values.
#' @param model_f A fitted \code{gls} or \code{glmmTMB} full model object.
#'   Required if not supplying raw deviance values.
#' @param deviance_r Numeric. Deviance for the reduced model. Only required
#'   when not supplying model objects (e.g. when using output from other
#'   software).
#' @param deviance_f Numeric. Deviance for the full model. Only required
#'   when not supplying model objects.
#' @param p_r Integer. Number of parameters in the reduced model. Only
#'   required when not supplying model objects.
#' @param p_f Integer. Number of parameters in the full model. Only required
#'   when not supplying model objects.
#'
#' @return Invisibly returns a \code{data.frame} with columns \code{deviance_r},
#'   \code{deviance_f}, \code{p_r}, \code{p_f}, \code{chi2}, \code{df},
#'   and \code{p}.
#'
#' @examples
#' \dontrun{
#' # Using model objects
#' deviance_diff_test(model_r = m_reduced, model_f = m_full)
#'
#' # Using raw values from other software
#' deviance_diff_test(deviance_r = 3689.3, deviance_f = 3674.1,
#'                    p_r = 10, p_f = 12)
#' }
#' @importFrom stats pchisq
#' @importFrom glmmTMB fixef
#' @export
deviance_diff_test <- function(model_r    = NULL,
                               model_f    = NULL,
                               deviance_r = NULL,
                               deviance_f = NULL,
                               p_r        = NULL,
                               p_f        = NULL) {

  chisquare_symbol <- paste0("\u03C7", "\u00B2")

  # -- helper: format p-value ------------------------------------------------
  fmt_p <- function(p) {
    if (p < .0001) "< .0001" else sprintf("= %.4f", p)
  }

  # -- guard against mixing model objects and raw values ---------------------
  if (length(deviance_r) == 1 & length(deviance_f) == 1 &
      (!is.null(model_r) | !is.null(model_f)))
    stop("Because you provided model objects as arguments, you do not need to ",
         "provide values for deviances.", call. = FALSE)

  # -- raw deviance values supplied ------------------------------------------
  if (length(deviance_r) == 1 & length(deviance_f) == 1) {
    deviance.reduced <- deviance_r
    deviance.full    <- deviance_f
    p.full           <- p_f
    p.reduced        <- p_r

    if (p.reduced >= p.full)
      stop("Number of parameters for full model must be greater than number ",
           "of parameters for reduced model.", call. = FALSE)

    chi <- deviance.reduced - deviance.full
    df  <- p.full - p.reduced
  }

  # -- guard against mixing model classes ------------------------------------
  if (!is.null(model_r) && !is.null(model_f) &&
      inherits(model_r, "gls") && inherits(model_f, "glmmTMB"))
    stop("Both model_r and model_f must be estimated using the same package ",
         "(i.e., both gls or both glmmTMB).", call. = FALSE)

  if (!is.null(model_r) && !is.null(model_f) &&
      inherits(model_r, "glmmTMB") && inherits(model_f, "gls"))
    stop("Both model_r and model_f must be estimated using the same package ",
         "(i.e., both gls or both glmmTMB).", call. = FALSE)

  # -- gls objects -----------------------------------------------------------
  if (!is.null(model_r) && !is.null(model_f) &&
      inherits(model_r, "gls") && inherits(model_f, "gls")) {

    fixeff_difference <- model_f$dims$p - model_r$dims$p

    if (fixeff_difference != 0) {
      if (eval(model_r$method) == "REML" || eval(model_f$method) == "REML")
        stop("Both model_r and model_f must be estimated using Maximum ",
             "Likelihood when the number of fixed effects differs across ",
             "the models.", call. = FALSE)
    }

    deviance.reduced <- -2 * model_r$logLik
    deviance.full    <- -2 * model_f$logLik
    p.full           <- model_f$dims$p + nrow(model_f$apVar)
    p.reduced        <- model_r$dims$p + nrow(model_r$apVar)

    if (p.reduced >= p.full)
      stop("Number of parameters for full model must be greater than number ",
           "of parameters for reduced model.", call. = FALSE)

    chi <- deviance.reduced - deviance.full
    df  <- p.full - p.reduced
  }

  # -- glmmTMB objects -------------------------------------------------------
  if (!is.null(model_r) && !is.null(model_f) &&
      inherits(model_r, "glmmTMB") && inherits(model_f, "glmmTMB")) {

    fixeff_difference <- length(fixef(model_f)$cond) -
      length(fixef(model_r)$cond)

    if (fixeff_difference != 0) {
      if (isTRUE(eval(model_r$call$REML)) || isTRUE(eval(model_f$call$REML)))
        stop("Both model_r and model_f must be estimated using Maximum ",
             "Likelihood when the number of fixed effects differs across ",
             "the models.", call. = FALSE)

      deviance.reduced <- -2 * (-(model_r$obj$fn()[[1]]))
      deviance.full    <- -2 * (-(model_f$obj$fn()[[1]]))
      p.full           <- length(model_f$obj$par)
      p.reduced        <- length(model_r$obj$par)

    } else {
      model_r_ml <- suppressWarnings(update(model_r, REML = FALSE))
      model_f_ml <- suppressWarnings(update(model_f, REML = FALSE))

      deviance.reduced <- -2 * (-(model_r_ml$obj$fn()[[1]]))
      deviance.full    <- -2 * (-(model_f_ml$obj$fn()[[1]]))
      p.full           <- length(model_f_ml$obj$par)
      p.reduced        <- length(model_r_ml$obj$par)
    }

    if (p.reduced >= p.full)
      stop("Number of parameters for full model must be greater than number ",
           "of parameters for reduced model.", call. = FALSE)

    chi <- deviance.reduced - deviance.full
    df  <- p.full - p.reduced
  }

  pvalue <- pchisq(q = chi, df = df, lower.tail = FALSE)

  # -- print -----------------------------------------------------------------
  cat("\n=== Deviance Difference Test ===\n\n")
  cat(sprintf("  Deviance (reduced) : %.3f  (df = %d)\n",
              deviance.reduced, p.reduced))
  cat(sprintf("  Deviance (full)    : %.3f  (df = %d)\n",
              deviance.full, p.full))
  cat(sprintf("  %s(%d) = %.3f, p %s\n\n",
              chisquare_symbol, df, chi, fmt_p(pvalue)))

  invisible(data.frame(
    deviance_r = round(deviance.reduced, 3),
    deviance_f = round(deviance.full,    3),
    p_r        = p.reduced,
    p_f        = p.full,
    chi2       = round(chi,    3),
    df         = df,
    p          = round(pvalue, 5),
    row.names  = NULL
  ))
}
