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
#' @return A named list with the following elements:
#'   \describe{
#'     \item{deviance_r}{Deviance for the reduced model.}
#'     \item{deviance_f}{Deviance for the full model.}
#'     \item{p_r}{Number of parameters in the reduced model.}
#'     \item{p_f}{Number of parameters in the full model.}
#'     \item{deviance_diff_test}{Chi-square statistic, degrees of freedom,
#'       and p-value for the deviance difference test.}
#'   }
#'
#' @examples
#' \dontrun{
#' # Using model objects
#' result <- deviance_diff_test(model_r = m_reduced, model_f = m_full)
#' result$deviance_diff_test
#'
#' # Using raw values from other software
#' result <- deviance_diff_test(deviance_r = 3689.3, deviance_f = 3674.1,
#'                              p_r = 10, p_f = 12)
#' result$deviance_diff_test
#' }
#' @importFrom stats pchisq
#' @importFrom glmmTMB fixef
deviance_diff_test <- function(model_r = NULL, model_f = NULL,
                               deviance_r = NULL, deviance_f = NULL,
                               p_r = NULL, p_f = NULL) {

  chisquare_symbol <- paste0("\u03C7", "\u00B2")

  # Guard against mixing model objects and raw deviance values
  if (length(deviance_r) == 1 & length(deviance_f) == 1 &
      (!is.null(model_r) | !is.null(model_f)))
    stop("Because you provided model objects as arguments, you do not need to ",
         "provide values for deviances.", call. = FALSE)

  # --- Raw deviance values supplied ---
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

  # --- Model objects supplied ---

  # Guard against mixing model classes
  if (!is.null(model_r) && !is.null(model_f) &&
      inherits(model_r, "gls") && inherits(model_f, "glmmTMB"))
    stop("Both model_r and model_f must be estimated using the same package ",
         "(i.e., both gls or both glmmTMB).", call. = FALSE)

  if (!is.null(model_r) && !is.null(model_f) &&
      inherits(model_r, "glmmTMB") && inherits(model_f, "gls"))
    stop("Both model_r and model_f must be estimated using the same package ",
         "(i.e., both gls or both glmmTMB).", call. = FALSE)

  # --- gls objects ---
  if (!is.null(model_r) && !is.null(model_f) &&
      inherits(model_r, "gls") && inherits(model_f, "gls")) {

    fixeff_difference <- model_f$dims$p - model_r$dims$p

    if (fixeff_difference != 0) {
      # Fixed effects differ — both must be ML
      if (eval(model_r$method) == "REML" || eval(model_f$method) == "REML")
        stop("Both model_r and model_f must be estimated using Maximum ",
             "Likelihood when the number of fixed effects differs across ",
             "the models.", call. = FALSE)

      deviance.reduced <- -2 * model_r$logLik
      deviance.full    <- -2 * model_f$logLik
      p.full           <- model_f$dims$p + nrow(model_f$apVar)
      p.reduced        <- model_r$dims$p + nrow(model_r$apVar)

    } else {
      # Fixed effects same — REML acceptable
      deviance.reduced <- -2 * model_r$logLik
      deviance.full    <- -2 * model_f$logLik
      p.full           <- model_f$dims$p + nrow(model_f$apVar)
      p.reduced        <- model_r$dims$p + nrow(model_r$apVar)
    }

    if (p.reduced >= p.full)
      stop("Number of parameters for full model must be greater than number ",
           "of parameters for reduced model.", call. = FALSE)

    chi <- deviance.reduced - deviance.full
    df  <- p.full - p.reduced
  }

  # --- glmmTMB objects ---
  if (!is.null(model_r) && !is.null(model_f) &&
      inherits(model_r, "glmmTMB") && inherits(model_f, "glmmTMB")) {

    fixeff_difference <- length(fixef(model_f)$cond) - length(fixef(model_r)$cond)

    if (fixeff_difference != 0) {
      # Fixed effects differ — both must be ML
      if (isTRUE(eval(model_r$call$REML)) || isTRUE(eval(model_f$call$REML)))
        stop("Both model_r and model_f must be estimated using Maximum ",
             "Likelihood when the number of fixed effects differs across ",
             "the models.", call. = FALSE)

      deviance.reduced <- -2 * (-(model_r$obj$fn()[[1]]))
      deviance.full    <- -2 * (-(model_f$obj$fn()[[1]]))
      p.full           <- length(model_f$obj$par)
      p.reduced        <- length(model_r$obj$par)

    } else {
      # Fixed effects same — REML acceptable; refit both with ML to get
      # comparable deviances before computing the test.
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

  text_results <- list(
    deviance_r        = paste0("Deviance for Reduced Model = ",
                                round(deviance.reduced, 3), "."),
    deviance_f        = paste0("Deviance for Full Model = ",
                                round(deviance.full, 3), "."),
    p_r               = paste0("Number of Parameters in Reduced Model = ",
                                p.reduced, "."),
    p_f               = paste0("Number of Parameters in Full Model = ",
                                p.full, "."),
    deviance_diff_test = paste0("Deviance Difference Test: ",
                                 chisquare_symbol, "(", df, ") = ",
                                 round(chi, 3), ", p = ",
                                 round(pvalue, 5), ".")
  )

  return(text_results)
}
