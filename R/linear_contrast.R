# =============================================================================
# INTERNAL: extract coefficients, vcov, and residual df
# =============================================================================

#' @keywords internal
.extract_coefs <- function(model) {
  if(inherits(model, "glmmTMB")) {
    list(coefs    = fixef(model)$cond,
         vcov_mat = as.matrix(vcov(model)$cond),
         df_resid = NULL)
  } else if(inherits(model, "lme")) {
    list(coefs    = fixef(model),
         vcov_mat = as.matrix(vcov(model)),
         df_resid = min(model$fixDF$X))  # most conservative df across all terms
  } else if(inherits(model, "gls")) {
    dims <- summary(model)$dims
    list(coefs    = coef(model),
         vcov_mat = as.matrix(vcov(model)),
         df_resid = dims$N - dims$p)
  } else {
    stop("model must be a glmmTMB, lme, or gls object.")
  }
}


#' @keywords internal
.find_int_term <- function(coef_names, pred, modx) {
  hits <- coef_names[grepl(pred,  coef_names, fixed = TRUE) &
                       grepl(modx, coef_names, fixed = TRUE)]
  if(length(hits) == 0)
    stop("No interaction term found for '", pred, "' and '", modx, "'.")
  hits[1]
}


# =============================================================================
# EXPORTED: linear_contrast
# =============================================================================

#' Test a Linear Contrast of Fixed Effects
#'
#' @description
#' Tests an arbitrary linear combination of fixed effects \eqn{L'\beta} for
#' models fitted with \code{glmmTMB}, \code{gls}, or \code{lme}. Three
#' methods are available for computing degrees of freedom, controlled by the
#' \code{df_method} argument:
#'
#' \describe{
#'   \item{\code{"z"} (default)}{Large-sample z-test for all model types.}
#'   \item{\code{"t"}}{t-test using residual df from the model. Available for
#'     \code{gls} and \code{lme} only. For \code{lme}, residual df are used
#'     for all terms rather than per-term df. Not available for
#'     \code{glmmTMB}; falls back to \code{"z"} with a warning.}
#'   \item{\code{"satterthwaite"}}{Satterthwaite df via numerical
#'     differentiation of the contrast variance. Available for \code{glmmTMB}
#'     with \code{dispformula = ~0} only. Falls back to \code{"t"} for
#'     \code{gls} and \code{lme} with a warning.}
#' }
#'
#' @param model A fitted model object of class \code{glmmTMB}, \code{gls},
#'   or \code{lme}.
#' @param L A named numeric vector of contrast weights. Names must match
#'   coefficient names in the model. Coefficients not named in \code{L} are
#'   assumed to have weight zero.
#' @param label A character string label for the contrast. Default is
#'   \code{"Contrast"}.
#' @param alpha Numeric. Significance level for confidence intervals.
#'   Default is \code{0.05}.
#' @param df_method Character string specifying the method for degrees of
#'   freedom. One of \code{"z"} (default), \code{"t"}, or
#'   \code{"satterthwaite"}. See Details.
#' @param eps Numeric. Step size for Jacobian in Satterthwaite computation.
#'   Only used when \code{df_method = "satterthwaite"} and model is
#'   \code{glmmTMB} with \code{dispformula = ~0}. Default is \code{1e-3}.
#' @param n_cores Integer. Number of parallel workers for Satterthwaite
#'   computation. Default is all physical cores.
#' @param verbose Logical. Print progress during Satterthwaite computation.
#'   Default is \code{FALSE}.
#'
#' @return Invisibly returns a \code{data.frame} with columns:
#'   \describe{
#'     \item{label}{Contrast label.}
#'     \item{estimate}{Estimated value of \eqn{L'\beta}.}
#'     \item{se}{Standard error.}
#'     \item{df}{Degrees of freedom.}
#'     \item{t}{t- or z-statistic.}
#'     \item{p}{Two-tailed p-value.}
#'     \item{ci_lower}{Lower confidence interval bound.}
#'     \item{ci_upper}{Upper confidence interval bound.}
#'   }
#'
#' @examples
#' \dontrun{
#' # simple slope of c_PosBehavior_A at low support (w = -1.31)
#' # glmmTMB with Satterthwaite df
#' linear_contrast(
#'   model     = ind_moderation_socA,
#'   L         = c(c_PosBehavior_A = 1,
#'                 "c_PosBehavior_A:c_Support_A" = -1.31),
#'   label     = "Actor simple slope at Low Support",
#'   df_method = "satterthwaite"
#' )
#'
#' # gls with residual df t-test
#' linear_contrast(
#'   model     = ind_moderation_bdyad_gls,
#'   L         = c(c_PosBehavior_A = 1,
#'                 "c_PosBehavior_A:c_Rellengthyrs" = -6.04),
#'   label     = "Actor simple slope at Low Rellengthyrs",
#'   df_method = "t"
#' )
#'
#' # lme with residual df t-test
#' linear_contrast(
#'   model     = ind_moderation_bdyad_lme,
#'   L         = c(c_PosBehavior_A = 1,
#'                 "c_PosBehavior_A:c_Rellengthyrs" = -6.04),
#'   label     = "Actor simple slope at Low Rellengthyrs",
#'   df_method = "t"
#' )
#' }
#'
#' @importFrom stats qt pt setNames
#' @export
linear_contrast <- function(model, L, label = "Contrast", alpha = .05,
                            df_method = c("z", "t", "satterthwaite"),
                            eps       = 1e-3,
                            n_cores   = parallel::detectCores(logical = FALSE),
                            verbose   = FALSE) {

  df_method   <- match.arg(df_method)
  is_glmmTMB  <- inherits(model, "glmmTMB")
  is_nlme     <- inherits(model, "gls") || inherits(model, "lme")

  if(!is_glmmTMB && !is_nlme)
    stop("model must be a glmmTMB, gls, or lme object.")

  # -- extract model components -----------------------------------------------
  ce       <- .extract_coefs(model)
  coefs    <- ce$coefs
  vcov_mat <- ce$vcov_mat
  df_resid <- ce$df_resid

  # -- align L to full coefficient vector -------------------------------------
  L_full <- setNames(rep(0, length(coefs)), names(coefs))
  for(nm in names(L)) {
    if(!nm %in% names(coefs))
      stop("Contrast term '", nm, "' not found in model coefficients.")
    L_full[nm] <- L[nm]
  }
  L_full <- as.numeric(L_full)

  # -- estimate and SE --------------------------------------------------------
  estimate <- sum(L_full * coefs)
  se       <- sqrt(as.numeric(t(L_full) %*% vcov_mat %*% L_full))

  # -- degrees of freedom -----------------------------------------------------
  df      <- Inf
  df_note <- "z-test"

  if(df_method == "t") {
    if(is_glmmTMB) {
      warning("df_method = 't' not available for glmmTMB; falling back to z-test.")
    } else {
      df      <- df_resid
      df_note <- if(inherits(model, "lme")) {
        df_range <- range(model$fixDF$X)
        sprintf("t-test with residual df = %d (minimum df across terms used; df range was %d to %d)",
                df_resid, df_range[1], df_range[2])
      } else {
        sprintf("t-test with residual df = %d", df_resid)
      }
    }
  }

  if(df_method == "satterthwaite") {
    if(is_nlme) {
      warning("df_method = 'satterthwaite' not available for gls/lme; ",
              "falling back to t-test with residual df.")
      df      <- df_resid
      df_note <- if(inherits(model, "lme")) {
        df_range <- range(model$fixDF$X)
        sprintf("t-test with residual df = %d (Satterthwaite not available; minimum df across terms used; df range was %d to %d)",
                df_resid, df_range[1], df_range[2])
      } else {
        sprintf("t-test with residual df = %d (Satterthwaite not available)", df_resid)
      }
    } else {
      has_resid_var <- !identical(deparse(model$call$dispformula), "~0")
      if(has_resid_var) {
        warning("Satterthwaite not available for glmmTMB with residual variance; ",
                "falling back to z-test.")
      } else {
        if(verbose) message("Computing Satterthwaite df for contrast...")
        df      <- .satterthwaite_contrast(model, L[L != 0], eps = eps,
                                           n_cores = n_cores, verbose = verbose)
        df_note <- "Satterthwaite df"
      }
    }
  }

  # -- test statistic and p-value ---------------------------------------------
  t_val    <- estimate / se
  p_val    <- 2 * pt(abs(t_val), df = df, lower.tail = FALSE)
  crit     <- qt(1 - alpha / 2, df = df)
  ci_lower <- estimate - crit * se
  ci_upper <- estimate + crit * se

  # -- print ------------------------------------------------------------------
  cat(sprintf("\n=== %s ===\n", label))
  cat(sprintf("Estimate : %8.4f\n", estimate))
  cat(sprintf("SE       : %8.4f\n", se))
  cat(sprintf("df       : %s\n",
              ifelse(is.infinite(df), sprintf("Inf (%s)", df_note),
                     sprintf("%.3f (%s)", df, df_note))))
  cat(sprintf("t        : %8.4f\n", t_val))
  cat(sprintf("p        : %8.4f\n", p_val))
  cat(sprintf("%d%% CI  : [%.4f, %.4f]\n",
              round((1 - alpha) * 100), ci_lower, ci_upper))

  invisible(data.frame(label, estimate, se, df, t = t_val, p = p_val,
                       ci_lower, ci_upper, row.names = NULL))
}
