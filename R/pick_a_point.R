# =============================================================================
# EXPORTED: pick_a_point
# =============================================================================

#' Pick-a-Point Analysis for Two-Way Interactions
#'
#' @description
#' Computes simple slopes of a focal predictor at specified values of a
#' moderator variable (the pick-a-point or spotlight approach). Supports
#' models fitted with \code{glmmTMB}, \code{lme}, or \code{gls}. Three
#' methods are available for degrees of freedom via \code{df_method}:
#'
#' \describe{
#'   \item{\code{"z"} (default)}{Large-sample z-test for all model types.}
#'   \item{\code{"t"}}{t-test using residual df. Available for \code{gls}
#'     and \code{lme} only. For \code{lme}, residual df are used for all
#'     simple slopes rather than per-term df. Not available for
#'     \code{glmmTMB}; falls back to \code{"z"} with a warning.}
#'   \item{\code{"satterthwaite"}}{Satterthwaite df computed separately for
#'     each simple slope via numerical differentiation of the contrast
#'     variance. Available for \code{glmmTMB} with \code{dispformula = ~0}
#'     only. Falls back to \code{"t"} for \code{gls} and \code{lme} with
#'     a warning.}
#' }
#'
#' @param model A fitted model object of class \code{glmmTMB}, \code{lme},
#'   or \code{gls} containing an interaction between \code{pred} and
#'   \code{modx}.
#' @param pred Character string. Name of the focal predictor.
#' @param modx Character string. Name of the moderator variable.
#' @param data A \code{data.frame} containing the variables used to fit
#'   \code{model}.
#' @param points Either a character vector of keywords (\code{"mean-sd"},
#'   \code{"mean"}, \code{"mean+sd"}) or a numeric vector of specific
#'   moderator values. Default is \code{c("mean-sd", "mean", "mean+sd")}.
#' @param alpha Numeric. Significance level. Default is \code{0.05}.
#' @param df_method Character string. One of \code{"z"} (default),
#'   \code{"t"}, or \code{"satterthwaite"}. See Description.
#' @param eps Numeric. Step size for Jacobian in Satterthwaite computation.
#'   Default is \code{1e-3}.
#' @param n_cores Integer. Number of parallel workers for Satterthwaite
#'   computation. Default is all physical cores.
#' @param verbose Logical. Print progress during Satterthwaite computation.
#'   Default is \code{FALSE}.
#'
#' @return Invisibly returns a \code{data.frame} with columns:
#'   \describe{
#'     \item{Label}{Description of the moderator value.}
#'     \item{Moderator_Value}{Numeric value of the moderator.}
#'     \item{Simple_Slope}{Estimated simple slope of \code{pred}.}
#'     \item{SE}{Standard error of the simple slope.}
#'     \item{df}{Degrees of freedom.}
#'     \item{t}{t- or z-statistic.}
#'     \item{p}{Two-tailed p-value.}
#'     \item{CI_lower}{Lower confidence interval bound.}
#'     \item{CI_upper}{Upper confidence interval bound.}
#'     \item{Sig}{Significance indicator (\code{*} if \code{p < alpha}).}
#'   }
#'
#' @examples
#' \dontrun{
#' # glmmTMB with Satterthwaite df
#' pick_a_point(
#'   model     = ind_moderation_socA,
#'   pred      = "c_PosBehavior_A",
#'   modx      = "c_Support_A",
#'   data      = pairwise_indisting,
#'   df_method = "satterthwaite"
#' )
#'
#' # gls with residual df t-test
#' pick_a_point(
#'   model     = ind_moderation_bdyad_gls,
#'   pred      = "c_PosBehavior_A",
#'   modx      = "c_Rellengthyrs",
#'   data      = pairwise_indisting,
#'   df_method = "t"
#' )
#'
#' # lme with residual df t-test
#' pick_a_point(
#'   model     = ind_moderation_bdyad_lme,
#'   pred      = "c_PosBehavior_A",
#'   modx      = "c_Rellengthyrs",
#'   data      = pairwise_indisting,
#'   df_method = "t"
#' )
#' }
#'
#' @importFrom stats sd qt pt setNames
#' @export
pick_a_point <- function(model, pred, modx, data,
                         points    = c("mean-sd", "mean", "mean+sd"),
                         alpha     = .05,
                         df_method = c("z", "t", "satterthwaite"),
                         eps       = 1e-3,
                         n_cores   = parallel::detectCores(logical = FALSE),
                         verbose   = FALSE) {

  df_method  <- match.arg(df_method)
  is_glmmTMB <- inherits(model, "glmmTMB")
  is_nlme    <- inherits(model, "gls") || inherits(model, "lme")

  # -- extract model components -----------------------------------------------
  ce       <- .extract_coefs(model)
  coefs    <- ce$coefs
  vcov_mat <- ce$vcov_mat
  df_resid <- ce$df_resid

  int_term <- .find_int_term(names(coefs), pred, modx)
  cat("Using interaction term:", int_term, "\n\n")

  b_pred   <- coefs[pred]
  b_int    <- coefs[int_term]
  var_pred <- vcov_mat[pred, pred]
  var_int  <- vcov_mat[int_term, int_term]
  cov_pi   <- vcov_mat[pred, int_term]

  # -- moderator values and labels --------------------------------------------
  mod_mean <- mean(data[[modx]], na.rm = TRUE)
  mod_sd   <- sd(data[[modx]],   na.rm = TRUE)

  if(is.numeric(points)) {
    mod_values <- points
    labels     <- paste("Value =", round(points, 3))
  } else {
    mod_values <- sapply(points, function(p) {
      switch(p,
             "mean-sd"  = mod_mean - mod_sd,
             "mean"     = mod_mean,
             "mean+sd"  = mod_mean + mod_sd,
             stop("Unknown point '", p,
                  "'. Use 'mean-sd', 'mean', 'mean+sd', or numeric values."))
    })
    labels <- sapply(points, function(p) {
      switch(p,
             "mean-sd"  = "Low (M - 1SD)",
             "mean"     = "Average (M)",
             "mean+sd"  = "High (M + 1SD)")
    })
  }

  # -- resolve df method ------------------------------------------------------
  use_sw   <- FALSE
  use_t    <- FALSE
  df_fixed <- Inf
  df_note  <- "z-test"

  if(df_method == "t") {
    if(is_glmmTMB) {
      warning("df_method = 't' not available for glmmTMB; falling back to z-test.")
    } else {
      use_t    <- TRUE
      df_fixed <- df_resid
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
      use_t    <- TRUE
      df_fixed <- df_resid
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
        use_sw  <- TRUE
        df_note <- "Satterthwaite df"
        cat("Satterthwaite df will be computed for each simple slope.\n\n")
      }
    }
  }

  # -- compute results --------------------------------------------------------
  results <- lapply(seq_along(mod_values), function(i) {
    w     <- mod_values[i]
    slope <- b_pred + b_int * w
    se    <- sqrt(var_pred + 2 * w * cov_pi + w^2 * var_int)

    if(use_sw) {
      if(verbose) cat("Computing Satterthwaite df for:", labels[i], "\n")
      L      <- setNames(c(1, w), c(pred, int_term))
      df_val <- .satterthwaite_contrast(model, L, eps = eps,
                                        n_cores = n_cores, verbose = verbose)
    } else {
      df_val <- df_fixed
    }

    t_val    <- slope / se
    p_val    <- 2 * pt(abs(t_val), df = df_val, lower.tail = FALSE)
    crit     <- qt(1 - alpha / 2, df = df_val)
    ci_lower <- slope - crit * se
    ci_upper <- slope + crit * se
    sig      <- ifelse(p_val < alpha, "*", "")

    data.frame(Label           = labels[i],
               Moderator_Value = round(w, 3),
               Simple_Slope    = round(slope, 4),
               SE              = round(se, 4),
               df              = ifelse(is.infinite(df_val), Inf,
                                        round(df_val, 2)),
               t               = round(t_val, 4),
               p               = round(p_val, 4),
               CI_lower        = round(ci_lower, 4),
               CI_upper        = round(ci_upper, 4),
               Sig             = sig)
  })

  results_df <- do.call(rbind, results)

  # -- print ------------------------------------------------------------------
  cat("Pick-a-Point: Simple slopes of", pred, "at values of", modx, "\n")
  cat(rep("=", 65), "\n", sep = "")
  print(results_df, row.names = FALSE)
  cat("\n* p <", alpha, "\n")
  cat("Note:", df_note, "\n")

  invisible(results_df)
}
