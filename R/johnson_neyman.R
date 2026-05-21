# =============================================================================
# EXPORTED: johnson_neyman
# =============================================================================

#' Johnson-Neyman Regions of Significance for Two-Way Interactions
#'
#' @description
#' Computes and plots Johnson-Neyman regions of significance for the effect
#' of a focal predictor across the observed range of a moderator variable.
#' Supports models fitted with \code{glmmTMB}, \code{lme}, or \code{gls}.
#' Three methods are available for degrees of freedom via \code{df_method}:
#'
#' \describe{
#'   \item{\code{"z"} (default)}{Large-sample z-test for all model types.}
#'   \item{\code{"t"}}{t-test using residual df from the model. Available
#'     for \code{gls} and \code{lme} only. For \code{lme}, residual df are
#'     used rather than per-term df. Not available for \code{glmmTMB};
#'     falls back to \code{"z"} with a warning.}
#'   \item{\code{"satterthwaite"}}{Satterthwaite df computed on a coarse
#'     grid and interpolated across the moderator range. Available for
#'     \code{glmmTMB} with \code{dispformula = ~0} only. Falls back to
#'     \code{"t"} for \code{gls} and \code{lme} with a warning.}
#' }
#'
#' @param model A fitted model object of class \code{glmmTMB}, \code{lme},
#'   or \code{gls} containing an interaction between \code{pred} and
#'   \code{modx}.
#' @param pred Character string. Name of the focal predictor.
#' @param modx Character string. Name of the moderator variable.
#' @param data A \code{data.frame} containing the variables used to fit
#'   \code{model}.
#' @param alpha Numeric. Significance level. Default is \code{0.05}.
#' @param n_mod Integer. Number of points in the moderator grid. Default
#'   is \code{1000}.
#' @param df_method Character string. One of \code{"z"} (default),
#'   \code{"t"}, or \code{"satterthwaite"}. See Description.
#' @param sw_grid Integer. Number of coarse grid points for Satterthwaite
#'   df computation when \code{df_method = "satterthwaite"}. Df values are
#'   interpolated to the full \code{n_mod}-point grid. Default is
#'   \code{20}.
#' @param eps Numeric. Step size for Jacobian in Satterthwaite computation.
#'   Default is \code{1e-3}.
#' @param n_cores Integer. Number of parallel workers for Satterthwaite
#'   computation. Default is all physical cores.
#' @param title Character string. Plot title. Auto-generated if \code{NULL}.
#' @param ylabel Character string. Y-axis label. Auto-generated if
#'   \code{NULL}.
#' @param verbose Logical. Print progress during Satterthwaite computation.
#'   Default is \code{FALSE}.
#'
#' @return Invisibly returns a \code{data.frame} with columns:
#'   \describe{
#'     \item{moderator}{Moderator values across the grid.}
#'     \item{simple_slope}{Simple slope of \code{pred}.}
#'     \item{se_slope}{Standard error of the simple slope.}
#'     \item{df}{Degrees of freedom.}
#'     \item{p_value}{Two-tailed p-value.}
#'   }
#'
#' @examples
#' \dontrun{
#' # glmmTMB with Satterthwaite df
#' johnson_neyman(
#'   model     = ind_moderation_socA,
#'   pred      = "c_PosBehavior_A",
#'   modx      = "c_Support_A",
#'   data      = pairwise_indisting,
#'   df_method = "satterthwaite",
#'   sw_grid   = 20
#' )
#'
#' # gls with residual df t-test
#' johnson_neyman(
#'   model     = ind_moderation_bdyad_gls,
#'   pred      = "c_PosBehavior_A",
#'   modx      = "c_Rellengthyrs",
#'   data      = pairwise_indisting,
#'   df_method = "t"
#' )
#'
#' # lme with residual df t-test
#' johnson_neyman(
#'   model     = ind_moderation_bdyad_lme,
#'   pred      = "c_PosBehavior_A",
#'   modx      = "c_Rellengthyrs",
#'   data      = pairwise_indisting,
#'   df_method = "t"
#' )
#' }
#'
#' @importFrom stats approx qt pt setNames
#' @importFrom graphics plot lines polygon abline mtext legend
#' @importFrom grDevices adjustcolor
#' @export
johnson_neyman <- function(model, pred, modx, data,
                           alpha     = .05,
                           n_mod     = 1000,
                           df_method = c("z", "t", "satterthwaite"),
                           sw_grid   = 20,
                           eps       = 1e-3,
                           n_cores   = parallel::detectCores(logical = FALSE),
                           title     = NULL,
                           ylabel    = NULL,
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
  cat("Using interaction term:", int_term, "\n")

  b_pred   <- coefs[pred]
  b_int    <- coefs[int_term]
  var_pred <- vcov_mat[pred, pred]
  var_int  <- vcov_mat[int_term, int_term]
  cov_pi   <- vcov_mat[pred, int_term]

  # -- moderator grid ---------------------------------------------------------
  mod_range    <- seq(min(data[[modx]], na.rm = TRUE),
                      max(data[[modx]], na.rm = TRUE),
                      length.out = n_mod)
  simple_slope <- b_pred + b_int * mod_range
  se_slope     <- sqrt(var_pred + 2 * mod_range * cov_pi + mod_range^2 * var_int)

  # -- resolve df method ------------------------------------------------------
  df_values <- rep(Inf, n_mod)
  df_note   <- "z-test"

  if(df_method == "t") {
    if(is_glmmTMB) {
      warning("df_method = 't' not available for glmmTMB; falling back to z-test.")
    } else {
      df_values <- rep(df_resid, n_mod)
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
      df_values <- rep(df_resid, n_mod)
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
        cat(sprintf(
          "Computing Satterthwaite df on %d-point grid (then interpolating)...\n",
          sw_grid))
        mod_grid <- seq(min(data[[modx]], na.rm = TRUE),
                        max(data[[modx]], na.rm = TRUE),
                        length.out = sw_grid)
        df_grid <- sapply(mod_grid, function(w) {
          if(verbose) cat("  df at moderator =", round(w, 3), "\n")
          L <- setNames(c(1, w), c(pred, int_term))
          .satterthwaite_contrast(model, L, eps = eps,
                                  n_cores = n_cores, verbose = FALSE)
        })
        df_values <- approx(mod_grid, df_grid, xout = mod_range)$y
        df_note   <- sprintf("Satterthwaite df interpolated from %d grid points",
                             sw_grid)
        cat("Done.\n")
      }
    }
  }

  # -- t/z statistics and p-values --------------------------------------------
  t_stat  <- simple_slope / se_slope
  p_value <- mapply(function(t, df) 2 * pt(abs(t), df = df, lower.tail = FALSE),
                    t_stat, df_values)
  crit    <- sapply(df_values, function(df) qt(1 - alpha / 2, df = df))

  # -- significance regions ---------------------------------------------------
  sig_idx     <- which(p_value <  alpha)
  nonsig_idx  <- which(p_value >= alpha)

  if(length(sig_idx) == 0) {
    cat("Effect of", pred, "is non-significant across entire range of", modx, "\n")
  } else if(length(sig_idx) == n_mod) {
    cat("Effect of", pred, "is significant across entire range of", modx, "\n")
  } else {
    transitions <- mod_range[sig_idx[c(TRUE, diff(sig_idx) > 1)]]
    transitions <- c(transitions,
                     mod_range[sig_idx[c(diff(sig_idx) > 1, TRUE)]])
    cat("Regions of significance for", pred, ":\n")
    for(k in seq(1, length(transitions), by = 2)) {
      cat(sprintf("  From %.3f to %.3f\n", transitions[k],
                  ifelse(k + 1 <= length(transitions),
                         transitions[k + 1], max(mod_range))))
    }
  }

  # -- plot -------------------------------------------------------------------
  ylim_range <- range(c(simple_slope + crit * se_slope,
                        simple_slope - crit * se_slope))

  plot(mod_range, simple_slope,
       type = "l", lwd = 2, col = "black",
       ylim = ylim_range,
       xlab = modx,
       ylab = ylabel %||% paste("Simple Slope of", pred),
       main = title  %||% paste("Johnson-Neyman:", pred))

  abline(h = 0, lty = 2, col = "grey50")

  if(length(nonsig_idx) > 0)
    polygon(c(mod_range[nonsig_idx], rev(mod_range[nonsig_idx])),
            c(simple_slope[nonsig_idx] + crit[nonsig_idx] * se_slope[nonsig_idx],
              rev(simple_slope[nonsig_idx] - crit[nonsig_idx] * se_slope[nonsig_idx])),
            col = adjustcolor("grey80", alpha.f = 0.3), border = NA)

  if(length(sig_idx) > 0)
    polygon(c(mod_range[sig_idx], rev(mod_range[sig_idx])),
            c(simple_slope[sig_idx] + crit[sig_idx] * se_slope[sig_idx],
              rev(simple_slope[sig_idx] - crit[sig_idx] * se_slope[sig_idx])),
            col = adjustcolor("grey50", alpha.f = 0.3), border = NA)

  lines(mod_range, simple_slope + crit * se_slope, lty = 2, col = "grey40")
  lines(mod_range, simple_slope - crit * se_slope, lty = 2, col = "grey40")
  lines(mod_range, simple_slope, lwd = 2, col = "black")

  legend("topleft", bty = "n",
         legend = c("p < alpha", "p >= alpha"),
         fill   = c(adjustcolor("grey50", alpha.f = 0.3),
                    adjustcolor("grey80", alpha.f = 0.3)),
         border = NA)

  mtext(df_note, side = 3, cex = 0.75, col = "grey40")

  invisible(data.frame(moderator    = mod_range,
                       simple_slope = simple_slope,
                       se_slope     = se_slope,
                       df           = df_values,
                       p_value      = p_value))
}
