# =============================================================================
# EXPORTED: run_dyadic_rsa
# =============================================================================

#' Fit a Dyadic Response Surface Analysis Model
#'
#' @description
#' Fits a dyadic response surface analysis (RSA) model using \code{glmmTMB},
#' computes RSA parameters a1--a5 for the overall surface and optionally for
#' subgroups defined by a moderator variable, and tests moderation of the
#' response surface. Validated against the dyadic RSA framework of Schönbrodt,
#' Humberg, & Nestler (2018).
#'
#' @details
#' The function automatically constructs the five polynomial terms required for
#' RSA (actor², partner², actor×partner) and builds the fixed-effects formula.
#' When a moderator is supplied, all five RSA terms are interacted with the
#' moderator, and surfaces are estimated at each value supplied in
#' \code{mod_values} as well as at the grand mean (moderator = 0), which under
#' effect coding equals the unweighted average of the group surfaces.
#'
#' RSA parameters are estimated as linear contrasts of the fixed effects using
#' \code{\link{linear_contrast}}:
#' \describe{
#'   \item{a1}{Slope along the line of congruence (actor + partner).}
#'   \item{a2}{Curvature along the line of congruence (actor² + actor×partner + partner²).}
#'   \item{a3}{Slope along the line of incongruence (actor - partner).}
#'   \item{a4}{Curvature along the line of incongruence (actor² - actor×partner + partner²).}
#'   \item{a5}{Rotation of the principal axis (actor² - partner²).}
#' }
#'
#' Results have been verified to be numerically identical (differences < 1e-5)
#' to the lavaan-based dyadic RSA implementation of Schönbrodt et al. (2018)
#' when both are estimated with maximum likelihood.
#'
#' @param data A data frame.
#' @param outcome Character string. Name of the outcome variable.
#' @param actor Character string. Name of the (centered) actor predictor.
#' @param partner Character string. Name of the (centered) partner predictor.
#' @param re_term Character string. The random effects term exactly as it
#'   should appear in the model formula, e.g.
#'   \code{"cs(0 + as.factor(Gender) | DyadID)"}. Passed verbatim.
#' @param moderator Character string or \code{NULL}. Name of the moderator
#'   variable (e.g. a distinguishing variable such as gender). When supplied,
#'   all five RSA terms are interacted with the moderator. Default \code{NULL}.
#' @param mod_values Named list or \code{NULL}. Values of the moderator at
#'   which to estimate group-specific surfaces, e.g.
#'   \code{list(Men = 1, Women = -1)}. Required when \code{moderator} is
#'   supplied and moderation tests are desired. Exactly two values triggers
#'   moderation tests. Default \code{NULL}.
#' @param alpha Numeric. Significance level for confidence intervals.
#'   Default \code{0.05}.
#' @param df_method Character string. Degrees-of-freedom method passed to
#'   \code{\link{linear_contrast}}. One of \code{"z"} (default),
#'   \code{"t"}, or \code{"satterthwaite"}.
#' @param reml Logical. Whether to use REML estimation. Default \code{TRUE},
#'   which is preferred for variance component estimation. Set to \code{FALSE}
#'   to match maximum likelihood estimators (e.g. for validation against
#'   lavaan or for likelihood ratio tests).
#'
#' @return Invisibly returns a named list with components:
#'   \describe{
#'     \item{model}{The fitted \code{glmmTMB} model object.}
#'     \item{results}{A named list of \code{data.frame}s with RSA contrast
#'       results. Always includes \code{$overall}. When \code{mod_values} is
#'       supplied, also includes one element per group and \code{$moderation}.}
#'     \item{coefs}{Named numeric vector of fixed-effect estimates.}
#'     \item{nms}{Character vector of fixed-effect names.}
#'   }
#'
#' @references
#' Schönbrodt, F. D., Humberg, S., & Nestler, S. (2018). Testing similarity
#' effects with dyadic response surface analysis. \emph{European Journal of
#' Personality}, 32(6), 627--641. \doi{10.1002/per.2169}
#'
#' @examples
#' \dontrun{
#' # Without moderator — single overall surface
#' rsa_out <- run_dyadic_rsa(
#'   data      = pairwise_data,
#'   outcome   = "RelSat_A",
#'   actor     = "c_PosEmotion_A",
#'   partner   = "c_PosEmotion_P",
#'   re_term   = "cs(0 + as.factor(ECGender_A) | DyadID)",
#'   df_method = "z"
#' )
#'
#' # With gender moderator — overall, group, and moderation surfaces
#' rsa_out <- run_dyadic_rsa(
#'   data       = pairwise_data,
#'   outcome    = "RelSat_A",
#'   actor      = "c_PosEmotion_A",
#'   partner    = "c_PosEmotion_P",
#'   re_term    = "cs(0 + as.factor(ECGender_A) | DyadID)",
#'   moderator  = "ECGender_A",
#'   mod_values = list(Men = 1, Women = -1),
#'   df_method  = "z"
#' )
#'
#' # Access results
#' rsa_out$model                 # glmmTMB model object
#' rsa_out$results$overall       # overall surface
#' rsa_out$results$Men           # men's surface
#' rsa_out$results$Women         # women's surface
#' rsa_out$results$moderation    # moderation tests
#' }
#'
#' @importFrom glmmTMB glmmTMB fixef
#' @importFrom stats as.formula update
#' @export
run_dyadic_rsa <- function(
    data,
    outcome,
    actor,
    partner,
    re_term,
    moderator  = NULL,
    mod_values = NULL,
    alpha      = 0.05,
    df_method  = "z",
    reml       = TRUE
) {

  # -- 1. Build polynomial terms -----------------------------------------------
  a2_nm <- paste0(actor,   "2")
  p2_nm <- paste0(partner, "2")
  ap_nm <- paste0(actor,   "_x_", partner)

  data[[a2_nm]] <- data[[actor]]^2
  data[[p2_nm]] <- data[[partner]]^2
  data[[ap_nm]] <- data[[actor]] * data[[partner]]

  # -- 2. Build formula --------------------------------------------------------
  rsa_terms <- paste(actor, partner, a2_nm, ap_nm, p2_nm, sep = " + ")

  if (is.null(moderator)) {
    fixed <- paste0(outcome, " ~ ", rsa_terms)
  } else {
    mod_interactions <- paste(
      paste0(actor,   ":", moderator),
      paste0(partner, ":", moderator),
      paste0(a2_nm,   ":", moderator),
      paste0(ap_nm,   ":", moderator),
      paste0(p2_nm,   ":", moderator),
      sep = " + "
    )
    fixed <- paste0(outcome, " ~ ", rsa_terms, " + ", moderator,
                    " + ", mod_interactions)
  }

  full_formula <- as.formula(paste(fixed, "+", re_term))
  cat("=== Model formula ===\n")
  print(full_formula)

  # -- 3. Fit model ------------------------------------------------------------
  cat("\nFitting model...\n")
  model <- glmmTMB(
    full_formula,
    dispformula = ~0,
    REML        = reml,
    data        = data
  )
  cat(sprintf("Done. (REML = %s)\n\n", reml))

  # -- 4. Push to global environment if Satterthwaite -------------------------
  # apim's Satterthwaite refits the model in parallel workers which need
  # to find the data and model in .GlobalEnv — fails inside a closure
  # In run_dyadic_rsa.R, update the Satterthwaite block (step 4):

  if (df_method == "satterthwaite") {
    # Force REML to a literal in the model call so workers don't need 'reml' symbol
    model$call$REML <- reml

    assign(".rsa_data_tmp", data,  envir = .GlobalEnv)
    model <- update(model, data = .rsa_data_tmp)
    assign(".rsa_model_tmp", model, envir = .GlobalEnv)

    on.exit({
      rm(list = c(".rsa_model_tmp", ".rsa_data_tmp"),
         envir    = .GlobalEnv,
         inherits = FALSE)
    }, add = TRUE)
  }

  # -- 5. Extract coefficients -------------------------------------------------
  coefs <- fixef(model)$cond
  nms   <- names(coefs)

  # -- 6. Helper: named contrast vector ----------------------------------------
  make_contrast <- function(weights_named) {
    v <- setNames(rep(0, length(coefs)), nms)
    for (nm in names(weights_named)) {
      if (!nm %in% nms) stop("Coefficient not found: ", nm)
      v[nm] <- weights_named[[nm]]
    }
    v
  }

  # -- 7. Helper: RSA contrasts for a moderator value --------------------------
  get_rsa_contrasts <- function(mod_val = NULL) {
    if (is.null(mod_val)) {
      c_a1 <- make_contrast(setNames(c(1,  1),     c(actor,  partner)))
      c_a2 <- make_contrast(setNames(c(1,  1,  1), c(a2_nm,  ap_nm,  p2_nm)))
      c_a3 <- make_contrast(setNames(c(1, -1),     c(actor,  partner)))
      c_a4 <- make_contrast(setNames(c(1, -1,  1), c(a2_nm,  ap_nm,  p2_nm)))
      c_a5 <- make_contrast(setNames(c(1, -1),     c(a2_nm,  p2_nm)))
    } else {
      a_int  <- paste0(actor,   ":", moderator)
      p_int  <- paste0(partner, ":", moderator)
      a2_int <- paste0(a2_nm,   ":", moderator)
      ap_int <- paste0(ap_nm,   ":", moderator)
      p2_int <- paste0(p2_nm,   ":", moderator)

      c_a1 <- make_contrast(setNames(
        c(1,  1,  mod_val,  mod_val),
        c(actor, partner, a_int, p_int)
      ))
      c_a2 <- make_contrast(setNames(
        c(1,  1,  1,  mod_val,  mod_val,  mod_val),
        c(a2_nm, ap_nm, p2_nm, a2_int, ap_int, p2_int)
      ))
      c_a3 <- make_contrast(setNames(
        c(1, -1,  mod_val, -mod_val),
        c(actor, partner, a_int, p_int)
      ))
      c_a4 <- make_contrast(setNames(
        c(1, -1,  1,  mod_val, -mod_val,  mod_val),
        c(a2_nm, ap_nm, p2_nm, a2_int, ap_int, p2_int)
      ))
      c_a5 <- make_contrast(setNames(
        c(1, -1,  mod_val, -mod_val),
        c(a2_nm, p2_nm, a2_int, p2_int)
      ))
    }
    list(a1 = c_a1, a2 = c_a2, a3 = c_a3, a4 = c_a4, a5 = c_a5)
  }

  # -- 8. Helper: run all 5 contrasts and print as horizontal table ----------
  run_rsa_contrasts <- function(contrasts, label_prefix) {

    model_for_contrast <- if (df_method == "satterthwaite") {
      .GlobalEnv$.rsa_model_tmp
    } else {
      model
    }

    # Run contrasts silently — suppress individual vertical printing
    results <- lapply(c("a1", "a2", "a3", "a4", "a5"), function(p) {
      out <- NULL
      capture.output(
        out <- linear_contrast(
          model_for_contrast,
          contrasts[[p]],
          label     = paste0(p, " \u2014 ", label_prefix),
          alpha     = alpha,
          df_method = df_method
        )
      )
      out
    })

    df <- do.call(rbind, results)

    # Significance stars
    stars <- ifelse(df$p < .001, "***",
                    ifelse(df$p < .01,  "**",
                           ifelse(df$p < .05,  "*",
                                  ifelse(df$p < .10,  ".",  ""))))

    # Stat label
    stat_label <- ifelse(df_method == "z", "z", "t")

    # Print horizontal table
    cat(sprintf("\n  %-4s  %9s  %9s  %9s  %9s  %5s  %10s  %10s\n",
                "par", "Estimate", "SE", stat_label, "p", "sig", "ci.lower", "ci.upper"))
    cat(" ", strrep("-", 80), "\n")

    for (i in 1:5) {
      cat(sprintf("  %-4s  %9.4f  %9.4f  %9.4f  %9.4f  %5s  %10.4f  %10.4f\n",
                  c("a1", "a2", "a3", "a4", "a5")[i],
                  df$estimate[i],
                  df$se[i],
                  df$t[i],
                  df$p[i],
                  stars[i],
                  df$ci_lower[i],
                  df$ci_upper[i]
      ))
    }
    cat("\n")

    invisible(df)
  }

  # -- 9. Run RSA --------------------------------------------------------------
  all_results <- list()

  if (is.null(moderator)) {

    cat("\n========================================\n")
    cat("RSA \u2014 Overall surface\n")
    cat("========================================\n")
    all_results$overall <- run_rsa_contrasts(
      get_rsa_contrasts(mod_val = NULL), "Overall"
    )

  } else {

    cat("\n========================================\n")
    cat("RSA \u2014 Overall surface (moderator = 0)\n")
    cat("========================================\n")
    all_results$overall <- run_rsa_contrasts(
      get_rsa_contrasts(mod_val = 0), "Overall"
    )

    if (!is.null(mod_values)) {
      for (grp in names(mod_values)) {
        val <- mod_values[[grp]]
        cat("\n========================================\n")
        cat("RSA \u2014", grp, "(moderator =", val, ")\n")
        cat("========================================\n")
        all_results[[grp]] <- run_rsa_contrasts(
          get_rsa_contrasts(mod_val = val), grp
        )
      }
    }

    if (!is.null(mod_values) && length(mod_values) == 2) {
      vals         <- unlist(mod_values)
      scale_factor <- vals[1] - vals[2]

      cat("\n========================================\n")
      cat("Moderation tests (surface difference)\n")
      cat("========================================\n")

      a_int  <- paste0(actor,   ":", moderator)
      p_int  <- paste0(partner, ":", moderator)
      a2_int <- paste0(a2_nm,   ":", moderator)
      ap_int <- paste0(ap_nm,   ":", moderator)
      p2_int <- paste0(p2_nm,   ":", moderator)

      mod_contrasts <- list(
        a1 = make_contrast(setNames(
          c(scale_factor,  scale_factor),
          c(a_int, p_int)
        )),
        a2 = make_contrast(setNames(
          c(scale_factor,  scale_factor,  scale_factor),
          c(a2_int, ap_int, p2_int)
        )),
        a3 = make_contrast(setNames(
          c(scale_factor, -scale_factor),
          c(a_int, p_int)
        )),
        a4 = make_contrast(setNames(
          c(scale_factor, -scale_factor,  scale_factor),
          c(a2_int, ap_int, p2_int)
        )),
        a5 = make_contrast(setNames(
          c(scale_factor, -scale_factor),
          c(a2_int, p2_int)
        ))
      )
      all_results$moderation <- run_rsa_contrasts(
        mod_contrasts, "Moderation"
      )
    }
  }

  invisible(list(model   = model,
                 results = all_results,
                 coefs   = coefs,
                 nms     = nms))
}
