# =============================================================================
# EXPORTED: check_rsa_requirements
# =============================================================================

#' Check Preliminary Requirements for Dyadic Response Surface Analysis
#'
#' @description
#' Evaluates whether the data meet the three requirements necessary for
#' meaningful response surface analysis (RSA) results: (1) adequate variance
#' in both actor and partner scores, (2) sufficient representation of
#' over-estimators and under-estimators (at least \code{min_pct}\% in each
#' category), and (3) a statistically significant R-squared from the
#' polynomial mixed model. Results are reported for the overall sample and,
#' optionally, separately for each level of a moderator variable.
#'
#' @details
#' Participants are classified into three categories based on the standardized
#' difference between actor and partner scores, matching the convention in
#' \code{RSA::summary.RSA} (Schönbrodt et al., 2018). Both variables are
#' first standardized to a common pooled mean and SD, then the difference
#' (partner_z - actor_z) is computed:
#' \describe{
#'   \item{Over-estimator}{Standardized difference > \code{congruence_threshold}
#'     (partner exceeds actor).}
#'   \item{Under-estimator}{Standardized difference < -\code{congruence_threshold}
#'     (actor exceeds partner).}
#'   \item{Congruent}{|standardized difference| <= \code{congruence_threshold}.}
#' }
#'
#' R-squared and its significance are estimated using \code{\link{estimate_Rsq}},
#' which fits null and polynomial \code{glmmTMB} models with the same random
#' effect structure specified in \code{re_term} and uses a deviance difference
#' test. This properly accounts for dyadic non-independence, unlike the OLS
#' approach used in \code{RSA::summary.RSA}. R-squared is always estimated on
#' the full sample. When a moderator is supplied, bias percentages are reported
#' separately per group but R-squared is not re-estimated per group.
#'
#' @param data A data frame.
#' @param outcome Character string. Name of the outcome variable.
#' @param actor Character string. Name of the (centered) actor predictor.
#' @param partner Character string. Name of the (centered) partner predictor.
#' @param re_term Character string. The random effects term exactly as it
#'   should appear in the model formula, e.g.
#'   \code{"cs(0 + as.factor(ECGender_A) | DyadID)"}. Passed verbatim to
#'   \code{glmmTMB}.
#' @param moderator Character string or \code{NULL}. Name of a grouping
#'   variable. When supplied, bias percentages are reported separately for
#'   each level. R-squared is always from the full sample. Default \code{NULL}.
#' @param mod_labels Named character vector or \code{NULL}. Labels for
#'   moderator levels, e.g. \code{c("-1" = "Women", "1" = "Men")}.
#'   Default \code{NULL} uses raw moderator values as labels.
#' @param congruence_threshold Numeric. Cutpoint applied to the standardized
#'   difference (partner_z - actor_z) for classifying cases as congruent.
#'   Default \code{0.5}, matching \code{RSA::summary.RSA}.
#' @param min_pct Numeric. Minimum percentage of cases required in each
#'   bias category. Default \code{10}.
#'
#' @return Invisibly returns a \code{data.frame} with one row per group.
#'
#' @references
#' Schönbrodt, F. D., Humberg, S., & Nestler, S. (2018). Testing similarity
#' effects with dyadic response surface analysis. \emph{European Journal of
#' Personality}, 32(6), 627--641. \doi{10.1002/per.2169}
#'
#' Edwards, J. R. (1994). The study of congruence in organizational behavior
#' research: Critique and a proposed alternative. \emph{Organizational
#' Behavior and Human Decision Processes}, 58(1), 51--100.
#'
#' @examples
#' \dontrun{
#' check_rsa_requirements(
#'   data       = pairwise_data,
#'   outcome    = "RelSat_A",
#'   actor      = "c_PosEmotion_A",
#'   partner    = "c_PosEmotion_P",
#'   re_term    = "cs(0 + as.factor(ECGender_A) | DyadID)",
#'   moderator  = "ECGender_A",
#'   mod_labels = c("-1" = "Women", "1" = "Men")
#' )
#' }
#'
#' @importFrom glmmTMB glmmTMB
#' @importFrom stats sd cor as.formula
#' @export
check_rsa_requirements <- function(
    data,
    outcome,
    actor,
    partner,
    re_term,
    moderator             = NULL,
    mod_labels            = NULL,
    congruence_threshold  = 0.5,
    min_pct               = 10
) {

  # -- build polynomial term names -------------------------------------------
  a2_nm <- paste0(actor,   "2")
  p2_nm <- paste0(partner, "2")
  ap_nm <- paste0(actor,   "_x_", partner)

  data[[a2_nm]] <- data[[actor]]^2
  data[[p2_nm]] <- data[[partner]]^2
  data[[ap_nm]] <- data[[actor]] * data[[partner]]

  # -- helper: bias classification for any data frame ------------------------
  get_bias <- function(df) {
    common_m  <- mean(c(df[[actor]], df[[partner]]), na.rm = TRUE)
    common_sd <- sd(c(df[[actor]], df[[partner]]),   na.rm = TRUE)
    actor_z   <- (df[[actor]]   - common_m) / common_sd
    partner_z <- (df[[partner]] - common_m) / common_sd
    z_diff    <- partner_z - actor_z

    bias_cat <- ifelse(
      z_diff >  congruence_threshold,  "Over-estimator",
      ifelse(z_diff < -congruence_threshold, "Under-estimator", "Congruent")
    )

    bias_tbl  <- prop.table(table(bias_cat)) * 100
    pct_over  <- round(unname(bias_tbl["Over-estimator"]),  1)
    pct_under <- round(unname(bias_tbl["Under-estimator"]), 1)
    pct_cong  <- round(unname(bias_tbl["Congruent"]),       1)

    list(
      pct_over  = ifelse(is.na(pct_over),  0, pct_over),
      pct_under = ifelse(is.na(pct_under), 0, pct_under),
      pct_cong  = ifelse(is.na(pct_cong),  0, pct_cong)
    )
  }

  # -- step 1: R2 from full sample via estimate_Rsq --------------------------
  assign(".rsa_req_data_tmp", data, envir = .GlobalEnv)
  on.exit(
    rm(list = ".rsa_req_data_tmp", envir = .GlobalEnv, inherits = FALSE),
    add = TRUE
  )

  null_formula <- as.formula(paste(outcome, "~ 1 +", re_term))
  poly_formula <- as.formula(paste(
    outcome, "~", actor, "+", partner, "+",
    a2_nm, "+", ap_nm, "+", p2_nm, "+", re_term
  ))

  model_poly <- suppressWarnings(
    glmmTMB::glmmTMB(poly_formula,
                     dispformula = ~0,
                     REML        = FALSE,
                     data        = .rsa_req_data_tmp)
  )
  model_poly$call$data <- quote(.rsa_req_data_tmp)

  r2_result <- estimate_Rsq(
    model_full        = model_poly,
    indistinguishable = TRUE
  )

  # parse R2 and p-value from returned text string
  r2_text <- as.character(r2_result)
  r2_val  <- as.numeric(
    regmatches(r2_text,
               regexpr("(?<=R\u00b2 = )[.0-9]+", r2_text, perl = TRUE))
  )
  p_val   <- as.numeric(
    regmatches(r2_text,
               regexpr("(?<=p = )[.0-9eE+-]+", r2_text, perl = TRUE))
  )

  # fallback if unicode matching fails
  if (length(r2_val) == 0 || is.na(r2_val)) {
    r2_val <- as.numeric(
      regmatches(r2_text,
                 regexpr("(?<=R2 = )[.0-9]+", r2_text, perl = TRUE))
    )
  }

  r2_sig <- !is.na(p_val) && p_val < .05
  stars  <- ifelse(is.na(p_val),   "",
                   ifelse(p_val < .001, "***",
                          ifelse(p_val < .01,  "**",
                                 ifelse(p_val < .05,  "*", ""))))
  r2_fmt <- if (!is.na(r2_val)) paste0(sprintf("%.3f", r2_val), stars) else "NA"

  # -- step 2: build output rows ---------------------------------------------
  make_row <- function(df, group_label) {
    b       <- get_bias(df)
    meets_bias <- b$pct_over >= min_pct & b$pct_under >= min_pct
    meets_all  <- meets_bias & r2_sig
    data.frame(
      Group           = group_label,
      N               = nrow(df),
      SD_actor        = round(sd(df[[actor]],   na.rm = TRUE), 2),
      SD_partner      = round(sd(df[[partner]], na.rm = TRUE), 2),
      r_actor_partner = round(cor(df[[actor]], df[[partner]],
                                  use = "complete.obs"), 2),
      R2              = r2_fmt,
      Pct_under       = b$pct_under,
      Pct_congruent   = b$pct_cong,
      Pct_over        = b$pct_over,
      Meets_bias_req  = ifelse(meets_bias, "YES", "NO"),
      Meets_R2_req    = ifelse(r2_sig,     "YES", "NO"),
      Meets_all       = ifelse(meets_all,  "YES", "NO"),
      stringsAsFactors = FALSE
    )
  }

  results              <- list()
  results[["Overall"]] <- make_row(data, "Overall")

  if (!is.null(moderator)) {
    grp_vals <- sort(unique(data[[moderator]]))
    for (val in grp_vals) {
      df_sub <- data[data[[moderator]] == val, ]
      lbl    <- if (!is.null(mod_labels) &&
                    as.character(val) %in% names(mod_labels)) {
        mod_labels[as.character(val)]
      } else {
        as.character(val)
      }
      results[[lbl]] <- make_row(df_sub, lbl)
    }
  }

  out <- do.call(rbind, results)
  rownames(out) <- NULL

  # -- print -----------------------------------------------------------------
  cat("\n=== RSA Requirements Check ===\n\n")
  cat("Outcome  :", outcome, "\n")
  cat("Actor    :", actor,   "\n")
  cat("Partner  :", partner, "\n")
  cat("RE term  :", re_term, "\n")
  cat(sprintf("Congruence threshold : |standardized A - P| <= %g\n",
              congruence_threshold))
  cat(sprintf("Minimum %% per bias category : %g%%\n\n", min_pct))

  cat(sprintf(
    "%-12s %5s %8s %8s %8s %10s %8s %10s %8s %12s %10s %10s\n",
    "Group", "N", "SD_A", "SD_P", "r(A,P)", "R2",
    "%Under", "%Congruent", "%Over", "Bias req?", "R2 req?", "Meets all?"
  ))
  cat(strrep("-", 110), "\n")

  for (i in seq_len(nrow(out))) {
    cat(sprintf(
      "%-12s %5d %8.2f %8.2f %8.2f %10s %8.1f %10.1f %8.1f %12s %10s %10s\n",
      out$Group[i], out$N[i],
      out$SD_actor[i], out$SD_partner[i], out$r_actor_partner[i],
      out$R2[i],
      out$Pct_under[i], out$Pct_congruent[i], out$Pct_over[i],
      out$Meets_bias_req[i], out$Meets_R2_req[i], out$Meets_all[i]
    ))
  }

  cat("\nNote. Over-estimator = Partner > Actor; ",
      "Under-estimator = Actor > Partner.\n")
  cat("Bias classification uses standardized difference (partner_z - actor_z),\n")
  cat("matching RSA::summary.RSA convention (Schonbrodt et al., 2018).\n")
  cat("R2 estimated from full sample via glmmTMB deviance difference test.\n")
  cat("R2 significance: * p < .05, ** p < .01, *** p < .001\n")

  invisible(out)
}
