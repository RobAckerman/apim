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
#' polynomial regression. Results are reported for the overall sample and,
#' optionally, separately for each level of a moderator variable.
#'
#' @details
#' Participants are classified into three categories based on the difference
#' between actor and partner scores:
#' \describe{
#'   \item{Over-estimator}{Actor score > Partner score by more than
#'     \code{congruence_threshold}.}
#'   \item{Under-estimator}{Actor score < Partner score by more than
#'     \code{congruence_threshold}.}
#'   \item{Congruent}{|Actor - Partner| <= \code{congruence_threshold}.}
#' }
#'
#' The R-squared and its significance are derived from an OLS polynomial
#' regression of the outcome on the five RSA terms (actor, partner, actor²,
#' actor×partner, partner²). This is used as a preliminary check only; final
#' inference should come from the mixed model fitted by
#' \code{\link{run_dyadic_rsa}}.
#'
#' @param data A data frame.
#' @param outcome Character string. Name of the outcome variable.
#' @param actor Character string. Name of the (centered) actor predictor.
#' @param partner Character string. Name of the (centered) partner predictor.
#' @param moderator Character string or \code{NULL}. Name of a grouping
#'   variable. When supplied, checks are run separately for each level in
#'   addition to the overall sample. Default \code{NULL}.
#' @param mod_labels Named character vector or \code{NULL}. Labels for
#'   moderator levels, e.g. \code{c("-1" = "Women", "1" = "Men")}.
#'   Default \code{NULL} uses raw moderator values as labels.
#' @param congruence_threshold Numeric. Maximum absolute difference
#'   |Actor - Partner| for a case to be classified as congruent.
#'   Default \code{0} (exact equality required for congruence).
#' @param min_pct Numeric. Minimum percentage of cases required in each
#'   bias category (over- and under-estimators) for the requirement to be
#'   met. Default \code{10}.
#'
#' @return Invisibly returns a \code{data.frame} with one row per group
#'   (overall plus any moderator levels) and columns:
#'   \describe{
#'     \item{Group}{Group label.}
#'     \item{N}{Sample size.}
#'     \item{SD_actor}{Standard deviation of actor scores.}
#'     \item{SD_partner}{Standard deviation of partner scores.}
#'     \item{r_actor_partner}{Correlation between actor and partner scores.}
#'     \item{R2}{R-squared from polynomial regression, with significance stars.}
#'     \item{Pct_under}{Percentage of under-estimators.}
#'     \item{Pct_congruent}{Percentage of congruent cases.}
#'     \item{Pct_over}{Percentage of over-estimators.}
#'     \item{Meets_bias_req}{Whether both bias categories meet \code{min_pct}.}
#'     \item{Meets_R2_req}{Whether the polynomial R-squared is significant.}
#'     \item{Meets_all}{Whether all requirements are met.}
#'   }
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
#'   moderator  = "ECGender_A",
#'   mod_labels = c("-1" = "Women", "1" = "Men"),
#'   min_pct    = 10
#' )
#' }
#'
#' @importFrom stats sd cor lm pf
#' @export
check_rsa_requirements <- function(
    data,
    outcome,
    actor,
    partner,
    moderator             = NULL,
    mod_labels            = NULL,
    congruence_threshold  = 0,
    min_pct               = 10
) {

  # -- build polynomial term names --------------------------------------------
  a2_nm <- paste0(actor,   "2")
  p2_nm <- paste0(partner, "2")
  ap_nm <- paste0(actor,   "_x_", partner)

  data[[a2_nm]] <- data[[actor]]^2
  data[[p2_nm]] <- data[[partner]]^2
  data[[ap_nm]] <- data[[actor]] * data[[partner]]

  # -- classify bias ----------------------------------------------------------
  data$diff_score <- data[[actor]] - data[[partner]]
  data$bias_cat   <- ifelse(
    data$diff_score >  congruence_threshold, "Over-estimator",
    ifelse(
      data$diff_score < -congruence_threshold, "Under-estimator",
      "Congruent"
    )
  )

  # -- core check function for one (sub)group ---------------------------------
  run_checks <- function(df, group_label) {

    n          <- nrow(df)
    sd_actor   <- sd(df[[actor]],   na.rm = TRUE)
    sd_partner <- sd(df[[partner]], na.rm = TRUE)
    r_ap       <- cor(df[[actor]], df[[partner]], use = "complete.obs")

    bias_tbl  <- prop.table(table(df$bias_cat)) * 100
    pct_over  <- round(unname(bias_tbl["Over-estimator"]),  1)
    pct_under <- round(unname(bias_tbl["Under-estimator"]), 1)
    pct_cong  <- round(unname(bias_tbl["Congruent"]),       1)

    pct_over  <- ifelse(is.na(pct_over),  0, pct_over)
    pct_under <- ifelse(is.na(pct_under), 0, pct_under)
    pct_cong  <- ifelse(is.na(pct_cong),  0, pct_cong)

    meets_bias <- pct_over >= min_pct & pct_under >= min_pct

    formula_poly <- as.formula(paste(
      outcome, "~", actor, "+", partner, "+",
      a2_nm, "+", ap_nm, "+", p2_nm
    ))
    lm_fit  <- lm(formula_poly, data = df)
    lm_sum  <- summary(lm_fit)
    r2      <- round(lm_sum$r.squared, 3)
    f_stat  <- lm_sum$fstatistic
    p_val   <- pf(f_stat[1], f_stat[2], f_stat[3], lower.tail = FALSE)
    r2_sig  <- p_val < .05

    stars  <- ifelse(p_val < .001, "***",
                     ifelse(p_val < .01,  "**",
                            ifelse(p_val < .05,  "*", "")))
    r2_fmt <- paste0(sprintf("%.3f", r2), stars)

    meets_all <- meets_bias & r2_sig

    data.frame(
      Group           = group_label,
      N               = n,
      SD_actor        = round(sd_actor,   2),
      SD_partner      = round(sd_partner, 2),
      r_actor_partner = round(r_ap,       2),
      R2              = r2_fmt,
      Pct_under       = pct_under,
      Pct_congruent   = pct_cong,
      Pct_over        = pct_over,
      Meets_bias_req  = ifelse(meets_bias, "YES", "NO"),
      Meets_R2_req    = ifelse(r2_sig,     "YES", "NO"),
      Meets_all       = ifelse(meets_all,  "YES", "NO"),
      stringsAsFactors = FALSE
    )
  }

  # -- run overall + by group -------------------------------------------------
  results        <- list()
  results[["Overall"]] <- run_checks(data, "Overall")

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
      results[[lbl]] <- run_checks(df_sub, lbl)
    }
  }

  out <- do.call(rbind, results)
  rownames(out) <- NULL

  # -- print ------------------------------------------------------------------
  cat("\n=== RSA Requirements Check ===\n\n")
  cat("Outcome  :", outcome,  "\n")
  cat("Actor    :", actor,    "\n")
  cat("Partner  :", partner,  "\n")
  cat(sprintf("Congruence threshold : |A - P| <= %g\n", congruence_threshold))
  cat(sprintf("Minimum %% per bias category : %g%%\n\n", min_pct))

  cat(sprintf("%-12s %5s %8s %8s %8s %10s %8s %10s %8s %12s %10s %10s\n",
              "Group", "N", "SD_A", "SD_P", "r(A,P)", "R2",
              "%Under", "%Congruent", "%Over", "Bias req?", "R2 req?", "Meets all?"))
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

  cat("\nNote. Over-estimator = Actor > Partner; ",
      "Under-estimator = Actor < Partner.\n")
  cat("R2 significance: * p < .05, ** p < .01, *** p < .001\n")

  invisible(out)
}
