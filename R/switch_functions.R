#' Switch a Distinguishable APIM model with two-intercept approach to the interaction approach
#'
#' @description
#' Reformulates a fitted \code{gls} or \code{glmmTMB} Distinguishable Actor-Partner
#' Interdependence Model (APIM) from the two-intercept parameterisation to the
#' equivalent interaction approach. The two parameterisations are statistically
#' equivalent; a numerical ML deviance check confirms this before returning.
#'
#' For \code{gls} models, non-formula arguments (\code{correlation},
#' \code{weights}, \code{na.action}) are preserved automatically by
#' \code{update()}.
#'
#' @param object A fitted \code{gls} or \code{glmmTMB} model object specified
#'   using the two-intercept approach.
#' @param disting_variable Character. Name of the effect-coded distinguishing
#'
#' @return A fitted model object of the same class as \code{object},
#'   reformulated using the interaction approach and estimated with the same
#'   method as the original (REML for \code{glmmTMB}; the original
#'   \code{method} for \code{gls}). Stops with an informative message if the
#'   ML deviance of the reformulated model differs from the original.
#'
#' @examples
#' \dontrun{
#' # gls
#' m_twoint <- gls(RelSat_A ~ 0 + man + woman +
#'                             man:c_Amity_A + woman:c_Amity_A +
#'                             man:c_Amity_P + woman:c_Amity_P,
#'                 correlation = corCompSymm(form = ~1 | DyadID),
#'                 weights     = varIdent(form = ~1 | ECGender_A),
#'                 data        = pairwise_disting)
#'
#' m_interact <- switch_to_interact(object         = m_twoint,
#'                                  disting_variable = "ECGender_A")
#'
#' # glmmTMB
#' m_twoint_tmb <- glmmTMB(RelSat_A ~ 0 + man + woman +
#'                                     man:c_Amity_A + woman:c_Amity_A +
#'                                     man:c_Amity_P + woman:c_Amity_P +
#'                                     cs(0 + man + woman | DyadID),
#'                          dispformula = ~0, REML = TRUE,
#'                          data = pairwise_disting)
#'
#' m_interact_tmb <- switch_to_interact(object           = m_twoint_tmb,
#'                                      disting_variable = "ECGender_A")
#' }
#' @importFrom stats formula update as.formula
switch_to_interact <- function(object = NULL, disting_variable = NULL) {
  # Function for gls objects
  if (!is.null(object) && inherits(object, "gls"))
  {
    fixeff_formula <- formula(object, fixed.only = TRUE)
    outcome_var <- as.character(fixeff_formula[2])
    extracted_vars <- all.vars(fixeff_formula)
    omit_vars_int <- c(":", "_A$", "_P$", outcome_var)
    vars_int <- extracted_vars[!grepl(paste(omit_vars_int, collapse = "|"), extracted_vars)]
    int_1 <- vars_int[[1]]
    int_2 <- vars_int[[2]]
    omit_vars_predictor <- c(":", int_1, int_2, outcome_var)
    vars_predictor <- extracted_vars[!grepl(paste(omit_vars_predictor, collapse = "|"), extracted_vars)]
    newpredictorlist <- paste0(disting_variable, "*", vars_predictor)
    interaction_formula <- paste(newpredictorlist, collapse = " + ")
    updatedmodel <- update(object, paste(outcome_var, " ~ ", interaction_formula))
    oldmodel <- suppressWarnings(update(object, method = "ML"))
    oldmodelML <- summary(oldmodel)
    oldmodelML_deviance <- (-2*oldmodelML$logLik)
    newmodel <- suppressWarnings(update(updatedmodel, method = "ML"))
    newmodelML <- summary(newmodel)
    newmodelML_deviance <- (-2*newmodelML$logLik)
    if (round(oldmodelML_deviance, 6) == round(newmodelML_deviance, 6)) {
      return(updatedmodel)
    } else { stop("The deviance values from both models run with ML estimation do not match.
           There was a mistake. Please specify interaction approach on your own as you cannot trust these results.")
    }
  }
  # Function for glmmTMB objects
  if (!is.null(object) && inherits(object, "glmmTMB"))
  {
    fixeff_formula <- formula(object, fixed.only = TRUE)
    outcome_var <- as.character(fixeff_formula[2])
    extracted_vars <- all.vars(fixeff_formula)
    omit_vars_int <- c(":", "_A$", "_P$", outcome_var)
    vars_int <- extracted_vars[!grepl(paste(omit_vars_int, collapse = "|"), extracted_vars)]
    int_1 <- vars_int[[1]]
    int_2 <- vars_int[[2]]
    omit_vars_predictor <- c(":", int_1, int_2, outcome_var)
    vars_predictor <- extracted_vars[!grepl(paste(omit_vars_predictor, collapse = "|"), extracted_vars)]
    newpredictorlist <- paste0(disting_variable, "*", vars_predictor)
    interaction_formula <- paste(newpredictorlist, collapse = " + ")
    formula_string <- gsub("\\s+", " ", paste(deparse(formula(object)), collapse = " "))
    rhs <- trimws(sub("^[^~]*~\\s*", "", formula_string))
    # Get the fixed-effects RHS as a string, then find where it ends in the
    # full RHS. Everything after the last fixed-effect term is the RE portion.
    fixeff_string <- gsub("\\s+", " ", paste(deparse(formula(object, fixed.only = TRUE)), collapse = " "))
    fixeff_rhs <- trimws(sub("^[^~]*~\\s*", "", fixeff_string))
    # The RE portion starts after the fixed-effects RHS in the full RHS
    re_start <- nchar(fixeff_rhs) + 1L
    random_effect_term <- trimws(sub("^\\s*\\+\\s*", "", substr(rhs, re_start, nchar(rhs))))
    newformula <- paste0(outcome_var, " ~ ", interaction_formula, " + ", random_effect_term)
    updatedmodel <- update(object, formula. = newformula)
    oldmodel <- suppressWarnings(update(object, REML = FALSE))
    oldmodelML <- summary(oldmodel)
    oldmodelML_deviance <- -2*(oldmodelML$logLik[1])
    newmodel <- suppressWarnings(update(updatedmodel, REML = FALSE))
    newmodelML <- summary(newmodel)
    newmodelML_deviance <- -2*(newmodelML$logLik[1])
    if (round(oldmodelML_deviance, 6) == round(newmodelML_deviance, 6)) {
      return(updatedmodel)
    } else { stop("The deviance values from both models run with ML estimation do not match.
           There was a mistake. Please specify interaction approach on your own as you cannot trust these results.")
    }
  }
}


#' Switch a Distinguishable APIM model with interaction approach to the two-intercept approach
#'
#' @description
#' Reformulates a fitted \code{gls} or \code{glmmTMB} Distinguishable Actor-Partner
#' Interdependence Model (APIM) from the interaction approach parameterisation to the
#' equivalent two-intercept approach. The two parameterisations are
#' statistically equivalent; a numerical ML deviance check confirms this before
#' returning.
#'
#' The two-intercept approach suppresses the global intercept and estimates a
#' separate intercept for each dyad member, with all predictors crossed with
#' their respective group dummy variable using \code{:}.
#'
#' For \code{gls} models, non-formula arguments (\code{correlation},
#' \code{weights}, \code{na.action}) are preserved automatically by
#' \code{update()}.
#'
#' @param object A fitted \code{gls} or \code{glmmTMB} model object specified
#'   using the two-intercept approach.
#' @param disting_variable Character. Name of the effect-coded distinguishing
#' @param disting_level_1 Character. Name of the dummy-coded variable for the
#'   first dyad member in the dataset (e.g. \code{"man"}).
#' @param disting_level_2 Character. Name of the dummy-coded variable for the
#'   second dyad member in the dataset (e.g. \code{"woman"}).
#'
#' @return A fitted model object of the same class as \code{object},
#'   reformulated using the two-intercept approach and estimated with the same
#'   method as the original (REML for \code{glmmTMB}; the original
#'   \code{method} for \code{gls}). Stops with an informative message if the
#'   ML deviance of the reformulated model differs from the original.
#'
#' @examples
#' \dontrun{
#' # gls
#' m_interact <- gls(RelSat_A ~ c_Amity_A * ECGender_A +
#'                               c_Amity_P * ECGender_A,
#'                   correlation = corCompSymm(form = ~1 | DyadID),
#'                   weights     = varIdent(form = ~1 | ECGender_A),
#'                   data        = pairwise_disting)
#'
#' m_twoint <- switch_to_twoint(object           = m_interact,
#'                              disting_variable  = "ECGender_A",
#'                              disting_level_1   = "man",
#'                              disting_level_2   = "woman")
#'
#' # glmmTMB
#' m_interact_tmb <- glmmTMB(RelSat_A ~ c_Amity_A * ECGender_A +
#'                                       c_Amity_P * ECGender_A +
#'                                       cs(0 + man + woman | DyadID),
#'                            dispformula = ~0, REML = TRUE,
#'                            data = pairwise_disting)
#'
#' m_twoint_tmb <- switch_to_twoint(object           = m_interact_tmb,
#'                                  disting_variable  = "ECGender_A",
#'                                  disting_level_1   = "man",
#'                                  disting_level_2   = "woman")
#' }
#' @importFrom stats formula update as.formula
switch_to_twoint <- function(object = NULL, disting_variable = NULL, disting_level_1 = NULL, disting_level_2 = NULL) {
  # Function for gls objects
  if (!is.null(object) && inherits(object, "gls"))
  {
    fixeff_formula <- formula(object, fixed.only = TRUE)
    outcome_var <- as.character(fixeff_formula[2])
    extracted_vars <- all.vars(fixeff_formula)
    omit_vars_levels <- c(disting_variable, outcome_var)
    # Bug fix: was omit_vars_level.1 (undefined variable) in the original
    vars_levels <- extracted_vars[!grepl(paste(omit_vars_levels, collapse = "|"), extracted_vars)]
    new_interceptlist <- paste0("0 + ", disting_level_1, " + ", disting_level_2)
    new_level_1_predictorlist <- paste0(disting_level_1, ":", vars_levels)
    new_level_2_predictorlist <- paste0(disting_level_2, ":", vars_levels)
    twoint_formula_level_1 <- paste(new_level_1_predictorlist, collapse = " + ")
    twoint_formula_level_2 <- paste(new_level_2_predictorlist, collapse = " + ")
    twoint_formula <- paste0(new_interceptlist, " + ", twoint_formula_level_1, " + ", twoint_formula_level_2)
    updatedmodel <- update(object, paste(outcome_var, " ~ ", twoint_formula))
    oldmodel <- suppressWarnings(update(object, method = "ML"))
    oldmodelML <- summary(oldmodel)
    oldmodelML_deviance <- (-2*oldmodelML$logLik)
    newmodel <- suppressWarnings(update(updatedmodel, method = "ML"))
    newmodelML <- summary(newmodel)
    newmodelML_deviance <- (-2*newmodelML$logLik)
    if (round(oldmodelML_deviance, 6) == round(newmodelML_deviance, 6)) {
      return(updatedmodel)
    } else { stop("The deviance values from both models run with ML estimation do not match.
           There was a mistake. Please specify two-intercept approach on your own as you cannot trust these results.")
    }
  }
  # Function for glmmTMB objects
  if (!is.null(object) && inherits(object, "glmmTMB"))
  {
    fixeff_formula <- formula(object, fixed.only = TRUE)
    outcome_var <- as.character(fixeff_formula[2])
    extracted_vars <- all.vars(fixeff_formula)
    omit_vars_levels <- c(disting_variable, outcome_var)
    vars_levels <- extracted_vars[!grepl(paste(omit_vars_levels, collapse = "|"), extracted_vars)]
    new_interceptlist <- paste0("0 + ", disting_level_1, " + ", disting_level_2)
    new_level_1_predictorlist <- paste0(disting_level_1, ":", vars_levels)
    new_level_2_predictorlist <- paste0(disting_level_2, ":", vars_levels)
    twoint_formula_level_1 <- paste(new_level_1_predictorlist, collapse = " + ")
    twoint_formula_level_2 <- paste(new_level_2_predictorlist, collapse = " + ")
    twoint_formula <- paste0(new_interceptlist, " + ", twoint_formula_level_1, " + ", twoint_formula_level_2)
    formula_string <- gsub("\\s+", " ", paste(deparse(formula(object)), collapse = " "))
    rhs <- trimws(sub("^[^~]*~\\s*", "", formula_string))
    # Get the fixed-effects RHS as a string, then find where it ends in the
    # full RHS. Everything after the last fixed-effect term is the RE portion.
    fixeff_string <- gsub("\\s+", " ", paste(deparse(formula(object, fixed.only = TRUE)), collapse = " "))
    fixeff_rhs <- trimws(sub("^[^~]*~\\s*", "", fixeff_string))
    # The RE portion starts after the fixed-effects RHS in the full RHS
    re_start <- nchar(fixeff_rhs) + 1L
    random_effect_term <- trimws(sub("^\\s*\\+\\s*", "", substr(rhs, re_start, nchar(rhs))))
    newformula_string <- paste0(outcome_var, " ~ ", twoint_formula, " + ", random_effect_term)
    newformula <- as.formula(newformula_string)
    updatedmodel <- update(object, formula. = newformula)
    oldmodel <- suppressWarnings(update(object, REML = FALSE))
    oldmodelML <- summary(oldmodel)
    oldmodelML_deviance <- -2*(oldmodelML$logLik[1])
    newmodel <- suppressWarnings(update(updatedmodel, REML = FALSE))
    newmodelML <- summary(newmodel)
    newmodelML_deviance <- -2*(newmodelML$logLik[1])
    if (round(oldmodelML_deviance, 6) == round(newmodelML_deviance, 6)) {
      return(updatedmodel)
    } else { stop("The deviance values from both models run with ML estimation do not match.
           There was a mistake. Please specify two-intercept approach on your own as you cannot trust these results.")
    }
  }
}
