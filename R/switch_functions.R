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
#'   variable (e.g. \code{"ECGender_A"}).
#'
#' @return A fitted model object of the same class as \code{object},
#'   reformulated using the interaction approach and estimated with the same
#'   method as the original (REML for \code{glmmTMB}; the original
#'   \code{method} for \code{gls}). Stops with an informative message if the
#'   ML deviance of the reformulated model differs from the original.
#'
#' @examples
#' \dontrun{
#' # gls with dummy codes
#' m_twoint <- gls(RelSat_A ~ 0 + man + woman +
#'                             man:c_Amity_A + woman:c_Amity_A +
#'                             man:c_Amity_P + woman:c_Amity_P,
#'                 correlation = corCompSymm(form = ~1 | DyadID),
#'                 weights     = varIdent(form = ~1 | ECGender_A),
#'                 data        = pairwise_disting)
#'
#' m_interact <- switch_to_interact(object           = m_twoint,
#'                                  disting_variable = "ECGender_A")
#'
#' # gls with as.factor() distinguishing variable
#' m_twoint_f <- gls(RelSat_A ~ 0 + as.factor(ECGender_A) +
#'                               as.factor(ECGender_A):c_Amity_A +
#'                               as.factor(ECGender_A):c_Amity_P,
#'                   correlation = corCompSymm(form = ~1 | DyadID),
#'                   weights     = varIdent(form = ~1 | ECGender_A),
#'                   data        = pairwise_disting)
#'
#' m_interact_f <- switch_to_interact(object           = m_twoint_f,
#'                                    disting_variable = "ECGender_A")
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
#' @importFrom stats formula update as.formula terms
#' @export
switch_to_interact <- function(object = NULL, disting_variable = NULL) {

  # ---------------------------------------------------------------------------
  # Helper: extract predictor names from interaction terms robustly.
  # Handles three cases:
  #   1. Effect-coded: ECGender_A*c_Amity_A  -> terms contain disting_variable
  #   2. as.factor():  as.factor(ECGender_A):c_Amity_A -> terms contain disting_variable
  #   3. Dummy codes:  man:c_Amity_A or c_Amity_A:man -> neither part is disting_variable
  #      In this case, person codes are identified as terms without ":"
  # ---------------------------------------------------------------------------
  .get_predictors <- function(all_terms, disting_variable) {
    vars_with_colon <- all_terms[grepl(":", all_terms)]
    if (length(vars_with_colon) == 0) return(character(0))

    has_disting <- any(grepl(disting_variable, vars_with_colon, fixed = TRUE))

    if (has_disting) {
      # effect-coded or as.factor() model
      # keep the part of each interaction that does NOT contain disting_variable
      preds <- unique(sapply(vars_with_colon, function(term) {
        parts       <- strsplit(term, ":")[[1]]
        non_disting <- parts[!grepl(disting_variable, parts, fixed = TRUE)]
        if (length(non_disting) == 1) non_disting else NA_character_
      }))
    } else {
      # dummy-code two-intercept model
      # person codes = terms without ":" (pure intercept terms e.g. man, woman)
      person_codes <- all_terms[!grepl(":", all_terms)]

      # keep the part of each interaction that is NOT a person code
      preds <- unique(sapply(vars_with_colon, function(term) {
        parts      <- strsplit(term, ":")[[1]]
        non_person <- parts[!parts %in% person_codes]
        if (length(non_person) == 1) non_person else NA_character_
      }))
    }

    preds[!is.na(preds)]
  }

  # ---------------------------------------------------------------------------
  # Helper: extract RE string using bracket-matching — robust to RE term
  # appearing anywhere in the formula.
  # ---------------------------------------------------------------------------
  .get_re_string <- function(object) {
    formula_string <- gsub("\\s+", " ",
                           paste(deparse(formula(object)), collapse = " "))
    rhs      <- trimws(sub("^[^~]*~\\s*", "", formula_string))
    pipe_pos <- regexpr("\\|", rhs)
    if (pipe_pos == -1) return("")
    chars <- strsplit(rhs, "")[[1]]

    depth <- 0
    left  <- pipe_pos
    for (i in pipe_pos:1) {
      if (chars[i] == ")") depth <- depth + 1
      if (chars[i] == "(") {
        if (depth == 0) { left <- i; break }
        depth <- depth - 1
      }
    }
    depth <- 0
    right <- pipe_pos
    for (i in pipe_pos:nchar(rhs)) {
      if (chars[i] == "(") depth <- depth + 1
      if (chars[i] == ")") {
        if (depth == 0) { right <- i; break }
        depth <- depth - 1
      }
    }
    func_start <- left
    for (i in (left - 1):1) {
      if (grepl("[a-zA-Z0-9_\\.]", chars[i])) {
        func_start <- i
      } else {
        break
      }
    }
    trimws(substr(rhs, func_start, right))
  }

  # -- gls -------------------------------------------------------------------
  if (!is.null(object) && inherits(object, "gls")) {
    fixeff_formula <- formula(object, fixed.only = TRUE)
    outcome_var    <- as.character(fixeff_formula[2])
    all_terms      <- attr(terms(fixeff_formula), "term.labels")

    vars_predictor      <- .get_predictors(all_terms, disting_variable)
    newpredictorlist    <- paste0(disting_variable, "*", vars_predictor)
    interaction_formula <- paste(newpredictorlist, collapse = " + ")

    updatedmodel <- update(object,
                           paste(outcome_var, " ~ ", interaction_formula))

    oldmodel       <- suppressWarnings(update(object, method = "ML"))
    oldmodelML_dev <- -2 * summary(oldmodel)$logLik
    newmodel       <- suppressWarnings(update(updatedmodel, method = "ML"))
    newmodelML_dev <- -2 * summary(newmodel)$logLik

    if (round(oldmodelML_dev, 6) == round(newmodelML_dev, 6)) {
      return(updatedmodel)
    } else {
      stop("The deviance values from both models run with ML estimation do ",
           "not match. There was a mistake. Please specify the interaction ",
           "approach on your own as you cannot trust these results.")
    }
  }

  # -- glmmTMB ---------------------------------------------------------------
  if (!is.null(object) && inherits(object, "glmmTMB")) {
    fixeff_formula <- formula(object, fixed.only = TRUE)
    outcome_var    <- as.character(fixeff_formula[2])
    all_terms      <- attr(terms(fixeff_formula), "term.labels")
    re_string      <- .get_re_string(object)

    vars_predictor      <- .get_predictors(all_terms, disting_variable)
    newpredictorlist    <- paste0(disting_variable, "*", vars_predictor)
    interaction_formula <- paste(newpredictorlist, collapse = " + ")

    newformula   <- paste0(outcome_var, " ~ ", interaction_formula,
                           " + ", re_string)
    updatedmodel <- update(object, formula. = as.formula(newformula))

    oldmodel       <- suppressWarnings(update(object, REML = FALSE))
    oldmodelML_dev <- -2 * summary(oldmodel)$logLik[1]
    newmodel       <- suppressWarnings(update(updatedmodel, REML = FALSE))
    newmodelML_dev <- -2 * summary(newmodel)$logLik[1]

    if (round(oldmodelML_dev, 6) == round(newmodelML_dev, 6)) {
      return(updatedmodel)
    } else {
      stop("The deviance values from both models run with ML estimation do ",
           "not match. There was a mistake. Please specify the interaction ",
           "approach on your own as you cannot trust these results.")
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
#'   using the interaction approach.
#' @param disting_variable Character. Name of the effect-coded distinguishing
#'   variable (e.g. \code{"ECGender_A"}).
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
#' @importFrom stats formula update as.formula terms
#' @export
switch_to_twoint <- function(object           = NULL,
                             disting_variable  = NULL,
                             disting_level_1   = NULL,
                             disting_level_2   = NULL) {

  # ---------------------------------------------------------------------------
  # Helper: extract predictor names from interaction terms robustly.
  # ---------------------------------------------------------------------------
  .get_predictors <- function(all_terms, disting_variable) {
    vars_with_colon <- all_terms[grepl(":", all_terms)]
    if (length(vars_with_colon) == 0) return(character(0))

    has_disting <- any(grepl(disting_variable, vars_with_colon, fixed = TRUE))

    if (has_disting) {
      preds <- unique(sapply(vars_with_colon, function(term) {
        parts       <- strsplit(term, ":")[[1]]
        non_disting <- parts[!grepl(disting_variable, parts, fixed = TRUE)]
        if (length(non_disting) == 1) non_disting else NA_character_
      }))
    } else {
      person_codes <- all_terms[!grepl(":", all_terms)]
      preds <- unique(sapply(vars_with_colon, function(term) {
        parts      <- strsplit(term, ":")[[1]]
        non_person <- parts[!parts %in% person_codes]
        if (length(non_person) == 1) non_person else NA_character_
      }))
    }

    preds[!is.na(preds)]
  }

  # ---------------------------------------------------------------------------
  # Helper: extract RE string using bracket-matching.
  # ---------------------------------------------------------------------------
  .get_re_string <- function(object) {
    formula_string <- gsub("\\s+", " ",
                           paste(deparse(formula(object)), collapse = " "))
    rhs      <- trimws(sub("^[^~]*~\\s*", "", formula_string))
    pipe_pos <- regexpr("\\|", rhs)
    if (pipe_pos == -1) return("")
    chars <- strsplit(rhs, "")[[1]]

    depth <- 0
    left  <- pipe_pos
    for (i in pipe_pos:1) {
      if (chars[i] == ")") depth <- depth + 1
      if (chars[i] == "(") {
        if (depth == 0) { left <- i; break }
        depth <- depth - 1
      }
    }
    depth <- 0
    right <- pipe_pos
    for (i in pipe_pos:nchar(rhs)) {
      if (chars[i] == "(") depth <- depth + 1
      if (chars[i] == ")") {
        if (depth == 0) { right <- i; break }
        depth <- depth - 1
      }
    }
    func_start <- left
    for (i in (left - 1):1) {
      if (grepl("[a-zA-Z0-9_\\.]", chars[i])) {
        func_start <- i
      } else {
        break
      }
    }
    trimws(substr(rhs, func_start, right))
  }

  # -- gls -------------------------------------------------------------------
  if (!is.null(object) && inherits(object, "gls")) {
    fixeff_formula <- formula(object, fixed.only = TRUE)
    outcome_var    <- as.character(fixeff_formula[2])
    all_terms      <- attr(terms(fixeff_formula), "term.labels")

    vars_predictor <- .get_predictors(all_terms, disting_variable)

    new_interceptlist         <- paste0("0 + ", disting_level_1, " + ",
                                        disting_level_2)
    new_level_1_predictorlist <- paste0(disting_level_1, ":", vars_predictor)
    new_level_2_predictorlist <- paste0(disting_level_2, ":", vars_predictor)
    twoint_formula <- paste0(
      new_interceptlist, " + ",
      paste(new_level_1_predictorlist, collapse = " + "), " + ",
      paste(new_level_2_predictorlist, collapse = " + ")
    )

    updatedmodel   <- update(object,
                             paste(outcome_var, " ~ ", twoint_formula))
    oldmodel       <- suppressWarnings(update(object, method = "ML"))
    oldmodelML_dev <- -2 * summary(oldmodel)$logLik
    newmodel       <- suppressWarnings(update(updatedmodel, method = "ML"))
    newmodelML_dev <- -2 * summary(newmodel)$logLik

    if (round(oldmodelML_dev, 6) == round(newmodelML_dev, 6)) {
      return(updatedmodel)
    } else {
      stop("The deviance values from both models run with ML estimation do ",
           "not match. There was a mistake. Please specify the two-intercept ",
           "approach on your own as you cannot trust these results.")
    }
  }

  # -- glmmTMB ---------------------------------------------------------------
  if (!is.null(object) && inherits(object, "glmmTMB")) {
    fixeff_formula <- formula(object, fixed.only = TRUE)
    outcome_var    <- as.character(fixeff_formula[2])
    all_terms      <- attr(terms(fixeff_formula), "term.labels")
    re_string      <- .get_re_string(object)

    vars_predictor <- .get_predictors(all_terms, disting_variable)

    new_interceptlist         <- paste0("0 + ", disting_level_1, " + ",
                                        disting_level_2)
    new_level_1_predictorlist <- paste0(disting_level_1, ":", vars_predictor)
    new_level_2_predictorlist <- paste0(disting_level_2, ":", vars_predictor)
    twoint_formula <- paste0(
      new_interceptlist, " + ",
      paste(new_level_1_predictorlist, collapse = " + "), " + ",
      paste(new_level_2_predictorlist, collapse = " + ")
    )

    newformula   <- as.formula(paste0(outcome_var, " ~ ", twoint_formula,
                                      " + ", re_string))
    updatedmodel <- update(object, formula. = newformula)

    oldmodel       <- suppressWarnings(update(object, REML = FALSE))
    oldmodelML_dev <- -2 * summary(oldmodel)$logLik[1]
    newmodel       <- suppressWarnings(update(updatedmodel, REML = FALSE))
    newmodelML_dev <- -2 * summary(newmodel)$logLik[1]

    if (round(oldmodelML_dev, 6) == round(newmodelML_dev, 6)) {
      return(updatedmodel)
    } else {
      stop("The deviance values from both models run with ML estimation do ",
           "not match. There was a mistake. Please specify the two-intercept ",
           "approach on your own as you cannot trust these results.")
    }
  }
}
