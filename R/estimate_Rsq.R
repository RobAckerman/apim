#' Estimate and Test R-squared for APIM Models
#'
#' @description
#' Reports the estimate(s) of R-squared and the corresponding deviance
#' difference test(s) for the indistinguishable and distinguishable
#' Actor-Partner Interdependence Model (APIM).
#'
#' If researchers use the \code{nlme} (with the \code{gls} function) or
#' \code{glmmTMB} packages, they need only specify the full model object with
#' the \code{model_full} argument. Either ML or REML estimation can be used
#' for the supplied model; the function handles refitting internally.
#'
#' By default, the function assumes an indistinguishable APIM. If the model
#' is distinguishable, researchers must specify \code{indistinguishable = FALSE}
#' and supply the names of the dummy-coded variables for the two dyad members
#' via \code{person_1} and \code{person_2}.
#'
#' If researchers use other software to estimate the models, they must provide
#' all individual quantities needed to estimate and test R-squared. When
#' supplying these values, researchers must ensure that the deviances are from
#' models estimated with ML and that the residual variances are from models
#' estimated with REML.
#'
#' @param model_full A fitted \code{gls} or \code{glmmTMB} full model object.
#' @param indistinguishable Logical. If \code{TRUE} (default), the function
#'   assumes an indistinguishable APIM. If \code{FALSE}, a distinguishable
#'   APIM is assumed.
#' @param deviance_null Numeric. ML deviance for the null model. Only required
#'   when not supplying model objects.
#' @param deviance_full Numeric. ML deviance for the full model. Only required
#'   when not supplying model objects.
#' @param deviance_p1 Numeric. ML deviance for the person 1 model (i.e., the
#'   model containing only person 1's actor and partner effects). Only required
#'   for the distinguishable case when not supplying model objects.
#' @param deviance_p2 Numeric. ML deviance for the person 2 model. Only
#'   required for the distinguishable case when not supplying model objects.
#' @param p_null Integer. Number of parameters in the null model. Only
#'   required when not supplying model objects.
#' @param p_full Integer. Number of parameters in the full model. Only
#'   required when not supplying model objects.
#' @param p_p1 Integer. Number of parameters in the person 1 model. Only
#'   required for the distinguishable case when not supplying model objects.
#' @param p_p2 Integer. Number of parameters in the person 2 model. Only
#'   required for the distinguishable case when not supplying model objects.
#' @param resvar_null Numeric. REML residual variance from the null model.
#'   Only required for the indistinguishable case when not supplying model
#'   objects.
#' @param resvar_full Numeric. REML residual variance from the full model.
#'   Only required for the indistinguishable case when not supplying model
#'   objects.
#' @param resvar_null_p1 Numeric. REML residual variance for person 1 from
#'   the null model. Only required for the distinguishable case when not
#'   supplying model objects.
#' @param resvar_null_p2 Numeric. REML residual variance for person 2 from
#'   the null model. Only required for the distinguishable case when not
#'   supplying model objects.
#' @param resvar_full_p1 Numeric. REML residual variance for person 1 from
#'   the full model. Only required for the distinguishable case when not
#'   supplying model objects.
#' @param resvar_full_p2 Numeric. REML residual variance for person 2 from
#'   the full model. Only required for the distinguishable case when not
#'   supplying model objects.
#' @param person_1 Character. Name of the dummy-coded variable for the first
#'   dyad member in the dataset (e.g. \code{"man"}). Used as a display label
#'   when supplying raw values, and as a column name when supplying model
#'   objects for the distinguishable case. Defaults to \code{"Person 1"}.
#' @param person_2 Character. Name of the dummy-coded variable for the second
#'   dyad member in the dataset (e.g. \code{"woman"}). Defaults to
#'   \code{"Person 2"}.
#'
#' @return For the indistinguishable case, a character string reporting R-squared
#'   and the deviance difference test. For the distinguishable case, a list of
#'   three character strings: the omnibus test, and individual tests for
#'   person 1 and person 2.
#'
#' @examples
#' \dontrun{
#' # Indistinguishable APIM — supply model object
#' estimate_Rsq(model_full = fullmodel_REML)
#'
#' # Indistinguishable APIM — supply raw values
#' estimate_Rsq(deviance_null = 1472.455, deviance_full = 1341.497,
#'              p_null = 3, p_full = 5,
#'              resvar_null = 0.722, resvar_full = 0.562)
#'
#' # Distinguishable APIM — supply model object
#' estimate_Rsq(model_full = fulldismodel_REML,
#'              indistinguishable = FALSE,
#'              person_1 = "man", person_2 = "woman")
#'
#' # Distinguishable APIM — supply raw values
#' estimate_Rsq(indistinguishable = FALSE,
#'              deviance_null = 1187.010, deviance_full = 1137.756,
#'              deviance_p1 = 1176.509, deviance_p2 = 1169.441,
#'              p_null = 5, p_full = 9, p_p1 = 7, p_p2 = 7,
#'              resvar_null_p1 = 0.429, resvar_null_p2 = 0.534,
#'              resvar_full_p1 = 0.389, resvar_full_p2 = 0.473,
#'              person_1 = "man", person_2 = "woman")
#' }
#' @importFrom stats pchisq formula update coef
#' @export
estimate_Rsq <- function(model_full        = NULL,
                         indistinguishable = TRUE,
                         deviance_null     = NULL,
                         deviance_full     = NULL,
                         deviance_p1       = NULL,
                         deviance_p2       = NULL,
                         p_null            = NULL,
                         p_full            = NULL,
                         p_p1              = NULL,
                         p_p2              = NULL,
                         resvar_null       = NULL,
                         resvar_full       = NULL,
                         resvar_null_p1    = NULL,
                         resvar_null_p2    = NULL,
                         resvar_full_p1    = NULL,
                         resvar_full_p2    = NULL,
                         person_1          = "Person 1",
                         person_2          = "Person 2") {

  chisquare_symbol <- paste0("\u03C7", "\u00B2")
  rsquare_symbol   <- paste0("R", "\u00B2")

  # ---------------------------------------------------------------------------
  # Helper: extract random-effects string from a glmmTMB model using the
  # fixeff_string subtraction approach, which handles complex structures
  # like cs(), homcs(), indisting(), and multiple RE terms correctly.
  # ---------------------------------------------------------------------------
  .get_re_string <- function(object) {
    formula_string <- gsub("\\s+", " ",
                           paste(deparse(formula(object)), collapse = " "))
    rhs <- trimws(sub("^[^~]*~\\s*", "", formula_string))
    fixeff_string <- gsub("\\s+", " ",
                          paste(deparse(formula(object, fixed.only = TRUE)),
                                collapse = " "))
    fixeff_rhs <- trimws(sub("^[^~]*~\\s*", "", fixeff_string))
    re_start <- nchar(fixeff_rhs) + 1L
    trimws(sub("^\\s*\\+\\s*", "", substr(rhs, re_start, nchar(rhs))))
  }

  # ===========================================================================
  # INDISTINGUISHABLE APIM
  # ===========================================================================
  if (indistinguishable == TRUE) {

    # --- Raw values supplied ---
    if (!is.null(deviance_null) && !is.null(deviance_full) &&
        !is.null(p_null) && !is.null(p_full) &&
        !is.null(resvar_null) && !is.null(resvar_full)) {

      if (p_null >= p_full)
        stop("Number of parameters for full model must be greater than ",
             "number of parameters for null model.", call. = FALSE)

      model_null_sigma2 <- resvar_null
      model_full_sigma2 <- resvar_full
      R2       <- 1 - (model_full_sigma2 / model_null_sigma2)
      chiR2    <- deviance_null - deviance_full
      dfR2     <- p_full - p_null
      pvalueR2 <- pchisq(q = chiR2, df = dfR2, lower.tail = FALSE)
    }

    # --- gls object supplied ---
    if (!is.null(model_full) && inherits(model_full, "gls")) {

      fixeff_formula <- formula(model_full, fixed.only = TRUE)
      outcome_var    <- as.character(fixeff_formula[2])
      model_null     <- suppressWarnings(
        update(model_full, paste(outcome_var, "~ 1"))
      )

      p_full_val <- model_full$dims$p + nrow(model_full$apVar)
      p_null_val <- model_null$dims$p + nrow(model_null$apVar)

      if (p_null_val >= p_full_val)
        stop("Number of parameters for full model must be greater than ",
             "number of parameters for null model.", call. = FALSE)

      model_null_REML <- suppressWarnings(update(model_null, method = "REML"))
      model_full_REML <- suppressWarnings(update(model_full, method = "REML"))
      model_null_ML   <- suppressWarnings(update(model_null, method = "ML"))
      model_full_ML   <- suppressWarnings(update(model_full, method = "ML"))

      deviance_null_val  <- -2 * model_null_ML$logLik
      model_null_sigma2  <- model_null_REML$sigma^2
      deviance_full_val  <- -2 * model_full_ML$logLik
      model_full_sigma2  <- model_full_REML$sigma^2

      R2       <- 1 - (model_full_sigma2 / model_null_sigma2)
      chiR2    <- deviance_null_val - deviance_full_val
      dfR2     <- p_full_val - p_null_val
      pvalueR2 <- pchisq(q = chiR2, df = dfR2, lower.tail = FALSE)
    }

    # --- glmmTMB object supplied ---
    if (!is.null(model_full) && inherits(model_full, "glmmTMB")) {

      grouping_variable <- model_full$modelInfo$grpVar
      fixeff_formula    <- formula(model_full, fixed.only = TRUE)
      outcome_var       <- as.character(fixeff_formula[2])
      re_string         <- .get_re_string(model_full)

      model_null <- suppressWarnings(
        update(model_full,
               as.formula(paste(outcome_var, "~ 1 +", re_string)))
      )

      model_null_REML <- suppressWarnings(update(model_null, REML = TRUE))
      model_full_REML <- suppressWarnings(update(model_full, REML = TRUE))
      model_null_ML   <- suppressWarnings(update(model_null, REML = FALSE))
      model_full_ML   <- suppressWarnings(update(model_full, REML = FALSE))

      p_full_val        <- length(model_full_ML$obj$par)
      p_null_val        <- length(model_null_ML$obj$par)
      deviance_null_val <- -2 * (-(model_null_ML$obj$fn()[[1]]))
      deviance_full_val <- -2 * (-(model_full_ML$obj$fn()[[1]]))
      model_null_sigma2 <- summary(model_null_REML)$varcor$cond[[grouping_variable]][1]
      model_full_sigma2 <- summary(model_full_REML)$varcor$cond[[grouping_variable]][1]

      if (p_null_val >= p_full_val)
        stop("Number of parameters for full model must be greater than ",
             "number of parameters for null model.", call. = FALSE)

      R2       <- 1 - (model_full_sigma2 / model_null_sigma2)
      chiR2    <- deviance_null_val - deviance_full_val
      dfR2     <- p_full_val - p_null_val
      pvalueR2 <- pchisq(q = chiR2, df = dfR2, lower.tail = FALSE)
    }

    text_results <- paste0(
      "Multiple Correlation Coefficient and Test: ",
      rsquare_symbol, " = ", round(R2, 3), ", ",
      chisquare_symbol, "(", dfR2, ") = ", round(chiR2, 3),
      ", p = ", round(pvalueR2, 5)
    )

  # ===========================================================================
  # DISTINGUISHABLE APIM
  # ===========================================================================
  } else {

    # --- Raw values supplied ---
    if (!is.null(deviance_null) && !is.null(deviance_full) &&
        !is.null(deviance_p1) && !is.null(deviance_p2) &&
        !is.null(p_null) && !is.null(p_full) &&
        !is.null(p_p1) && !is.null(p_p2) &&
        !is.null(resvar_null_p1) && !is.null(resvar_full_p1) &&
        !is.null(resvar_null_p2) && !is.null(resvar_full_p2)) {

      if (p_null >= p_full)
        stop("Number of parameters for full model must be greater than ",
             "number of parameters for null model.", call. = FALSE)

      omnibus_chiR2    <- deviance_null - deviance_full
      omnibus_dfR2     <- p_full - p_null
      omnibus_pvalueR2 <- pchisq(q = omnibus_chiR2, df = omnibus_dfR2,
                                 lower.tail = FALSE)
      p1_R2       <- 1 - (resvar_full_p1 / resvar_null_p1)
      p1_chiR2    <- deviance_p2 - deviance_full
      p1_dfR2     <- p_full - p_p2
      p1_pvalueR2 <- pchisq(q = p1_chiR2, df = p1_dfR2, lower.tail = FALSE)

      p2_R2       <- 1 - (resvar_full_p2 / resvar_null_p2)
      p2_chiR2    <- deviance_p1 - deviance_full
      p2_dfR2     <- p_full - p_p1
      p2_pvalueR2 <- pchisq(q = p2_chiR2, df = p2_dfR2, lower.tail = FALSE)
    }

    # --- gls object supplied ---
    if (!is.null(model_full) && inherits(model_full, "gls")) {

      fixeff_formula <- formula(model_full, fixed.only = TRUE)
      outcome_var    <- as.character(fixeff_formula[2])
      extracted_vars <- all.vars(fixeff_formula)

      # Identify person dummy codes and predictors using the same two-step
      # filter as the switch functions: first isolate the two level variables
      # (person_1, person_2), then extract everything else as predictors.
      omit_vars_person    <- c(":", "_A$", "_P$", outcome_var)
      vars_person         <- extracted_vars[!grepl(
        paste(omit_vars_person, collapse = "|"), extracted_vars)]
      person_1_col        <- vars_person[[1]]
      person_2_col        <- vars_person[[2]]
      omit_vars_predictor <- c(":", person_1_col, person_2_col, outcome_var)
      vars_predictor      <- extracted_vars[!grepl(
        paste(omit_vars_predictor, collapse = "|"), extracted_vars)]

      # Null, person 1, and person 2 models
      model_null <- suppressWarnings(update(model_full,
        paste(outcome_var, "~ 0 +", person_1_col, "+", person_2_col)))

      newpredictorlist_p1 <- paste0(person_1_col, ":", vars_predictor)
      model_p1 <- suppressWarnings(update(model_null,
        paste(". ~ . +", paste(newpredictorlist_p1, collapse = " + "))))

      newpredictorlist_p2 <- paste0(person_2_col, ":", vars_predictor)
      model_p2 <- suppressWarnings(update(model_null,
        paste(". ~ . +", paste(newpredictorlist_p2, collapse = " + "))))

      # REML and ML versions
      model_null_REML <- suppressWarnings(update(model_null, method = "REML"))
      model_full_REML <- suppressWarnings(update(model_full, method = "REML"))
      model_null_ML   <- suppressWarnings(update(model_null, method = "ML"))
      model_full_ML   <- suppressWarnings(update(model_full, method = "ML"))
      model_p1_ML     <- suppressWarnings(update(model_p1,   method = "ML"))
      model_p2_ML     <- suppressWarnings(update(model_p2,   method = "ML"))

      p_full_val <- model_full$dims$p + nrow(model_full$apVar)
      p_null_val <- model_null$dims$p + nrow(model_null$apVar)
      p_p1_val   <- model_p1$dims$p   + nrow(model_p1$apVar)
      p_p2_val   <- model_p2$dims$p   + nrow(model_p2$apVar)

      p1_model_null_sigma2 <- (summary(model_null_REML)$sigma * 1.00)^2
      p2_model_null_sigma2 <- as.numeric(
        (summary(model_null_REML)$sigma *
           coef(model_null_REML$modelStruct$varStruct, uncons = FALSE))^2)
      p1_model_full_sigma2 <- (summary(model_full_REML)$sigma * 1.00)^2
      p2_model_full_sigma2 <- as.numeric(
        (summary(model_full_REML)$sigma *
           coef(model_full_REML$modelStruct$varStruct, uncons = FALSE))^2)

      p1_R2    <- 1 - (p1_model_full_sigma2 / p1_model_null_sigma2)
      p2_R2    <- 1 - (p2_model_full_sigma2 / p2_model_null_sigma2)

      omnibus_chiR2    <- (-2 * model_null_ML$logLik) - (-2 * model_full_ML$logLik)
      omnibus_dfR2     <- p_full_val - p_null_val
      omnibus_pvalueR2 <- pchisq(q = omnibus_chiR2, df = omnibus_dfR2,
                                 lower.tail = FALSE)
      p1_chiR2    <- (-2 * model_p2_ML$logLik) - (-2 * model_full_ML$logLik)
      p2_chiR2    <- (-2 * model_p1_ML$logLik) - (-2 * model_full_ML$logLik)
      p1_dfR2     <- p_full_val - p_p2_val
      p2_dfR2     <- p_full_val - p_p1_val
      p1_pvalueR2 <- pchisq(q = p1_chiR2, df = p1_dfR2, lower.tail = FALSE)
      p2_pvalueR2 <- pchisq(q = p2_chiR2, df = p2_dfR2, lower.tail = FALSE)

      person_1 <- person_1_col
      person_2 <- person_2_col
    }

    # --- glmmTMB object supplied ---
    if (!is.null(model_full) && inherits(model_full, "glmmTMB")) {

      grouping_variable <- model_full$modelInfo$grpVar
      fixeff_formula    <- formula(model_full, fixed.only = TRUE)
      outcome_var       <- as.character(fixeff_formula[2])
      extracted_vars    <- all.vars(fixeff_formula)
      re_string         <- .get_re_string(model_full)

      # Identify person dummy codes and predictors
      omit_vars_person    <- c(":", "_A$", "_P$", outcome_var)
      vars_person         <- extracted_vars[!grepl(
        paste(omit_vars_person, collapse = "|"), extracted_vars)]
      person_1_col        <- vars_person[[1]]
      person_2_col        <- vars_person[[2]]
      omit_vars_predictor <- c(":", person_1_col, person_2_col, outcome_var)
      vars_predictor      <- extracted_vars[!grepl(
        paste(omit_vars_predictor, collapse = "|"), extracted_vars)]

      # Null model
      model_null <- suppressWarnings(update(model_full,
        as.formula(paste(outcome_var, "~ 0 +",
                         person_1_col, "+", person_2_col, "+", re_string))))

      # Person 1 model — built from null using . ~ . + to avoid RE duplication
      newpredictorlist_p1 <- paste0(person_1_col, ":", vars_predictor)
      model_p1 <- suppressWarnings(update(model_null,
        paste(". ~ . +", paste(newpredictorlist_p1, collapse = " + "))))

      # Person 2 model
      newpredictorlist_p2 <- paste0(person_2_col, ":", vars_predictor)
      model_p2 <- suppressWarnings(update(model_null,
        paste(". ~ . +", paste(newpredictorlist_p2, collapse = " + "))))

      # REML and ML versions
      model_null_REML <- suppressWarnings(update(model_null, REML = TRUE))
      model_full_REML <- suppressWarnings(update(model_full, REML = TRUE))
      model_null_ML   <- suppressWarnings(update(model_null, REML = FALSE))
      model_full_ML   <- suppressWarnings(update(model_full, REML = FALSE))
      model_p1_ML     <- suppressWarnings(update(model_p1,   REML = FALSE))
      model_p2_ML     <- suppressWarnings(update(model_p2,   REML = FALSE))

      p_full_val <- length(model_full_ML$obj$par)
      p_null_val <- length(model_null_ML$obj$par)
      p_p1_val   <- length(model_p1_ML$obj$par)
      p_p2_val   <- length(model_p2_ML$obj$par)

      p1_model_null_sigma2 <- summary(model_null_REML)$varcor$cond[[grouping_variable]][person_1_col, person_1_col]
      p2_model_null_sigma2 <- summary(model_null_REML)$varcor$cond[[grouping_variable]][person_2_col, person_2_col]
      p1_model_full_sigma2 <- summary(model_full_REML)$varcor$cond[[grouping_variable]][person_1_col, person_1_col]
      p2_model_full_sigma2 <- summary(model_full_REML)$varcor$cond[[grouping_variable]][person_2_col, person_2_col]

      p1_R2 <- 1 - (p1_model_full_sigma2 / p1_model_null_sigma2)
      p2_R2 <- 1 - (p2_model_full_sigma2 / p2_model_null_sigma2)

      omnibus_chiR2    <- (-2 * (summary(model_null_ML)$logLik[1])) -
                          (-2 * (summary(model_full_ML)$logLik[1]))
      omnibus_dfR2     <- p_full_val - p_null_val
      omnibus_pvalueR2 <- pchisq(q = omnibus_chiR2, df = omnibus_dfR2,
                                 lower.tail = FALSE)
      p1_chiR2    <- (-2 * (summary(model_p2_ML)$logLik[1])) -
                     (-2 * (summary(model_full_ML)$logLik[1]))
      p2_chiR2    <- (-2 * (summary(model_p1_ML)$logLik[1])) -
                     (-2 * (summary(model_full_ML)$logLik[1]))
      p1_dfR2     <- p_full_val - p_p2_val
      p2_dfR2     <- p_full_val - p_p1_val
      p1_pvalueR2 <- pchisq(q = p1_chiR2, df = p1_dfR2, lower.tail = FALSE)
      p2_pvalueR2 <- pchisq(q = p2_chiR2, df = p2_dfR2, lower.tail = FALSE)

      person_1 <- person_1_col
      person_2 <- person_2_col
    }

    omnibus_text <- paste0(
      "Omnibus Test of Multiple Correlation Coefficients for ",
      person_1, " and ", person_2, ": ",
      chisquare_symbol, "(", omnibus_dfR2, ") = ",
      round(omnibus_chiR2, 3), ", p = ", round(omnibus_pvalueR2, 5))

    p1_text <- paste0(
      "Multiple Correlation Coefficient and Test for ", person_1, ": ",
      rsquare_symbol, " = ", round(p1_R2, 3), ", ",
      chisquare_symbol, "(", p1_dfR2, ") = ",
      round(p1_chiR2, 3), ", p = ", round(p1_pvalueR2, 5))

    p2_text <- paste0(
      "Multiple Correlation Coefficient and Test for ", person_2, ": ",
      rsquare_symbol, " = ", round(p2_R2, 3), ", ",
      chisquare_symbol, "(", p2_dfR2, ") = ",
      round(p2_chiR2, 3), ", p = ", round(p2_pvalueR2, 5))

    text_results <- list(omnibus_text, p1_text, p2_text)
  }

  return(text_results)
}
