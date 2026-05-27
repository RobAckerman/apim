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
#' For \code{gls} models with \code{indistinguishable = FALSE}, the model
#' must be in the two-intercept form (e.g.
#' \code{RelSat_A ~ 0 + man + woman + man:c_Amity_A + woman:c_Amity_A}).
#' If your model uses the interaction approach, use \code{switch_to_twoint()}
#' first. The \code{person_1} argument is assumed to correspond to the
#' reference level of \code{varIdent} (the level with ratio = 1.0, determined
#' by the first value of the grouping variable encountered in the data) and
#' \code{person_2} to the non-reference level. A note is printed at runtime
#' showing which factor levels correspond to \code{person_1} and
#' \code{person_2} so you can verify the mapping is correct.
#'
#' For \code{glmmTMB} models with \code{indistinguishable = FALSE}, the
#' random effects term must use explicit dummy-coded variables
#' (e.g. \code{cs(0 + man + woman | DyadID)}) rather than \code{as.factor()}.
#' The names supplied to \code{person_1} and \code{person_2} must match the
#' dummy variable names exactly.
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
#' @param deviance_p1 Numeric. ML deviance for the person 1 model. Only
#'   required for the distinguishable case when not supplying model objects.
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
#' @param person_1 Character. Display label and variable name for the first
#'   dyad member. For \code{glmmTMB} models, must match the dummy-coded
#'   variable name in the random effects term exactly (e.g. \code{"man"}).
#'   For \code{gls} models, assumed to correspond to the reference level of
#'   \code{varIdent} (ratio = 1.0). Defaults to \code{"Person 1"}.
#' @param person_2 Character. Display label and variable name for the second
#'   dyad member. For \code{glmmTMB} models, must match the dummy-coded
#'   variable name exactly. For \code{gls} models, assumed to correspond to
#'   the non-reference level of \code{varIdent}. Defaults to \code{"Person 2"}.
#'
#' @return Invisibly returns a \code{data.frame} with columns \code{person},
#'   \code{R2}, \code{chi2}, \code{df}, and \code{p}. For the
#'   indistinguishable case the data frame has one row; for the distinguishable
#'   case it has three rows (omnibus, person 1, person 2).
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
#' # Distinguishable APIM — gls, must be in two-intercept form
#' m_twoint <- switch_to_twoint(m_interact,
#'                              disting_variable = "ECGender_A",
#'                              disting_level_1  = "man",
#'                              disting_level_2  = "woman")
#' estimate_Rsq(m_twoint,
#'              indistinguishable = FALSE,
#'              person_1 = "man", person_2 = "woman")
#'
#' # Distinguishable APIM — glmmTMB, RE must use explicit dummy codes
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
#' @importFrom stats pchisq formula update coef terms
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

  # ---------------------------------------------------------------------------
  # Helper: format p-value
  # ---------------------------------------------------------------------------
  fmt_p <- function(p) {
    if (p < .0001) "< .0001" else sprintf("= %.4f", p)
  }

  # ---------------------------------------------------------------------------
  # Helper: extract RE string using bracket-matching
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

  # ---------------------------------------------------------------------------
  # Helper: push original data to .GlobalEnv
  # ---------------------------------------------------------------------------
  .push_data <- function(model) {
    original_data <- tryCatch(
      eval(model$call$data, envir = .GlobalEnv),
      error = function(e) model$frame
    )
    assign(".estimate_rsq_data_tmp", original_data, envir = .GlobalEnv)
    model$call$data <- quote(.estimate_rsq_data_tmp)
    model
  }

  # ---------------------------------------------------------------------------
  # Helper: silence stderr messages from glmmTMB
  # ---------------------------------------------------------------------------
  quietly <- function(expr) {
    tmp <- textConnection("msg_sink", open = "w", local = TRUE)
    sink(tmp, type = "message")
    on.exit({
      sink(type = "message")
      close(tmp)
    })
    suppressWarnings(suppressMessages(expr))
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
        update(model_full, paste(outcome_var, "~ 1")))

      p_full_val <- model_full$dims$p + nrow(model_full$apVar)
      p_null_val <- model_null$dims$p + nrow(model_null$apVar)

      if (p_null_val >= p_full_val)
        stop("Number of parameters for full model must be greater than ",
             "number of parameters for null model.", call. = FALSE)

      model_null_REML <- suppressWarnings(update(model_null, method = "REML"))
      model_full_REML <- suppressWarnings(update(model_full, method = "REML"))
      model_null_ML   <- suppressWarnings(update(model_null, method = "ML"))
      model_full_ML   <- suppressWarnings(update(model_full, method = "ML"))

      deviance_null_val <- -2 * model_null_ML$logLik
      model_null_sigma2 <- model_null_REML$sigma^2
      deviance_full_val <- -2 * model_full_ML$logLik
      model_full_sigma2 <- model_full_REML$sigma^2

      R2       <- 1 - (model_full_sigma2 / model_null_sigma2)
      chiR2    <- deviance_null_val - deviance_full_val
      dfR2     <- p_full_val - p_null_val
      pvalueR2 <- pchisq(q = chiR2, df = dfR2, lower.tail = FALSE)
    }

    # --- glmmTMB object supplied ---
    if (!is.null(model_full) && inherits(model_full, "glmmTMB")) {

      model_full <- .push_data(model_full)
      on.exit(
        try(rm(".estimate_rsq_data_tmp", envir = .GlobalEnv), silent = TRUE),
        add = TRUE
      )

      grouping_variable <- model_full$modelInfo$grpVar
      fixeff_formula    <- formula(model_full, fixed.only = TRUE)
      outcome_var       <- as.character(fixeff_formula[2])
      re_string         <- .get_re_string(model_full)

      model_null <- suppressWarnings(
        update(model_full,
               as.formula(paste(outcome_var, "~ 1 +", re_string))))

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

    # -- print -----------------------------------------------------------------
    cat("\n=== R-squared and Model Test ===\n")
    cat(sprintf("  R\u00b2 = %.3f\n", R2))
    cat(sprintf("  \u03c7\u00b2(%d) = %.3f, p %s\n\n", dfR2, chiR2,
                fmt_p(pvalueR2)))

    invisible(data.frame(
      person = "Overall",
      R2     = round(R2, 3),
      chi2   = round(chiR2, 3),
      df     = dfR2,
      p      = round(pvalueR2, 5),
      row.names = NULL
    ))

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

      # -- check model is in two-intercept form --------------------------------
      if (!person_1 %in% extracted_vars || !person_2 %in% extracted_vars) {
        stop(
          "For the distinguishable gls case, estimate_Rsq() requires the ",
          "model to be in two-intercept form\n",
          "  (e.g. RelSat_A ~ 0 + man + woman + man:c_Amity_A + ",
          "woman:c_Amity_A + ...).\n",
          "  '", person_1, "' and/or '", person_2, "' were not found among ",
          "the model's fixed-effect terms.",
          call. = FALSE
        )
      }

      # -- note about varIdent reference level ---------------------------------
      model_data <- tryCatch(
        eval(model_full$call$data, envir = .GlobalEnv),
        error = function(e) model_full$data
      )
      var_struct_coefs <- coef(model_full$modelStruct$varStruct,
                               uncons = FALSE)
      weights_var <- tryCatch({
        wf <- as.character(model_full$call$weights[[2]])
        gsub("[[:space:]~]", "", sub(".*\\|", "", wf[2]))
      }, error = function(e) "distinguishing variable")

      all_levels <- tryCatch(
        as.character(sort(unique(model_data[[weights_var]]))),
        error = function(e) NULL
      )
      ref_level    <- if (!is.null(all_levels))
        all_levels[!all_levels %in% names(var_struct_coefs)] else "unknown"
      nonref_level <- names(var_struct_coefs)
      first_val    <- tryCatch(
        as.character(model_data[[weights_var]][1]),
        error = function(e) "unknown"
      )

      message(
        "Note: For gls models, person_1 ('", person_1, "') is assumed to ",
        "correspond to the reference level of varIdent\n",
        "  (", weights_var, " = ", paste(ref_level, collapse = "/"), ") and ",
        "person_2 ('", person_2, "') is assumed to correspond to the ",
        "non-reference level\n",
        "  (", weights_var, " = ", paste(nonref_level, collapse = "/"), ").\n",
        "  The reference level is the first value of ", weights_var,
        " encountered in the data (row 1 = ", first_val, ").\n",
        "  Verify this matches your model specification before interpreting ",
        "the person-specific R\u00b2 values."
      )

      # -- use person_1 and person_2 directly ----------------------------------
      person_1_col <- person_1
      person_2_col <- person_2

      # extract predictors using term labels — robust to any variable naming
      all_terms       <- attr(terms(fixeff_formula), "term.labels")
      vars_with_colon <- all_terms[grepl(":", all_terms)]
      vars_predictor  <- unique(sapply(vars_with_colon, function(term) {
        parts      <- strsplit(term, ":")[[1]]
        non_person <- parts[!parts %in% c(person_1_col, person_2_col)]
        if (length(non_person) == 1) non_person else NA_character_
      }))
      vars_predictor <- vars_predictor[!is.na(vars_predictor)]

      model_null <- suppressWarnings(update(model_full,
                                            paste(outcome_var, "~ 0 +", person_1_col, "+", person_2_col)))
      newpredictorlist_p1 <- paste0(person_1_col, ":", vars_predictor)
      model_p1 <- suppressWarnings(update(model_null,
                                          paste(". ~ . +", paste(newpredictorlist_p1, collapse = " + "))))
      newpredictorlist_p2 <- paste0(person_2_col, ":", vars_predictor)
      model_p2 <- suppressWarnings(update(model_null,
                                          paste(". ~ . +", paste(newpredictorlist_p2, collapse = " + "))))

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

      p1_R2 <- 1 - (p1_model_full_sigma2 / p1_model_null_sigma2)
      p2_R2 <- 1 - (p2_model_full_sigma2 / p2_model_null_sigma2)

      omnibus_chiR2    <- (-2 * model_null_ML$logLik) -
        (-2 * model_full_ML$logLik)
      omnibus_dfR2     <- p_full_val - p_null_val
      omnibus_pvalueR2 <- pchisq(q = omnibus_chiR2, df = omnibus_dfR2,
                                 lower.tail = FALSE)
      p1_chiR2    <- (-2 * model_p2_ML$logLik) - (-2 * model_full_ML$logLik)
      p2_chiR2    <- (-2 * model_p1_ML$logLik) - (-2 * model_full_ML$logLik)
      p1_dfR2     <- p_full_val - p_p2_val
      p2_dfR2     <- p_full_val - p_p1_val
      p1_pvalueR2 <- pchisq(q = p1_chiR2, df = p1_dfR2, lower.tail = FALSE)
      p2_pvalueR2 <- pchisq(q = p2_chiR2, df = p2_dfR2, lower.tail = FALSE)
    }

    # --- glmmTMB object supplied ---
    if (!is.null(model_full) && inherits(model_full, "glmmTMB")) {

      model_full <- .push_data(model_full)
      on.exit(
        try(rm(".estimate_rsq_data_tmp", envir = .GlobalEnv), silent = TRUE),
        add = TRUE
      )

      grouping_variable <- model_full$modelInfo$grpVar
      fixeff_formula    <- formula(model_full, fixed.only = TRUE)
      outcome_var       <- as.character(fixeff_formula[2])
      extracted_vars    <- all.vars(fixeff_formula)
      re_string         <- .get_re_string(model_full)
      dispformula       <- if (!is.null(model_full$call$dispformula)) {
        eval(model_full$call$dispformula)
      } else { ~1 }

      # use person_1 and person_2 directly
      person_1_col <- person_1
      person_2_col <- person_2

      omit_vars_predictor <- c(":", person_1_col, person_2_col, outcome_var)
      vars_predictor      <- extracted_vars[!grepl(
        paste(omit_vars_predictor, collapse = "|"), extracted_vars)]

      # build formulas explicitly
      f_null <- as.formula(paste(
        outcome_var, "~ 0 +", person_1_col, "+", person_2_col, "+", re_string
      ))
      f_full <- as.formula(paste(
        outcome_var, "~ 0 +", person_1_col, "+", person_2_col, "+", re_string,
        "+", paste(paste0(person_1_col, ":", vars_predictor), collapse = " + "),
        "+", paste(paste0(person_2_col, ":", vars_predictor), collapse = " + ")
      ))
      f_p1 <- as.formula(paste(
        outcome_var, "~ 0 +", person_1_col, "+", person_2_col, "+", re_string,
        "+", paste(paste0(person_1_col, ":", vars_predictor), collapse = " + ")
      ))
      f_p2 <- as.formula(paste(
        outcome_var, "~ 0 +", person_1_col, "+", person_2_col, "+", re_string,
        "+", paste(paste0(person_2_col, ":", vars_predictor), collapse = " + ")
      ))

      # fit all models directly
      model_null_REML <- quietly(glmmTMB::glmmTMB(
        f_null, dispformula = dispformula, REML = TRUE,
        data = .estimate_rsq_data_tmp))
      model_full_REML <- quietly(glmmTMB::glmmTMB(
        f_full, dispformula = dispformula, REML = TRUE,
        data = .estimate_rsq_data_tmp))
      model_null_ML <- quietly(glmmTMB::glmmTMB(
        f_null, dispformula = dispformula, REML = FALSE,
        data = .estimate_rsq_data_tmp))
      model_full_ML <- quietly(glmmTMB::glmmTMB(
        f_full, dispformula = dispformula, REML = FALSE,
        data = .estimate_rsq_data_tmp))
      model_p1_ML <- quietly(glmmTMB::glmmTMB(
        f_p1, dispformula = dispformula, REML = FALSE,
        data = .estimate_rsq_data_tmp))
      model_p2_ML <- quietly(glmmTMB::glmmTMB(
        f_p2, dispformula = dispformula, REML = FALSE,
        data = .estimate_rsq_data_tmp))

      p_full_val <- length(model_full_ML$obj$par)
      p_null_val <- length(model_null_ML$obj$par)
      p_p1_val   <- length(model_p1_ML$obj$par)
      p_p2_val   <- length(model_p2_ML$obj$par)

      # check person names exist in RE matrix
      re_names <- rownames(
        summary(model_null_REML)$varcor$cond[[grouping_variable]]
      )
      if (!person_1_col %in% re_names || !person_2_col %in% re_names) {
        stop(
          "Could not find '", person_1_col, "' and/or '", person_2_col,
          "' in the random effects variance-covariance matrix.\n",
          "  Found: ", paste(re_names, collapse = ", "), "\n",
          "  For the distinguishable case, the model must use explicit ",
          "dummy-coded variables in the random effects term\n",
          "  (e.g. cs(0 + man + woman | DyadID)) rather than as.factor().\n",
          "  The names supplied to person_1 and person_2 must match the ",
          "dummy variable names exactly.",
          call. = FALSE
        )
      }

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
    }

    # -- print -----------------------------------------------------------------
    cat("\n=== R-squared and Model Tests ===\n\n")
    cat(sprintf("  Omnibus (%s and %s)\n", person_1, person_2))
    cat(sprintf("    \u03c7\u00b2(%d) = %.3f, p %s\n\n",
                omnibus_dfR2, omnibus_chiR2, fmt_p(omnibus_pvalueR2)))
    cat(sprintf("  %s\n", person_1))
    cat(sprintf("    R\u00b2 = %.3f, \u03c7\u00b2(%d) = %.3f, p %s\n\n",
                p1_R2, p1_dfR2, p1_chiR2, fmt_p(p1_pvalueR2)))
    cat(sprintf("  %s\n", person_2))
    cat(sprintf("    R\u00b2 = %.3f, \u03c7\u00b2(%d) = %.3f, p %s\n\n",
                p2_R2, p2_dfR2, p2_chiR2, fmt_p(p2_pvalueR2)))

    invisible(data.frame(
      person = c("Omnibus", person_1, person_2),
      R2     = c(NA, round(p1_R2, 3), round(p2_R2, 3)),
      chi2   = round(c(omnibus_chiR2, p1_chiR2, p2_chiR2), 3),
      df     = c(omnibus_dfR2, p1_dfR2, p2_dfR2),
      p      = round(c(omnibus_pvalueR2, p1_pvalueR2, p2_pvalueR2), 5),
      row.names = NULL
    ))
  }
}
