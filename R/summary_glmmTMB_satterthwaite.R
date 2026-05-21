#' @keywords internal
.satterthwaite_contrast <- function(model, L, eps = 1e-3,
                                    n_cores = parallel::detectCores(logical = FALSE),
                                    verbose = FALSE) {

  coef_names <- names(fixef(model)$cond)

  # build full contrast vector aligned to all fixed effects
  L_full <- setNames(rep(0, length(coef_names)), coef_names)
  for(nm in names(L)) {
    if(!nm %in% coef_names)
      stop("Contrast term '", nm, "' not found in model coefficients.")
    L_full[nm] <- L[nm]
  }
  L_full <- as.numeric(L_full)

  # extract model quantities
  X         <- model.matrix(model)
  Z_full    <- glmmTMB::getME(model, "Z")
  theta_opt <- glmmTMB::getME(model, "theta")
  n_obs     <- nrow(X)
  p_theta   <- length(theta_opt)

  re_struc  <- model$modelInfo$reStruc$cond
  z_widths  <- vapply(re_struc, function(x) x$blockSize * x$blockReps, numeric(1))
  gp_manual <- c(0L, cumsum(z_widths))

  # asymptotic covariance of theta
  V_full         <- vcov(model, full = TRUE)
  nuisance_idx   <- grep("^theta", colnames(V_full))
  Sigma_nuisance <- V_full[nuisance_idx, nuisance_idx, drop = FALSE]
  Sigma_nuisance[is.na(Sigma_nuisance)] <- 0

  n_cores <- max(1L, min(as.integer(n_cores), p_theta))

  # refit function -- returns log variance of contrast L'beta
  make_refit_fn <- function(model_obj, X_, Z_full_, gp_manual_, n_obs_, L_) {
    function(t_new) {
      tmp_mod <- tryCatch(
        withCallingHandlers(
          update(model_obj,
                 start = list(theta = t_new),
                 map   = list(theta = factor(rep(NA_integer_, length(t_new))))),
          warning = function(w) invokeRestart("muffleWarning")
        ),
        error = function(e) stop("Refit failed: ", conditionMessage(e))
      )

      vc_new <- VarCorr(tmp_mod)$cond
      V_tmp  <- Matrix::Matrix(0, n_obs_, n_obs_, sparse = TRUE)

      for(i in seq_along(vc_new)) {
        Sigma_i <- vc_new[[i]]
        Z_i     <- Z_full_[, (gp_manual_[i] + 1L):gp_manual_[i + 1L], drop = FALSE]
        n_grps  <- ncol(Z_i) / ncol(Sigma_i)
        V_tmp   <- V_tmp +
          Z_i %*% Matrix::kronecker(Matrix::Diagonal(n_grps), Sigma_i) %*% Matrix::t(Z_i)
      }

      diag_scale <- mean(Matrix::diag(V_tmp))
      if(!is.finite(diag_scale) || diag_scale <= 0) diag_scale <- 1
      V_tmp <- V_tmp + Matrix::Diagonal(n_obs_, diag_scale * 1e-8)

      XtVinvX <- tryCatch(
        as.matrix(Matrix::t(X_) %*% Matrix::solve(V_tmp, X_)),
        error = function(e) stop("Sparse solve failed: ", conditionMessage(e))
      )

      var_contrast <- as.numeric(t(L_) %*% solve(XtVinvX) %*% L_)
      log(var_contrast)
    }
  }

  refit_fn <- make_refit_fn(model, X, Z_full, gp_manual, n_obs, L_full)
  h        <- eps * pmax(abs(theta_opt), 1e-4)

  # central difference Jacobian
  if(n_cores == 1L) {
    cols <- lapply(seq_len(p_theta), function(j) {
      x_fwd <- theta_opt; x_fwd[j] <- theta_opt[j] + h[j]
      x_bck <- theta_opt; x_bck[j] <- theta_opt[j] - h[j]
      (refit_fn(x_fwd) - refit_fn(x_bck)) / (2 * h[j])
    })
  } else {
    cl <- parallel::makeCluster(n_cores)
    tryCatch({
      parallel::clusterEvalQ(cl, { library(glmmTMB); library(Matrix) })
      parallel::clusterExport(
        cl,
        c("model", "X", "Z_full", "gp_manual", "n_obs", "L_full",
          "theta_opt", "h", "make_refit_fn"),
        envir = environment()
      )
      parallel::clusterEvalQ(cl,
        refit_fn <- make_refit_fn(model, X, Z_full, gp_manual, n_obs, L_full)
      )
      data_name <- tryCatch(as.character(model$call$data), error = function(e) NULL)
      if(!is.null(data_name) && length(data_name) == 1 &&
         exists(data_name, envir = .GlobalEnv))
        parallel::clusterExport(cl, data_name, envir = .GlobalEnv)

      cols <- parallel::parLapply(cl, seq_len(p_theta), function(j) {
        x_fwd <- theta_opt; x_fwd[j] <- theta_opt[j] + h[j]
        x_bck <- theta_opt; x_bck[j] <- theta_opt[j] - h[j]
        (refit_fn(x_fwd) - refit_fn(x_bck)) / (2 * h[j])
      })
    }, finally = { parallel::stopCluster(cl) })
  }

  g_log <- as.numeric(unlist(cols))
  v_den <- sum(g_log * (Sigma_nuisance %*% g_log))

  if(v_den <= 0) {
    warning("Non-positive df denominator -- check model convergence or increase 'eps'.")
    return(Inf)
  }

  df <- 2 / v_den
  if(verbose) cat("Satterthwaite df for contrast:", round(df, 3), "\n")
  return(df)
}
# Internal function -- not exported
# Satterthwaite degrees of freedom for fixed effects in a glmmTMB model.
#
# Two approaches depending on model specification:
#
#   dispformula = ~0 (no residual variance):
#     Uses central-difference numerical differentiation of the REML
#     variance function via constrained conditional refits.
#
#   Otherwise (residual variance estimated):
#     Routes to glmmTMB's built-in Satterthwaite df via
#     summary(model, ddf = "Satterthwaite").
#
# Arguments:
#   model   : a fitted glmmTMB model object
#   eps     : step size for central differences in the Jacobian
#             (only used when dispformula = ~0)
#   n_cores : number of parallel workers (default: all physical cores)
#             Set to 1 for sequential execution.
#             (only used when dispformula = ~0)
#   verbose : print progress messages and results
#
# Returns:
#   data.frame with columns: Term, Estimate, SE, df, t, p
.satterthwaite_glmmTMB_general <- function(model, eps = 1e-3,
                                           n_cores = parallel::detectCores(logical = FALSE),
                                           verbose = TRUE) {

  for (pkg in c("glmmTMB", "Matrix", "parallel")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop(sprintf("Package '%s' is required but not installed.", pkg))
  }

  # -- Route to glmmTMB built-in for models with residual variance -----------
  has_resid_var <- !identical(deparse(model$call$dispformula), "~0")

  if (has_resid_var) {
    if (verbose)
      message("Model has residual variance -- using glmmTMB's built-in Satterthwaite df.")
    s <- summary(model, ddf = "satterthwaite")$coefficients$cond
    result <- data.frame(
      Term      = rownames(s),
      Estimate  = s[, "Estimate"],
      SE        = s[, "Std. Error"],
      df        = s[, "ddf"],
      t         = s[, "t value"],
      p         = s[, "Pr(>|t|)"],
      row.names = NULL
    )
    if (verbose) {
      cat("\n=== Satterthwaite Results (glmmTMB built-in) ===\n")
      print(round(result[, c("Estimate", "SE", "df", "t", "p")], 4))
    }
    return(result)
  }

  # -- Quantities from the model ---------------------------------------------
  beta_vcov_orig <- vcov(model)$cond
  X              <- model.matrix(model)
  Z_full         <- glmmTMB::getME(model, "Z")
  theta_opt      <- glmmTMB::getME(model, "theta")
  n_obs          <- nrow(X)
  p_theta        <- length(theta_opt)

  re_struc <- model$modelInfo$reStruc$cond
  z_widths <- vapply(re_struc, function(x) x$blockSize * x$blockReps, numeric(1))

  if (any(is.na(z_widths)) || any(z_widths <= 0))
    stop("Could not extract valid Z block widths from model$modelInfo$reStruc$cond.")
  if (sum(z_widths) != ncol(Z_full))
    stop(sprintf(
      "Z block widths sum to %d but ncol(Z_full) = %d -- mapping is misaligned.",
      sum(z_widths), ncol(Z_full)
    ))

  gp_manual <- c(0L, cumsum(z_widths))

  # -- Asymptotic covariance of theta ----------------------------------------
  V_full       <- vcov(model, full = TRUE)
  nuisance_idx <- grep("^theta", colnames(V_full))
  if (length(nuisance_idx) == 0)
    stop("No '^theta' columns found in vcov(model, full=TRUE).")

  Sigma_nuisance <- V_full[nuisance_idx, nuisance_idx, drop = FALSE]
  Sigma_nuisance[is.na(Sigma_nuisance)] <- 0

  n_cores <- max(1L, min(as.integer(n_cores), p_theta))

  if (verbose)
    message(sprintf(
      "Computing Jacobian for %d theta parameter(s) using %d core(s).%s",
      p_theta, n_cores,
      if (p_theta > 6)
        "\n  Note: runtime scales with length(theta)."
      else ""
    ))

  # -- Conditional refit helper --------------------------------------------------
  make_refit_fn <- function(model_obj, X_, Z_full_, gp_manual_, n_obs_) {
    function(t_new) {
      tmp_mod <- tryCatch(
        withCallingHandlers(
          update(
            model_obj,
            start = list(theta = t_new),
            map   = list(theta = factor(rep(NA_integer_, length(t_new))))
          ),
          warning = function(w) invokeRestart("muffleWarning")
        ),
        error = function(e)
          stop("Model refit failed during Jacobian evaluation: ", conditionMessage(e))
      )

      vc_new <- VarCorr(tmp_mod)$cond
      V_tmp  <- Matrix::Matrix(0, n_obs_, n_obs_, sparse = TRUE)

      for (i in seq_along(vc_new)) {
        Sigma_i <- vc_new[[i]]
        Z_i     <- Z_full_[, (gp_manual_[i] + 1L):gp_manual_[i + 1L], drop = FALSE]
        n_grps  <- ncol(Z_i) / ncol(Sigma_i)
        V_tmp   <- V_tmp +
          Z_i %*% Matrix::kronecker(Matrix::Diagonal(n_grps), Sigma_i) %*% Matrix::t(Z_i)
      }

      diag_scale <- mean(Matrix::diag(V_tmp))
      if (!is.finite(diag_scale) || diag_scale <= 0) diag_scale <- 1
      V_tmp <- V_tmp + Matrix::Diagonal(n_obs_, diag_scale * 1e-8)

      XtVinvX <- tryCatch(
        as.matrix(Matrix::t(X_) %*% Matrix::solve(V_tmp, X_)),
        error = function(e) stop("Sparse solve failed: ", conditionMessage(e))
      )
df_vec <- sapply(coef_names, function(nm) {
  L <- setNames(1, nm)
  .satterthwaite_contrast(model, L, eps = eps, n_cores = n_cores, verbose = FALSE)
})
    }
  }

  refit_fn <- make_refit_fn(model, X, Z_full, gp_manual, n_obs)

  # -- Parallel central-difference Jacobian ----------------------------------
  h <- eps * pmax(abs(theta_opt), 1e-4)

  compute_jacobian <- function() {
    if (n_cores == 1L) {
      cols <- lapply(seq_len(p_theta), function(j) {
        x_fwd <- theta_opt; x_fwd[j] <- theta_opt[j] + h[j]
        x_bck <- theta_opt; x_bck[j] <- theta_opt[j] - h[j]
        (refit_fn(x_fwd) - refit_fn(x_bck)) / (2 * h[j])
      })
    } else {
      cl <- parallel::makeCluster(n_cores)
      tryCatch({
        parallel::clusterEvalQ(cl, {
          requireNamespace("glmmTMB", quietly = TRUE)
          requireNamespace("Matrix",  quietly = TRUE)
          library(glmmTMB)
          library(Matrix)
        })
        parallel::clusterExport(
          cl,
          c("model", "X", "Z_full", "gp_manual", "n_obs",
            "theta_opt", "h", "make_refit_fn"),
          envir = environment()
        )
        parallel::clusterEvalQ(cl,
                               refit_fn <- make_refit_fn(model, X, Z_full, gp_manual, n_obs)
        )
        data_name <- tryCatch(
          as.character(model$call$data),
          error = function(e) NULL
        )
        if (!is.null(data_name) && length(data_name) == 1 &&
            exists(data_name, envir = .GlobalEnv)) {
          parallel::clusterExport(cl, data_name, envir = .GlobalEnv)
        }
        cols <- parallel::parLapply(cl, seq_len(p_theta), function(j) {
          x_fwd <- theta_opt; x_fwd[j] <- theta_opt[j] + h[j]
          x_bck <- theta_opt; x_bck[j] <- theta_opt[j] - h[j]
          (refit_fn(x_fwd) - refit_fn(x_bck)) / (2 * h[j])
        })
      }, finally = {
        parallel::stopCluster(cl)
      })
    }
    do.call(cbind, cols)
  }

  g_log <- tryCatch(
    compute_jacobian(),
    error = function(e) stop("Jacobian computation failed: ", conditionMessage(e))
  )

  v_den_log <- rowSums((g_log %*% Sigma_nuisance) * g_log)

  n_bad <- sum(v_den_log <= 0)
  if (n_bad > 0)
    warning(sprintf(
      "%d term(s) have non-positive df denominator -- check model convergence or increase 'eps'.",
      n_bad
    ))

  # -- Results ---------------------------------------------------------------
  result <- data.frame(
    Term      = names(fixef(model)$cond),
    Estimate  = fixef(model)$cond,
    SE        = sqrt(diag(beta_vcov_orig)),
    df        = 2 / pmax(v_den_log, 1e-12),
    row.names = NULL
  )
  result$t <- result$Estimate / result$SE
  result$p  <- 2 * pt(abs(result$t), df = result$df, lower.tail = FALSE)

  if (verbose) {
    cat("\n=== Satterthwaite Results ===\n")
    print(round(result[, c("Estimate", "SE", "df", "t", "p")], 4))
  }

  return(result)
}
#' Satterthwaite Summary for glmmTMB Models
#'
#' @description
#' Produces a formatted summary for a fitted \code{glmmTMB} model using
#' Satterthwaite degrees of freedom for fixed effects, styled after
#' \code{glmmTMB}'s summary output.
#'
#' Two paths are used depending on the model specification:
#'
#' \strong{dispformula = ~0 (no residual variance):}
#' Full custom formatting -- family, formula, fit statistics, random effects,
#' and a fixed effects table with 5 decimal places. Satterthwaite df are
#' computed using a numerical differentiation approach with a parallel central-difference
#' Jacobian.
#'
#' \strong{Otherwise (residual variance estimated):}
#' Prints everything from \code{glmmTMB}'s native summary (including residual
#' variance and dispersion estimate) but replaces the fixed effects table with
#' the same 5 decimal place formatting using \code{glmmTMB}'s built-in
#' Satterthwaite df.
#'
#' @param model A fitted \code{glmmTMB} model object.
#' @param eps Numeric. Step size for central differences in the Jacobian.
#'   Only used when \code{dispformula = ~0}. Default is \code{1e-3}.
#' @param digits Integer. Number of digits for rounding in the coefficient
#'   table. Default is \code{5}.
#' @param verbose Logical. If \code{TRUE}, prints progress messages during
#'   Jacobian computation. Default is \code{FALSE}.
#'
#' @return Invisibly returns a list with elements:
#'   \describe{
#'     \item{satterthwaite}{A data frame with columns Term, Estimate, SE, df,
#'       t, and p from the Satterthwaite computation.}
#'     \item{vc}{VarCorr output from the model.}
#'     \item{fit_stats}{Named vector of AIC, BIC, logLik, -2*log(L), and
#'       df.resid.}
#'   }
#'
#' @examples
#' \dontrun{
#' m <- glmmTMB(RelSat_A ~ c_Amity_A * ECGender_A + c_Amity_P * ECGender_A +
#'                cs(0 + man + woman | DyadID),
#'              dispformula = ~0, REML = TRUE, data = pairwise_disting)
#'
#' summary_glmmTMB_satterthwaite(m)
#' }
#' @importFrom glmmTMB getME fixef VarCorr
#' @importFrom Matrix Matrix Diagonal kronecker
#' @importFrom stats vcov model.matrix pt logLik AIC BIC model.frame cov2cor update formula coef family symnum
#' @importFrom utils capture.output
#' @export
summary_glmmTMB_satterthwaite <- function(model, eps = 1e-3, digits = 5,
                                          verbose = FALSE) {

  # -- Helper: format and print fixed effects table --------------------------
  print_coef_table <- function(sw) {
    p_fmt <- ifelse(
      sw$p < 5e-7,
      formatC(sw$p, format = "e", digits = 3),
      formatC(sw$p, format = "f", digits = 6)
    )
    p_fmt[sw$p < 2e-16] <- "< 2e-16"

    stars <- symnum(
      sw$p,
      corr = FALSE, na = FALSE,
      cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
      symbols   = c("***", "**", "*", ".", " ")
    )

    print_df <- data.frame(
      Estimate     = formatC(sw$Estimate, format = "f", digits = 5),
      `Std. Error` = formatC(sw$SE,       format = "f", digits = 5),
      `t value`    = formatC(sw$t,        format = "f", digits = 5),
      ddf          = formatC(sw$df,       format = "f", digits = 3),
      `Pr(>|t|)`   = p_fmt,
      ` `          = as.character(stars),
      check.names  = FALSE,
      row.names    = sw$Term
    )

    print(print_df)
    cat("---\n")
    cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  }

  # -- 0. Satterthwaite df ---------------------------------------------------
  sw <- .satterthwaite_glmmTMB_general(model, eps = eps, verbose = verbose)

  # -- Detect model type -----------------------------------------------------
  has_resid_var <- !identical(deparse(model$call$dispformula), "~0")

  # ==========================================================================
  # PATH A: Model without dispformula = ~0
  # ==========================================================================
  if (has_resid_var) {

    s_full   <- summary(model, ddf = "satterthwaite")
    full_out <- capture.output(print(s_full))

    # Split at the second "Conditional model:" -- first is random effects,
    # second is fixed effects which we replace with custom formatting
    coef_start <- grep("^Conditional model:", full_out)
    if (length(coef_start) >= 2) {
      cat(paste(full_out[1:(coef_start[2] - 1)], collapse = "\n"), "\n")
    } else if (length(coef_start) == 1) {
      cat(paste(full_out[1:(coef_start[1] - 1)], collapse = "\n"), "\n")
    } else {
      cat(paste(full_out, collapse = "\n"), "\n")
    }

    cat("\nConditional model:\n")
    print_coef_table(sw)

    invisible(list(
      satterthwaite = sw,
      vc            = s_full$varcor,
      fit_stats     = s_full$AICtab
    ))

  } else {

    # ==========================================================================
    # PATH B: Model with dispformula = ~0
    # ==========================================================================

    # -- 1. Family / formula / data header ----------------------------------
    fam     <- family(model)
    fam_str <- paste0("Family: ", fam$family, "  ( ", fam$link, " )")

    raw_formula <- paste(deparse(formula(model), width.cutoff = 500L), collapse = " ")
    data_name   <- tryCatch(deparse(model$call$data), error = function(e) "<unknown>")

    cat(fam_str, "\n")
    cat(paste0("Formula:          ", raw_formula), "\n")
    cat("Dispersion:              ~0\n")
    cat("Data:", data_name, "\n\n")

    # -- 2. Fit statistics ---------------------------------------------------
    ll      <- logLik(model)
    n_obs   <- nrow(model.frame(model))
    p_fixed <- length(fixef(model)$cond)
    df_res  <- n_obs - p_fixed

    fit_stats <- c(
      AIC          = AIC(model),
      BIC          = BIC(model),
      logLik       = as.numeric(ll),
      `-2*log(L)`  = -2 * as.numeric(ll),
      df.resid     = df_res
    )

    stat_names <- names(fit_stats)
    stat_vals  <- sprintf("%.2f", fit_stats)
    col_widths <- pmax(nchar(stat_names), nchar(stat_vals)) + 1L

    name_line <- paste(mapply(function(nm, w) formatC(nm, width = w, flag = " "),
                              stat_names, col_widths), collapse = " ")
    val_line  <- paste(mapply(function(v,  w) formatC(v,  width = w, flag = " "),
                              stat_vals,  col_widths), collapse = " ")
    cat(name_line, "\n")
    cat(val_line,  "\n\n")

    # -- 3. Random effects ---------------------------------------------------
    cat("Random effects:\n\n")
    cat("Conditional model:\n")

    vc        <- VarCorr(model)$cond
    re_struc  <- model$modelInfo$reStruc$cond
    grp_names <- names(vc)

    re_rows <- lapply(grp_names, function(grp) {
      Sig   <- vc[[grp]]
      pars  <- rownames(Sig)
      n_par <- length(pars)
      vars  <- diag(Sig)
      sds   <- sqrt(vars)

      corr_mat <- if (n_par > 1L) cov2cor(Sig) else matrix(1)

      struc_entry <- re_struc[[grp]]
      cov_type    <- if (!is.null(struc_entry$covtype)) {
        tolower(as.character(struc_entry$covtype))
      } else "us"
      is_cs <- cov_type %in% c("cs", "hcs", "homcs")

      corr_strs <- character(n_par)
      if (n_par > 1L) {
        if (is_cs) {
          corr_strs[2L] <- paste0(
            formatC(corr_mat[2L, 1L], format = "f", digits = 3), " (cs)"
          )
        } else {
          for (i in 2L:n_par) {
            vals <- corr_mat[i, 1L:(i - 1L)]
            corr_strs[i] <- paste(
              formatC(vals, format = "f", digits = 3), collapse = "  "
            )
          }
        }
      }

      data.frame(
        Groups   = c(grp, rep("", n_par - 1L)),
        Name     = pars,
        Variance = formatC(vars, format = "f", digits = 5),
        Std.Dev. = formatC(sds,  format = "f", digits = 5),
        Corr     = corr_strs,
        stringsAsFactors = FALSE
      )
    })

    re_df <- do.call(rbind, re_rows)

    corr_width   <- 7L
    corr_tokens  <- strsplit(trimws(re_df$Corr), "\\s+")
    plain_tokens <- lapply(corr_tokens, function(x) x[!grepl("\\(cs\\)", x)])
    max_corrs    <- max(lengths(plain_tokens))

    re_df$Corr <- sapply(seq_len(nrow(re_df)), function(i) {
      raw    <- re_df$Corr[i]
      cs_tag <- if (grepl("\\(cs\\)", raw)) " (cs)" else ""
      toks   <- plain_tokens[[i]]
      if (length(toks) == 0L || all(toks == "")) {
        return(strrep(" ", corr_width * max_corrs))
      }
      padded   <- formatC(toks, width = corr_width, flag = " ")
      trailing <- strrep(" ", corr_width * (max_corrs - length(toks)))
      paste0(paste(padded, collapse = ""), trailing, cs_tag)
    })

    # Build table manually with leading space to match native glmmTMB style
    w_groups <- max(nchar(c("Groups",   re_df$Groups)))   + 1L
    w_name   <- max(nchar(c("Name",     re_df$Name)))     + 1L
    w_var    <- max(nchar(c("Variance", re_df$Variance))) + 1L
    w_sd     <- max(nchar(c("Std.Dev.", re_df$Std.Dev.))) + 1L

    header <- paste0(
      " ",
      formatC("Groups",   width = w_groups, flag = "-"),
      formatC("Name",     width = w_name,   flag = "-"),
      formatC("Variance", width = w_var,    flag = "-"),
      formatC("Std.Dev.", width = w_sd,     flag = "-"),
      " Corr"
    )
    cat(header, "\n")

    for (i in seq_len(nrow(re_df))) {
      corr_str <- if (nchar(trimws(re_df$Corr[i])) > 0) re_df$Corr[i] else ""
      row <- paste0(
        " ",
        formatC(re_df$Groups[i],   width = w_groups, flag = "-"),
        formatC(re_df$Name[i],     width = w_name,   flag = "-"),
        formatC(re_df$Variance[i], width = w_var,    flag = "-"),
        formatC(re_df$Std.Dev.[i], width = w_sd,     flag = "-"),
        corr_str
      )
      cat(row, "\n")
    }

    # -- Number of obs / groups ----------------------------------------------
    mf <- model.frame(model)
    n_grps <- sapply(grp_names, function(g) {
      if (g %in% names(mf)) {
        length(unique(mf[[g]]))
      } else if (grepl(":", g)) {
        parts <- strsplit(g, ":")[[1]]
        if (all(parts %in% names(mf))) {
          length(unique(interaction(mf[parts], drop = TRUE)))
        } else NA_integer_
      } else NA_integer_
    })

    grp_parts <- mapply(
      function(g, n) if (!is.na(n)) paste0(g, ", ", n) else g,
      grp_names, n_grps
    )
    cat(sprintf("\nNumber of obs: %d, groups:  %s\n", n_obs,
                paste(grp_parts, collapse = ";  ")))

    # -- 4. Fixed effects ----------------------------------------------------
    cat("\nConditional model:\n")
    print_coef_table(sw)

    invisible(list(
      satterthwaite = sw,
      vc            = vc,
      fit_stats     = fit_stats
    ))
  }
}
