# =============================================================================
# EXPORTED: plot_dyadic_rsa
# =============================================================================

#' Plot Dyadic Response Surface Analysis Surfaces
#'
#' @description
#' Produces interactive 3D response surface plots from the output of
#' \code{\link{run_dyadic_rsa}} using \code{plotly}. One panel is produced
#' for each surface (overall, and one per group in \code{mod_values}).
#' All panels share a common z-axis range for direct visual comparability.
#' The line of congruence (where actor = partner) is overlaid on each surface
#' as a dashed black line.
#'
#' @param rsa_out A list returned by \code{\link{run_dyadic_rsa}}.
#' @param actor Character string. Name of the actor predictor — must match
#'   the value passed to \code{\link{run_dyadic_rsa}}.
#' @param partner Character string. Name of the partner predictor — must
#'   match the value passed to \code{\link{run_dyadic_rsa}}.
#' @param moderator Character string or \code{NULL}. Name of the moderator
#'   variable — must match the value passed to \code{\link{run_dyadic_rsa}}.
#'   Default \code{NULL}.
#' @param mod_values Named list or \code{NULL}. Moderator values for
#'   group-specific surfaces, e.g. \code{list(Men = 1, Women = -1)}.
#'   Must match the value passed to \code{\link{run_dyadic_rsa}}.
#'   Default \code{NULL}.
#' @param xlim Numeric vector of length 2. Range of the actor axis.
#'   Default \code{c(-3, 3)}.
#' @param ylim Numeric vector of length 2. Range of the partner axis.
#'   Default \code{c(-3, 3)}.
#' @param n Integer. Number of grid points along each axis for surface
#'   rendering. Default \code{60}.
#' @param xlab Character string. Label for the actor axis.
#'   Default \code{"Actor"}.
#' @param ylab Character string. Label for the partner axis.
#'   Default \code{"Partner"}.
#' @param zlab Character string. Label for the outcome axis.
#'   Default \code{"Outcome"}.
#'
#' @return Invisibly returns a named list of \code{plotly} figure objects,
#'   one per panel. Each figure is also printed to the viewer.
#'
#' @references
#' Schönbrodt, F. D., Humberg, S., & Nestler, S. (2018). Testing similarity
#' effects with dyadic response surface analysis. \emph{European Journal of
#' Personality}, 32(6), 627--641. \doi{10.1002/per.2169}
#'
#' @examples
#' \dontrun{
#' rsa_out <- run_dyadic_rsa(
#'   data       = pairwise_data,
#'   outcome    = "RelSat_A",
#'   actor      = "c_PosEmotion_A",
#'   partner    = "c_PosEmotion_P",
#'   re_term    = "cs(0 + as.factor(ECGender_A) | DyadID)",
#'   moderator  = "ECGender_A",
#'   mod_values = list(Men = 1, Women = -1)
#' )
#'
#' plot_dyadic_rsa(
#'   rsa_out    = rsa_out,
#'   actor      = "c_PosEmotion_A",
#'   partner    = "c_PosEmotion_P",
#'   moderator  = "ECGender_A",
#'   mod_values = list(Men = 1, Women = -1),
#'   xlim       = c(-3, 3),
#'   ylim       = c(-3, 3),
#'   xlab       = "Actor positive emotion",
#'   ylab       = "Partner positive emotion",
#'   zlab       = "Relationship satisfaction"
#' )
#' }
#'
#' @importFrom stats setNames
#' @export
plot_dyadic_rsa <- function(
    rsa_out,
    actor,
    partner,
    moderator  = NULL,
    mod_values = NULL,
    xlim       = c(-3, 3),
    ylim       = c(-3, 3),
    n          = 60,
    xlab       = "Actor",
    ylab       = "Partner",
    zlab       = "Outcome"
) {

  if (!requireNamespace("plotly", quietly = TRUE))
    stop("Package 'plotly' is required. Install it with install.packages('plotly').")

  coefs <- rsa_out$coefs
  a2_nm <- paste0(actor,   "2")
  p2_nm <- paste0(partner, "2")
  ap_nm <- paste0(actor,   "_x_", partner)

  # -- get b-values for a moderator value -------------------------------------
  get_bvals <- function(mod_val = NULL) {
    b0 <- unname(coefs["(Intercept)"])
    b1 <- unname(coefs[actor])
    b2 <- unname(coefs[partner])
    b3 <- unname(coefs[a2_nm])
    b4 <- unname(coefs[ap_nm])
    b5 <- unname(coefs[p2_nm])

    if (!is.null(mod_val) && !is.null(moderator)) {
      b0 <- b0 + unname(coefs[moderator])                       * mod_val
      b1 <- b1 + unname(coefs[paste0(actor,   ":", moderator)]) * mod_val
      b2 <- b2 + unname(coefs[paste0(partner, ":", moderator)]) * mod_val
      b3 <- b3 + unname(coefs[paste0(a2_nm,   ":", moderator)]) * mod_val
      b4 <- b4 + unname(coefs[paste0(ap_nm,   ":", moderator)]) * mod_val
      b5 <- b5 + unname(coefs[paste0(p2_nm,   ":", moderator)]) * mod_val
    }
    list(b0 = b0, b1 = b1, b2 = b2, b3 = b3, b4 = b4, b5 = b5)
  }

  # -- surface matrix ---------------------------------------------------------
  make_surface <- function(bv) {
    x <- seq(xlim[1], xlim[2], length.out = n)
    y <- seq(ylim[1], ylim[2], length.out = n)
    z <- outer(x, y, function(a, p)
      bv$b0 + bv$b1*a + bv$b2*p +
        bv$b3*a^2 + bv$b4*a*p + bv$b5*p^2)
    list(x = x, y = y, z = z)
  }

  # -- line of congruence (A = P) ---------------------------------------------
  make_loc <- function(bv) {
    x <- seq(xlim[1], xlim[2], length.out = 50)
    z <- bv$b0 + (bv$b1 + bv$b2)*x + (bv$b3 + bv$b4 + bv$b5)*x^2
    list(x = x, z = z)
  }

  # -- single plotly panel ----------------------------------------------------
  make_plot <- function(bv, title_text, zmin, zmax) {
    s   <- make_surface(bv)
    loc <- make_loc(bv)

    plotly::plot_ly() |>
      plotly::add_surface(
        x            = s$x, y = s$y, z = s$z,
        cmin         = zmin, cmax = zmax,
        colorscale   = "RdBu",
        reversescale = TRUE,
        showscale    = TRUE,
        opacity      = 0.88
      ) |>
      plotly::add_trace(
        type       = "scatter3d",
        mode       = "lines",
        x          = loc$x,
        y          = loc$x,
        z          = loc$z,
        line       = list(color = "black", width = 5, dash = "dash"),
        name       = "Line of congruence",
        showlegend = TRUE
      ) |>
      plotly::layout(
        title = list(text = title_text, font = list(size = 14)),
        scene = list(
          xaxis      = list(title = xlab, range = xlim),
          yaxis      = list(title = ylab, range = ylim),
          zaxis      = list(title = zlab, range = c(zmin, zmax)),
          camera     = list(eye = list(x = 1.5, y = -1.5, z = 0.8)),
          aspectmode = "cube"
        )
      )
  }

  # -- collect panels ---------------------------------------------------------
  panel_bvals        <- list()
  panel_bvals[["Overall"]] <- get_bvals(
    mod_val = if (!is.null(moderator)) 0 else NULL
  )
  if (!is.null(mod_values)) {
    for (grp in names(mod_values)) {
      panel_bvals[[grp]] <- get_bvals(mod_val = mod_values[[grp]])
    }
  }

  # shared z range
  all_z <- lapply(panel_bvals, function(bv) make_surface(bv)$z)
  zmin  <- floor(min(sapply(all_z, min))   * 10) / 10
  zmax  <- ceiling(max(sapply(all_z, max)) * 10) / 10

  # -- print each panel separately --------------------------------------------
  figs <- list()
  for (nm in names(panel_bvals)) {
    fig       <- make_plot(panel_bvals[[nm]], nm, zmin, zmax)
    figs[[nm]] <- fig
    print(fig)
  }

  invisible(figs)
}
