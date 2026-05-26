# =============================================================================
# EXPORTED: plot_dyadic_rsa
# =============================================================================

#' Plot Dyadic Response Surface Analysis Surfaces
#'
#' @description
#' Produces response surface plots from the output of
#' \code{\link{run_dyadic_rsa}}. Two plot types are available: an interactive
#' 3D surface using \code{plotly} (default), and a static black-and-white
#' wireframe plot using \code{lattice} matching the style of
#' \code{RSA::plotRSA}. One panel is produced per surface (overall, and one
#' per group in \code{mod_values}). The actor axis is reversed by default
#' to match \code{RSA::plotRSA} orientation. The line of congruence
#' (actor = partner) and line of incongruence (actor = -partner) are overlaid
#' on each surface.
#'
#' @param rsa_out A list returned by \code{\link{run_dyadic_rsa}}.
#' @param actor Character string. Name of the actor predictor — must match
#'   the value passed to \code{\link{run_dyadic_rsa}}.
#' @param partner Character string. Name of the partner predictor — must
#'   match the value passed to \code{\link{run_dyadic_rsa}}.
#' @param moderator Character string or \code{NULL}. Name of the moderator
#'   variable. Default \code{NULL}.
#' @param mod_values Named list or \code{NULL}. Moderator values for
#'   group-specific surfaces. Default \code{NULL}.
#' @param type Character string. Plot type: \code{"interactive"} for plotly
#'   (default) or \code{"static"} for lattice wireframe.
#' @param xlim Numeric vector of length 2. Range of the actor axis.
#'   Default \code{c(-3, 3)}. Note: the actor axis is displayed reversed
#'   (high values at front-left) to match \code{RSA::plotRSA} orientation.
#' @param ylim Numeric vector of length 2. Range of the partner axis.
#'   Default \code{c(-3, 3)}.
#' @param n Integer. Number of grid points along each axis.
#'   Default \code{60} for interactive, \code{21} for static.
#' @param xlab Character string. Label for the actor axis.
#'   Default \code{"Actor"}.
#' @param ylab Character string. Label for the partner axis.
#'   Default \code{"Partner"}.
#' @param zlab Character string. Label for the outcome axis.
#'   Default \code{"Outcome"}.
#' @param rotation Named list. Rotation angles for static plot.
#'   Default \code{list(x = -63, y = 32, z = 15)} matching
#'   \code{RSA::plotRSA}.
#'
#' @return Invisibly returns a named list of plot objects, one per panel.
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
#' # Interactive plotly
#' plot_dyadic_rsa(
#'   rsa_out    = rsa_out,
#'   actor      = "c_PosEmotion_A",
#'   partner    = "c_PosEmotion_P",
#'   moderator  = "ECGender_A",
#'   mod_values = list(Men = 1, Women = -1),
#'   type       = "interactive"
#' )
#'
#' # Static B&W wireframe
#' plot_dyadic_rsa(
#'   rsa_out    = rsa_out,
#'   actor      = "c_PosEmotion_A",
#'   partner    = "c_PosEmotion_P",
#'   moderator  = "ECGender_A",
#'   mod_values = list(Men = 1, Women = -1),
#'   type       = "static"
#' )
#' }
#'
#' @importFrom lattice wireframe panel.3dwire panel.3dscatter
#' @importFrom grDevices colorRampPalette
#' @importFrom stats setNames
#' @export
plot_dyadic_rsa <- function(
    rsa_out,
    actor,
    partner,
    moderator  = NULL,
    mod_values = NULL,
    type       = "interactive",
    xlim       = c(-3, 3),
    ylim       = c(-3, 3),
    n          = NULL,
    xlab       = "Actor",
    ylab       = "Partner",
    zlab       = "Outcome",
    rotation   = list(x = -63, y = 32, z = 15)
) {

  type <- match.arg(type, c("interactive", "static"))
  if (is.null(n)) n <- if (type == "static") 21L else 60L

  coefs <- rsa_out$coefs
  a2_nm <- paste0(actor,   "2")
  p2_nm <- paste0(partner, "2")
  ap_nm <- paste0(actor,   "_x_", partner)

  # -- get b-values for a moderator value ------------------------------------
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
    list(b0=b0, b1=b1, b2=b2, b3=b3, b4=b4, b5=b5)
  }

  # -- surface grid ----------------------------------------------------------
  make_grid <- function(bv, xlim, ylim, n) {
    x <- seq(xlim[1], xlim[2], length.out = n)
    y <- seq(ylim[1], ylim[2], length.out = n)
    grid <- expand.grid(x = x, y = y)
    grid$z <- with(grid,
                   bv$b0 + bv$b1*x + bv$b2*y +
                     bv$b3*x^2 + bv$b4*x*y + bv$b5*y^2)
    list(x=x, y=y, grid=grid)
  }

  # -- line of congruence (A = P) and incongruence (A = -P) -----------------
  make_loc <- function(bv, xlim) {
    x <- seq(xlim[1], xlim[2], length.out = 100)
    z <- bv$b0 + (bv$b1+bv$b2)*x + (bv$b3+bv$b4+bv$b5)*x^2
    data.frame(x=x, y=x, z=z)
  }

  make_loic <- function(bv, xlim) {
    x <- seq(xlim[1], xlim[2], length.out = 100)
    z <- bv$b0 + (bv$b1-bv$b2)*x + (bv$b3-bv$b4+bv$b5)*x^2
    data.frame(x=x, y=-x, z=z)
  }

  # -- collect panels --------------------------------------------------------
  panel_bvals        <- list()
  panel_bvals[["Overall"]] <- get_bvals(
    mod_val = if (!is.null(moderator)) 0 else NULL
  )
  if (!is.null(mod_values)) {
    for (grp in names(mod_values))
      panel_bvals[[grp]] <- get_bvals(mod_val = mod_values[[grp]])
  }

  # shared z range
  all_z  <- lapply(panel_bvals, function(bv) make_grid(bv, xlim, ylim, n)$grid$z)
  zmin   <- min(sapply(all_z, min))
  zmax   <- max(sapply(all_z, max))

  # -- B&W palette matching plotRSA(bw=TRUE) ---------------------------------
  bw_pal <- colorRampPalette(c("#FFFFFF", "#AAAAAA", "#030303"), bias=2)(11)

  # ==========================================================================
  # STATIC LATTICE PLOT
  # ==========================================================================
  make_static <- function(bv, title_text, zmin, zmax) {

    g   <- make_grid(bv, xlim, ylim, n)
    loc  <- make_loc(bv,  xlim)
    loic <- make_loic(bv, xlim)

    # clip LOC and LOIC to plotting region
    clip_lines <- function(df) {
      df[df$y >= min(ylim) & df$y <= max(ylim) &
           df$x >= min(xlim) & df$x <= max(xlim), ]
    }
    loc  <- clip_lines(loc)
    loic <- clip_lines(loic)

    at <- seq(zmin, zmax, length.out = length(bw_pal) - 1)

    # custom panel: draw contours on floor, then LOC/LOIC on surface
    my_panel <- function(x, y, z, xlim, ylim, zlim,
                         xlim.scaled, ylim.scaled, zlim.scaled, ...) {

      # rescale helpers
      rz <- function(z1) zlim.scaled[1] +
        diff(zlim.scaled) * (z1 - zlim[1]) / diff(zlim)
      rx <- function(x1) xlim.scaled[1] +
        diff(xlim.scaled) * (x1 - xlim[1]) / diff(xlim)
      ry <- function(y1) ylim.scaled[1] +
        diff(ylim.scaled) * (y1 - ylim[1]) / diff(ylim)

      floor_z <- rz(zmin - 0.01 * diff(range(zmin, zmax)))

      # contour lines projected on floor
      cs <- contourLines(
        x = seq(xlim[1], xlim[2], length.out = n),
        y = seq(ylim[1], ylim[2], length.out = n),
        z = matrix(g$grid$z, nrow = n, ncol = n),
        nlevels = 10
      )
      for (cl in cs) {
        panel.3dscatter(
          x = rx(cl$x), y = ry(cl$y),
          z = rep(floor_z, length(cl$x)),
          xlim=xlim, ylim=ylim, zlim=zlim,
          xlim.scaled=xlim.scaled, ylim.scaled=ylim.scaled,
          zlim.scaled=zlim.scaled,
          type="l", col.line="grey40", lty="solid", lwd=0.8, ...
        )
      }

      # wireframe surface
      panel.3dwire(x=x, y=y, z=z,
                   xlim=xlim, ylim=ylim, zlim=zlim,
                   xlim.scaled=xlim.scaled, ylim.scaled=ylim.scaled,
                   zlim.scaled=zlim.scaled,
                   col="grey10", lwd=0.4, ...)

      # LOC on surface
      if (nrow(loc) > 1) {
        panel.3dscatter(
          x=rx(loc$x), y=ry(loc$y), z=rz(loc$z),
          xlim=xlim, ylim=ylim, zlim=zlim,
          xlim.scaled=xlim.scaled, ylim.scaled=ylim.scaled,
          zlim.scaled=zlim.scaled,
          type="l", col.line="black", lty="solid", lwd=2, ...
        )
        # LOC projected on floor
        panel.3dscatter(
          x=rx(loc$x), y=ry(loc$y),
          z=rep(floor_z, nrow(loc)),
          xlim=xlim, ylim=ylim, zlim=zlim,
          xlim.scaled=xlim.scaled, ylim.scaled=ylim.scaled,
          zlim.scaled=zlim.scaled,
          type="l", col.line="black", lty="solid", lwd=2, ...
        )
      }

      # LOIC on surface
      if (nrow(loic) > 1) {
        panel.3dscatter(
          x=rx(loic$x), y=ry(loic$y), z=rz(loic$z),
          xlim=xlim, ylim=ylim, zlim=zlim,
          xlim.scaled=xlim.scaled, ylim.scaled=ylim.scaled,
          zlim.scaled=zlim.scaled,
          type="l", col.line="black", lty="solid", lwd=2, ...
        )
        # LOIC projected on floor
        panel.3dscatter(
          x=rx(loic$x), y=ry(loic$y),
          z=rep(floor_z, nrow(loic)),
          xlim=xlim, ylim=ylim, zlim=zlim,
          xlim.scaled=xlim.scaled, ylim.scaled=ylim.scaled,
          zlim.scaled=zlim.scaled,
          type="l", col.line="black", lty="solid", lwd=2, ...
        )
      }
    }

    # actor axis reversed: pass rev(xlim) to wireframe
    p <- wireframe(
      z ~ x * y, data = g$grid,
      drape        = TRUE,
      xlim         = rev(xlim),   # reversed actor axis
      ylim         = ylim,
      zlim         = c(zmin, zmax),
      xlab         = list(label = xlab, rot = 19),
      ylab         = list(label = ylab, rot = -40),
      zlab         = list(label = zlab, rot = 92),
      main         = title_text,
      screen       = rotation,
      at           = at,
      col.regions  = bw_pal,
      colorkey     = list(labels = list(cex = 0.8)),
      scales       = list(arrows = FALSE, col = "black",
                          font = 1, tck = c(1.5,1.5,1.5),
                          distance = c(1.3,1.3,1.4)),
      par.settings = list(
        axis.line  = list(col = "transparent"),
        box.3d     = list(col = "black")
      ),
      panel.3d.wireframe = my_panel
    )
    print(p)
    p
  }

  # ==========================================================================
  # INTERACTIVE PLOTLY PLOT
  # ==========================================================================
  make_interactive <- function(bv, title_text, zmin, zmax) {
    if (!requireNamespace("plotly", quietly = TRUE))
      stop("Package 'plotly' required. Install with install.packages('plotly').")

    g    <- make_grid(bv, xlim, ylim, n)
    loc  <- make_loc(bv,  xlim)
    loic <- make_loic(bv, xlim)
    zmat <- matrix(g$grid$z, nrow=n, ncol=n)

    p <- plotly::plot_ly() |>
      plotly::add_surface(
        x            = g$x,
        y            = g$y,
        z            = t(zmat),
        cmin         = zmin, cmax = zmax,
        colorscale   = "Greys",
        reversescale = TRUE,
        showscale    = TRUE,
        opacity      = 0.88
      ) |>
      plotly::add_trace(
        type="scatter3d", mode="lines",
        x=loc$x, y=loc$y, z=loc$z,
        line=list(color="black", width=5),
        name="Line of congruence"
      ) |>
      plotly::add_trace(
        type="scatter3d", mode="lines",
        x=loic$x, y=loic$y, z=loic$z,
        line=list(color="black", width=5, dash="dash"),
        name="Line of incongruence"
      ) |>
      plotly::layout(
        title = list(text=title_text, font=list(size=14)),
        scene = list(
          xaxis = list(title=xlab, range=rev(xlim),  # reversed
                       autorange="reversed"),
          yaxis = list(title=ylab, range=ylim),
          zaxis = list(title=zlab, range=c(zmin,zmax)),
          camera     = list(eye=list(x=1.5, y=-1.5, z=0.8)),
          aspectmode = "cube"
        )
      )
    print(p)
    p
  }

  # -- render all panels -----------------------------------------------------
  figs <- list()
  for (nm in names(panel_bvals)) {
    figs[[nm]] <- if (type == "static") {
      make_static(panel_bvals[[nm]], nm, zmin, zmax)
    } else {
      make_interactive(panel_bvals[[nm]], nm, zmin, zmax)
    }
  }

  invisible(figs)
}
