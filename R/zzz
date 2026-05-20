.onAttach <- function(libname, pkgname) {
  vc <- tryCatch(
    getFromNamespace(".valid_covstruct", "glmmTMB"),
    error = function(e) NULL
  )
  if (is.null(vc) || !"indisting" %in% names(vc)) {
    packageStartupMessage(
      "Note: apim requires a patched version of glmmTMB with indisting() support.\n",
      "Please run: remotes::install_github(\"RobAckerman/glmmTMB\", subdir = \"glmmTMB\")\n",
      "then restart R."
    )
  }
}
