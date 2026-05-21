#' Null coalescing operator
#'
#' Returns \code{a} if it is not \code{NULL}, otherwise returns \code{b}.
#'
#' @param a Left-hand side value.
#' @param b Right-hand side default value.
#' @return \code{a} if not \code{NULL}, else \code{b}.
#' @keywords internal
`%||%` <- function(a, b) if(!is.null(a)) a else b
