# ---------------------------------------------------------------------------
# utils-pipe.R
# Re-export the magrittr pipe so it can be used inside package code and by
# package users without attaching magrittr explicitly.
# ---------------------------------------------------------------------------

#' Pipe operator
#'
#' Re-exports \code{\link[magrittr:pipe]{\%>\%}} from \pkg{magrittr} so that
#' the pipe can be used with \code{admetshiny::\%>\%} or after attaching the
#' package.
#'
#' @param lhs A value or the magrittr placeholder.
#' @param rhs A function call using the magrittr semantics.
#' @export
#' @name %>%
#' @rdname pipe
#' @usage lhs \%>\% rhs
#' @return The result of calling \code{rhs(lhs, ...)}.
#' @keywords internal
#' @importFrom magrittr %>%
NULL
