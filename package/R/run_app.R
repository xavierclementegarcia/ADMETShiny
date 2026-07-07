# ---------------------------------------------------------------------------
# run_app.R
# Public entry point to launch the ADMETShiny application.
# ---------------------------------------------------------------------------

#' Launch the ADMETShiny application
#'
#' Starts the interactive ADMETShiny Shiny application. The function registers
#' the package's static web assets (the hex sticker) and returns a Shiny app
#' object that can be printed or run directly.
#'
#' @param onStart An optional function called once when the app starts
#'   (forwarded to \code{\link[shiny:shinyApp]{shiny::shinyApp}}).
#' @param ... Additional arguments forwarded to
#'   \code{\link[shiny:shinyApp]{shiny::shinyApp}} (e.g. \code{options}).
#' @return A \code{shiny.appobj} object (invisibly); called for the side-effect
#'   of launching a Shiny application when printed.
#' @export
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   admetshiny::run_app()
#' }
#' }
run_app <- function(onStart = NULL, ...) {

  www_dir <- system.file("www", package = "admetshiny", mustWork = FALSE)
  if (dir.exists(www_dir) && nchar(www_dir) > 0) {
    shiny::addResourcePath("admetshiny", www_dir)
  }

  shiny::shinyApp(ui = app_ui(), server = app_server, onStart = onStart, ...)
}
