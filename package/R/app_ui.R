# ---------------------------------------------------------------------------
# app_ui.R
# Main UI builder: assembles all tabs into a navbarPage and attaches the
# dark-mode module.
# ---------------------------------------------------------------------------

#' ADMETShiny application UI
#'
#' Builds the full Shiny UI for the ADMETShiny application: a collapsible
#' \code{navbarPage} with the Home, CDK & webchem, ADMET Master Manager,
#' Report, Documentation, Tutorial and About tabs, plus the floating
#' dark-mode toggle button.
#'
#' This function is intended for internal use; end users should launch the app
#' with \code{\link{run_app}}.
#'
#' @return A \code{shiny.tag.list} suitable for \code{shiny::shinyApp()}.
#' @keywords internal
app_ui <- function() {
  tagList(
    navbarPage(
      id = "main_nav",
      title = "admetshiny",
      collapsible = TRUE,
      home_tab(),
      cdk_tab(),
      admet_master_tab(),
      report_tab(),
      docs_tab(),
      tutorial_tab(),
      about_tab()
    ),
    dark_mode_ui()
  )
}
