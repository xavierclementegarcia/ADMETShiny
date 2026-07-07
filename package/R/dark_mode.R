# ---------------------------------------------------------------------------
# dark_mode.R
# Modular dark-mode logic: CSS, custom message handler and toggle button.
# ---------------------------------------------------------------------------

#' Dark mode UI module
#'
#' Returns the UI elements required for the dark-mode toggle: the CSS rules
#' for \code{body.dark-mode}, the JavaScript custom message handler and the
#' floating toggle button.
#'
#' @return A \code{shiny.tag.list} with the CSS, script and action button.
#' @keywords internal
dark_mode_ui <- function() {
  tagList(
    tags$head(tags$style(HTML("
      body.dark-mode {
        background-color: #1e1e1e !important;
        color: #e0e0e0 !important;
      }
      body.dark-mode .well {
        background-color: #2d2d2d !important;
        border: 1px solid #444 !important;
        color: #e0e0e0 !important;
      }
      body.dark-mode .navbar-default {
        background-color: #252525 !important;
        border-color: #333 !important;
      }
      body.dark-mode .navbar-default .navbar-nav > li > a {
        color: #ccc !important;
      }
      body.dark-mode .navbar-default .navbar-nav > .active > a,
      body.dark-mode .navbar-default .navbar-nav > .active > a:focus,
      body.dark-mode .navbar-default .navbar-nav > .active > a:hover {
        color: #fff !important;
        background-color: #375a7f !important;
      }
      body.dark-mode .nav-tabs {
        border-bottom: 1px solid #444 !important;
      }
      body.dark-mode .nav-tabs > li > a {
        color: #ccc !important;
        border: 1px solid transparent !important;
      }
      body.dark-mode .nav-tabs > li.active > a,
      body.dark-mode .nav-tabs > li.active > a:hover,
      body.dark-mode .nav-tabs > li.active > a:focus {
        color: #fff !important;
        background-color: #2d2d2d !important;
        border: 1px solid #444 !important;
        border-bottom-color: transparent !important;
      }
      body.dark-mode input, body.dark-mode select, body.dark-mode textarea {
        background-color: #444 !important;
        color: #fff !important;
        border-color: #666 !important;
      }
      body.dark-mode .btn-default {
        background-color: #444 !important;
        color: #fff !important;
        border-color: #666 !important;
      }
      body.dark-mode .dataTables_wrapper {
        color: #e0e0e0 !important;
      }
      body.dark-mode table {
        color: #e0e0e0 !important;
      }
      body.dark-mode .table-striped > tbody > tr:nth-of-type(odd) {
        background-color: #353535 !important;
      }
      body.dark-mode .table-hover > tbody > tr:hover {
        background-color: #404040 !important;
      }
    "))),
    tags$script(HTML("
      Shiny.addCustomMessageHandler('toggle-dark-mode', function(isDark) {
        if(isDark) {
          document.body.classList.add('dark-mode');
        } else {
          document.body.classList.remove('dark-mode');
        }
      });
    ")),
    actionButton(
      "toggle_darkmode", label = "", icon = icon("moon"),
      style = paste(
        "position: fixed; bottom: 25px; right: 25px; z-index: 1000;",
        "width: 50px; height: 50px; border-radius: 50%;",
        "box-shadow: 0 4px 8px rgba(0,0,0,0.3);"
      )
    )
  )
}
