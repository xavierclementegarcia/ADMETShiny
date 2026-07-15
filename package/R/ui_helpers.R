# ---------------------------------------------------------------------------
# ui_helpers.R
# Shared CSS and small UI helpers used across the information tabs (About,
# Tutorial, Documentation) and the Home tab.
# ---------------------------------------------------------------------------

#' Juan Felipe, i love you, thank you ":?"
#' Shared info-card CSS
#'
#' Returns a \code{tags$head} element with the CSS rules for the
#' \code{.info-card} family of classes used across the application, including
#' dark-mode variants.
#'
#' @return A \code{shiny.tag}.
#' @keywords internal
info_card_css <- function() {
  tags$head(tags$style(HTML("
    .info-card {
      background-color: #ffffff;
      border: 1px solid #e3e6f0;
      border-radius: 0.5rem;
      box-shadow: 0 0.15rem 1.75rem 0 rgba(58, 59, 69, 0.1);
      margin-bottom: 1.5rem;
      height: 96%;
      display: flex;
      flex-direction: column;
    }
    .info-card-header {
      background-color: #f8f9fc;
      border-bottom: 1px solid #e3e6f0;
      padding: 1rem 1.35rem;
      border-radius: 0.5rem 0.5rem 0 0;
      font-weight: 700;
      font-size: 1.1rem;
      display: flex;
      align-items: center;
    }
    .info-card-body {
      padding: 1.35rem;
      flex-grow: 1;
    }
    .accent-blue .info-card-header { color: #3a5169; border-left: 4px solid #4e73df; }
    .accent-green .info-card-header { color: #3a5169; border-left: 4px solid #1cc88a; }
    .accent-orange .info-card-header { color: #3a5169; border-left: 4px solid #f6c23e; }
    .accent-red .info-card-header { color: #3a5169; border-left: 4px solid #e74a3b; }
    .accent-cyan .info-card-header { color: #3a5169; border-left: 4px solid #36b9cc; }

    body.dark-mode .info-card {
      background-color: #2d2d2d !important;
      border: 1px solid #444 !important;
      box-shadow: 0 0.15rem 1.75rem 0 rgba(0,0,0,0.4) !important;
    }
    body.dark-mode .info-card-header {
      background-color: #353535 !important;
      border-bottom: 1px solid #444 !important;
      color: #e0e0e0 !important;
    }
    body.dark-mode .accent-blue .info-card-header { border-left: 4px solid #4e73df; }
    body.dark-mode .accent-green .info-card-header { border-left: 4px solid #1cc88a; }
    body.dark-mode .accent-orange .info-card-header { border-left: 4px solid #f6c23e; }
    body.dark-mode .accent-red .info-card-header { border-left: 4px solid #e74a3b; }
    body.dark-mode .accent-cyan .info-card-header { border-left: 4px solid #36b9cc; }

    body.dark-mode .text-muted { color: #aaa !important; }
    body.dark-mode hr { border-top: 1px solid #555; }
  ")))
}
