# ---------------------------------------------------------------------------
# ui_docs.R
# Documentation tab.
# ---------------------------------------------------------------------------

#' Documentation tab UI
#'
#' Builds the Documentation tab of the ADMETShiny application, with reference
#' tables for the drug-likeness rules, ADMET models and data schemas.
#'
#' @return A \code{shiny.tabPanel}.
#' @keywords internal
docs_tab <- function() {
  tabPanel(
    title = tagList(icon("book-open"), "Documentation"),
    value = "docs",

    tags$head(tags$style(HTML("
      #docs_tabs > li > a { font-weight: 600; color: #555; }
      #docs_tabs > li.active > a { background-color: #4e73df !important; color: white !important; border: 1px solid #4e73df !important; }
      body.dark-mode #docs_tabs > li > a { background-color: #353535 !important; color: #ccc !important; border: 1px solid #444 !important; }
      body.dark-mode #docs_tabs > li.active > a { background-color: #4e73df !important; color: white !important; }

      .doc-table { width: 100%; font-size: 14px; margin-top: 10px; }
      .doc-table th { text-align: left; padding: 8px; border-bottom: 2px solid #ddd; }
      .doc-table td { padding: 8px; border-bottom: 1px solid #eee; }
      body.dark-mode .doc-table th { border-bottom: 2px solid #555; }
      body.dark-mode .doc-table td { border-bottom: 1px solid #444; }

      .math-block { background-color: #f8f9fc; border: 1px solid #e3e6f0; padding: 15px; border-radius: 8px; font-family: monospace; white-space: pre-wrap; color: #333; }
      body.dark-mode .math-block { background-color: #353535 !important; border: 1px solid #444 !important; color: #e0e0e0 !important; }

      .alert-info-dark { background-color: #e7f0fe; border: 1px solid #b8d4fc; border-left: 4px solid #3b82f6; color: #31708f; padding: 15px; border-radius: 8px; margin-top: 15px; }
      .alert-warning-dark { background-color: #fff8e1; border: 1px solid #ffe082; border-left: 4px solid #f39c12; color: #8a6d3b; padding: 15px; border-radius: 8px; margin-top: 15px; }
      body.dark-mode .alert-info-dark { background-color: #1e3a5f; border-color: #2a5a8f; color: #d1e7ff; }
      body.dark-mode .alert-warning-dark { background-color: #4d4000; border-color: #666200; color: #fff3cd; }
    "))),

    fluidRow(
      column(10, offset = 1,
             br(),
             tags$div(class = "info-card", style = "border-radius: 12px; text-align: center; margin-bottom: 25px;",
                      tags$div(class = "info-card-body", style = "padding: 30px;",
                               tags$h2("Package Documentation", style = "font-weight: 800; color: #2c3e50;"),
                               tags$p("Detailed reference for the algorithms, rules, and data mappings used under the hood.", style = "font-size: 1.2rem; color: #666;")
                      )
             ),

             tags$div(class = "nav-justified",
                      tabsetPanel(
                        id = "docs_tabs",

                        tabPanel(
                          title = tagList(icon("ruler"), " Drug-likeness Rules"),
                          br(),
                          tags$div(class = "info-card accent-orange",
                                   tags$div(class = "info-card-header", icon("ruler"), " Implemented Thresholds"),
                                   tags$div(class = "info-card-body",
                                            tags$p("ADMETShiny evaluates drug-likeness based on the classical literature rules. The app calculates the number of violations for each rule and allows you to filter compounds that fall within your desired tolerance."),
                                            tags$table(class = "doc-table",
                                                       tags$thead(tags$tr(tags$th("Rule"), tags$th("Property"), tags$th("Threshold / Range"))),
                                                       tags$tbody(
                                                         tags$tr(tags$td(tags$b("Lipinski")), tags$td("MW / LogP / HBA / HBD"), tags$td("<= 500 / <= 5 / <= 10 / <= 5")),
                                                         tags$tr(tags$td(tags$b("Veber")), tags$td("Rotatable Bonds / TPSA"), tags$td("<= 10 / <= 140 (or HBA+HBD <= 12)")),
                                                         tags$tr(tags$td(tags$b("Ghose")), tags$td("MW / MR / LogP / Heavy Atoms"), tags$td("160-480 / 40-130 / -0.4 to 5.6 / 20-70")),
                                                         tags$tr(tags$td(tags$b("Egan")), tags$td("TPSA / LogP"), tags$td("<= 131.6 / <= 5.88")),
                                                         tags$tr(tags$td(tags$b("Muegge")), tags$td("MW / LogP / HBA / HBD / RB / TPSA"), tags$td("200-600 / -2 to 5 / <= 10 / <= 5 / <= 15 / <= 150"))
                                                       )
                                            )
                                   )
                          )
                        ),

                        tabPanel(
                          title = tagList(icon("egg"), " ADMET Models"),
                          br(),
                          fluidRow(
                            column(6,
                                   tags$div(class = "info-card accent-red",
                                            tags$div(class = "info-card-header", icon("egg"), " BOILED-Egg Model"),
                                            tags$div(class = "info-card-body",
                                                     tags$p("The BOILED-Egg model (Daina & Zoete, 2016) is computed using the physicochemical properties ", tags$b("TPSA"), " (x-axis) and ", tags$b("LogP"), " (y-axis). The model was originally calibrated with WLOGP; the app uses the generic LogP column as an approximation."),
                                                     tags$p(tags$b("GI Absorption (White Region):")),
                                                     tags$p("Classification uses the official polygon coordinates from the BOILED-Egg supplementary data (Data S3), via point-in-polygon testing. This is more accurate than the ellipse approximation."),
                                                     tags$p(tags$b("BBB Permeability (Yellow Yolk):")),
                                                     tags$p("Same point-in-polygon approach using the official BBB polygon coordinates.")
                                            )
                                   )
                            ),
                            column(6,
                                   tags$div(class = "info-card accent-blue",
                                            tags$div(class = "info-card-header", icon("shield-virus"), " P-glycoprotein (P-gp)"),
                                            tags$div(class = "info-card-body",
                                                     tags$p("admetshiny ships a Random Forest classifier that predicts whether a small molecule is a substrate of P-glycoprotein (P-gp/ABCB1). The model is stores as an R list of 100 decision trees."),
                                                     tags$hr(),
                                                     tags$b("Training data & feature importance"),
                                                     tags$p("The 882 compound training set comes from Metrabase. The 9 CDK descriptors and their normalized Gini feature importances."),
                                                     tags$hr(),
                                                     tags$b("The accuracy:"),
                                                     tags$p("The Random Forest was built taking the reference of heuristics calculus based on Seeling (1998) and Didzpetris (2003). This Random Forest model accuracy reaches to 0.694 and sensitivity to 0.700.")
                                            )
                                   )
                            )
                          )
                        ),


                      )
             ),
             br(), br()
      )
    )
  )
}
