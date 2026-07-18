# ---------------------------------------------------------------------------
# ui_tutorial.R
# First Steps / interactive tutorial tab.
# ---------------------------------------------------------------------------

#' Tutorial tab UI
#'
#' Builds the "First Steps" interactive tutorial tab of the ADMETShiny
#' application.
#'
#' @return A \code{shiny.tabPanel}.
#' @keywords internal
tutorial_tab <- function() {
  tabPanel(
    title = tagList(icon("graduation-cap"), "First Steps"),
    value = "tutorial",

    tags$head(tags$style(HTML("
      #tutorial_tabs > li > a { font-weight: 600; color: #555; }
      #tutorial_tabs > li.active > a { background-color: #4e73df !important; color: white !important; border: 1px solid #4e73df !important; }
      body.dark-mode #tutorial_tabs > li > a { background-color: #353535 !important; color: #ccc !important; border: 1px solid #444 !important; }
      body.dark-mode #tutorial_tabs > li.active > a { background-color: #4e73df !important; color: white !important; }

      .step-icon { font-size: 3.5rem; color: #4e73df; margin-bottom: 1rem; text-align: center; }
      body.dark-mode .step-icon { color: #8ab4f8; }

      .alert-info-dark { background-color: #e7f0fe; border: 1px solid #b8d4fc; border-left: 4px solid #3b82f6; color: #31708f; padding: 15px; border-radius: 8px; margin-top: 15px; }
      .alert-warning-dark { background-color: #fff8e1; border: 1px solid #ffe082; border-left: 4px solid #f39c12; color: #8a6d3b; padding: 15px; border-radius: 8px; margin-top: 15px; }
      .alert-success-dark { background-color: #e8f5e9; border: 1px solid #c8e6c9; border-left: 4px solid #2ecc71; color: #3c763d; padding: 15px; border-radius: 8px; margin-top: 15px; }

      body.dark-mode .alert-info-dark { background-color: #1e3a5f; border-color: #2a5a8f; color: #d1e7ff; }
      body.dark-mode .alert-warning-dark { background-color: #4d4000; border-color: #666200; color: #fff3cd; }
      body.dark-mode .alert-success-dark { background-color: #144d2e; border-color: #1e6e42; color: #d1e7dd; }
    "))),

    fluidRow(
      column(10, offset = 1,
             br(),
             tags$div(class = "info-card", style = "border-radius: 12px; text-align: center; margin-bottom: 25px;",
                      tags$div(class = "info-card-body", style = "padding: 30px;",
                               tags$h2("Interactive Tutorial", style = "font-weight: 800; color: #2c3e50;"),
                               tags$p("Learn how to get the most out of ADMETShiny in 4 simple steps.", style = "font-size: 1.2rem; color: #666;")
                      )
             ),

             tags$div(class = "nav-justified",
                      tabsetPanel(
                        id = "tutorial_tabs",

                        tabPanel(
                          title = tagList(icon("filter"), "ADMET Master Manager"),
                          br(),
                          fluidRow(
                            column(2, class = "step-icon", icon("upload")),
                            column(10,
                                   tags$div(class = "info-card accent-blue",
                                            tags$div(class = "info-card-header", icon("info-circle"), " The ADMET Master Manager Module"),
                                            tags$div(class = "info-card-body",
                                                     tags$p("This module is a four step wizard that brings any tabular dataset into the standardised scheam so the same filters and plots."),
                                                     tags$hr(),
                                                     tags$b("Workflow:"),
                                                     tags$ol(
                                                       tags$li(tags$b("Upload Dataset:"), " Use the file input in the sidebar to load your CSV / Excel"),
                                                       tags$li(tags$b("Map Colums"), " Manual column mapping to standard schema"),
                                                       tags$li(tags$b("Filter"), " Apply drug-likeness filters"),
                                                       tags$li(tags$b("Plots"), " Render any of 14 plots types")
                                                     ),
                                                     tags$div(class = "alert-info-dark",
                                                              tags$b(icon("lightbulb"), " Pro Tip: "), "You can download the filtered dataset at any time using the 'Download filtered dataset' button in the sidebar."
                                                     )
                                            )
                                   )
                            )
                          )
                        ),

                        tabPanel(
                          title = tagList(icon("flask"), " CDK + webchem"),
                          br(),
                          fluidRow(
                            column(2, class = "step-icon", icon("atom")),
                            column(10,
                                   tags$div(class = "info-card accent-green",
                                            tags$div(class = "info-card-header", icon("flask"), " Step 2: CDK & webchem Integration"),
                                            tags$div(class = "info-card-body",
                                                     tags$p("Don't have a ADMET dataset? No problem. You can generate descriptors entirely locally using the CDK & webchem module."),
                                                     tags$hr(),
                                                     tags$b("Workflow:"),
                                                     tags$ol(
                                                       tags$li(tags$b("Fetch SMILES:"), " Enter common names or CIDs (e.g., 'aspirin', 'ibuprofen') and retrieve their canonical SMILES from PubChem. Alternatively, paste SMILES directly or upload a CSV file with a SMILES column."),
                                                       tags$li(tags$b("Visualize:"), " Use the built-in 2D molecule viewer to inspect the retrieved structures."),
                                                       tags$li(tags$b("Calculate Descriptors:"), " Select the physicochemical properties you need (MW, TPSA, LogP, etc.) and calculate them locally with CDK."),
                                                       tags$li(tags$b("Filter & Plot:"), " The calculated descriptors are automatically mapped to the app's standards, allowing you to apply drug-likeness filters and generate plots.")
                                                     ),
                                                     tags$div(class = "alert-warning-dark",
                                                              tags$b(icon("exclamation-triangle"), " Note: "), "Because this module calculates descriptors locally, it requires the 'rcdk' package and a valid Java JDK installation."
                                                     )
                                            )
                                   )
                            )
                          )
                        ),

                        tabPanel(
                          title = tagList(icon("vials"), " Plots & Palette"),
                          br(),
                          fluidRow(
                            column(2, class = "step-icon", icon("brain")),
                            column(10,
                                   tags$div(class = "info-card accent-orange",
                                            tags$div(class = "info-card-header", icon("vials"), " The Plots & palette option"),
                                            tags$div(class = "info-card-body",
                                                     tags$p("The package exports 14 plots functions. All are ggplot2 based, except Radar, Tanimoto and Cluster Heatmap plots"),
                                                     tags$hr(),
                                                     tags$b("What's a palette?"),
                                                     tags$ul(
                                                       tags$li(tags$b("Apply_palette()"), " Detects whether the colour column is numeric (continuous) or not (discrte) and applies the appropiate ggplot2 scale"),
                                                       tags$li(tags$b("About the themes"), " Different themes could are apply, for example: viridis, magma, inferno, set1,set2, dark2 & spectral")

                                                     )
                                            )
                                   )
                            )
                          )
                        ),

                        tabPanel(
                          title = tagList(icon("palette"), " Report & Export"),
                          br(),
                          fluidRow(
                            column(2, class = "step-icon", icon("chart-pie")),
                            column(10,
                                   tags$div(class = "info-card accent-red",
                                            tags$div(class = "info-card-header", icon("palette"), " Report & Export"),
                                            tags$div(class = "info-card-body",
                                                     tags$p("The report tab collects the state of all loaded datasets and renders a comprehensive analysis report in HTML, PDF & Word format via rmarkdown."),
                                                     tags$hr(),
                                                     tags$b("Key Features:"),
                                                     tags$ul(
                                                       tags$li(tags$b("Render Report:"), " The Shiny app uses three download handlers for the PDF, Word and HTML formats. PDF export requires a LaTex installation & Word requires pandoc"),
                                                       tags$li(tags$b("Statistics Computation:"), " Several property stats (mean, median, sd, min & max), filters stats (per filter pass rate, mean violations), admet stats and a score composite drug-likeness score"),
                                                       tags$li(tags$b("Composite drug-likeness score:"), " Configure the format (PNG, PDF, SVG, TIFF), resolution (up to 1200 DPI), and dimensions using the 'Export Figure' panel, then download with one click.")
                                                     ),
                                                     tags$div(class = "alert-success-dark",
                                                              tags$b(icon("check-circle"), " Best Practice: "), "For publication-ready images, use vector formats like PDF or SVG, and a minimum of 600 DPI for raster formats."
                                                     )
                                            )
                                   )
                            )
                          )
                        )
                      )
             ),
             br(), br()
      )
    )
  )
}
