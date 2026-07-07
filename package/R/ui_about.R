# ---------------------------------------------------------------------------
# ui_about.R
# About tab.
# ---------------------------------------------------------------------------

#' About tab UI
#'
#' Builds the About tab of the ADMETShiny application.
#'
#' @return A \code{shiny.tabPanel}.
#' @keywords internal
about_tab <- function() {
  tabPanel(
    title = tagList(icon("circle-info"), "About"),
    value = "about",

    info_card_css(),

    fluidRow(
      column(12,
             tags$div(class = "info-card",
                      tags$div(class = "info-card-body", style = "text-align: center; padding: 40px 20px;",
                               tags$img(src = "admetshiny/hexicon.png", width = "140px",
                                        style = "margin-bottom: 20px;"),
                               tags$h2("ADMETShiny",
                                       style = "font-weight: 800; color: #2c3e50;"),
                               tags$h4(style = "color: #666; font-weight: 400;",
                                       "An open-source web application for management, calculation, filtering and selection of drug-likeness potential molecules in R with Shiny"),
                               tags$span(class = "badge",
                                         style = "background-color: #4ca1af; color: white; font-size: 14px; padding: 8px 15px; margin-top: 15px;",
                                         tags$b(paste0("Version ", ADMETSHINY_VERSION, " ")),
                                         tags$em(paste0("- ", ADMETSHINY_CODENAME)))
                      )
             )
      )
    ),

    fluidRow(
      column(6,
             tags$div(class = "info-card accent-blue",
                      tags$div(class = "info-card-header", icon("circle-info"), " Overview"),
                      tags$div(class = "info-card-body",
                               tags$p(style = "color:#555; font-size: 15px; line-height: 1.6;",
                                      "ADMETShiny is an R package that provides an interactive environment for the calculation, filtering, visualization, and exploratory analysis of molecular descriptors and ADMET properties. It integrates cheminformatics and bioinformatics workflows with intuitive dashboards to support the prioritization of compounds in early stage drug discovery.")
                      )
             )
      ),
      column(6,
             tags$div(class = "info-card accent-green",
                      tags$div(class = "info-card-header", icon("star"), " Main Features"),
                      tags$div(class = "info-card-body",
                               tags$ul(style = "color:#555; font-size: 15px; line-height: 1.8; padding-left: 20px;",
                                       tags$li("Interactive compound filtering (Lipinski, Veber, Ghose, Egan, Muegge)"),
                                       tags$li("BOILED-Egg visualization & Tanimoto/AGNES clustering"),
                                       tags$li("Chemical space exploration (PCA, t-SNE, Parallel Coordinates)"),
                                       tags$li("Local descriptor calculation via CDK & webchem"),
                                       tags$li("Integration with SwissADME, ADMETlab 3.0 & Deep-PK datasets")
                               )
                      )
             )
      )
    ),

    fluidRow(
      column(4,
             tags$div(class = "info-card accent-cyan",
                      tags$div(class = "info-card-header", icon("user"), " Developer"),
                      tags$div(class = "info-card-body",
                               tags$b("Xavier Clemente Garcia Cevallos"), br(),
                               "Biologist", br(),
                               "Universidad del Cauca", br(),
                               "Popayan, Colombia", br(),
                               tags$a(href = "mailto:xgarcia@unicauca.edu.co", "xgarcia@unicauca.edu.co", style = "color: #36b9cc;"), br(),
                               tags$a(href = "https://github.com/xavierclementegarcia", target = "_blank", icon("github"), " GitHub Profile", style = "color: #36b9cc;"), br(),
                               tags$a(href = "https://xavierclementegarcia.github.io/", target = "_blank", icon("globe"), " My Website", style = "color: #36b9cc;")
                      )
             )
      ),
      column(4,
             tags$div(class = "info-card accent-orange",
                      tags$div(class = "info-card-header", icon("box"), " Package Information"),
                      tags$div(class = "info-card-body",
                               tags$table(class = "table table-borderless", style = "margin-bottom: 0;",
                                          tags$tr(tags$td(tags$b("Package:")), tags$td("admetshiny")),
                                          tags$tr(tags$td(tags$b("Version:")), tags$td(ADMETSHINY_VERSION)),
                                          tags$tr(tags$td(tags$b("Codename:")), tags$td(tags$em(ADMETSHINY_CODENAME))),
                                          tags$tr(tags$td(tags$b("License:")), tags$td("CC BY 4.0"))
                               )
                      )
             )
      ),
      column(4,
             tags$div(class = "info-card accent-red",
                      tags$div(class = "info-card-header", icon("shield-alt"), " Disclaimer"),
                      tags$div(class = "info-card-body",
                               tags$p(style = "color:#555; font-size: 14px;",
                                      "Predictions and visualizations generated by this package are intended for research purposes only and should not be interpreted as experimental or clinical evidence.")
                      )
             )
      )
    ),

    fluidRow(
      column(6,
             tags$div(class = "info-card accent-blue",
                      tags$div(class = "info-card-header", icon("heart"), " Acknowledgements"),
                      tags$div(class = "info-card-body",
                               tags$p(style = "color:#555;", "This work would not be possible without the following open-source projects:"),
                               tags$div(style = "display: flex; flex-wrap: wrap; gap: 10px;",
                                        tags$span(class = "badge bg-light text-dark", "SwissADME"),
                                        tags$span(class = "badge bg-light text-dark", "PubChem"),
                                        tags$span(class = "badge bg-light text-dark", "Chemistry Development Kit (CDK)"),
                                        tags$span(class = "badge bg-light text-dark", "webchem"),
                                        tags$span(class = "badge bg-light text-dark", "rcdk"),
                                        tags$span(class = "badge bg-light text-dark", "Shiny"),
                                        tags$span(class = "badge bg-light text-dark", "ggplot2"),
                                        tags$span(class = "badge bg-light text-dark", "DT")
                               )
                      )
             )
      ),
      column(6,
             tags$div(class = "info-card accent-green",
                      tags$div(class = "info-card-header", icon("book"), " Citation & Codename"),
                      tags$div(class = "info-card-body",
                               tags$p(style = "color:#555;", "If you use admetshiny in your research, please cite the corresponding publication once available."),
                               tags$hr(),
                               tags$h5(tags$b(paste0("Version Codename: ", ADMETSHINY_CODENAME))),
                               tags$p(style = "color:#555; font-size: 14px; text-align: justify;",
                                      tags$em(ADMETSHINY_CODENAME),
                                      " is named after one of Colombia's iconic frailejon species, native to the high-altitude paramo ecosystems of the northern Andes. This codename reflects the country's unique botanical diversity, the ecological importance of paramo conservation, and the potential of native plants as a source of bioactive compounds for drug discovery."
                               )
                      )
             )
      )
    ),

    br()
  )
}
