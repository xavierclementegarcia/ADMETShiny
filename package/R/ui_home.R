# ---------------------------------------------------------------------------
# ui_home.R
# Home / landing tab.
# ---------------------------------------------------------------------------

#' Home tab UI
#'
#' Builds the Home (landing) tab of the ADMETShiny application.
#'
#' @return A \code{shiny.tabPanel}.
#' @keywords internal
home_tab <- function() {
  tabPanel(
    title = tagList(icon("home"), "Home"),
    value = "home",

    tags$head(tags$style(HTML("
      .nav-card-btn {
        display: flex; flex-direction: column; align-items: center; justify-content: center;
        text-align: center; height: 140px; border-radius: 12px; border: 1px solid #e3e6f0;
        background-color: #ffffff; box-shadow: 0 0.15rem 1.75rem 0 rgba(58, 59, 69, 0.1);
        transition: all 0.3s ease; cursor: pointer; color: #4a4a4a; width: 100%; margin-bottom: 20px;
      }
      .nav-card-btn:hover {
        transform: translateY(-4px); text-decoration: none; color: #4e73df;
        border-color: #4e73df; box-shadow: 0 0.5rem 2rem 0 rgba(78, 115, 223, 0.2);
      }
      .nav-card-btn i { font-size: 32px; margin-bottom: 12px; }
      body.dark-mode .nav-card-btn { background-color: #2d2d2d; border-color: #444; color: #e0e0e0; }
      body.dark-mode .nav-card-btn:hover { color: #fff; border-color: #4e73df; box-shadow: 0 0.5rem 2rem 0 rgba(0,0,0,0.4); }
      body.dark-mode pre { background-color: #353535 !important; color: #e0e0e0 !important; border: 1px solid #444 !important; }
    "))),

    fluidRow(
      column(12,
             tags$div(class = "info-card", style = "border-radius: 15px; margin-top: 20px;",
                      tags$div(class = "info-card-body", style = "text-align: center; padding: 40px 20px;",
                               tags$img(src = "admetshiny/hexicon.png", width = "160px",
                                        style = "filter: drop-shadow(0 4px 8px rgba(0,0,0,0.2));"),
                               tags$h2("ADMETShiny",
                                       style = "font-weight: 800; color: #2c3e50; margin-top: 20px;"),
                               tags$h5(style = "color: #666; font-weight: 400;",
                                       "An open-source web application for management, calculation, filtering and selection of drug-likeness potential molecules in R with Shiny"),
                               tags$span(class = "badge",
                                         style = "background-color: #4ca1af; color: white; font-size: 14px; padding: 8px 15px; margin-top: 15px;",
                                         tags$b(paste0("Version ", ADMETSHINY_VERSION, " ")))
                      )
             )
      )
    ),

    fluidRow(
      column(6,
             tags$div(class = "info-card accent-blue",
                      tags$div(class = "info-card-header", icon("info-circle"), " Description"),
                      tags$div(class = "info-card-body",
                               tags$p(style = "color:#555; line-height: 1.6;",
                                      "ADMETShiny is an R package that provides an interactive environment for the calculation, filtering, visualization, and exploratory analysis of molecular descriptors and ADMET properties. It integrates cheminformatics and bioinformatics workflows with intuitive dashboards to support the prioritization of compounds in early stage drug discovery.")
                      )
             )
      ),
      column(6,
             tags$div(class = "info-card accent-green",
                      tags$div(class = "info-card-header", icon("quote-right"), " How to cite?"),
                      tags$div(class = "info-card-body",
                               tags$pre(style = "background:#f8f9fc; padding:15px; border-radius:8px; font-size:12px; white-space: pre-wrap; border: 1px solid #e3e6f0; color: #333;",
                                        paste0(
                                          "Garcia Cevallos, X. C. (", format(Sys.Date(), "%Y"), "). ADMETShiny: ",
                                          "A open-source web application for management, calculation, filtering and selection of drug-likeness.\n",
                                          "Package & App: https://github.com/xavierclementegarcia/admetshiny"
                                        )
                               )
                      )
             )
      )
    ),

    fluidRow(
      column(12,
             tags$div(class = "info-card accent-orange",
                      tags$div(class = "info-card-header", icon("book"), " References"),
                      tags$div(class = "info-card-body",
                               tags$ul(style = "color:#555; line-height: 1.8; padding-left: 20px;",
                                       tags$li("Daina, A., Michielin, O., & Zoete, V. (2017). SwissADME: a free web tool to evaluate pharmacokinetics, drug-likeness and medicinal chemistry friendliness of small molecules. ", tags$i("Scientific reports"), ", 7(1), 42717."),
                                       tags$li("Myung, Y., de Sa, A. G., & Ascher, D. B. (2024). Deep-PK: deep learning for small molecule pharmacokinetic and toxicity prediction. ", tags$i("Nucleic acids research"), ", 52(W1), W469-W475."),
                                       tags$li("Fu, L., Shi, S., Yi, J., Wang, N., He, Y., Wu, Z., ... & Cao, D. (2024). ADMETlab 3.0: an updated comprehensive online ADMET prediction platform enhanced with broader coverage, improved performance, API functionality and decision support. ", tags$i("Nucleic acids research"), ", 52(W1), W422-W431.")
                               )
                      )
             )
      )
    ),

    fluidRow(
      column(12, align = "center",
             tags$h3(style = "margin-bottom: 25px; color: #2c3e50;", "Quick Access")
      ),
      column(4, actionButton("go_cdk", label = tagList(icon("flask"), br(), "CDK & webchem"), class = "nav-card-btn")),
      column(4, actionButton("go_master", label = tagList(icon("table"), br(), "ADMET Master Manager"), class = "nav-card-btn")),
      column(4, actionButton("go_report", label = tagList(icon("file-lines"), br(), "Generate Report"), class = "nav-card-btn",
                                 style = "background-color: #4e73df; color: white; border-color: #4e73df;"))
    ),

    br(), br()
  )
}
