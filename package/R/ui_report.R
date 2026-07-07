# ---------------------------------------------------------------------------
# ui_report.R
# Report tab UI: generate and download comprehensive analysis reports.
# ---------------------------------------------------------------------------

#' Report tab UI
#'
#' Builds the Report tab of the ADMETShiny application, which generates a
#' comprehensive report of all analyses performed across the four data source
#' modules (SwissADME Manager, CDK & webchem, ADMETlab 3.0, Deep-PK Manager).
#'
#' @return A \code{shiny.tabPanel}.
#' @keywords internal
report_tab <- function() {
  tabPanel(
    title = tagList(icon("file-lines"), "Report"),
    value = "report",

    info_card_css(),

    sidebarLayout(

      sidebarPanel(

        width = 3,

        tags$div(class = "info-card accent-blue",
                 tags$div(class = "info-card-header", icon("file-lines"), " Report Generation"),
                 tags$div(class = "info-card-body",
                          tags$p(style = "color:#555; font-size: 14px;",
                                 "Generate a comprehensive report of all analyses performed across the four data source modules.")
                 )
        ),

        tags$hr(),

        h5("Generate"),
        actionButton("generate_report_btn", "Generate Preview",
                     icon = icon("eye"), class = "btn-primary",
                     width = "100%"),
        br(), br(),

        tags$hr(),

        h5("Download Report"),

        downloadButton("download_report_pdf", "Download PDF",
                       class = "btn-success", style = "width: 100%;"),
        br(), br(),
        downloadButton("download_report_doc", "Download Word (.docx)",
                       class = "btn-success", style = "width: 100%;"),
        br(), br(),
        downloadButton("download_report_html", "Download HTML",
                       class = "btn-success", style = "width: 100%;"),

        tags$hr(),

        tags$div(
          class = "alert alert-warning",
          icon("circle-info"),
          tags$b(" Note: "),
          "PDF export requires a LaTeX installation (e.g., ",
          tags$code("tinytex::install_tinytex()"),
          "). Word export requires pandoc (bundled with RStudio)."
        ),

        tags$hr(),

        h5("Report Contents"),
        tags$ul(
          style = "color:#555; font-size: 13px; padding-left: 18px;",
          tags$li("Executive summary with all datasets"),
          tags$li("Per-dataset statistics (MW, LogP, TPSA, etc.)"),
          tags$li("Drug-likeness filter pass rates"),
          tags$li("BOILED-Egg ADMET classification"),
          tags$li("Additional metrics:"),
          tags$ul(
            style = "padding-left: 18px;",
            tags$li("Pfizer 3/75 rule (Hughes 2008)"),
            tags$li("GSK 4/400 rule (Gleeson 2011)"),
            tags$li("Lead-likeness (Teague 1999)"),
            tags$li("Egan bioavailability proxy"),
            tags$li("Lipinski Ro5 (original)"),
            tags$li("Veber oral bioavailability")
          ),
          tags$li("Composite drug-likeness score (0-100)"),
          tags$li("Cross-dataset comparison"),
          tags$li("Visualizations (BOILED-Egg, distributions, heatmaps)"),
          tags$li("Literature references")
        )
      ),

      mainPanel(

        width = 9,

        uiOutput("report_preview_container")
      )
    )
  )
}
