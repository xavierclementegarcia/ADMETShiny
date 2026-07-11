# ---------------------------------------------------------------------------
# ui_swissadme.R
# SwissADME Manager tab: upload, filter, table and plots.
# ---------------------------------------------------------------------------

#' SwissADME Manager tab UI
#'
#' Builds the SwissADME Manager tab of the ADMETShiny application.
#'
#' @return A \code{shiny.tabPanel}.
#' @keywords internal
swissadme_tab <- function() {
  tabPanel(
    title = tagList(icon("flask"), "SwissADME Manager"),
    value = "swissadme",

    sidebarLayout(

      sidebarPanel(

        fileInput("archivo", "Upload SwissADME dataset", accept = ".csv"),

        tags$hr(),

        checkboxGroupInput(
          "filters", "Drug-likeness Filters",
          choices = c("Lipinski", "Veber", "Ghose", "Egan", "Muegge"),
          selected = "Lipinski"
        ),

        conditionalPanel(
          condition = "input.filters.includes('Lipinski')",
          tags$b("Lipinski"),
          numericInput("mw", "Max. Molecular Weight", 500),
          numericInput("logp", "Max. LogP", 5),
          numericInput("hba", "Max. H-bond acceptors", 10),
          numericInput("hbd", "Max. H-bond donors", 5),
          numericInput("violations", "Maximum Lipinski Violations", 0)
        ),

        conditionalPanel(
          condition = "input.filters.includes('Veber')",
          tags$b("Veber"),
          numericInput("v_rb", "Max. Rotatable Bonds", 10),
          numericInput("v_tpsa", "Max. TPSA", 140),
          numericInput("v_hb_sum", "Max. HBA + HBD", 12),
          numericInput("veber_violations", "Maximum Veber Violations", 0)
        ),

        conditionalPanel(
          condition = "input.filters.includes('Ghose')",
          tags$b("Ghose"),
          sliderInput("g_mw", "Molecular Weight range", min = 100, max = 700, value = c(160, 480)),
          sliderInput("g_mr", "Molar Refractivity range", min = 0, max = 200, value = c(40, 130)),
          sliderInput("g_logp", "LogP range", min = -5, max = 10, value = c(-0.4, 5.6)),
          sliderInput("g_ha", "Heavy Atoms range", min = 0, max = 100, value = c(20, 70)),
          numericInput("ghose_violations", "Maximum Ghose Violations", 0)
        ),

        conditionalPanel(
          condition = "input.filters.includes('Egan')",
          tags$b("Egan"),
          numericInput("e_tpsa", "Max. TPSA", 131.6),
          numericInput("e_logp", "Max. LogP", 5.88),
          numericInput("egan_violations", "Maximum Egan Violations", 0)
        ),

        conditionalPanel(
          condition = "input.filters.includes('Muegge')",
          tags$b("Muegge"),
          sliderInput("m_mw", "Molecular Weight range", min = 50, max = 800, value = c(200, 600)),
          sliderInput("m_logp", "LogP range", min = -10, max = 10, value = c(-2, 5)),
          numericInput("m_hba", "Max. H-bond acceptors", 10),
          numericInput("m_hbd", "Max. H-bond donors", 5),
          numericInput("m_rb", "Max. Rotatable Bonds", 15),
          numericInput("m_tpsa", "Max. TPSA", 150),
          numericInput("muegge_violations", "Maximum Muegge Violations", 0)
        ),

        tags$hr(),
        actionButton("run", "Apply Filter", class = "btn-primary"),
        br(), br(),
        downloadButton("download", "Download filtered dataset")
      ),

      mainPanel(

        tabsetPanel(

          tabPanel("Data Preview", br(), h4("SwissADME dataset"), DTOutput("preview")),

          tabPanel("Filtered Molecules", br(), h4("Filtered Results"), DTOutput("tabla")),

          tabPanel(
            "Plots", br(),

            selectInput("plot_type", "Select plot", choices = c(
              "Boiled Egg", "Molecular Weight", "TPSA", "LogP",
              "Radar plot (Chemical profile)",
              "Tanimoto / AGNES (Structural similarity)",
              "Correlation Heatmap",
              "Principal Component Analysys (PCA - Chemical space)",
              "t-SNE (Chemical space)",
              "Cluster Heatmap (Dendrogram)",
              "Parallel Coordinates",
              "Violin Plot"
            )),

            conditionalPanel(
              condition = "input.plot_type == 'Radar plot (Chemical profile)'",
              tags$div(class = "alert-primary", icon("book-open"), tags$b(" Function: "),
                       "Simultaneously compares multiple physicochemical properties of one or more compounds."),
              tags$div(class = "alert-primary", icon("book-open"), tags$b(" Analysis: "),
                       "Each axis corresponds to a descriptor. Similar profiles indicate similar properties; deviations highlight differences between molecules."),
              tags$hr(),
              uiOutput("radar_id_selector"),
              helpText("Select an identifier column (e.g., name or SMILES). Max 5 molecules in the radar plot")
            ),

            conditionalPanel(
              condition = "input.plot_type == 'Tanimoto / AGNES (Structural similarity)'",
              tags$div(class = "alert-primary", icon("book-open"), tags$b(" Function: "),
                       "Groups compounds according to their structural or physicochemical similarity."),
              tags$div(class = "alert-primary", icon("book-open"), tags$b(" Analysis: "),
                       "The dendrogram shows the relationship between molecules; branches close together indicate greater similarity."),
              tags$div(class = "alert-primary", icon("book-open"), tags$b(" Function: "),
                       "Quantifies the structural similarity between molecules using molecular fingerprints."),
              tags$div(class = "alert-danger", icon("binoculars"), tags$b(" Analysis: "),
                       "Values close to 1 indicate high structural similarity, while values close to 0 reflect different structures."),
              tags$hr(),
              uiOutput("smiles_col_selector"),
              uiOutput("label_col_selector"),
              numericInput("tanimoto_max_n", "Max molecules comparition", 40, min = 5, max = 200),
              selectInput("agnes_method", "Linking method (AGNES)",
                          choices = c("average", "complete", "single", "ward"), selected = "average"),
              helpText("Select the column that contains the SMILES for the similarity calculation and the column you want to use as dendrogram labels.")
            ),

            conditionalPanel(
              condition = "input.plot_type =='Parallel Coordinates'",
              tags$div(class = "alert-primary", icon("book-open"), tags$b(" Function: "),
                       "Groups compounds according to their structural or physicochemical similarity."),
              tags$div(class = "alert-primary", icon("book-open"), tags$b(" Analysis: "),
                       "The dendrogram shows the relationship between molecules; branches close together indicate greater similarity."),
              tags$hr(),
              uiOutput("parallel_controls")
            ),

            conditionalPanel(
              condition = "input.plot_type == 'Principal Component Analysys (PCA - Chemical space)'",
              tags$div(class = "alert-primary", icon("book-open"), tags$b(" Function: "),
                       "Reduces the dimensionality of data to identify global patterns."),
              tags$div(class = "alert-primary", icon("book-open"), tags$b(" Analysis: "),
                       "Compounds located close together exhibit similar physicochemical profiles, while those located far apart show different characteristics."),
              tags$hr(),
              uiOutput("pca_controls")
            ),

            conditionalPanel(
              condition = "input.plot_type == 'Violin Plot'",
              tags$div(class = "alert-primary", icon("book-open"), tags$b(" Function: "),
                       "Displays the distribution and density of a descriptor."),
              tags$div(class = "alert-primary", icon("book-open"), tags$b(" Analysis: "),
                       "Wider areas indicate where more observations are concentrated; allows for comparison of variability between groups."),
              tags$hr(),
              uiOutput("violin_controls")
            ),

            conditionalPanel(
              condition = "input.plot_type == 't-SNE (Chemical space)'",
              tags$div(class = "alert-primary", icon("book-open"), tags$b(" Function: "),
                       "Groups compounds according to the similarity of their properties while preserving local relationships."),
              tags$div(class = "alert-primary", icon("book-open"), tags$b(" Analysis: "),
                       "Tight clusters represent similar molecules; the distance between groups reflects differences in their profiles."),
              tags$hr(),
              uiOutput("tsne_controls")
            ),

            conditionalPanel(
              condition = "input.plot_type == 'Cluster Heatmap (Dendrogram)'",
              tags$div(class = "alert-primary", icon("book-open"), tags$b(" Function: "),
                       "Combines hierarchical clustering with a heatmap to reveal groups of similar molecules."),
              tags$div(class = "alert-primary", icon("book-open"), tags$b(" Analysis: "),
                       "Branches close together indicate similar compounds; the color scale shows z-scored property values."),
              tags$hr(),
              uiOutput("cluster_heatmap_controls")
            ),

            palette_selector_ui("swiss_palette"),
            br(),

            plotOutput("plot", height = "650px"),
            tags$hr(),
            h4("Export Figure"),
            fluidRow(
              column(4, selectInput("format", "Format",
                                    choices = c("png", "pdf", "svg", "jpeg", "tiff"), selected = "png")),
              column(4, selectInput("dpi", "Resolution",
                                    choices = c("300 dpi" = 300, "600 dpi" = 600, "1200 dpi" = 1200), selected = 600))
            ),
            fluidRow(
              column(4, sliderInput("width", "Width (in)", min = 4, max = 12, value = 7, step = 0.5)),
              column(4, sliderInput("height", "Height (in)", min = 4, max = 12, value = 6, step = 0.5))
            ),
            br(),
            downloadButton("downloadPlot", "Download Figure", class = "btn-success")
          )
        )
      )
    )
  )
}
