# ---------------------------------------------------------------------------
# ui_admetlab.R
# ADMETlab 3.0 Manager tab.
# ---------------------------------------------------------------------------

#' ADMETlab 3.0 Manager tab UI
#'
#' Builds the ADMETlab 3.0 Manager tab of the ADMETShiny application.
#'
#' @return A \code{shiny.tabPanel}.
#' @keywords internal
admetlab_tab <- function() {
  tabPanel(
    title = tagList(icon("syringe"), "ADMETlab 3.0 Manager"),
    value = "admetlab",

    sidebarLayout(

      sidebarPanel(
        fileInput("admetlab_archivo", "Upload ADMETlab CSV", accept = ".csv"),
        tags$hr(),

        checkboxGroupInput("admetlab_filters", "Drug-likeness Filters",
                           choices = c("Lipinski", "Veber", "Ghose", "Egan", "Muegge"),
                           selected = "Lipinski"),

        conditionalPanel(condition = "input.admetlab_filters.includes('Lipinski')",
                         tags$b("Lipinski"),
                         numericInput("admetlab_mw", "Max. Molecular Weight", 500),
                         numericInput("admetlab_logp", "Max. LogP", 5),
                         numericInput("admetlab_hba", "Max. H-bond acceptors", 10),
                         numericInput("admetlab_hbd", "Max. H-bond donors", 5),
                         numericInput("admetlab_violations", "Maximum Lipinski Violations", 0)),

        conditionalPanel(condition = "input.admetlab_filters.includes('Veber')",
                         tags$b("Veber"),
                         numericInput("admetlab_v_rb", "Max. Rotatable Bonds", 10),
                         numericInput("admetlab_v_tpsa", "Max. TPSA", 140),
                         numericInput("admetlab_v_hb_sum", "Max. HBA + HBD", 12),
                         numericInput("admetlab_veber_violations", "Maximum Veber Violations", 0)),

        conditionalPanel(condition = "input.admetlab_filters.includes('Ghose')",
                         tags$b("Ghose"),
                         sliderInput("admetlab_g_mw", "Molecular Weight range", min = 100, max = 700, value = c(160, 480)),
                         sliderInput("admetlab_g_mr", "Molar Refractivity range", min = 0, max = 200, value = c(40, 130)),
                         sliderInput("admetlab_g_logp", "LogP range", min = -5, max = 10, value = c(-0.4, 5.6)),
                         sliderInput("admetlab_g_ha", "Heavy Atoms range", min = 0, max = 100, value = c(20, 70)),
                         numericInput("admetlab_ghose_violations", "Maximum Ghose Violations", 0)),

        conditionalPanel(condition = "input.admetlab_filters.includes('Egan')",
                         tags$b("Egan"),
                         numericInput("admetlab_e_tpsa", "Max. TPSA", 131.6),
                         numericInput("admetlab_e_logp", "Max. LogP", 5.88),
                         numericInput("admetlab_egan_violations", "Maximum Egan Violations", 0)),

        conditionalPanel(condition = "input.admetlab_filters.includes('Muegge')",
                         tags$b("Muegge"),
                         sliderInput("admetlab_m_mw", "Molecular Weight range", min = 50, max = 800, value = c(200, 600)),
                         sliderInput("admetlab_m_logp", "LogP range", min = -10, max = 10, value = c(-2, 5)),
                         numericInput("admetlab_m_hba", "Max. H-bond acceptors", 10),
                         numericInput("admetlab_m_hbd", "Max. H-bond donors", 5),
                         numericInput("admetlab_m_rb", "Max. Rotatable Bonds", 15),
                         numericInput("admetlab_m_tpsa", "Max. TPSA", 150),
                         numericInput("admetlab_muegge_violations", "Maximum Muegge Violations", 0)),

        tags$hr(),
        actionButton("admetlab_run", "Apply Filter", class = "btn-primary"),
        br(), br(),
        downloadButton("admetlab_download", "Download filtered dataset")
      ),

      mainPanel(
        tabsetPanel(
          tabPanel("Data Preview", br(), h4("ADMETlab dataset"), DTOutput("admetlab_preview")),
          tabPanel("Filtered Molecules", br(), h4("Filtered Results"), DTOutput("admetlab_tabla")),
          tabPanel("Plots", br(),

                   selectInput("admetlab_plot_type", "Select plot", choices = c(
                     "Boiled Egg", "TPSA", "LogP",
                     "Radar plot (Chemical profile)",
                     "Tanimoto / AGNES (Structural similarity)",
                     "Correlation Heatmap",
                     "Principal Component Analysys (PCA - Chemical space)",
                     "t-SNE (Chemical space)",
                     "Cluster Heatmap (Dendrogram)",
                     "Parallel Coordinates",
                     "Violin Plot"
                   )),

                   conditionalPanel(condition = "input.admetlab_plot_type == 'Radar plot (Chemical profile)'",
                                    tags$div(class = "alert-primary", icon("book-open"), tags$b(" Function: "),
                                             "Simultaneously compares multiple physicochemical properties of one or more compounds."),
                                    tags$div(class = "alert-primary", icon("book-open"), tags$b(" Analysis: "),
                                             "Each axis corresponds to a descriptor. Similar profiles indicate similar properties; deviations highlight differences between molecules."),
                                    tags$hr(),
                                    uiOutput("admetlab_radar_id_selector")),

                   conditionalPanel(condition = "input.admetlab_plot_type == 'Tanimoto / AGNES (Structural similarity)'",
                                    tags$div(class = "alert-primary", icon("book-open"), tags$b(" Function: "),
                                             "Groups compounds according to their structural or physicochemical similarity."),
                                    tags$div(class = "alert-primary", icon("book-open"), tags$b(" Analysis: "),
                                             "The dendrogram shows the relationship between molecules; branches close together indicate greater similarity."),
                                    tags$div(class = "alert-primary", icon("book-open"), tags$b(" Function: "),
                                             "Quantifies the structural similarity between molecules using molecular fingerprints."),
                                    tags$div(class = "alert-danger", icon("binoculars"), tags$b(" Analysis: "),
                                             "Values close to 1 indicate high structural similarity, while values close to 0 reflect different structures."),
                                    tags$hr(),
                                    uiOutput("admetlab_smiles_col_selector"),
                                    uiOutput("admetlab_label_col_selector"),
                                    numericInput("admetlab_tanimoto_max_n", "Max molecules", 40, min = 5, max = 200),
                                    selectInput("admetlab_agnes_method", "AGNES method",
                                                choices = c("average", "complete", "single", "ward"), selected = "average")),

                   conditionalPanel(condition = "input.admetlab_plot_type =='Parallel Coordinates'",
                                    tags$div(class = "alert-primary", icon("book-open"), tags$b(" Function: "),
                                             "Groups compounds according to their structural or physicochemical similarity."),
                                    tags$div(class = "alert-primary", icon("book-open"), tags$b(" Analysis: "),
                                             "The dendrogram shows the relationship between molecules; branches close together indicate greater similarity."),
                                    tags$hr(),
                                    uiOutput("admetlab_parallel_controls")),

                   conditionalPanel(condition = "input.admetlab_plot_type == 'Principal Component Analysys (PCA - Chemical space)'",
                                    tags$div(class = "alert-primary", icon("book-open"), tags$b(" Function: "),
                                             "Reduces the dimensionality of data to identify global patterns."),
                                    tags$div(class = "alert-primary", icon("book-open"), tags$b(" Analysis: "),
                                             "Compounds located close together exhibit similar physicochemical profiles, while those located far apart show different characteristics."),
                                    tags$hr(),
                                    uiOutput("admetlab_pca_controls")),

                   conditionalPanel(condition = "input.admetlab_plot_type == 'Violin Plot'",
                                    tags$div(class = "alert-primary", icon("book-open"), tags$b(" Function: "),
                                             "Displays the distribution and density of a descriptor."),
                                    tags$div(class = "alert-primary", icon("book-open"), tags$b(" Analysis: "),
                                             "Wider areas indicate where more observations are concentrated; allows for comparison of variability between groups."),
                                    tags$hr(),
                                    uiOutput("admetlab_violin_controls")),

                   conditionalPanel(condition = "input.admetlab_plot_type == 't-SNE (Chemical space)'",
                                    tags$div(class = "alert-primary", icon("book-open"), tags$b(" Function: "),
                                             "Groups compounds according to the similarity of their properties while preserving local relationships."),
                                    tags$div(class = "alert-primary", icon("book-open"), tags$b(" Analysis: "),
                                             "Tight clusters represent similar molecules; the distance between groups reflects differences in their profiles."),
                                    tags$hr(),
                                    uiOutput("admetlab_tsne_controls")),

                   conditionalPanel(condition = "input.admetlab_plot_type == 'Cluster Heatmap (Dendrogram)'",
                                    tags$div(class = "alert-primary", icon("book-open"), tags$b(" Function: "),
                                             "Combines hierarchical clustering with a heatmap to reveal groups of similar molecules."),
                                    tags$div(class = "alert-primary", icon("book-open"), tags$b(" Analysis: "),
                                             "Branches close together indicate similar compounds; the color scale shows z-scored property values."),
                                    tags$hr(),
                                    uiOutput("admetlab_cluster_heatmap_controls")),

                   palette_selector_ui("admetlab_palette"),
                   br(),
                   plotOutput("admetlab_plot", height = "650px"),
                   tags$hr(),
                   h4("Export Figure"),
                   fluidRow(
                     column(4, selectInput("admetlab_format", "Format", choices = c("png", "pdf", "svg", "jpeg", "tiff"), selected = "png")),
                     column(4, selectInput("admetlab_dpi", "Resolution", choices = c("300 dpi" = 300, "600 dpi" = 600, "1200 dpi" = 1200), selected = 600))
                   ),
                   fluidRow(
                     column(4, sliderInput("admetlab_width", "Width (in)", min = 4, max = 12, value = 7, step = 0.5)),
                     column(4, sliderInput("admetlab_height", "Height (in)", min = 4, max = 12, value = 6, step = 0.5))
                   ),
                   br(),
                   downloadButton("admetlab_downloadPlot", "Download Figure", class = "btn-success")
          )
        )
      )
    )
  )
}
