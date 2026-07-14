# ---------------------------------------------------------------------------
# ui_deeppk.R
# Deep-PK Manager tab.
# ---------------------------------------------------------------------------

#' Deep-PK Manager tab UI
#'
#' Builds the Deep-PK Manager tab of the ADMETShiny application.
#'
#' @return A \code{shiny.tabPanel}.
#' @keywords internal
deeppk_tab <- function() {
  tabPanel(
    title = tagList(icon("pills"), "Deep-PK Manager"),
    value = "deeppk",

    sidebarLayout(
      sidebarPanel(
        fileInput("deeppk_archivo", "Upload Deep-PK CSV", accept = ".csv"),
        tags$hr(),

        checkboxGroupInput("deeppk_filters", "Drug-likeness Filters",
                           choices = c("Lipinski", "Veber", "Ghose", "Egan", "Muegge"),
                           selected = "Lipinski"),

        conditionalPanel(condition = "input.deeppk_filters.includes('Lipinski')",
                         tags$b("Lipinski"),
                         numericInput("deeppk_mw", "Max. MW", 500),
                         numericInput("deeppk_logp", "Max. LogP", 5),
                         numericInput("deeppk_hba", "Max. HBA", 10),
                         numericInput("deeppk_hbd", "Max. HBD", 5),
                         numericInput("deeppk_violations", "Max Violations", 0)),

        conditionalPanel(condition = "input.deeppk_filters.includes('Veber')",
                         tags$b("Veber"),
                         numericInput("deeppk_v_rb", "Max. Rotatable Bonds", 10),
                         numericInput("deeppk_v_tpsa", "Max. TPSA", 140),
                         numericInput("deeppk_v_hb_sum", "Max. HBA + HBD", 12),
                         numericInput("deeppk_veber_violations", "Max Violations", 0)),

        conditionalPanel(condition = "input.deeppk_filters.includes('Ghose')",
                         tags$b("Ghose"),
                         sliderInput("deeppk_g_mw", "MW range", min = 100, max = 700, value = c(160, 480)),
                         sliderInput("deeppk_g_mr", "MR range", min = 0, max = 200, value = c(40, 130)),
                         sliderInput("deeppk_g_logp", "LogP range", min = -5, max = 10, value = c(-0.4, 5.6)),
                         sliderInput("deeppk_g_ha", "Heavy Atoms range", min = 0, max = 100, value = c(20, 70)),
                         numericInput("deeppk_ghose_violations", "Max Violations", 0)),

        conditionalPanel(condition = "input.deeppk_filters.includes('Egan')",
                         tags$b("Egan"),
                         numericInput("deeppk_e_tpsa", "Max. TPSA", 131.6),
                         numericInput("deeppk_e_logp", "Max. LogP", 5.88),
                         numericInput("deeppk_egan_violations", "Max Violations", 0)),

        conditionalPanel(condition = "input.deeppk_filters.includes('Muegge')",
                         tags$b("Muegge"),
                         sliderInput("deeppk_m_mw", "MW range", min = 50, max = 800, value = c(200, 600)),
                         sliderInput("deeppk_m_logp", "LogP range", min = -10, max = 10, value = c(-2, 5)),
                         numericInput("deeppk_m_hba", "Max. HBA", 10),
                         numericInput("deeppk_m_hbd", "Max. HBD", 5),
                         numericInput("deeppk_m_rb", "Max. Rotatable Bonds", 15),
                         numericInput("deeppk_m_tpsa", "Max. TPSA", 150),
                         numericInput("deeppk_muegge_violations", "Max Violations", 0)),

        tags$hr(),
        actionButton("deeppk_run", "Apply Filter", class = "btn-primary"),
        br(), br(),
        downloadButton("deeppk_download", "Download filtered dataset")
      ),

      mainPanel(tabsetPanel(
        tabPanel("Data Preview", br(), h4("Deep-PK dataset"), DTOutput("deeppk_preview")),
        tabPanel("Filtered Molecules", br(), h4("Filtered Results"), DTOutput("deeppk_tabla")),
        tabPanel("Plots", br(),
                 selectInput("deeppk_plot_type", "Select plot", choices = c(
                   "Boiled Egg", "TPSA", "LogP",
                   "Radar plot (Chemical profile)",
                   "Tanimoto / AGNES (Structural similarity)",
                   "Correlation Heatmap",
                   "Principal Component Analysys (PCA - Chemical space)",
                   "t-SNE (Chemical space)",
                   "UMAP (Chemical space)",
                   "Cluster Heatmap (Dendrogram)",
                   "Parallel Coordinates",
                   "Violin Plot"
                 )),

                 conditionalPanel(condition = "input.deeppk_plot_type == 'Radar plot (Chemical profile)'",
                                  tags$div(class = "alert-primary", icon("book-open"), tags$b(" Function: "),
                                           "Simultaneously compares multiple physicochemical properties of one or more compounds."),
                                  tags$div(class = "alert-primary", icon("book-open"), tags$b(" Analysis: "),
                                           "Each axis corresponds to a descriptor. Similar profiles indicate similar properties; deviations highlight differences between molecules."),
                                  tags$hr(),
                                  uiOutput("deeppk_radar_id_selector")),

                 conditionalPanel(condition = "input.deeppk_plot_type == 'Tanimoto / AGNES (Structural similarity)'",
                                  tags$div(class = "alert-primary", icon("book-open"), tags$b(" Function: "),
                                           "Groups compounds according to their structural or physicochemical similarity."),
                                  tags$div(class = "alert-primary", icon("book-open"), tags$b(" Analysis: "),
                                           "The dendrogram shows the relationship between molecules; branches close together indicate greater similarity."),
                                  tags$div(class = "alert-primary", icon("book-open"), tags$b(" Function: "),
                                           "Quantifies the structural similarity between molecules using molecular fingerprints."),
                                  tags$div(class = "alert-danger", icon("binoculars"), tags$b(" Analysis: "),
                                           "Values close to 1 indicate high structural similarity, while values close to 0 reflect different structures."),
                                  tags$hr(),
                                  uiOutput("deeppk_smiles_col_selector"),
                                  uiOutput("deeppk_label_col_selector"),
                                  numericInput("deeppk_tanimoto_max_n", "Max molecules", 40, min = 5, max = 200),
                                  selectInput("deeppk_agnes_method", "AGNES method",
                                              choices = c("average", "complete", "single", "ward"), selected = "average")),

                 conditionalPanel(condition = "input.deeppk_plot_type =='Parallel Coordinates'",
                                  tags$div(class = "alert-primary", icon("book-open"), tags$b(" Function: "),
                                           "Groups compounds according to their structural or physicochemical similarity."),
                                  tags$div(class = "alert-primary", icon("book-open"), tags$b(" Analysis: "),
                                           "The dendrogram shows the relationship between molecules; branches close together indicate greater similarity."),
                                  tags$hr(),
                                  uiOutput("deeppk_parallel_controls")),

                 conditionalPanel(condition = "input.deeppk_plot_type == 'Principal Component Analysys (PCA - Chemical space)'",
                                  tags$div(class = "alert-primary", icon("book-open"), tags$b(" Function: "),
                                           "Reduces the dimensionality of data to identify global patterns."),
                                  tags$div(class = "alert-primary", icon("book-open"), tags$b(" Analysis: "),
                                           "Compounds located close together exhibit similar physicochemical profiles, while those located far apart show different characteristics."),
                                  tags$hr(),
                                  uiOutput("deeppk_pca_controls")),

                 conditionalPanel(condition = "input.deeppk_plot_type == 'Violin Plot'",
                                  tags$div(class = "alert-primary", icon("book-open"), tags$b(" Function: "),
                                           "Displays the distribution and density of a descriptor."),
                                  tags$div(class = "alert-primary", icon("book-open"), tags$b(" Analysis: "),
                                           "Wider areas indicate where more observations are concentrated; allows for comparison of variability between groups."),
                                  tags$hr(),
                                  uiOutput("deeppk_violin_controls")),

                 conditionalPanel(condition = "input.deeppk_plot_type == 't-SNE (Chemical space)'",
                                  tags$div(class = "alert-primary", icon("book-open"), tags$b(" Function: "),
                                           "Groups compounds according to the similarity of their properties while preserving local relationships."),
                                  tags$div(class = "alert-primary", icon("book-open"), tags$b(" Analysis: "),
                                           "Tight clusters represent similar molecules; the distance between groups reflects differences in their profiles."),
                                  tags$hr(),
                                  uiOutput("deeppk_tsne_controls")),

                 conditionalPanel(condition = "input.deeppk_plot_type == 'UMAP (Chemical space)'",
                                  tags$div(class = "alert-primary", icon("book-open"), tags$b(" Function: "),
                                           "Non-linear dimensionality reduction that preserves both local and global structure of the chemical space."),
                                  tags$div(class = "alert-primary", icon("book-open"), tags$b(" Analysis: "),
                                           "Nearby points correspond to molecules with similar physicochemical profiles; well-separated clusters indicate distinct chemical series."),
                                  tags$hr(),
                                  uiOutput("deeppk_umap_controls")),

                 conditionalPanel(condition = "input.deeppk_plot_type == 'Cluster Heatmap (Dendrogram)'",
                                  tags$div(class = "alert-primary", icon("book-open"), tags$b(" Function: "),
                                           "Combines hierarchical clustering with a heatmap to reveal groups of similar molecules."),
                                  tags$div(class = "alert-primary", icon("book-open"), tags$b(" Analysis: "),
                                           "Branches close together indicate similar compounds; the color scale shows z-scored property values."),
                                  tags$hr(),
                                  uiOutput("deeppk_cluster_heatmap_controls")),

                 palette_selector_ui("deeppk_palette"),
                 br(),
                 plotOutput("deeppk_plot", height = "650px"),
                 tags$hr(),
                 h4("Export Figure"),
                 fluidRow(
                   column(4, selectInput("deeppk_format", "Format", choices = c("png", "pdf", "svg", "jpeg", "tiff"), selected = "png")),
                   column(4, selectInput("deeppk_dpi", "Resolution", choices = c("300 dpi" = 300, "600 dpi" = 600, "1200 dpi" = 1200), selected = 600))
                 ),
                 fluidRow(
                   column(4, sliderInput("deeppk_width", "Width (in)", min = 4, max = 12, value = 7, step = 0.5)),
                   column(4, sliderInput("deeppk_height", "Height (in)", min = 4, max = 12, value = 6, step = 0.5))
                 ),
                 br(),
                 downloadButton("deeppk_downloadPlot", "Download Figure", class = "btn-success")
        )
      ))
    )
  )
}
