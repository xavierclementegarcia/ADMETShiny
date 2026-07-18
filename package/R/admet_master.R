# ---------------------------------------------------------------------------
# admet_master.R
# ADMET Master Manager tab: a wizard-style UI that lets the user upload any
# CSV/XLSX dataset, map its columns to the application's standard schema,
# optionally calculate missing descriptors with CDK, run the standard
# drug-likeness filters and produce the same plot catalogue as the other
# modules.
# ---------------------------------------------------------------------------

#' ADMET Master Manager tab UI
#'
#' Builds the "ADMET Master Manager" tab of the ADMETShiny application. This
#' module is a four-step wizard (Upload & Preview -> Map Columns -> Filter ->
#' Plots) that lets users bring in any tabular dataset (CSV or Excel) and map
#' its columns to the application's standard schema, so that the same
#' drug-likeness filters and plots as the other modules can be applied to data
#' coming from any source.
#'
#' @return A \code{shiny.tabPanel}.
#' @keywords internal
#' @seealso \code{\link{mapADMETColumns}}, \code{\link{detectColumnTypes}},
#'   \code{\link{detectSMILESColumn}}.
admet_master_tab <- function() {
  tabPanel(
    title = tagList(icon("wand-magic-sparkles"), "ADMET Master Manager"),
    value = "admet_master",

    fluidPage(

      h2("ADMET Master Manager"),
      p("Upload any molecular dataset (CSV or Excel), map its columns to the",
        " application's standard schema, optionally calculate missing descriptors",
        " with CDK, then apply the standard drug-likeness filters and produce",
        " the same rich set of plots as the other modules."),

      tags$div(
        class = "alert alert-info",
        icon("circle-info"),
        tags$b(" How it works: "),
        "Step 1 uploads your data and inspects it. Step 2 maps your columns to",
        " the standard schema (SMILES, MW, LogP, TPSA, HBD, HBA, ...). Step 3",
        " applies the drug-likeness filters. Step 4 visualises the filtered",
        " results with the same plots available in the SwissADME / CDK /",
        " ADMETlab / Deep-PK modules."
      ),

      tags$hr(),

      tabsetPanel(

        ## ===================================================================
        ## Tab 1: Upload & Preview
        ## ===================================================================
        tabPanel(
          title = tagList(icon("upload"), " 1. Upload & Preview"),
          value = "master_step1",

          br(),

          fluidRow(
            column(6,
                   fileInput("master_file", "Upload dataset",
                             accept = c(".csv", ".xlsx", ".xls",
                                        "text/csv",
                                        "text/comma-separated-values",
                                        "application/vnd.ms-excel",
                                        "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"))
            ),
            column(6,
                   tags$div(
                     class = "alert alert-warning",
                     icon("triangle-exclamation"),
                     tags$b(" Note: "),
                     "Accepted formats are CSV and Excel (.xlsx / .xls).",
                     " For Excel files the first sheet is read by default."
                   )
            )
          ),

          tags$hr(),

          h4("Detected SMILES column"),
          uiOutput("master_smiles_detection"),

          tags$hr(),

          h4("Column types"),
          p("The table below lists every column detected in the uploaded file,",
            " its inferred type (string / numeric), the number of unique",
            " values and a small sample."),
          DTOutput("master_column_info"),

          tags$hr(),

          h4("Raw data preview"),
          DTOutput("master_preview")
        ),

        ## ===================================================================
        ## Tab 2: Map Columns
        ## ===================================================================
        tabPanel(
          title = tagList(icon("arrows-left-right-to-line"), " 2. Map Columns"),
          value = "master_step2",

          br(),

          tags$div(
            class = "alert alert-info",
            icon("circle-info"),
            tags$b(" Instructions: "),
            "For each column in your dataset, pick the standard ADMET field",
            " it corresponds to. Columns mapped to \"None (skip)\" are kept in",
            " the dataset but ignored by the filters and plots. Required",
            " fields for the drug-likeness filters are: MW, LogP, TPSA,",
            " #H-bond donors, #H-bond acceptors, #Rotatable bonds,",
            " Molar Refractivity, #Heavy atoms."
          ),

          h4("Column mapping"),
          uiOutput("master_mapping_ui"),

          tags$hr(),

          fluidRow(
            column(12,
                   checkboxInput("master_calculate_cdk",
                                 "Calculate missing descriptors with CDK (requires SMILES + Java)",
                                 value = TRUE),
                   helpText("When checked, any missing MW, LogP (ALogP), TPSA,",
                            " HBD, HBA, rotatable bonds, MR, #Heavy atoms and",
                            " #Aromatic heavy atoms are calculated from the",
                            " mapped SMILES column using the Chemistry",
                            " Development Kit (via the rcdk package).")
            )
          ),

          br(),
          actionButton("master_calculate", "Calculate & Standardize",
                       icon = icon("calculator"), class = "btn-primary"),

          tags$hr(),

          h4("Mapping summary"),
          DTOutput("master_mapping_summary"),

          tags$hr(),

          h4("Filter availability"),
          uiOutput("master_validation")
        ),

        ## ===================================================================
        ## Tab 3: Filter
        ## ===================================================================
        tabPanel(
          title = tagList(icon("filter"), " 3. Filter"),
          value = "master_step3",

          br(),

          sidebarLayout(

            sidebarPanel(

              checkboxGroupInput(
                "master_filters", "Drug-likeness Filters",
                choices = c("Lipinski", "Veber", "Ghose", "Egan", "Muegge"),
                selected = "Lipinski"
              ),

              conditionalPanel(
                condition = "input.master_filters.includes('Lipinski')",
                tags$b("Lipinski"),
                numericInput("master_mw", "Max. Molecular Weight", 500),
                numericInput("master_logp", "Max. LogP", 5),
                numericInput("master_hba", "Max. H-bond acceptors", 10),
                numericInput("master_hbd", "Max. H-bond donors", 5),
                numericInput("master_violations", "Maximum Lipinski Violations", 0)
              ),

              conditionalPanel(
                condition = "input.master_filters.includes('Veber')",
                tags$b("Veber"),
                numericInput("master_v_rb", "Max. Rotatable Bonds", 10),
                numericInput("master_v_tpsa", "Max. TPSA", 140),
                numericInput("master_v_hb_sum", "Max. HBA + HBD", 12),
                numericInput("master_veber_violations", "Maximum Veber Violations", 0)
              ),

              conditionalPanel(
                condition = "input.master_filters.includes('Ghose')",
                tags$b("Ghose"),
                sliderInput("master_g_mw", "Molecular Weight range",
                            min = 100, max = 700, value = c(160, 480)),
                sliderInput("master_g_mr", "Molar Refractivity range",
                            min = 0, max = 200, value = c(40, 130)),
                sliderInput("master_g_logp", "LogP range",
                            min = -5, max = 10, value = c(-0.4, 5.6)),
                sliderInput("master_g_ha", "Heavy Atoms range",
                            min = 0, max = 100, value = c(20, 70)),
                numericInput("master_ghose_violations", "Maximum Ghose Violations", 0)
              ),

              conditionalPanel(
                condition = "input.master_filters.includes('Egan')",
                tags$b("Egan"),
                numericInput("master_e_tpsa", "Max. TPSA", 131.6),
                numericInput("master_e_logp", "Max. LogP", 5.88),
                numericInput("master_egan_violations", "Maximum Egan Violations", 0)
              ),

              conditionalPanel(
                condition = "input.master_filters.includes('Muegge')",
                tags$b("Muegge"),
                sliderInput("master_m_mw", "Molecular Weight range",
                            min = 50, max = 800, value = c(200, 600)),
                sliderInput("master_m_logp", "LogP range",
                            min = -10, max = 10, value = c(-2, 5)),
                numericInput("master_m_hba", "Max. H-bond acceptors", 10),
                numericInput("master_m_hbd", "Max. H-bond donors", 5),
                numericInput("master_m_rb", "Max. Rotatable Bonds", 15),
                numericInput("master_m_tpsa", "Max. TPSA", 150),
                numericInput("master_muegge_violations", "Maximum Muegge Violations", 0)
              ),

              tags$hr(),
              actionButton("master_run", "Apply Filter", class = "btn-primary"),
              br(), br(),
              downloadButton("master_download", "Download filtered dataset")
            ),

            mainPanel(
              br(),
              h4("Filtered Results"),
              DTOutput("master_tabla")
            )
          )
        ),

        ## ===================================================================
        ## Tab 4: Plots
        ## ===================================================================
        tabPanel(
          title = tagList(icon("chart-line"), " 4. Plots"),
          value = "master_step4",

          br(),

          selectInput("master_plot_type", "Select plot", choices = c(
            "Boiled Egg", "Molecular Weight", "TPSA", "LogP",
            "Radar plot (Chemical profile)",
            "Tanimoto / AGNES (Structural similarity)",
            "Correlation Heatmap",
            "Principal Component Analysys (PCA - Chemical space)",
            "t-SNE (Chemical space)",
            "UMAP (Chemical space)",
            "Cluster Heatmap (Dendrogram)",
            "Parallel Coordinates",
            "Violin Plot",
            "Custom Histogram"
          )),

          conditionalPanel(
            condition = "input.master_plot_type == 'Boiled Egg'",
            tags$div(class = "alert-primary", icon("book-open"), tags$b(" Function: "),
                     "Predicts gastrointestinal absorption and blood-brain barrier permeability using the BOILED-Egg model."),
            tags$div(class = "alert-primary", icon("book-open"), tags$b(" Analysis: "),
                     "Points inside the white region (yolk) are predicted to be absorbed by the GI tract; points inside the yellow region (white) are predicted to cross the BBB."),
            tags$hr(),
            uiOutput("master_boiled_egg_logp_ui"),
            helpText("Select which LogP variant to use on the y-axis. WLOGP uses the official BOILED-Egg polygons (Daina & Zoete 2016). Other LogP variants use ALogP-trained polygons as approximation.")
          ),

          conditionalPanel(
            condition = "input.master_plot_type == 'Radar plot (Chemical profile)'",
            tags$div(class = "alert-primary", icon("book-open"), tags$b(" Function: "),
                     "Simultaneously compares multiple physicochemical properties of one or more compounds."),
            tags$div(class = "alert-primary", icon("book-open"), tags$b(" Analysis: "),
                     "Each axis corresponds to a descriptor. Similar profiles indicate similar properties; deviations highlight differences between molecules."),
            tags$hr(),
            uiOutput("master_radar_id_selector"),
            helpText("Select an identifier column (e.g., name or SMILES). Max 5 molecules in the radar plot")
          ),

          conditionalPanel(
            condition = "input.master_plot_type == 'Tanimoto / AGNES (Structural similarity)'",
            tags$div(class = "alert-primary", icon("book-open"), tags$b(" Function: "),
                     "Groups compounds according to their structural or physicochemical similarity."),
            tags$div(class = "alert-primary", icon("book-open"), tags$b(" Analysis: "),
                     "The dendrogram shows the relationship between molecules; branches close together indicate greater similarity."),
            tags$div(class = "alert-primary", icon("book-open"), tags$b(" Function: "),
                     "Quantifies the structural similarity between molecules using molecular fingerprints."),
            tags$div(class = "alert-danger", icon("binoculars"), tags$b(" Analysis: "),
                     "Values close to 1 indicate high structural similarity, while values close to 0 reflect different structures."),
            tags$hr(),
            uiOutput("master_smiles_col_selector"),
            uiOutput("master_label_col_selector"),
            numericInput("master_tanimoto_max_n", "Max molecules comparition", 40, min = 5, max = 200),
            selectInput("master_agnes_method", "Linking method (AGNES)",
                        choices = c("average", "complete", "single", "ward"), selected = "average"),
            helpText("Select the column that contains the SMILES for the similarity calculation and the column you want to use as dendrogram labels.")
          ),

          conditionalPanel(
            condition = "input.master_plot_type =='Parallel Coordinates'",
            tags$div(class = "alert-primary", icon("book-open"), tags$b(" Function: "),
                     "Groups compounds according to their structural or physicochemical similarity."),
            tags$div(class = "alert-primary", icon("book-open"), tags$b(" Analysis: "),
                     "The dendrogram shows the relationship between molecules; branches close together indicate greater similarity."),
            tags$hr(),
            uiOutput("master_parallel_controls")
          ),

          conditionalPanel(
            condition = "input.master_plot_type == 'Principal Component Analysys (PCA - Chemical space)'",
            tags$div(class = "alert-primary", icon("book-open"), tags$b(" Function: "),
                     "Reduces the dimensionality of data to identify global patterns."),
            tags$div(class = "alert-primary", icon("book-open"), tags$b(" Analysis: "),
                     "Compounds located close together exhibit similar physicochemical profiles, while those located far apart show different characteristics."),
            tags$hr(),
            uiOutput("master_pca_controls")
          ),

          conditionalPanel(
            condition = "input.master_plot_type == 'Violin Plot'",
            tags$div(class = "alert-primary", icon("book-open"), tags$b(" Function: "),
                     "Displays the distribution and density of a descriptor."),
            tags$div(class = "alert-primary", icon("book-open"), tags$b(" Analysis: "),
                     "Wider areas indicate where more observations are concentrated; allows for comparison of variability between groups."),
            tags$hr(),
            uiOutput("master_violin_controls")
          ),

          conditionalPanel(
            condition = "input.master_plot_type == 't-SNE (Chemical space)'",
            tags$div(class = "alert-primary", icon("book-open"), tags$b(" Function: "),
                     "Groups compounds according to the similarity of their properties while preserving local relationships."),
            tags$div(class = "alert-primary", icon("book-open"), tags$b(" Analysis: "),
                     "Tight clusters represent similar molecules; the distance between groups reflects differences in their profiles."),
            tags$hr(),
            uiOutput("master_tsne_controls")
          ),

          conditionalPanel(
            condition = "input.master_plot_type == 'UMAP (Chemical space)'",
            tags$div(class = "alert-primary", icon("book-open"), tags$b(" Function: "),
                     "Non-linear dimensionality reduction that preserves both local and global structure of the chemical space."),
            tags$div(class = "alert-primary", icon("book-open"), tags$b(" Analysis: "),
                     "Nearby points correspond to molecules with similar physicochemical profiles; well-separated clusters indicate distinct chemical series."),
            tags$hr(),
            uiOutput("master_umap_controls")
          ),

          conditionalPanel(
            condition = "input.master_plot_type == 'Cluster Heatmap (Dendrogram)'",
            tags$div(class = "alert-primary", icon("book-open"), tags$b(" Function: "),
                     "Combines hierarchical clustering with a heatmap to reveal groups of similar molecules."),
            tags$div(class = "alert-primary", icon("book-open"), tags$b(" Analysis: "),
                     "Branches close together indicate similar compounds; the color scale shows z-scored property values."),
            tags$hr(),
            uiOutput("master_cluster_heatmap_controls")
          ),

          conditionalPanel(
            condition = "input.master_plot_type == 'Custom Histogram'",
            tags$div(class = "alert-primary", icon("book-open"), tags$b(" Function: "),
                     "Highly customizable histogram of any numeric column in the dataset."),
            tags$div(class = "alert-primary", icon("book-open"), tags$b(" Analysis: "),
                     "Reveals the distribution shape, modality and outliers of any descriptor or ADMET probability."),
            tags$hr(),
            uiOutput("master_histogram_custom_controls")
          ),

          palette_selector_ui("master_palette"),
          br(),

          tags$div(style = "max-width: 900px; margin: 0 auto;",
                   plotOutput("master_plot", height = "600px")),
          tags$hr(),
          h4("Export Figure"),
          fluidRow(
            column(4, selectInput("master_format", "Format",
                                  choices = c("png", "pdf", "svg", "jpeg", "tiff"), selected = "png")),
            column(4, selectInput("master_dpi", "Resolution",
                                  choices = c("300 dpi" = 300, "600 dpi" = 600, "1200 dpi" = 1200), selected = 600))
          ),
          fluidRow(
            column(4, sliderInput("master_width", "Width (in)", min = 4, max = 12, value = 7, step = 0.5)),
            column(4, sliderInput("master_height", "Height (in)", min = 4, max = 12, value = 6, step = 0.5))
          ),
          br(),
          downloadButton("master_downloadPlot", "Download Figure", class = "btn-success")
        )
      )
    )
  )
}
