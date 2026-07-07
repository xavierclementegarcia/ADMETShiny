# ---------------------------------------------------------------------------
# ui_cdk.R
# CDK & webchem tab: SMILES retrieval, descriptor calculation, filtering.
# ---------------------------------------------------------------------------

#' CDK & webchem tab UI
#'
#' Builds the CDK & webchem tab of the ADMETShiny application.
#'
#' @return A \code{shiny.tabPanel}.
#' @keywords internal
cdk_tab <- function() {
  tabPanel(
    title = tagList(icon("vial"), "CDK & webchem"),
    value = "cdk",

    fluidPage(

      h2("Obtaining SMILES & molecular descriptors "),

      tags$div(
        class = "alert alert-warning",
        icon("triangle-exclamation"),
        tags$b(" Warning: "),
        "This section requires an internet connection (PubChem access) and the rcdk package with rJava and a Java Development Kit (JDK) installed (a JRE alone is insufficient). The LogP (ALogP) and MR (AMR) values calculated by CDK are approximations and may differ slightly from values computed by other platforms."
      ),

      tags$hr(),

      h3("Step 1: Obtain SMILES"),

      p("Choose one of three methods to provide SMILES for your molecules: ",
        tags$b("PubChem identifiers"), ", ",
        tags$b("direct SMILES entry"), ", or ",
        tags$b("CSV upload"), "."),

      tabsetPanel(
        id = "cdk_input_method",

        ## ---- Method 1: PubChem identifiers ----
        tabPanel(
          title = tagList(icon("magnifying-glass"), " PubChem"),
          value = "pubchem",
          br(),
          fluidRow(
            column(6,
                   textAreaInput("cdk_ids", "Place the identifier on each line",
                                 rows = 6, placeholder = "aspirin\nibuprofen\n50-78-2")
            ),
            column(6,
                   selectInput("cdk_id_type", "Identifier type",
                               choices = c("Name" = "name", "CAS" = "cas",
                                           "InChIKey" = "inchikey", "PubChem CID" = "cid")),
                   br(),
                   actionButton("fetch_smiles", "Obtain SMILES",
                                icon = icon("magnifying-glass"), class = "btn-primary")
            )
          ),
          tags$div(
            class = "alert alert-info",
            icon("circle-info"),
            tags$b(" Tip: "),
            "Queries PubChem via the webchem package. Requires an internet connection."
          )
        ),

        ## ---- Method 2: Direct SMILES entry ----
        tabPanel(
          title = tagList(icon("keyboard"), " Manual SMILES"),
          value = "manual",
          br(),
          textAreaInput("cdk_manual_smiles",
                        "Paste SMILES (one per line)",
                        rows = 8,
                        placeholder = "CCO\nCC(=O)Oc1ccccc1C(=O)O\nc1ccccc1\nCC(C)CC(N)C(=O)O"),
          fluidRow(
            column(6,
                   textInput("cdk_manual_names", "Optional: molecule names (comma-separated)",
                             placeholder = "ethanol,aspirin,benzene,valine")
            ),
            column(6, br(),
                   actionButton("use_manual_smiles", "Use these SMILES",
                                icon = icon("check"), class = "btn-primary")
            )
          ),
          tags$div(
            class = "alert alert-info",
            icon("circle-info"),
            tags$b(" Tip: "),
            "Enter one canonical SMILES per line. Optionally, provide names separated by commas ",
            "in the same order as the SMILES. No internet connection required."
          )
        ),

        ## ---- Method 3: CSV upload ----
        tabPanel(
          title = tagList(icon("upload"), " Upload CSV"),
          value = "csv",
          br(),
          fileInput("cdk_smiles_csv", "Upload a CSV file with SMILES",
                    accept = c(".csv", "text/csv", "text/comma-separated-values")),
          fluidRow(
            column(6,
                   selectInput("cdk_csv_smiles_col", "Column containing SMILES",
                               choices = NULL),
                   helpText("Select the column that contains the SMILES strings.")
            ),
            column(6,
                   checkboxInput("cdk_csv_has_header", "File has a header row", TRUE),
                   br(),
                   actionButton("use_csv_smiles", "Load SMILES from CSV",
                                icon = icon("file-import"), class = "btn-primary")
            )
          ),
          tags$div(
            class = "alert alert-info",
            icon("circle-info"),
            tags$b(" Tip: "),
            "The CSV can have any number of columns; you will be asked to select which ",
            "one contains the SMILES. All other columns are preserved in the dataset."
          )
        )
      ),

      br(),
      DTOutput("cdk_smiles_table"),

      tags$hr(),
      fluidRow(
        column(4, radioButtons("smiles_export_format", "Export format",
                               choices = c("CSV" = "csv", "Excel (.xlsx)" = "xlsx"), inline = TRUE)),
        column(3, br(), downloadButton("download_smiles", "Dataset download"))
      ),

      tags$hr(),

      h3("Molecular viewer"),
      p("Select a molecule from the list obtained in Step 1 to view its 2D structure (image from PubChem)."),

      fluidRow(
        column(width = 4,
               wellPanel(style = "min-height: 400px;",
                         selectInput("cdk_mol_selector", "Select molecule:", choices = NULL, width = "100%"),
                         tags$hr(),
                         uiOutput("cdk_mol_info")
               )
        ),
        column(width = 8,
               wellPanel(style = "text-align: center; min-height: 400px; display: flex; align-items: center; justify-content: center; background: #fafafa;",
                         uiOutput("cdk_molecule_image")
               )
        )
      ),

      tags$hr(),

      h3("Step 2: Molecular descriptors"),
      fluidRow(
        column(6,
               checkboxGroupInput("cdk_descriptors", "Molecular descriptors",
                                  choices = c(
                                    "Molecular weight (MW)"                 = "mw",
                                    "LogP (ALogP)"                          = "alogp",
                                    "TPSA"                                  = "tpsa",
                                    "H bond donors (HBD)"                   = "hbd",
                                    "H bond acceptors (HBA)"                = "hba",
                                    "Rotatable bonds"                       = "rotb",
                                    "Heavy atoms"                           = "heavy",
                                    "Aromatic heavy atoms"                  = "aroma",
                                    "Molar refractivity (aprox., AMR)"      = "mr"
                                  ),
                                  selected = c("mw", "alogp", "tpsa", "hbd", "hba", "rotb", "heavy", "aroma", "mr"))
        ),
        column(6, br(),
               actionButton("calc_cdk", "Calculate molecular descriptors",
                            icon = icon("calculator"), class = "btn-primary")
        )
      ),

      br(),
      DTOutput("cdk_results_table"),

      tags$hr(),

      h3("Step 3: Visualization and Drug-likeness filter"),
      p("It applies the same drug-likeness filters (Lipinski, Veber, Ghose, Egan, Muegge) to the descriptors calculated with CDK. The columns are automatically mapped to the SwissADME schema, and the #violations columns required by the filters are calculated."),

      sidebarLayout(

        sidebarPanel(

          checkboxGroupInput("cdk_filters", "Drug-likeness Filters",
                             choices = c("Lipinski", "Veber", "Ghose", "Egan", "Muegge"),
                             selected = "Lipinski"),

          conditionalPanel(condition = "input.cdk_filters.includes('Lipinski')",
                           tags$b("Lipinski"),
                           numericInput("cdk_mw", "Max. Molecular Weight", 500),
                           numericInput("cdk_logp", "Max. LogP", 5),
                           numericInput("cdk_hba", "Max. H-bond acceptors", 10),
                           numericInput("cdk_hbd", "Max. H-bond donors", 5),
                           numericInput("cdk_violations", "Maximum Lipinski Violations", 0)),

          conditionalPanel(condition = "input.cdk_filters.includes('Veber')",
                           tags$b("Veber"),
                           numericInput("cdk_v_rb", "Max. Rotatable Bonds", 10),
                           numericInput("cdk_v_tpsa", "Max. TPSA", 140),
                           numericInput("cdk_v_hb_sum", "Max. HBA + HBD", 12),
                           numericInput("cdk_veber_violations", "Maximum Veber Violations", 0)),

          conditionalPanel(condition = "input.cdk_filters.includes('Ghose')",
                           tags$b("Ghose"),
                           sliderInput("cdk_g_mw", "Molecular Weight range", min = 100, max = 700, value = c(160, 480)),
                           sliderInput("cdk_g_mr", "Molar Refractivity range", min = 0, max = 200, value = c(40, 130)),
                           sliderInput("cdk_g_logp", "LogP range", min = -5, max = 10, value = c(-0.4, 5.6)),
                           sliderInput("cdk_g_ha", "Heavy Atoms range", min = 0, max = 100, value = c(20, 70)),
                           numericInput("cdk_ghose_violations", "Maximum Ghose Violations", 0)),

          conditionalPanel(condition = "input.cdk_filters.includes('Egan')",
                           tags$b("Egan"),
                           numericInput("cdk_e_tpsa", "Max. TPSA", 131.6),
                           numericInput("cdk_e_logp", "Max. LogP", 5.88),
                           numericInput("cdk_egan_violations", "Maximum Egan Violations", 0)),

          conditionalPanel(condition = "input.cdk_filters.includes('Muegge')",
                           tags$b("Muegge"),
                           sliderInput("cdk_m_mw", "Molecular Weight range", min = 50, max = 800, value = c(200, 600)),
                           sliderInput("cdk_m_logp", "LogP range", min = -10, max = 10, value = c(-2, 5)),
                           numericInput("cdk_m_hba", "Max. H-bond acceptors", 10),
                           numericInput("cdk_m_hbd", "Max. H-bond donors", 5),
                           numericInput("cdk_m_rb", "Max. Rotatable Bonds", 15),
                           numericInput("cdk_m_tpsa", "Max. TPSA", 150),
                           numericInput("cdk_muegge_violations", "Maximum Muegge Violations", 0)),

          tags$hr(),
          actionButton("cdk_run", "Apply Filter", class = "btn-primary"),
          br(), br(),
          downloadButton("cdk_download", "Download filtered dataset")
        ),

        mainPanel(

          tabsetPanel(

            tabPanel("Filtered Molecules", br(), h4("CDK Filtered Results"), DTOutput("cdk_tabla")),

            tabPanel(
              "Plots", br(),

              selectInput("cdk_plot_type", "Select plot", choices = c(
                "Boiled Egg", "Molecular Weight", "TPSA", "LogP",
                "Radar plot (Chemical profile)",
                "Tanimoto / AGNES (Structural similarity)",
                "Correlation Heatmap",
                "Principal Component Analysys (PCA - Chemical space)",
                "t-SNE (Chemical space)",
                "Parallel Coordinates",
                "Violin Plot"
              )),

              conditionalPanel(
                condition = "input.cdk_plot_type == 'Radar plot (Chemical profile)'",
                tags$div(class = "alert-primary", icon("book-open"), tags$b(" Function: "),
                         "Simultaneously compares multiple physicochemical properties of one or more compounds."),
                tags$div(class = "alert-primary", icon("book-open"), tags$b(" Analysis: "),
                         "Each axis corresponds to a descriptor. Similar profiles indicate similar properties; deviations highlight differences between molecules."),
                tags$hr(),
                uiOutput("cdk_radar_id_selector"),
                helpText("Select an identifier column (e.g. name or SMILES) and up to 5 molecules to overlay on the radar.")
              ),

              conditionalPanel(
                condition = "input.cdk_plot_type == 'Tanimoto / AGNES (Structural similarity)'",
                tags$div(class = "alert-primary", icon("book-open"), tags$b(" Function: "),
                         "Groups compounds according to their structural or physicochemical similarity."),
                tags$div(class = "alert-primary", icon("book-open"), tags$b(" Analysis: "),
                         "The dendrogram shows the relationship between molecules; branches close together indicate greater similarity."),
                tags$div(class = "alert-primary", icon("book-open"), tags$b(" Function: "),
                         "Quantifies the structural similarity between molecules using molecular fingerprints."),
                tags$div(class = "alert-danger", icon("binoculars"), tags$b(" Analysis: "),
                         "Values close to 1 indicate high structural similarity, while values close to 0 reflect different structures."),
                tags$hr(),
                uiOutput("cdk_smiles_col_selector"),
                uiOutput("cdk_label_col_selector"),
                numericInput("cdk_tanimoto_max_n", "Max molecules comparition", 40, min = 5, max = 200),
                selectInput("cdk_agnes_method", "Linking method (AGNES)",
                            choices = c("average", "complete", "single", "ward"), selected = "average"),
                helpText("Select the column that contains the SMILES and the column for the dendrogram labels.")
              ),

              conditionalPanel(
                condition = "input.cdk_plot_type =='Parallel Coordinates'",
                tags$div(class = "alert-primary", icon("book-open"), tags$b(" Function: "),
                         "Groups compounds according to their structural or physicochemical similarity."),
                tags$div(class = "alert-primary", icon("book-open"), tags$b(" Analysis: "),
                         "The dendrogram shows the relationship between molecules; branches close together indicate greater similarity."),
                tags$hr(),
                uiOutput("cdk_parallel_controls")
              ),

              conditionalPanel(
                condition = "input.cdk_plot_type == 'Principal Component Analysys (PCA - Chemical space)'",
                tags$div(class = "alert-primary", icon("book-open"), tags$b(" Function: "),
                         "Reduces the dimensionality of data to identify global patterns."),
                tags$div(class = "alert-primary", icon("book-open"), tags$b(" Analysis: "),
                         "Compounds located close together exhibit similar physicochemical profiles, while those located far apart show different characteristics."),
                tags$hr(),
                uiOutput("cdk_pca_controls")
              ),

              conditionalPanel(
                condition = "input.cdk_plot_type == 'Violin Plot'",
                tags$div(class = "alert-primary", icon("book-open"), tags$b(" Function: "),
                         "Displays the distribution and density of a descriptor."),
                tags$div(class = "alert-primary", icon("book-open"), tags$b(" Analysis: "),
                         "Wider areas indicate where more observations are concentrated; allows for comparison of variability between groups."),
                tags$hr(),
                uiOutput("cdk_violin_controls")
              ),

              conditionalPanel(
                condition = "input.cdk_plot_type == 't-SNE (Chemical space)'",
                tags$div(class = "alert-primary", icon("book-open"), tags$b(" Function: "),
                         "Groups compounds according to the similarity of their properties while preserving local relationships."),
                tags$div(class = "alert-primary", icon("book-open"), tags$b(" Analysis: "),
                         "Tight clusters represent similar molecules; the distance between groups reflects differences in their profiles."),
                tags$hr(),
                uiOutput("cdk_tsne_controls")
              ),

              palette_selector_ui("cdk_palette"),
              br(),
              plotOutput("cdk_plot", height = "650px"),
              tags$hr(),
              h4("Export Figure"),
              fluidRow(
                column(4, selectInput("cdk_format", "Format",
                                      choices = c("png", "pdf", "svg", "jpeg", "tiff"), selected = "png")),
                column(4, selectInput("cdk_dpi", "Resolution",
                                      choices = c("300 dpi" = 300, "600 dpi" = 600, "1200 dpi" = 1200), selected = 600))
              ),
              fluidRow(
                column(4, sliderInput("cdk_width", "Width (in)", min = 4, max = 12, value = 7, step = 0.5)),
                column(4, sliderInput("cdk_height", "Height (in)", min = 4, max = 12, value = 6, step = 0.5))
              ),
              br(),
              downloadButton("cdk_downloadPlot", "Download Figure", class = "btn-success")
            )
          )
        )
      ),

      h3("Step 4: Use this dataset in drug-likeness filtering (BETA)"),
      p("The columns are automatically mapped to the same schema that produces",
        tags$code("fixSwissADME()"),
        "(MW, TPSA, #H-bond donors/acceptors, etc.) and calculate the #violations columns that the filters need. Continuing, you can apply Lipinski, Veber, Ghose, Egan, and Muegge and view the same existing tables and graphs, including BOILED-Egg, Radar plot, and Tanimoto/AGNES (the latter using the newly obtained SMILES column)."
      ),
      actionButton("send_to_filter", "Send to SwissADME Filtering",
                   icon = icon("arrow-right"), class = "btn-success"),
      br(), br(),
      actionButton("back_home_cdk", "Back to Home", icon = icon("arrow-left"))
    )
  )
}
