# ---------------------------------------------------------------------------
# ui_tutorial.R
# Tutorial / First Steps tab.
# ---------------------------------------------------------------------------

#' Tutorial tab UI
#'
#' Builds the Tutorial tab of the ADMETShiny application.
#'
#' @return A \code{shiny.tabPanel}.
#' @keywords internal
tutorial_tab <- function() {
  tabPanel(
    title = tagList(icon("graduation-cap"), "First Steps"),
    value = "tutorial",

    info_card_css(),

    fluidRow(
      column(12,
             tags$div(class = "info-card",
                      tags$div(class = "info-card-header", icon("graduation-cap"), " First Steps"),
                      tags$div(class = "info-card-body",
                               tags$h4("Getting Started with ADMETShiny"),
                               tags$hr(),

                               tags$h5(tags$b("Step 1: CDK & webchem Module")),
                               tags$p(style = "color:#555;",
                                      "Start by obtaining SMILES strings for your molecules. You can:"),
                               tags$ul(style = "color:#555; line-height: 1.8; padding-left: 20px;",
                                       tags$li(tags$b("Search by name:"), " Enter compound names (e.g., aspirin, ibuprofen) and ADMETShiny will query PubChem to retrieve canonical SMILES via the webchem package."),
                                       tags$li(tags$b("Enter SMILES manually:"), " Paste SMILES strings directly into the text area."),
                                       tags$li(tags$b("Upload a CSV:"), " Upload a CSV file with a SMILES column. You can also load the built-in example dataset (31 common drugs).")
                               ),
                               tags$p(style = "color:#555;",
                                      "Once SMILES are loaded, select the descriptors to calculate (MW, LogP, TPSA, HBD, HBA, RB, Heavy Atoms, Aromatic Atoms, MR) and click 'Calculate molecular descriptors'. CDK will process each molecule locally."),
                               tags$p(style = "color:#555;",
                                      "The calculated descriptors are automatically mapped to the standard schema, drug-likeness violation columns are computed, and BOILED-Egg ADMET properties (GI absorption, BBB permeant, Pgp substrate) are predicted."),

                               tags$hr(),
                               tags$h5(tags$b("Step 2: ADMET Master Manager")),
                               tags$p(style = "color:#555;",
                                      "This module accepts any ADMET dataset from any source. The workflow has 4 steps:"),
                               tags$ol(style = "color:#555; line-height: 1.8; padding-left: 20px;",
                                       tags$li(tags$b("Upload & Preview:"), " Upload your CSV or Excel file. The system auto-detects column types (numeric/string) and identifies the SMILES column."),
                                       tags$li(tags$b("Map Columns:"), " For each column in your dataset, select which standard field it corresponds to (e.g., map 'mol_weight' to 'MW', 'logp' to 'LogP'). Check 'Calculate missing descriptors with CDK' if you want the system to compute missing properties from SMILES. Click 'Calculate & Standardize'."),
                                       tags$li(tags$b("Filter:"), " Apply drug-likeness filters (Lipinski, Veber, Ghose, Egan, Muegge) with the same interactive interface. Download the filtered dataset."),
                                       tags$li(tags$b("Plots:"), " Visualize your data with 14 plot types including BOILED-Egg, PCA, t-SNE, UMAP, cluster heatmap, violin, and custom histograms.")
                               ),

                               tags$hr(),
                               tags$h5(tags$b("Step 3: Generate Report")),
                               tags$p(style = "color:#555;",
                                      "Navigate to the Report tab to generate a comprehensive report. The report includes:"),
                               tags$ul(style = "color:#555; line-height: 1.8; padding-left: 20px;",
                                       tags$li("Dataset statistics and drug-likeness scores"),
                                       tags$li("BOILED-Egg ADMET classification"),
                                       tags$li("Additional metrics (Pfizer 3/75, GSK 4/400, Golden Triangle)"),
                                       tags$li("Available in HTML preview, PDF, Word, and HTML download formats")
                               ),

                               tags$hr(),
                               tags$h5(tags$b("Tips")),
                               tags$ul(style = "color:#555; line-height: 1.8; padding-left: 20px;",
                                       tags$li(tags$b("BOILED-Egg LogP source:"), " When using the BOILED-Egg plot, select which LogP variant to use. WLOGP uses the official polygons; other LogP variants use the ALogP-trained polygons (approximation)."),
                                       tags$li(tags$b("CDK descriptors:"), " If your dataset lacks MR, Heavy Atoms, or Aromatic Atoms (needed for Ghose filter), enable CDK calculation in the ADMET Master module."),
                                       tags$li(tags$b("Plot export:"), " All plots can be exported as PNG, PDF, SVG, JPEG, or TIFF at 300/600/1200 DPI with customizable dimensions."),
                                       tags$li(tags$b("Custom Histogram:"), " Use the Custom Histogram plot to visualize any numeric column in your dataset, including ADMET probabilities and user-supplied data.")
                               )
                      )
             )
      )
    ),

    br()
  )
}
