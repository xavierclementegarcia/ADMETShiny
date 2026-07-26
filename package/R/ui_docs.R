# ---------------------------------------------------------------------------
# ui_docs.R
# Documentation tab.
# ---------------------------------------------------------------------------

#' Documentation tab UI
#'
#' Builds the Documentation tab of the ADMETShiny application.
#'
#' @return A \code{shiny.tabPanel}.
#' @keywords internal
docs_tab <- function() {
  tabPanel(
    title = tagList(icon("book"), "Documentation"),
    value = "docs",

    info_card_css(),

    fluidRow(
      column(12,
             tags$div(class = "info-card",
                      tags$div(class = "info-card-header", icon("book"), " Documentation"),
                      tags$div(class = "info-card-body",
                               tags$h4("ADMETShiny Documentation"),
                               tags$hr(),
                               tags$h5(tags$b("Modules")),
                               tags$p(style = "color:#555;",
                                      "ADMETShiny consists of two main modules:"),
                               tags$ul(style = "color:#555; line-height: 1.8; padding-left: 20px;",
                                       tags$li(tags$b("CDK & webchem:"), " Retrieve SMILES from PubChem via compound names, CAS numbers, or InChIKeys, or enter SMILES manually. Calculate 9 physicochemical descriptors locally using the Chemistry Development Kit (CDK) through rcdk. Apply drug-likeness filters and visualize results."),
                                       tags$li(tags$b("ADMET Master Manager:"), " Upload any ADMET dataset (CSV or Excel) from any source. Manually map your columns to the application's standard schema. Optionally calculate missing descriptors with CDK. Apply filters, visualize, and generate reports.")
                               ),

                               tags$hr(),
                               tags$h5(tags$b("Drug-Likeness Filters")),
                               tags$p(style = "color:#555;",
                                      "The package implements five drug-likeness filters with thresholds from the original publications:"),
                               tags$ul(style = "color:#555; line-height: 1.8; padding-left: 20px;",
                                       tags$li(tags$b("Lipinski (1997):"), " MW <= 500, LogP <= 5, HBA <= 10, HBD <= 5"),
                                       tags$li(tags$b("Veber (2002):"), " RB <= 10 AND (TPSA <= 140 OR HBA+HBD <= 12)"),
                                       tags$li(tags$b("Ghose (1999):"), " MW 160-480, MR 40-130, LogP -0.4 to 5.6, heavy atoms 20-70"),
                                       tags$li(tags$b("Egan (2000):"), " TPSA <= 131.6, LogP <= 5.88"),
                                       tags$li(tags$b("Muegge (2001):"), " MW 200-600, LogP -2 to 5, HBA <= 10, HBD <= 5, RB <= 15, TPSA <= 150, pharmacophore points (HBA+HBD) >= 4")
                               ),

                               tags$hr(),
                               tags$h5(tags$b("BOILED-Egg Model")),
                               tags$p(style = "color:#555;",
                                      "The BOILED-Egg model (Daina & Zoete, 2016) predicts gastrointestinal absorption (HIA) and blood-brain barrier permeability (BBB) from TPSA and LogP."),
                               tags$p(style = "color:#555;",
                                      "ADMETShiny uses dual polygon systems:"),
                               tags$ul(style = "color:#555; line-height: 1.8; padding-left: 20px;",
                                       tags$li(tags$b("Official WLOGP polygons:"), " Used when the data contains a WLOGP column (e.g., SwissADME exports). These are the original coordinates from Daina & Zoete (2016) Data S3."),
                                       tags$li(tags$b("ALogP-trained polygons:"), " Used when WLOGP is not available (CDK, ADMET Master). Trained from 439 molecules using Monte Carlo optimization. Independently validated against 632 HIA and 240 BBB experimental compounds from Daina's dataset (HIA: 90.5% accuracy, MCC=0.525; BBB: 82.9% accuracy, MCC=0.643).")
                               ),

                               tags$hr(),
                               tags$h5(tags$b("P-glycoprotein Substrate Prediction")),
                               tags$p(style = "color:#555;",
                                      "P-glycoprotein (P-gp/ABCB1) substrate prediction uses a Random Forest model (100 trees, max depth 10) trained from 882 experimental compounds from Metrabase (J. Cheminformatics 2015, 7:21)."),
                               tags$p(style = "color:#555;",
                                      "The model uses 9 CDK descriptors (MW, LogP, TPSA, HBD, HBA, RB, HeavyAtoms, AromAtoms, MR) and achieves 69.4% cross-validated accuracy (MCC=0.383), significantly improving over the previous heuristic (55% accuracy, MCC=0.145)."),
                               tags$p(style = "color:#555;",
                                      "The model is stored as pure R code (100 decision trees as parallel arrays) and requires no Python or external ML packages at runtime."),

                               tags$hr(),
                               tags$h5(tags$b("Column Mapping (ADMET Master)")),
                               tags$p(style = "color:#555;",
                                      "The ADMET Master Manager allows users to upload any ADMET dataset and manually map columns to the standard schema. The system supports 20 standard fields including SMILES, MW, LogP, WLOGP, TPSA, HBD, HBA, RB, MR, Heavy Atoms, Aromatic Heavy Atoms, GI Absorption (categorical or numeric 0-1), BBB Permeant, Pgp Substrate, LogS, and LogD."),
                               tags$p(style = "color:#555;",
                                      "When a SMILES column is available, missing descriptors can be calculated on-the-fly using CDK. ADMET numeric probabilities (0-1) are automatically converted to categorical labels (High/Low, Yes/No) using a 0.5 threshold."),

                               tags$hr(),
                               tags$h5(tags$b("Visualization")),
                               tags$p(style = "color:#555;",
                                      "Available plot types (14 total):"),
                               tags$ul(style = "color:#555; line-height: 1.8; padding-left: 20px;",
                                       tags$li("BOILED-Egg, MW/TPSA/LogP histograms"),
                                       tags$li("Radar plot, Tanimoto/AGNES clustering"),
                                       tags$li("Correlation heatmap, Cluster heatmap (dendrogram)"),
                                       tags$li("PCA, t-SNE, UMAP (chemical space)"),
                                       tags$li("Parallel coordinates, Violin plot"),
                                       tags$li("Custom histogram (any numeric column)")
                               ),

                               tags$hr(),
                               tags$h5(tags$b("Report Generation")),
                               tags$p(style = "color:#555;",
                                      "The Report module generates comprehensive reports (HTML preview, PDF, Word, HTML download) including statistics, drug-likeness scores, BOILED-Egg classification, and additional metrics (Pfizer 3/75, GSK 4/400, Golden Triangle).")
                      )
             )
      )
    ),

    br()
  )
}
