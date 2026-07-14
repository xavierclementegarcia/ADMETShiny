#' admetshiny: Interactive ADMET and Drug-Likeness Analysis
#'
#' The \pkg{admetshiny} package provides an interactive Shiny application and a
#' collection of R functions for the management, calculation, filtering,
#' visualization and exploratory analysis of molecular descriptors and ADMET
#' properties of small molecules. It integrates data exported from SwissADME,
#' ADMETlab 3.0 and Deep-PK, computes descriptors locally via the Chemistry
#' Development Kit (CDK) through \pkg{rcdk}, and offers drug-likeness filters
#' (Lipinski, Veber, Ghose, Egan, Muegge), the BOILED-Egg model, PCA, t-SNE,
#' radar plots and Tanimoto/AGNES structural similarity clustering.
#'
#' To launch the interactive application simply run:
#'
#' \preformatted{
#' admetshiny::run_app()
#' }
#'
#' The exported functions (filters, plots, data fixers and CDK helpers) can
#' also be used programmatically in plain R scripts.
#' By Xavier Clemente Garcia Cevallos :3 "Hecho con amor para la ciencia"
#'
#' @section Drug-likeness filters:
#'
#' \itemize{
#'   \item \code{\link{lipinskiFilter}} - Lipinski Rule-of-Five
#'   \item \code{\link{veberFilter}} - Veber oral bioavailability
#'   \item \code{\link{ghoseFilter}} - Ghose qualifying range
#'   \item \code{\link{eganFilter}} - Egan PSA-LogP rule
#'   \item \code{\link{mueggeFilter}} - Muegge pharmacophore filter
#'   \item \code{\link{applyFilters}} - Apply multiple filters in sequence
#'   \item \code{\link{computeViolationColumns}} - Compute violation columns
#' }
#'
#' @section Data normalization:
#'
#' \itemize{
#'   \item \code{\link{fixSwissADME}} - Normalize SwissADME data
#'   \item \code{\link{fixADMETlab}} - Normalize ADMETlab 3.0 data
#'   \item \code{\link{fixDeepPK}} - Normalize Deep-PK data
#'   \item \code{\link{mapCDKDescriptors}} - Map CDK descriptors to standard schema
#' }
#'
#' @section CDK & webchem:
#'
#' \itemize{
#'   \item \code{\link{getSmilesFromIdentifiers}} - Retrieve SMILES from PubChem
#'   \item \code{\link{calcCDKDescriptors}} - Calculate descriptors with CDK
#'   \item \code{\link{computeADMETProperties}} - Compute BOILED-Egg ADMET properties
#' }
#'
#' @section Visualizations:
#'
#' \itemize{
#'   \item \code{\link{plotBoiledEgg}} - BOILED-Egg model
#'   \item \code{\link{plotMW}}, \code{\link{plotTPSA}}, \code{\link{plotLogP}} - Distributions
#'   \item \code{\link{plotPCA}} - PCA chemical space
#'   \item \code{\link{plotTSNE}} - t-SNE chemical space
#'   \item \code{\link{plotParallel}} - Parallel coordinates
#'   \item \code{\link{plotViolin}} - Violin plot
#'   \item \code{\link{plotRadar}} - Radar chart
#'   \item \code{\link{plotTanimoto}} - Tanimoto/AGNES clustering
#'   \item \code{\link{plotCorrHeatmap}} - Correlation heatmap
#'   \item \code{\link{palette_selector_ui}}, \code{\link{apply_palette}} - Colour customization
#' }
#'
#' @section Application:
#'
#' \itemize{
#'   \item \code{\link{run_app}} - Launch the Shiny application
#'   \item \code{\link{generateReport}} - Generate a report from the R console
#' }
#'
#' @aliases admetshiny-package
#' @name admetshiny-package
#' @docType package
#' @import shiny
#' @import ggplot2
#' @importFrom dplyr filter mutate select
#' @importFrom magrittr %>%
#' @importFrom DT DTOutput renderDT
#' @importFrom stats prcomp cor complete.cases as.dist median sd setNames var dist hclust as.dendrogram heatmap
#' @importFrom utils read.csv write.csv
#' @importFrom grDevices rainbow colorRampPalette
#' @importFrom graphics legend par layout image axis mtext
#' @importFrom grid unit
"_PACKAGE"

# ---------------------------------------------------------------------------
# globalVariables
#
# Column names that are referenced as bare symbols inside dplyr::filter /
# mutate and ggplot2::aes calls. These are data columns, not true global
# variables, so we declare them here to silence R CMD check "no visible
# binding" NOTEs. This is the standard, CRAN-accepted approach.
# ---------------------------------------------------------------------------
utils::globalVariables(c(
  # ---- filters.R (dplyr::filter / mutate bare column names) ----
  "MW", "LogP", "TPSA", "MR",
  "#H-bond acceptors", "#H-bond donors", "#Rotatable bonds",
  "#Heavy atoms", "#Aromatic heavy atoms",
  "Lipinski #violations", "Ghose #violations", "Veber #violations",
  "Egan #violations", "Muegge #violations",
  ".hb_sum",
  # ---- plots.R (ggplot2::aes bare column names) ----
  "HIA", "BBB", "PGP", "x", "y", "Var1", "Var2", "Correlation",
  # ---- boiled_egg_data.R (internal data objects) ----
  ".hia_polygon", ".bbb_polygon",
  # ---- pca.R (ggplot2::aes bare column names) ----
  "PC1", "PC2", "Colour", "Label", "Variable",
  # ---- t_sne.R (ggplot2::aes bare column names) ----
  "TSNE1", "TSNE2", "Color",
  # ---- umap (ggplot2::aes bare column names) ----
  "UMAP1", "UMAP2",
  # ---- parallel.R (ggplot2::aes / ggparcoord) ----
  "Group",
  # ---- report.R (plotViolationsSummary / plotDruglikenessScore) ----
  "Rule", "Violations", "Freq", "Score", "Category", "BinMid"
))
