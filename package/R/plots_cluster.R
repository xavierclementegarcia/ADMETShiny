# ---------------------------------------------------------------------------
# plots_cluster.R
# Cluster heatmap with dendrogram for chemical space analysis.
#
# Combines hierarchical clustering (on rows/compounds) with a heatmap of
# physicochemical properties, producing a publication-quality visualization
# that reveals clusters of similar molecules.
# ---------------------------------------------------------------------------

#' Cluster heatmap with dendrogram
#'
#' Produces a heatmap of physicochemical properties with a dendrogram
#' coupled to the left side (compounds) showing hierarchical clustering.
#' The properties (columns) are also clustered and the column dendrogram
#' appears on top. Values are z-score standardized per property.
#'
#' The plot is rendered using base R graphics so it works without additional
#' dependencies beyond what the package already requires.
#'
#' @param data A data.frame with physicochemical property columns.
#' @param variables Character vector of numeric column names to use. If
#'   \code{NULL}, a default set of common descriptors is used.
#' @param id_col Character. Name of the column to use as row labels. If
#'   \code{NULL} or \code{"None"}, row numbers are used.
#' @param method Character. Agglomeration method for hierarchical clustering:
#'   one of \code{"ward.D"}, \code{"ward.D2"}, \code{"single"}, \code{"complete"},
#'   \code{"average"} (= UPGMA), \code{"mcquitty"}, \code{"median"}, \code{"centroid"}.
#'   Default \code{"ward.D2"}.
#' @param scale_data Logical. Whether to z-score standardize each property
#'   before clustering. Default \code{TRUE}.
#' @param palette Character. Colour palette name: one of \code{"default"},
#'   \code{"viridis"}, \code{"magma"}, \code{"inferno"}, \code{"RdBu"},
#'   \code{"Set1"}, \code{"Set2"}. Default \code{"default"} (red-blue).
#' @return Invisible \code{NULL}; called for the side-effect of drawing the
#'   plot on the current graphics device.
#' @export
#' @references Murtagh, F., & Contreras, P. (2012). Algorithms for
#'   hierarchical clustering: an overview. \emph{Wiley Interdisciplinary
#'   Reviews: Data Mining and Knowledge Discovery}, 2(1), 86-97.
#' @examples
#' \dontrun{
#' d <- data.frame(
#'   Name = c("A", "B", "C", "D"),
#'   MW = c(300, 350, 600, 320),
#'   LogP = c(2, 3, 5, 1),
#'   TPSA = c(50, 60, 120, 40),
#'   check.names = FALSE
#' )
#' plotClusterHeatmap(d, id_col = "Name")
#' }
plotClusterHeatmap <- function(data, variables = NULL, id_col = NULL,
                               method = "ward.D2", scale_data = TRUE,
                               palette = "default") {

  default_vars <- c(
    "MW", "LogP", "TPSA", "MR",
    "#H-bond acceptors", "#H-bond donors",
    "#Rotatable bonds", "#Heavy atoms"
  )

  if (is.null(variables) || length(variables) == 0) {
    variables <- intersect(default_vars, names(data))
  }
  variables <- intersect(variables, names(data))
  if (length(variables) < 2) {
    stop("At least two numeric variables are required.", call. = FALSE)
  }

  ## Coerce to numeric explicitly
  mat_data <- as.data.frame(
    lapply(data[, variables, drop = FALSE],
           function(x) suppressWarnings(as.numeric(as.character(x)))),
    stringsAsFactors = FALSE
  )
  names(mat_data) <- variables

  ## Remove rows with NA
  cc <- complete.cases(mat_data)
  mat_data <- mat_data[cc, , drop = FALSE]
  if (nrow(mat_data) < 3) {
    stop("At least 3 complete observations are required.", call. = FALSE)
  }

  ## Remove zero-variance columns
  col_vars <- sapply(mat_data, var, na.rm = TRUE)
  zero_var <- names(mat_data)[is.na(col_vars) | col_vars == 0]
  if (length(zero_var) > 0) {
    mat_data <- mat_data[, !(names(mat_data) %in% zero_var), drop = FALSE]
    variables <- setdiff(variables, zero_var)
  }
  if (ncol(mat_data) < 2) {
    stop("After removing zero-variance columns, fewer than 2 variables remain.",
         call. = FALSE)
  }

  ## Row labels -- safely handle NULL, "None", character(0), NA, and ""
  ## The input can be:
  ##   - NULL (not selected)
  ##   - character(0) (empty vector from selectInput reset)
  ##   - NA (missing value)
  ##   - "None" (the default choice)
  ##   - "" (empty string)
  ##   - a valid column name
  use_id_col <- !is.null(id_col) &&
                is.character(id_col) &&
                length(id_col) == 1 &&
                !is.na(id_col[1]) &&
                id_col[1] != "None" &&
                id_col[1] != "" &&
                id_col[1] %in% names(data)

  if (use_id_col) {
    row_labels <- as.character(data[[id_col[1]]][cc])
    row_labels[is.na(row_labels)] <- "NA"
    row_labels <- make.unique(row_labels)
  } else {
    row_labels <- paste0("C", seq_len(nrow(mat_data)))
  }

  ## Convert to matrix
  mat <- as.matrix(mat_data)
  ## CRITICAL: rownames must be set BEFORE calling heatmap() because
  ## heatmap() uses them for labRow when labRow is not explicitly passed.
  ## We set them here and also pass labRow explicitly for safety.
  rownames(mat) <- row_labels

  ## Scale per column (z-score)
  if (scale_data) {
    mat <- scale(mat)
    mat[is.nan(mat) | is.infinite(mat)] <- 0
    ## scale() preserves rownames, but remove attributes that can
    ## interfere with heatmap()
    attr(mat, "scaled:center") <- NULL
    attr(mat, "scaled:scale") <- NULL
  }

  ## Ensure rownames are still set (safety net)
  rownames(mat) <- row_labels

  ## Get colour palette
  colors <- .get_heatmap_palette(palette)

  ## Use the built-in heatmap() function which couples the dendrogram
  ## directly to the heatmap. This produces a single integrated figure
  ## with the row dendrogram on the left and the column dendrogram on top,
  ## both aligned with the heatmap cells.
  ##
  ## heatmap() handles:
  ##   - Row and column clustering
  ##   - Reordering of the matrix
  ##   - Drawing the dendrograms coupled to the heatmap
  ##   - The color scale
  ##
  ## We compute the dendrograms ourselves to have control over the method.
  row_dist <- dist(mat, method = "euclidean")
  row_hclust <- hclust(row_dist, method = method)
  col_dist <- dist(t(mat), method = "euclidean")
  col_hclust <- hclust(col_dist, method = method)

  ## Truncate row labels if too many (to avoid overlap)
  n_rows <- nrow(mat)
  if (n_rows > 50) {
    lab_row <- rep("", n_rows)
  } else if (n_rows > 30) {
    lab_row <- substr(rownames(mat), 1, 10)
  } else {
    lab_row <- rownames(mat)
  }

  ## Calculate appropriate cexRow based on number of rows
  cex_row <- if (n_rows > 40) 0.5 else if (n_rows > 25) 0.6 else
             if (n_rows > 15) 0.7 else 0.8

  ## Calculate left margin based on longest label
  max_label_len <- max(nchar(lab_row))
  left_margin <- max(8, min(max_label_len, 20))

  ## Draw the coupled heatmap + dendrogram
  heatmap(mat,
          Rowv = as.dendrogram(row_hclust),
          Colv = as.dendrogram(col_hclust),
          scale = "none",
          col = colors,
          margins = c(8, left_margin),
          labRow = lab_row,
          labCol = colnames(mat),
          main = "Cluster Heatmap (Dendrogram)",
          xlab = "", ylab = "",
          cexCol = 0.9,
          cexRow = cex_row)

  invisible(NULL)
}

#' Get heatmap colour palette
#'
#' Internal helper that returns a colour vector for the cluster heatmap
#' based on the selected palette name.
#'
#' @param palette_name Character. Palette name.
#' @return A character vector of colour hex codes.
#' @keywords internal
.get_heatmap_palette <- function(palette_name) {

  n_colors <- 100

  if (palette_name == "viridis") {
    if (requireNamespace("viridisLite", quietly = TRUE)) {
      return(viridisLite::viridis(n_colors))
    }
    return(colorRampPalette(c("#440154", "#3b528b", "#21918c", "#5ec962",
                               "#fde725"))(n_colors))
  }

  if (palette_name == "magma") {
    if (requireNamespace("viridisLite", quietly = TRUE)) {
      return(viridisLite::magma(n_colors))
    }
    return(colorRampPalette(c("#000004", "#3b0f70", "#8c2981", "#de4968",
                               "#fe9f6d", "#fcfdbf"))(n_colors))
  }

  if (palette_name == "inferno") {
    if (requireNamespace("viridisLite", quietly = TRUE)) {
      return(viridisLite::inferno(n_colors))
    }
    return(colorRampPalette(c("#000004", "#420a68", "#932667", "#dd513a",
                               "#fca50a", "#fcffa4"))(n_colors))
  }

  if (palette_name == "RdBu") {
    return(colorRampPalette(c("#2166AC", "#67A9CF", "white",
                               "#F4A582", "#B2182B"))(n_colors))
  }

  if (palette_name == "Set1") {
    return(colorRampPalette(c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3",
                               "#ff7f00", "#ffff33", "#a65628"))(n_colors))
  }

  if (palette_name == "Set2") {
    return(colorRampPalette(c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3",
                               "#a6d854", "#ffd92f", "#e5c494"))(n_colors))
  }

  ## Default: blue-white-red (good for z-scores)
  colorRampPalette(c("#2166AC", "#67A9CF", "white",
                      "#F4A582", "#B2182B"))(n_colors)
}
