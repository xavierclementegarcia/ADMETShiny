# ---------------------------------------------------------------------------
# plots.R
# Plotting functions: BOILED-Egg, descriptor histograms, correlation heatmap,
# radar chart, Tanimoto/AGNES clustering, PCA, t-SNE, UMAP, parallel
# coordinates and violin plot.
# ---------------------------------------------------------------------------

## --------------------------- BOILED-Egg (official polygons) ---------------
## Uses the official polygon coordinates from the BOILED-Egg supplementary
## data (Daina & Zoete, 2016, Data S3). The polygon coordinates are defined
## in R/boiled_egg_data.R as internal objects (.hia_polygon, .bbb_polygon).
##
## Axis convention (matching SwissADME):
##   X-axis = TPSA (0 to 200)
##   Y-axis = LogP (-2 to 7)
##
## The generic "LogP" column is used as a proxy for WLOGP.
## Point-in-polygon classification uses sp::point.in.polygon() (preferred)
## or a ray-casting fallback.

## Internal helper: point-in-polygon test
## Uses sp::point.in.polygon if available, otherwise a ray-casting fallback
..point_in_polygon <- function(x, y, poly) {
  if (requireNamespace("sp", quietly = TRUE)) {
    result <- sp::point.in.polygon(x, y, poly[, 1], poly[, 2])
    return(result > 0)
  }
  ## Fallback: ray-casting algorithm
  n <- nrow(poly)
  inside <- rep(FALSE, length(x))
  for (k in seq_along(x)) {
    xi <- x[k]; yi <- y[k]
    is_inside <- FALSE
    j <- n
    for (i in seq_len(n)) {
      if (((poly[i, 2] > yi) != (poly[j, 2] > yi)) &&
          (xi < (poly[j, 1] - poly[i, 1]) * (yi - poly[i, 2]) /
           (poly[j, 2] - poly[i, 2]) + poly[i, 1])) {
        is_inside <- !is_inside
      }
      j <- i
    }
    inside[k] <- is_inside
  }
  inside
}

#' BOILED-Egg plot
#'
#' Produces the BOILED-Egg model plot (Daina & Zoete, 2016) using the
#' official polygon coordinates from the supplementary data (Data S3).
#' Points are coloured by P-gp substrate status when a
#' \code{"Pgp substrate"} column is available.
#'
#' The BOILED-Egg model was originally calibrated with WLOGP; the application
#' uses the generic \code{LogP} column as an approximation. The axes follow
#' the SwissADME convention: TPSA on the x-axis, LogP on the y-axis.
#'
#' @param data A data.frame that must contain the columns \code{LogP} and
#'   \code{TPSA}. A \code{"Pgp substrate"} column (with values "Yes"/"No") is
#'   optional and used for point colouring.
#' @return A \pkg{ggplot2} object.
#' @export
#' @references Daina, A., & Zoete, V. (2016). A boiled egg to predict
#'   gastrointestinal absorption and brain penetration of small molecules.
#'   \emph{ChemMedChem}, 11(11), 1117-1121.
#' @examples
#' \dontrun{
#' d <- data.frame(LogP = c(2, 5), TPSA = c(50, 150))
#' plotBoiledEgg(d)
#' }
plotBoiledEgg <- function(data) {

  req_cols <- c("LogP", "TPSA")
  miss <- setdiff(req_cols, names(data))
  if (length(miss) > 0) {
    stop("Missing required columns: ", paste(miss, collapse = ", "),
         call. = FALSE)
  }

  ## Classify points using official polygons (point-in-polygon)
  ## Polygons are stored as (TPSA, LogP) = (x, y)
  logp <- as.numeric(data$LogP)
  tpsa <- as.numeric(data$TPSA)
  data$HIA <- ..point_in_polygon(tpsa, logp, .hia_polygon)
  data$BBB <- ..point_in_polygon(tpsa, logp, .bbb_polygon)

  ## P-gp substrate colours (SwissADME style)
  if ("Pgp substrate" %in% names(data)) {
    pgp_raw <- tolower(trimws(as.character(data[["Pgp substrate"]])))
    data$PGP <- ifelse(
      is.na(data[["Pgp substrate"]]) | pgp_raw == "",
      "NA",
      ifelse(pgp_raw %in% c("yes", "true", "1"),
             "PGP+",
             "PGP-")
    )
  } else {
    data$PGP <- "NA"
  }

  colores <- c(
    "PGP+" = "#2C7FB8",
    "PGP-" = "#D7301F",
    "NA"   = "grey40"
  )

  ## Convert polygon matrices to data.frames for geom_polygon
  ## Polygons are already in (TPSA, LogP) = (x, y) format
  hia_df <- data.frame(x = .hia_polygon[, 1], y = .hia_polygon[, 2])
  bbb_df <- data.frame(x = .bbb_polygon[, 1], y = .bbb_polygon[, 2])

  ggplot(data, aes(x = TPSA, y = LogP)) +
    annotate("rect", xmin = 0, xmax = 200, ymin = -2, ymax = 7,
             fill = "#BDBDBD") +
    geom_polygon(data = hia_df, aes(x = x, y = y), inherit.aes = FALSE,
                 fill = "white", colour = NA) +
    geom_polygon(data = bbb_df, aes(x = x, y = y), inherit.aes = FALSE,
                 fill = "#F5D548", colour = NA) +
    geom_point(aes(fill = PGP), shape = 21, colour = "black",
               stroke = .25, size = 2.8) +
    scale_fill_manual(values = colores, name = "P-gp") +
    coord_cartesian(xlim = c(0, 200), ylim = c(-2, 7)) +
    scale_x_continuous(breaks = seq(0, 200, 20)) +
    scale_y_continuous(breaks = -2:7) +
    labs(title = "BOILED-Egg",
         x = expression(TPSA~(ring(A)^2)),
         y = "LogP") +
    theme_classic(base_size = 13) +
    theme(panel.grid = element_blank(),
          legend.position = "right",
          plot.title = element_text(hjust = .5, face = "bold"),
          axis.title = element_text(face = "bold"),
          axis.line = element_line(colour = "black"),
          legend.title = element_text(face = "bold"))
}

## ---------------------- Descriptor histograms -----------------------------

#' Molecular weight distribution histogram
#'
#' @param data A data.frame containing a numeric \code{MW} column.
#' @return A \pkg{ggplot2} object.
#' @examples
#' \dontrun{
#' d <- data.frame(MW = c(300, 400, 500))
#' plotMW(d)
#' }
#' @export
plotMW <- function(data) {
  ggplot(data, aes(x = MW)) +
    geom_histogram(bins = 30, fill = "#3366CC", color = "white", na.rm = TRUE) +
    labs(title = "Distribution of Molecular Weight",
         x = "MW (g/mol)", y = "Frequency") +
    theme_minimal(base_size = 13)
}

#' TPSA distribution histogram
#'
#' @param data A data.frame containing a numeric \code{TPSA} column.
#' @return A \pkg{ggplot2} object.
#' @examples
#' \dontrun{
#' d <- data.frame(TPSA = c(40, 80, 120))
#' plotTPSA(d)
#' }
#' @export
plotTPSA <- function(data) {
  ggplot(data, aes(x = TPSA)) +
    geom_histogram(bins = 30, fill = "#33A02C", color = "white", na.rm = TRUE) +
    labs(title = "Distribution of TPSA",
         x = expression(TPSA~(ring(A)^2)), y = "Frequency") +
    theme_minimal(base_size = 13)
}

#' LogP distribution histogram
#'
#' @param data A data.frame containing a numeric \code{LogP} column.
#' @return A \pkg{ggplot2} object.
#' @examples
#' \dontrun{
#' d <- data.frame(LogP = c(1, 2, 3))
#' plotLogP(d)
#' }
#' @export
plotLogP <- function(data) {
  ggplot(data, aes(x = LogP)) +
    geom_histogram(bins = 30, fill = "#E31A1C", color = "white", na.rm = TRUE) +
    labs(title = "Distribution of LogP",
         x = "LogP", y = "Frequency") +
    theme_minimal(base_size = 13)
}

## ---------------------- Radar (physicochemical profile) -------------------
## Compares up to 5 molecules simultaneously on 6 key properties, normalized
## 0-1 against the range observed in the (filtered) dataset. Requires the
## 'fmsb' package.

#' Radar plot of physicochemical profile
#'
#' Compares up to 5 molecules simultaneously on a radar chart using 6 key
#' physicochemical properties, normalized 0-1 against the range observed in
#' the full dataset. Requires the suggested package \pkg{fmsb}.
#'
#' @param data A data.frame with the physicochemical properties.
#' @param id_col Character. Name of the identifier column.
#' @param ids Character vector of molecule identifiers to compare (max 5
#'   recommended).
#' @param props Character vector of property column names to display. Defaults
#'   to \code{MW, LogP, TPSA, #H-bond acceptors, #H-bond donors,
#'   #Rotatable bonds}.
#' @return Invisibly \code{NULL}; called for the side-effect of drawing the
#'   radar chart on the current graphics device.
#' @examples
#' \dontrun{
#' d <- data.frame(
#'   Name = c("mol1", "mol2"),
#'   MW = c(300, 400), LogP = c(2, 3), TPSA = c(60, 90),
#'   "#H-bond acceptors" = c(4, 5), "#H-bond donors" = c(2, 1),
#'   "#Rotatable bonds" = c(3, 5), check.names = FALSE)
#' plotRadar(d, id_col = "Name", ids = c("mol1", "mol2"))
#' }
#' @export
#' @references Nakazawa, M. (2019). \emph{fmsb: Functions for Medical
#'   Statistics Book with some Demographic Data}. R package.
plotRadar <- function(data, id_col, ids, props = c(
  "MW", "LogP", "TPSA", "#H-bond acceptors", "#H-bond donors", "#Rotatable bonds"
)) {

  if (!requireNamespace("fmsb", quietly = TRUE)) {
    stop(
      "The 'fmsb' package is required for the radar plot. Install it with: ",
      "install.packages('fmsb')",
      call. = FALSE
    )
  }

  props <- intersect(props, names(data))
  if (length(props) < 3) {
    stop("At least 3 numeric properties are required for the radar plot.",
         call. = FALSE)
  }

  stopifnot(id_col %in% names(data))

  sub <- data[data[[id_col]] %in% ids, c(id_col, props), drop = FALSE]
  sub <- sub[!duplicated(sub[[id_col]]), ]

  if (nrow(sub) == 0) {
    stop("No molecules were found with the selected identifiers.",
         call. = FALSE)
  }

  ## Normalize 0-1 against the range of the full dataset (not just the
  ## selection) so the radar is comparable across selections.
  ranges <- vapply(data[props], function(x) range(x, na.rm = TRUE),
                   numeric(2))

  norm <- sub
  for (p in props) {
    rmin <- ranges[1, p]
    rmax <- ranges[2, p]
    norm[[p]] <- if (rmax > rmin) {
      (sub[[p]] - rmin) / (rmax - rmin)
    } else {
      0.5
    }
  }

  radar_df <- rbind(
    setNames(as.data.frame(matrix(1, nrow = 1, ncol = length(props))), props),
    setNames(as.data.frame(matrix(0, nrow = 1, ncol = length(props))), props),
    norm[props]
  )
  rownames(radar_df) <- c("max", "min", as.character(norm[[id_col]]))

  colors <- rainbow(nrow(norm), s = 0.7, v = 0.85)

  fmsb::radarchart(
    radar_df,
    pcol = colors,
    plwd = 2,
    cglcol = "grey80",
    cglty = 1,
    axislabcol = "grey40",
    vlcex = 0.9,
    title = "Physicochemical profile (radar)"
  )

  legend(
    "topright",
    legend = as.character(norm[[id_col]]),
    col = colors,
    lty = 1, lwd = 2, bty = "n", cex = 0.8
  )

  invisible(NULL)
}

## --------------------------- Tanimoto / AGNES -----------------------------
## Structural similarity clustering using extended fingerprints and the
## AGNES hierarchical clustering algorithm. Requires rcdk, fingerprint and
## cluster.

#' Tanimoto / AGNES structural similarity dendrogram
#'
#' Parses SMILES strings, computes extended fingerprints, builds a Tanimoto
#' similarity matrix and clusters the molecules with AGNES. Requires the
#' suggested packages \pkg{rcdk}, \pkg{fingerprint} and \pkg{cluster}, plus a
#' working Java JDK (for \pkg{rcdk}).
#'
#' @param data A data.frame containing a SMILES column.
#' @param smiles_col Character. Name of the SMILES column.
#' @param label_col Character. Name of the column to use as dendrogram labels,
#'   or \code{NULL} to use the SMILES themselves.
#' @param max_n Integer. Maximum number of molecules to include. Default 40.
#' @param method Character. AGNES linking method: one of \code{"average"},
#'   \code{"complete"}, \code{"single"}, \code{"ward"}. Default
#'   \code{"average"}.
#' @return Invisibly the Tanimoto similarity matrix.
#' @examples
#' \dontrun{
#' d <- data.frame(
#'   Name = c("ethanol", "aspirin", "benzene", "toluene", "phenol"),
#'   SMILES = c("CCO", "CC(=O)Oc1ccccc1C(=O)O", "c1ccccc1",
#'             "Cc1ccccc1", "Oc1ccccc1"),
#'   stringsAsFactors = FALSE)
#' plotTanimoto(d, smiles_col = "SMILES", label_col = "Name")
#' }
#' @export
#' @references Willett, P., Barnard, J. M., & Downs, G. M. (1998). Chemical
#'   similarity searching. \emph{Journal of Chemical Information and Computer
#'   Sciences}, 38(6), 983-996.
plotTanimoto <- function(data, smiles_col, label_col = NULL,
                         max_n = 40, method = "average") {

  if (!requireNamespace("rcdk", quietly = TRUE)) {
    stop("This plot requires the 'rcdk' package:\n",
         "  install.packages('rcdk')", call. = FALSE)
  }
  if (!requireNamespace("fingerprint", quietly = TRUE)) {
    stop("This plot requires the 'fingerprint' package:\n",
         "  install.packages('fingerprint')", call. = FALSE)
  }
  if (!requireNamespace("cluster", quietly = TRUE)) {
    stop("This plot requires the 'cluster' package.", call. = FALSE)
  }

  stopifnot(smiles_col %in% names(data))

  smiles <- trimws(as.character(data[[smiles_col]]))
  if (!is.null(label_col) && label_col %in% names(data)) {
    labels <- trimws(as.character(data[[label_col]]))
  } else {
    labels <- smiles
  }

  valid <- !is.na(smiles) & nzchar(smiles)
  smiles <- smiles[valid]
  labels <- labels[valid]

  if (length(smiles) > max_n) {
    idx <- seq_len(max_n)
    smiles <- smiles[idx]
    labels <- labels[idx]
  }

  mols <- rcdk::parse.smiles(smiles)
  ok <- !vapply(mols, is.null, logical(1))
  if (any(!ok)) {
    warning(sum(!ok), " SMILES could not be parsed and were omitted.")
    smiles <- smiles[ok]
    labels <- labels[ok]
    mols   <- mols[ok]
  }

  if (length(mols) < 3) {
    stop("At least 3 valid SMILES are required for Tanimoto/AGNES.",
         call. = FALSE)
  }

  labels <- make.unique(labels)

  fps <- lapply(mols, rcdk::get.fingerprint, type = "extended")

  sim_mat <- fingerprint::fp.sim.matrix(fps, method = "tanimoto")
  rownames(sim_mat) <- labels
  colnames(sim_mat) <- labels

  dist_mat <- as.dist(1 - sim_mat)

  agnes_fit <- cluster::agnes(dist_mat, method = method)

  plot(
    agnes_fit,
    which.plots = 2,
    labels = labels,
    main = paste0("AGNES (", method, ") on Tanimoto similarity"),
    xlab = "Molecules",
    sub = ""
  )

  invisible(sim_mat)
}

## --------------------------- Correlation Heatmap --------------------------
## Heatmap of correlations between key physicochemical properties; useful to
## detect collinearity before SAR/QSAR modelling.

#' Correlation heatmap of physicochemical properties
#'
#' Produces a heatmap of the Pearson correlation between key physicochemical
#' properties.
#'
#' @param data A data.frame with numeric property columns.
#' @param props Character vector of property column names. Defaults to a set
#'   of common descriptors.
#' @return A \pkg{ggplot2} object.
#' @examples
#' \dontrun{
#' d <- data.frame(MW = c(300, 400, 500), LogP = c(2, 3, 4),
#'   TPSA = c(60, 80, 100))
#' plotCorrHeatmap(d)
#' }
#' @export
plotCorrHeatmap <- function(data, props = c(
  "MW", "LogP", "TPSA", "MR",
  "#H-bond acceptors", "#H-bond donors", "#Rotatable bonds", "#Heavy atoms"
)) {

  props <- intersect(props, names(data))
  if (length(props) < 2) {
    stop("At least 2 numeric properties are required for the heatmap.",
         call. = FALSE)
  }

  ## Coerce to numeric explicitly and remove zero-variance columns
  num_data <- as.data.frame(
    lapply(data[, props, drop = FALSE],
           function(x) suppressWarnings(as.numeric(as.character(x)))),
    stringsAsFactors = FALSE
  )
  names(num_data) <- props

  col_vars <- sapply(num_data, var, na.rm = TRUE)
  zero_var <- names(num_data)[is.na(col_vars) | col_vars == 0]
  if (length(zero_var) > 0) {
    num_data <- num_data[, !(names(num_data) %in% zero_var), drop = FALSE]
  }
  if (ncol(num_data) < 2) {
    stop("After removing zero-variance columns, fewer than 2 variables remain.",
         call. = FALSE)
  }

  cor_mat <- cor(num_data, use = "pairwise.complete.obs")

  cor_df <- as.data.frame(as.table(cor_mat))
  names(cor_df) <- c("Var1", "Var2", "Correlation")

  ## Replace NaN with NA for cleaner display
  cor_df$Correlation[is.nan(cor_df$Correlation)] <- NA

  ## Format label: show "N/A" for missing correlations
  cor_df$Label <- ifelse(is.na(cor_df$Correlation), "N/A",
                         sprintf("%.2f", cor_df$Correlation))

  ggplot(cor_df, aes(x = Var1, y = Var2, fill = Correlation)) +
    geom_tile(color = "white") +
    geom_text(aes(label = Label), size = 3, na.rm = TRUE) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, limits = c(-1, 1), na.value = "grey90") +
    labs(title = "Correlation between physicochemical properties",
         x = NULL, y = NULL) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

## --------------------------- PCA (chemical space) -------------------------

#' Principal Component Analysis (PCA) chemical space plot
#'
#' Computes a PCA on the selected numeric physicochemical variables and
#' produces a scatterplot of the first two principal components, optionally
#' coloured by a grouping variable, with 95% confidence ellipses, observation
#' labels and variable loading arrows. Labels require the suggested package
#' \pkg{ggrepel}.
#'
#' @param data A data.frame with numeric property columns.
#' @param variables Character vector of numeric column names to use in the
#'   PCA. If \code{NULL}, a default set of common descriptors is used.
#' @param color_by Character. Name of the column to colour points by, or
#'   \code{"None"} for a single colour.
#' @param label_by Character. Name of the column to use as point labels, or
#'   \code{"None"} for no labels.
#' @param scale_data Logical. Whether to scale variables to unit variance
#'   before the PCA. Default \code{TRUE}.
#' @param ellipse Logical. Whether to draw 95% confidence ellipses per group.
#'   Default \code{TRUE}.
#' @return A \pkg{ggplot2} object.
#' @examples
#' \dontrun{
#' d <- data.frame(MW = c(300, 400, 500), LogP = c(2, 3, 4),
#'   TPSA = c(60, 80, 100), group = c("A", "A", "B"))
#' plotPCA(d, variables = c("MW", "LogP", "TPSA"), color_by = "group")
#' }
#' @export
plotPCA <- function(data, variables = NULL, color_by = "None",
                    label_by = "None", scale_data = TRUE, ellipse = TRUE) {

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

  ## Re-check that variables are actually numeric
  numeric_vars <- variables[sapply(data[variables], is.numeric)]
  if (length(numeric_vars) < 2) {
    stop("Not enough numeric variables available.", call. = FALSE)
  }

  ## Extract PCA data: explicitly coerce to numeric to avoid issues with
  ## factors or mixed-type columns that report is.numeric = TRUE
  pca_data <- as.data.frame(
    lapply(data[, numeric_vars, drop = FALSE],
           function(x) suppressWarnings(as.numeric(as.character(x)))),
    stringsAsFactors = FALSE
  )
  names(pca_data) <- numeric_vars

  ## Remove rows with NA values
  cc <- complete.cases(pca_data)
  pca_data <- pca_data[cc, , drop = FALSE]

  if (nrow(pca_data) < 3) {
    stop("Too few observations to compute a PCA.", call. = FALSE)
  }

  ## Remove zero-variance columns -- these cause scale() to divide by 0,
  ## which triggers the "length of 'scale'" error in prcomp
  col_vars <- sapply(pca_data, var, na.rm = TRUE)
  zero_var_cols <- names(pca_data)[is.na(col_vars) | col_vars == 0]
  if (length(zero_var_cols) > 0) {
    pca_data <- pca_data[, !(names(pca_data) %in% zero_var_cols), drop = FALSE]
    numeric_vars <- setdiff(numeric_vars, zero_var_cols)
  }
  if (length(numeric_vars) < 2) {
    stop("After removing zero-variance columns, fewer than 2 variables remain. ",
         "Please select different variables.", call. = FALSE)
  }

  ## Also keep the color/label columns aligned with the filtered rows
  df_aux <- data[cc, , drop = FALSE]

  ## Compute PCA -- wrap in tryCatch as a safety net
  pca <- tryCatch(
    prcomp(pca_data, center = TRUE, scale. = scale_data),
    error = function(e) {
      ## If scaling fails, retry without scaling
      prcomp(pca_data, center = TRUE, scale. = FALSE)
    }
  )

  scores <- as.data.frame(pca$x)

  ## Use isTRUE() to safely handle NA values in color_by / label_by
  ## (input selects can return NA before a choice is made)
  if (!is.null(color_by) && !is.na(color_by) && color_by != "None" &&
      color_by %in% names(df_aux)) {
    scores$Colour <- as.factor(df_aux[[color_by]])
  } else {
    scores$Colour <- factor("Compounds")
  }

  if (!is.null(label_by) && !is.na(label_by) && label_by != "None" &&
      label_by %in% names(df_aux)) {
    scores$Label <- as.character(df_aux[[label_by]])
  } else {
    scores$Label <- ""
  }

  loadings <- as.data.frame(pca$rotation)
  loadings$Variable <- rownames(loadings)

  ## Scale the loading arrows to fit the plot; handle edge cases where
  ## the range is zero (all loadings identical for a PC)
  range_pc1_scores <- max(scores$PC1) - min(scores$PC1)
  range_pc2_scores <- max(scores$PC2) - min(scores$PC2)
  range_pc1_load <- max(loadings$PC1) - min(loadings$PC1)
  range_pc2_load <- max(loadings$PC2) - min(loadings$PC2)

  mult <- if (range_pc1_load > 0 && range_pc2_load > 0) {
    min(range_pc1_scores / range_pc1_load, range_pc2_scores / range_pc2_load)
  } else if (range_pc1_load > 0) {
    range_pc1_scores / range_pc1_load
  } else if (range_pc2_load > 0) {
    range_pc2_scores / range_pc2_load
  } else {
    1
  }
  mult <- mult * 0.35
  loadings$PC1 <- loadings$PC1 * mult
  loadings$PC2 <- loadings$PC2 * mult

  var_exp <- summary(pca)$importance[2, ]
  xlabel <- paste0("PC1 (", round(var_exp[1] * 100, 1), "%)")
  ylabel <- paste0("PC2 (", round(var_exp[2] * 100, 1), "%)")

  g <- ggplot(scores, aes(PC1, PC2, colour = Colour)) +
    geom_point(size = 3, alpha = .9) +
    labs(title = "Chemical Space (PCA)", x = xlabel, y = ylabel, colour = NULL) +
    theme_classic(base_size = 13) +
    theme(plot.title = element_text(face = "bold", hjust = .5))

  ## Only draw ellipses if there are multiple groups AND each group has
  ## at least 3 points (otherwise stat_ellipse fails with "Too few points")
  if (ellipse && length(unique(scores$Colour)) > 1) {
    group_sizes <- table(scores$Colour)
    if (all(group_sizes >= 3)) {
      g <- g + stat_ellipse(level = .95, linewidth = .7)
    }
  }

  if (any(scores$Label != "")) {
    if (!requireNamespace("ggrepel", quietly = TRUE)) {
      warning("The 'ggrepel' package is required for non-overlapping labels; ",
              "install it with install.packages('ggrepel'). Labels will be ",
              "drawn with geom_text instead.")
      g <- g + geom_text(aes(label = Label), size = 3, show.legend = FALSE)
    } else {
      g <- g + ggrepel::geom_text_repel(
        aes(label = Label), size = 3, max.overlaps = 50, show.legend = FALSE
      )
    }
  }

  ## Loading arrows -- remove any rows with NA to avoid "Removed rows" warnings
  loadings_clean <- loadings[complete.cases(loadings[, c("PC1", "PC2")]), ,
                             drop = FALSE]

  g <- g +
    geom_segment(
      data = loadings_clean,
      aes(x = 0, y = 0, xend = PC1, yend = PC2),
      inherit.aes = FALSE,
      arrow = arrow(length = unit(0.25, "cm")),
      colour = "grey35"
    ) +
    geom_text(
      data = loadings_clean,
      aes(PC1, PC2, label = Variable),
      inherit.aes = FALSE,
      colour = "black", fontface = "bold", size = 4
    )

  g
}

## --------------------------- t-SNE (chemical space) -----------------------

#' t-SNE chemical space plot
#'
#' Computes a 2-dimensional t-SNE embedding of the selected numeric
#' physicochemical variables and produces a scatterplot, optionally coloured
#' and labelled. Requires the suggested packages \pkg{Rtsne} and
#' \pkg{ggrepel}.
#'
#' @param data A data.frame with numeric property columns.
#' @param variables Character vector of numeric column names to use.
#' @param color_by Character. Column to colour points by, or \code{"None"}.
#' @param label_by Character. Column to use as labels, or \code{"None"}.
#' @param perplexity Numeric. t-SNE perplexity. Default 30.
#' @param max_iter Integer. Number of iterations. Default 1000.
#' @param scale_data Logical. Whether to scale variables. Default \code{TRUE}.
#' @return A \pkg{ggplot2} object.
#' @examples
#' \dontrun{
#' d <- data.frame(MW = c(300, 400, 500), LogP = c(2, 3, 4),
#'   TPSA = c(60, 80, 100), group = c("A", "A", "B"))
#' plotTSNE(d, variables = c("MW", "LogP", "TPSA"), color_by = "group")
#' }
#' @export
plotTSNE <- function(data, variables, color_by = "None", label_by = "None",
                     perplexity = 30, max_iter = 1000, scale_data = TRUE) {

  if (!requireNamespace("Rtsne", quietly = TRUE)) {
    stop("Install the 'Rtsne' package: install.packages('Rtsne')",
         call. = FALSE)
  }
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop("Install the 'ggrepel' package: install.packages('ggrepel')",
         call. = FALSE)
  }

  if (is.null(variables) || length(variables) < 2) {
    stop("Select at least 2 variables for t-SNE.", call. = FALSE)
  }

  keep <- complete.cases(data[, variables, drop = FALSE])
  X <- data[keep, variables, drop = FALSE]
  original <- data[keep, , drop = FALSE]

  if (nrow(X) < 3) {
    stop("At least 3 observations are required for t-SNE.", call. = FALSE)
  }

  X <- as.data.frame(lapply(X, function(x) as.numeric(as.character(x))))
  keep2 <- complete.cases(X)
  X <- X[keep2, , drop = FALSE]
  original <- original[keep2, , drop = FALSE]

  if (nrow(X) < 3) {
    stop("After cleaning NAs there are not enough observations.",
         call. = FALSE)
  }

  ## Remove zero-variance columns -- scale() produces NaN for these
  col_vars <- sapply(X, var, na.rm = TRUE)
  zero_var_cols <- names(X)[is.na(col_vars) | col_vars == 0]
  if (length(zero_var_cols) > 0) {
    X <- X[, !(names(X) %in% zero_var_cols), drop = FALSE]
  }
  if (ncol(X) < 2) {
    stop("After removing zero-variance columns, fewer than 2 variables remain.",
         call. = FALSE)
  }

  if (scale_data) {
    X <- scale(X)
    ## Replace any NaN/Inf produced by scaling (shouldn't happen after
    ## zero-variance removal, but as a safety net)
    X[is.nan(X) | is.infinite(X)] <- 0
  }

  max_perplexity <- floor((nrow(X) - 1) / 3)
  perplexity <- min(perplexity, max_perplexity)
  if (perplexity < 2) {
    stop("Perplexity is too low for the number of observations.",
         call. = FALSE)
  }

  tsne <- Rtsne::Rtsne(
    X, dims = 2, perplexity = perplexity, max_iter = max_iter,
    check_duplicates = FALSE
  )

  plot_df <- data.frame(TSNE1 = tsne$Y[, 1], TSNE2 = tsne$Y[, 2])

  if (!is.null(color_by) && !is.na(color_by) && color_by != "None" &&
      color_by %in% names(original)) {
    plot_df$Color <- original[[color_by]]
  }
  if (!is.null(label_by) && !is.na(label_by) && label_by != "None" &&
      label_by %in% names(original)) {
    plot_df$Label <- original[[label_by]]
  }

  p <- ggplot(plot_df, aes(TSNE1, TSNE2)) + geom_point(size = 3)

  if ("Color" %in% names(plot_df)) {
    p <- p + aes(color = Color)
  }
  if ("Label" %in% names(plot_df)) {
    p <- p + ggrepel::geom_text_repel(
      aes(label = Label), size = 3, max.overlaps = Inf
    )
  }

  p + labs(title = "t-SNE Chemical Space", x = "t-SNE 1", y = "t-SNE 2") +
    theme_classic(base_size = 13)
}

## --------------------------- UMAP ----------------------------------------

#' UMAP chemical space plot
#'
#' Produces a 2D UMAP projection of the selected numeric variables for
#' exploratory analysis of the chemical space. Points can be coloured and
#' labelled. Requires the suggested packages \pkg{uwot} and \pkg{ggrepel}.
#'
#' @param data A data.frame with numeric property columns.
#' @param variables Character vector of numeric column names to use.
#' @param color_by Character. Column to colour points by, or \code{"None"}.
#' @param label_by Character. Column to use as labels, or \code{"None"}.
#' @param n_neighbors Integer. Number of nearest neighbours used by UMAP.
#'   Default 15.
#' @param min_dist Numeric. Minimum distance between points in the embedding.
#'   Default 0.1.
#' @param scale_data Logical. Whether to scale variables. Default \code{TRUE}.
#' @return A \pkg{ggplot2} object.
#' @examples
#' \dontrun{
#' d <- data.frame(MW = c(300, 400, 500), LogP = c(2, 3, 4),
#'   TPSA = c(60, 80, 100), group = c("A", "A", "B"))
#' plotUMAP(d, variables = c("MW", "LogP", "TPSA"), color_by = "group")
#' }
#' @export
plotUMAP <- function(data, variables, color_by = "None", label_by = "None",
                     n_neighbors = 15, min_dist = 0.1, scale_data = TRUE) {

  if (!requireNamespace("uwot", quietly = TRUE)) {
    stop("Install the 'uwot' package: install.packages('uwot')",
         call. = FALSE)
  }
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop("Install the 'ggrepel' package: install.packages('ggrepel')",
         call. = FALSE)
  }

  if (is.null(variables) || length(variables) < 2) {
    stop("Select at least 2 variables for UMAP.", call. = FALSE)
  }

  keep <- complete.cases(data[, variables, drop = FALSE])
  X <- data[keep, variables, drop = FALSE]
  original <- data[keep, , drop = FALSE]

  if (nrow(X) < 4) {
    stop("At least 4 observations are required for UMAP.", call. = FALSE)
  }

  X <- as.data.frame(lapply(X, function(x) as.numeric(as.character(x))))
  keep2 <- complete.cases(X)
  X <- X[keep2, , drop = FALSE]
  original <- original[keep2, , drop = FALSE]

  if (nrow(X) < 4) {
    stop("After cleaning NAs there are not enough observations.",
         call. = FALSE)
  }

  ## Remove zero-variance columns -- scale() produces NaN for these
  col_vars <- sapply(X, var, na.rm = TRUE)
  zero_var_cols <- names(X)[is.na(col_vars) | col_vars == 0]
  if (length(zero_var_cols) > 0) {
    X <- X[, !(names(X) %in% zero_var_cols), drop = FALSE]
  }
  if (ncol(X) < 2) {
    stop("After removing zero-variance columns, fewer than 2 variables remain.",
         call. = FALSE)
  }

  if (scale_data) {
    X <- scale(X)
    ## Replace any NaN/Inf produced by scaling
    X[is.nan(X) | is.infinite(X)] <- 0
  }

  ## Clamp n_neighbors to the valid range for uwot::umap
  n_obs <- nrow(X)
  n_neighbors <- as.integer(min(max(n_neighbors, 2), n_obs - 1))
  if (n_neighbors < 2) {
    stop("n_neighbors is too small for the number of observations.",
         call. = FALSE)
  }

  ## uwot::umap returns a matrix with `n_components` columns
  set.seed(42)
  emb <- uwot::umap(
    X, n_neighbors = n_neighbors, min_dist = min_dist,
    n_components = 2, init = "spectral", verbose = FALSE
  )

  plot_df <- data.frame(UMAP1 = emb[, 1], UMAP2 = emb[, 2])

  if (!is.null(color_by) && !is.na(color_by) && color_by != "None" &&
      color_by %in% names(original)) {
    plot_df$Color <- original[[color_by]]
  }
  if (!is.null(label_by) && !is.na(label_by) && label_by != "None" &&
      label_by %in% names(original)) {
    plot_df$Label <- original[[label_by]]
  }

  p <- ggplot(plot_df, aes(UMAP1, UMAP2)) + geom_point(size = 3)

  if ("Color" %in% names(plot_df)) {
    p <- p + aes(color = Color)
  }
  if ("Label" %in% names(plot_df)) {
    p <- p + ggrepel::geom_text_repel(
      aes(label = Label), size = 3, max.overlaps = Inf
    )
  }

  p + labs(title = "UMAP Chemical Space", x = "UMAP 1", y = "UMAP 2") +
    theme_classic(base_size = 13)
}

## --------------------------- Parallel coordinates -------------------------

#' Parallel coordinates plot
#'
#' Produces a parallel coordinates plot of the selected numeric variables,
#' optionally coloured by a grouping variable. Requires the suggested package
#' \pkg{GGally}.
#'
#' @param data A data.frame with numeric property columns.
#' @param variables Character vector of numeric column names to display.
#' @param color_by Character. Column to colour lines by, or \code{NULL}/
#'   \code{"None"} for a single colour.
#' @param scale_data Logical. Whether to standardize variables. Default
#'   \code{TRUE}.
#' @return A \pkg{ggplot2} object.
#' @examples
#' \dontrun{
#' d <- data.frame(MW = c(300, 400, 500), LogP = c(2, 3, 4),
#'   TPSA = c(60, 80, 100), group = c("A", "A", "B"))
#' plotParallel(d, variables = c("MW", "LogP", "TPSA"), color_by = "group")
#' }
#' @export
plotParallel <- function(data, variables, color_by = NULL,
                         scale_data = TRUE) {

  if (!requireNamespace("GGally", quietly = TRUE)) {
    stop("This plot requires the 'GGally' package:\n",
         "  install.packages('GGally')", call. = FALSE)
  }

  if (length(variables) < 2) {
    stop("Select at least two variables.", call. = FALSE)
  }

  ## Coerce to numeric and remove zero-variance columns
  df <- as.data.frame(
    lapply(data[, variables, drop = FALSE],
           function(x) suppressWarnings(as.numeric(as.character(x)))),
    stringsAsFactors = FALSE
  )
  names(df) <- variables

  col_vars <- sapply(df, var, na.rm = TRUE)
  zero_var <- names(df)[is.na(col_vars) | col_vars == 0]
  if (length(zero_var) > 0) {
    df <- df[, !(names(df) %in% zero_var), drop = FALSE]
    variables <- setdiff(variables, zero_var)
  }
  if (length(variables) < 2) {
    stop("After removing zero-variance columns, fewer than 2 variables remain.",
         call. = FALSE)
  }

  ## Remove rows with NA in the numeric columns
  cc <- complete.cases(df)
  df <- df[cc, , drop = FALSE]

  if (scale_data) {
    scale_type <- "std"
  } else {
    scale_type <- "uniminmax"
  }

  if (!is.null(color_by) && !is.na(color_by) && color_by != "None" &&
      color_by %in% names(data)) {
    df$Group <- as.factor(data[[color_by]][cc])
    p <- GGally::ggparcoord(
      data = df, columns = seq_along(variables),
      groupColumn = "Group", scale = scale_type, alphaLines = 0.5
    ) +
      theme_bw() +
      labs(x = "Descriptors", y = "Scaled value", color = color_by)
  } else {
    p <- GGally::ggparcoord(
      data = df, columns = seq_along(variables),
      scale = scale_type, alphaLines = 0.4
    ) +
      theme_bw() +
      labs(x = "Descriptors", y = "Scaled value")
  }

  p
}

## --------------------------- Violin plot ----------------------------------

#' Violin plot
#'
#' Produces a violin plot of a numeric variable grouped by a categorical
#' variable, with an optional inner boxplot and jittered points.
#'
#' @param data A data.frame.
#' @param variable Character. Name of the numeric variable to plot.
#' @param group_by Character. Name of the categorical grouping variable.
#' @param show_box Logical. Whether to overlay a boxplot. Default \code{TRUE}.
#' @param show_points Logical. Whether to overlay jittered points. Default
#'   \code{FALSE}.
#' @return A \pkg{ggplot2} object.
#' @examples
#' \dontrun{
#' d <- data.frame(MW = c(300, 400, 500),
#'   group = c("A", "A", "B"))
#' plotViolin(d, variable = "MW", group_by = "group")
#' }
#' @export
plotViolin <- function(data, variable, group_by,
                       show_box = TRUE, show_points = FALSE) {

  stopifnot(variable %in% names(data))
  stopifnot(group_by %in% names(data))

  p <- ggplot(data, aes(x = .data[[group_by]], y = .data[[variable]],
                        fill = .data[[group_by]])) +
    geom_violin(trim = FALSE, alpha = 0.7)

  if (show_box) {
    p <- p + geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.5)
  }
  if (show_points) {
    p <- p + geom_jitter(width = 0.1, alpha = 0.4, size = 1.3)
  }

  p + theme_bw() +
    labs(x = group_by, y = variable, fill = group_by) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))
}
