# ---------------------------------------------------------------------------
# report.R
# Comprehensive report generation for ADMETShiny.
#
# This module provides:
#   - Per-dataset statistics (molecules uploaded, filtered, property summaries)
#   - Drug-likeness filter pass rates
#   - BOILED-Egg ADMET classification
#   - Additional literature-supported metrics (Pfizer 3/75, GSK 4/400, etc.)
#   - Composite drug-likeness score (0-100)
#   - Report-specific plots (violations summary, score distribution)
#   - Cross-dataset comparison
#   - Multi-format rendering (HTML, PDF, Word) via rmarkdown
# ---------------------------------------------------------------------------

## ===================== Internal helpers ===================================

#' Safely extract a numeric column
#' @keywords internal
.get_col <- function(data, col) {
  if (col %in% names(data)) as.numeric(data[[col]]) else NULL
}

#' Summary statistics for a numeric vector
#' @keywords internal
.numeric_summary <- function(x) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) {
    return(c(Mean = NA, Median = NA, SD = NA, Min = NA, Max = NA))
  }
  c(
    Mean   = round(mean(x, na.rm = TRUE), 2),
    Median = round(median(x, na.rm = TRUE), 2),
    SD     = round(sd(x, na.rm = TRUE), 2),
    Min    = round(min(x, na.rm = TRUE), 2),
    Max    = round(max(x, na.rm = TRUE), 2)
  )
}

#' Safe rbind that filters out NULL entries
#'
#' \code{do.call(rbind, ...)} fails with "'dimnames' applied to non-array"
#' when the list contains NULL entries. This helper removes them first.
#' @keywords internal
.safe_rbind <- function(list_of_dfs) {
  list_of_dfs <- Filter(function(x) !is.null(x) && is.data.frame(x), list_of_dfs)
  if (length(list_of_dfs) == 0) return(NULL)
  do.call(rbind, list_of_dfs)
}

#' Build a Markdown pipe table from a data.frame or matrix
#' @keywords internal
.md_table <- function(df) {
  if (is.null(df) || (is.data.frame(df) && nrow(df) == 0)) return("(no data)")
  if (!is.data.frame(df)) df <- as.data.frame(df, stringsAsFactors = FALSE)
  if (nrow(df) == 0) return("(no data)")
  header <- paste0("| ", paste(names(df), collapse = " | "), " |")
  separator <- paste0("|", paste(rep("---|", ncol(df)), collapse = ""))
  rows <- apply(df, 1, function(row) {
    row[is.na(row)] <- ""
    paste0("| ", paste(row, collapse = " | "), " |")
  })
  paste(c(header, separator, rows), collapse = "\n")
}

## ===================== Statistics functions ==============================

#' Compute per-dataset statistics for the report
#'
#' @param data_raw data.frame of raw uploaded data (NULL if not loaded).
#' @param data_filtered data.frame of filtered data (NULL if not filtered).
#' @param filters_applied Character vector of filter names applied.
#' @param source_name Character; human-readable dataset source name.
#' @return A list with all computed statistics.
#' @keywords internal
computeDatasetStats <- function(data_raw, data_filtered, filters_applied,
                                source_name) {

  if (is.null(data_raw) || nrow(data_raw) == 0) return(NULL)

  filtered <- if (!is.null(data_filtered) && nrow(data_filtered) > 0) data_filtered else data_raw

  ## Key physicochemical properties for the summary table
  prop_cols <- c(
    "MW", "LogP", "TPSA", "MR",
    "#H-bond acceptors", "#H-bond donors",
    "#Rotatable bonds", "#Heavy atoms", "#Aromatic heavy atoms"
  )

  property_stats <- .safe_rbind(lapply(prop_cols, function(col) {
    vals <- .get_col(filtered, col)
    s <- .numeric_summary(vals)
    data.frame(Property = col, Mean = s["Mean"], Median = s["Median"],
               SD = s["SD"], Min = s["Min"], Max = s["Max"],
               stringsAsFactors = FALSE, row.names = NULL)
  }))
  if (!is.null(property_stats)) rownames(property_stats) <- NULL

  ## Drug-likeness filter pass rates
  rules <- c("Lipinski", "Ghose", "Veber", "Egan", "Muegge")
  filter_stats <- .safe_rbind(lapply(rules, function(rule) {
    col <- paste0(rule, " #violations")
    if (col %in% names(filtered)) {
      v <- as.numeric(filtered[[col]])
      n_total <- sum(!is.na(v))
      n_pass <- sum(v == 0, na.rm = TRUE)
      mean_v <- if (n_total > 0) round(mean(v, na.rm = TRUE), 2) else NA
      pct_pass <- if (n_total > 0) round(n_pass / n_total * 100, 1) else NA
    } else {
      n_total <- NA
      n_pass <- NA
      mean_v <- NA
      pct_pass <- NA
    }
    data.frame(
      Filter = rule,
      `Molecules Pass` = n_pass,
      `Total` = n_total,
      `% Pass` = pct_pass,
      `Mean Violations` = mean_v,
      check.names = FALSE, stringsAsFactors = FALSE, row.names = NULL
    )
  }))
  if (!is.null(filter_stats)) rownames(filter_stats) <- NULL

  ## BOILED-Egg ADMET classification
  admet_classes <- c("GI absorption", "BBB permeant", "Pgp substrate")
  admet_stats <- .safe_rbind(lapply(admet_classes, function(col) {
    if (col %in% names(filtered)) {
      vals <- filtered[[col]]
      vals <- vals[!is.na(vals)]
      n_total <- length(vals)
      ## Use named vector instead of table() to avoid 'dimnames' issues
      ## with empty/all-NA vectors
      counts <- c("High" = 0, "Low" = 0, "Yes" = 0, "No" = 0)
      if (n_total > 0) {
        tab <- table(vals)
        for (cat in names(counts)) {
          if (cat %in% names(tab)) counts[cat] <- as.integer(tab[cat])
        }
      }
    } else {
      n_total <- NA
      counts <- c("High" = NA, "Low" = NA, "Yes" = NA, "No" = NA)
    }
    if (col == "GI absorption") {
      data.frame(
        Property = col,
        `High` = counts["High"],
        `Low` = counts["Low"],
        `Yes` = NA,
        `No` = NA,
        check.names = FALSE, stringsAsFactors = FALSE, row.names = NULL
      )
    } else {
      data.frame(
        Property = col,
        `High` = NA,
        `Low` = NA,
        `Yes` = counts["Yes"],
        `No` = counts["No"],
        check.names = FALSE, stringsAsFactors = FALSE, row.names = NULL
      )
    }
  }))
  if (!is.null(admet_stats)) rownames(admet_stats) <- NULL

  ## Additional metrics
  additional <- computeAdditionalMetrics(filtered)

  ## Composite drug-likeness score
  score <- computeDruglikenessScore(filtered)

  list(
    source_name = source_name,
    n_raw = nrow(data_raw),
    n_filtered = nrow(filtered),
    pct_filtered = round(nrow(filtered) / nrow(data_raw) * 100, 1),
    filters_applied = filters_applied,
    property_stats = property_stats,
    filter_stats = filter_stats,
    admet_stats = admet_stats,
    additional = additional,
    score = score
  )
}

#' Compute a composite drug-likeness score (0-100)
#'
#' For each compound, counts how many of the five drug-likeness rules
#' (Lipinski, Ghose, Veber, Egan, Muegge) have zero violations, then scales
#' to 0-100. A compound passing all five rules scores 100; one passing none
#' scores 0.
#'
#' @param data A data.frame with violation columns.
#' @return A list with per-compound scores and summary statistics.
#' @keywords internal
computeDruglikenessScore <- function(data) {

  rules <- c("Lipinski", "Ghose", "Veber", "Egan", "Muegge")
  pass_cols <- paste0(rules, " #violations")

  n <- nrow(data)
  if (is.null(n) || n == 0) {
    return(list(
      per_compound = numeric(0),
      mean   = NA, median = NA, sd = NA,
      excellent = NA, acceptable = NA, poor = NA
    ))
  }

  ## Build a matrix of pass/fail (1/0) per compound per rule
  passes <- matrix(NA_integer_, nrow = n, ncol = length(rules))
  colnames(passes) <- rules
  for (i in seq_along(rules)) {
    col <- pass_cols[i]
    if (col %in% names(data)) {
      passes[, i] <- as.integer(as.numeric(data[[col]]) == 0)
    }
  }

  n_rules <- rowSums(!is.na(passes))
  score <- rowSums(passes, na.rm = TRUE) / pmax(n_rules, 1) * 100
  score[n_rules == 0] <- NA

  list(
    per_compound = score,
    mean   = round(mean(score, na.rm = TRUE), 1),
    median = round(median(score, na.rm = TRUE), 1),
    sd     = round(sd(score, na.rm = TRUE), 1),
    excellent  = round(mean(score >= 80, na.rm = TRUE) * 100, 1),
    acceptable = round(mean(score >= 60 & score < 80, na.rm = TRUE) * 100, 1),
    poor       = round(mean(score < 60, na.rm = TRUE) * 100, 1)
  )
}

#' Compute additional literature-supported drug-likeness metrics
#'
#' @param data A data.frame with physicochemical properties.
#' @return A data.frame with metric name, \% compliance, threshold, and reference.
#' @keywords internal
computeAdditionalMetrics <- function(data) {

  mw   <- .get_col(data, "MW")
  logp <- .get_col(data, "LogP")
  tpsa <- .get_col(data, "TPSA")
  hba  <- .get_col(data, "#H-bond acceptors")
  hbd  <- .get_col(data, "#H-bond donors")

  metrics <- list()

  ## Helper: safe percentage (handles n = 0)
  .pct <- function(pass, n) {
    if (is.null(n) || n == 0) return(NA_real_)
    round(pass / n * 100, 1)
  }

  ## 1. Pfizer 3/75 rule (Hughes et al., 2008)
  ## Compounds with TPSA >= 75 and LogP <= 3 have lower in vivo toxicity risk.
  if (!is.null(tpsa) && !is.null(logp)) {
    n <- sum(!is.na(tpsa) & !is.na(logp))
    pass <- sum(tpsa >= 75 & logp <= 3, na.rm = TRUE)
    metrics[[length(metrics) + 1]] <- data.frame(
      Metric = "Pfizer 3/75 rule",
      `Compliance (%)` = .pct(pass, n),
      Threshold = "TPSA >= 75 AND LogP <= 3",
      Reference = "Hughes et al., 2008",
      check.names = FALSE, stringsAsFactors = FALSE, row.names = NULL
    )
  }

  ## 2. GSK 4/400 rule (Gleeson et al., 2011)
  ## Compounds with LogP <= 4 and MW <= 400 have lower promiscuity.
  if (!is.null(mw) && !is.null(logp)) {
    n <- sum(!is.na(mw) & !is.na(logp))
    pass <- sum(mw <= 400 & logp <= 4, na.rm = TRUE)
    metrics[[length(metrics) + 1]] <- data.frame(
      Metric = "GSK 4/400 rule",
      `Compliance (%)` = .pct(pass, n),
      Threshold = "MW <= 400 AND LogP <= 4",
      Reference = "Gleeson et al., 2011",
      check.names = FALSE, stringsAsFactors = FALSE, row.names = NULL
    )
  }

  ## 3. Lead-likeness (Teague et al., 1999)
  ## Lead compounds: MW <= 450, LogP <= 4.5
  if (!is.null(mw) && !is.null(logp)) {
    n <- sum(!is.na(mw) & !is.na(logp))
    pass <- sum(mw <= 450 & logp <= 4.5, na.rm = TRUE)
    metrics[[length(metrics) + 1]] <- data.frame(
      Metric = "Lead-likeness",
      `Compliance (%)` = .pct(pass, n),
      Threshold = "MW <= 450 AND LogP <= 4.5",
      Reference = "Teague et al., 1999",
      check.names = FALSE, stringsAsFactors = FALSE, row.names = NULL
    )
  }

  ## 4. Abb oral bioavailability proxy (Egan et al., 2000)
  ## TPSA <= 131.6 and LogP <= 5.88 (same as Egan filter)
  if (!is.null(tpsa) && !is.null(logp)) {
    n <- sum(!is.na(tpsa) & !is.na(logp))
    pass <- sum(tpsa <= 131.6 & logp <= 5.88, na.rm = TRUE)
    metrics[[length(metrics) + 1]] <- data.frame(
      Metric = "Egan bioavailability proxy",
      `Compliance (%)` = .pct(pass, n),
      Threshold = "TPSA <= 131.6 AND LogP <= 5.88",
      Reference = "Egan et al., 2000",
      check.names = FALSE, stringsAsFactors = FALSE, row.names = NULL
    )
  }

  ## 5. Lipinski Rule-of-Five (original, Lipinski et al., 1997)
  if (!is.null(mw) && !is.null(logp) && !is.null(hba) && !is.null(hbd)) {
    n <- sum(!is.na(mw) & !is.na(logp) & !is.na(hba) & !is.na(hbd))
    pass <- sum(mw <= 500 & logp <= 5 & hba <= 10 & hbd <= 5, na.rm = TRUE)
    metrics[[length(metrics) + 1]] <- data.frame(
      Metric = "Lipinski Ro5 (original)",
      `Compliance (%)` = .pct(pass, n),
      Threshold = "MW <= 500, LogP <= 5, HBA <= 10, HBD <= 5",
      Reference = "Lipinski et al., 1997",
      check.names = FALSE, stringsAsFactors = FALSE, row.names = NULL
    )
  }

  ## 6. Veber oral bioavailability (Veber et al., 2002)
  if (!is.null(tpsa)) {
    rb <- .get_col(data, "#Rotatable bonds")
    if (!is.null(rb)) {
      n <- sum(!is.na(tpsa) & !is.na(rb))
      pass <- sum(rb <= 10 & tpsa <= 140, na.rm = TRUE)
      metrics[[length(metrics) + 1]] <- data.frame(
        Metric = "Veber oral bioavailability",
        `Compliance (%)` = .pct(pass, n),
        Threshold = "RB <= 10 AND TPSA <= 140",
        Reference = "Veber et al., 2002",
        check.names = FALSE, stringsAsFactors = FALSE, row.names = NULL
      )
    }
  }

  if (length(metrics) == 0) return(NULL)
  .safe_rbind(metrics)
}

## ===================== Report-specific plots =============================

#' Violations summary bar chart
#'
#' Produces a stacked bar chart showing, for each drug-likeness rule, the
#' distribution of compounds by number of violations (0, 1, 2, 3+).
#'
#' @param data A data.frame with violation columns.
#' @return A ggplot2 object.
#' @keywords internal
plotViolationsSummary <- function(data) {

  rules <- c("Lipinski", "Ghose", "Veber", "Egan", "Muegge")
  viol_cols <- paste0(rules, " #violations")

  existing <- viol_cols[viol_cols %in% names(data)]
  if (length(existing) == 0) {
    return(ggplot() + annotate("text", x = 1, y = 1,
                               label = "No violation data available") +
             theme_void() +
             labs(title = "Drug-likeness Violations Summary"))
  }

  plot_df <- .safe_rbind(lapply(existing, function(col) {
    v <- as.numeric(data[[col]])
    ## Remove NA values before cutting -- they would produce NA factors
    v <- v[!is.na(v)]
    if (length(v) == 0) return(NULL)
    v_cat <- cut(v, breaks = c(-0.5, 0.5, 1.5, 2.5, Inf),
                 labels = c("0", "1", "2", "3+"))
    tab <- as.data.frame(table(Violations = v_cat), stringsAsFactors = FALSE)
    tab$Rule <- gsub(" #violations", "", col)
    tab
  }))

  if (is.null(plot_df) || nrow(plot_df) == 0) {
    return(ggplot() + annotate("text", x = 1, y = 1,
                               label = "No violation data available") +
             theme_void() +
             labs(title = "Drug-likeness Violations Summary"))
  }

  ## Ensure Freq is numeric (table() can return integer)
  plot_df$Freq <- as.numeric(plot_df$Freq)

  ## Ensure all violation levels are present (fill = 0 if missing)
  all_levels <- c("0", "1", "2", "3+")
  plot_df$Violations <- factor(plot_df$Violations, levels = all_levels)
  plot_df$Rule <- factor(plot_df$Rule, levels = rules)

  colors <- c("0" = "#2ecc71", "1" = "#f1c40f", "2" = "#e67e22", "3+" = "#e74c3c")

  ggplot(plot_df, aes(x = Rule, y = Freq, fill = Violations)) +
    geom_bar(stat = "identity", position = "stack", na.rm = TRUE) +
    scale_fill_manual(values = colors, name = "Violations") +
    labs(title = "Drug-likeness Violations Summary",
         x = NULL, y = "Number of Compounds") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          axis.text.x = element_text(angle = 0, hjust = 0.5))
}

#' Drug-likeness composite score distribution
#'
#' Produces a histogram of the composite drug-likeness score (0-100) with
#' color-coded ranges: poor (<60), acceptable (60-79), excellent (>=80).
#'
#' @param data A data.frame with violation columns.
#' @return A ggplot2 object.
#' @keywords internal
plotDruglikenessScore <- function(data) {

  score_info <- computeDruglikenessScore(data)
  scores <- score_info$per_compound

  if (all(is.na(scores))) {
    return(ggplot() + annotate("text", x = 50, y = 1,
                               label = "No score data available") +
             theme_void() +
             labs(title = "Composite Drug-likeness Score Distribution"))
  }

  ## Remove NA scores
  scores <- scores[!is.na(scores)]
  if (length(scores) == 0) {
    return(ggplot() + annotate("text", x = 50, y = 1,
                               label = "No score data available") +
             theme_void() +
             labs(title = "Composite Drug-likeness Score Distribution"))
  }

  ## Manual binning into 5-point bins (0-5, 5-10, ..., 95-100).
  ## This is more robust than geom_histogram with fill = Category, which
  ## groups by the fill variable and misaligns bins.
  bin_breaks <- seq(0, 100, by = 5)
  scores_cut <- cut(scores, breaks = bin_breaks, include.lowest = TRUE)

  bin_counts <- as.data.frame(table(Bin = scores_cut), stringsAsFactors = FALSE)
  bin_counts$BinMid <- (bin_breaks[-length(bin_breaks)] + bin_breaks[-1]) / 2
  bin_counts$Freq <- as.numeric(bin_counts$Freq)

  ## Assign each bin a category based on its midpoint
  bin_counts$Category <- cut(bin_counts$BinMid,
                             breaks = c(-1, 60, 80, 101),
                             labels = c("Poor (<60)", "Acceptable (60-79)", "Excellent (>=80)"))

  ## Remove empty bins for cleaner rendering
  bin_counts <- bin_counts[bin_counts$Freq > 0, , drop = FALSE]

  if (nrow(bin_counts) == 0) {
    return(ggplot() + annotate("text", x = 50, y = 1,
                               label = "No score data available") +
             theme_void() +
             labs(title = "Composite Drug-likeness Score Distribution"))
  }

  colors <- c("Poor (<60)" = "#e74c3c",
              "Acceptable (60-79)" = "#f1c40f",
              "Excellent (>=80)" = "#2ecc71")

  ggplot(bin_counts, aes(x = BinMid, y = Freq, fill = Category)) +
    geom_col(color = "white", width = 4.5) +
    scale_fill_manual(values = colors, name = "Category") +
    scale_x_continuous(breaks = seq(0, 100, by = 10), limits = c(-2.5, 102.5)) +
    labs(title = "Composite Drug-likeness Score Distribution",
         x = "Score (0-100)", y = "Number of Compounds") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
}

## ===================== Plot generation for report ========================

#' Generate and save all plots for a dataset as PNG files
#'
#' @param data A data.frame (filtered data).
#' @param prefix Character; prefix for plot filenames.
#' @param plot_dir Directory where plots will be saved.
#' @return A named list of plot file paths.
#' @keywords internal
generateReportPlots <- function(data, prefix, plot_dir) {

  paths <- list()

  save_plot <- function(name, plot_obj, w = 7, h = 5) {
    path <- file.path(plot_dir, paste0(prefix, "_", name, ".png"))
    tryCatch({
      if (is.function(plot_obj)) {
        ## Base R graphics plot (e.g., plotClusterHeatmap, plotRadar,
        ## plotTanimoto) -- use grDevices::png instead of ggsave
        grDevices::png(path, width = w, height = h, units = "in",
                       res = 150, bg = "white")
        tryCatch(plot_obj(), error = function(e) NULL,
                 finally = grDevices::dev.off())
      } else {
        ## ggplot object
        ggplot2::ggsave(path, plot_obj, width = w, height = h,
                        dpi = 150, bg = "white")
      }
      paths[[name]] <<- path
    }, error = function(e) {
      warning(sprintf("Could not save plot '%s': %s", name, e$message))
    })
  }

  ## BOILED-Egg
  if (all(c("LogP", "TPSA") %in% names(data))) {
    save_plot("boiled_egg", plotBoiledEgg(data))
  }

  ## MW distribution
  if ("MW" %in% names(data)) {
    save_plot("mw", plotMW(data), w = 6, h = 4)
  }

  ## TPSA distribution
  if ("TPSA" %in% names(data)) {
    save_plot("tpsa", plotTPSA(data), w = 6, h = 4)
  }

  ## LogP distribution
  if ("LogP" %in% names(data)) {
    save_plot("logp", plotLogP(data), w = 6, h = 4)
  }

  ## Correlation heatmap
  numeric_props <- intersect(
    c("MW", "LogP", "TPSA", "MR", "#H-bond acceptors",
      "#H-bond donors", "#Rotatable bonds", "#Heavy atoms"),
    names(data)
  )
  if (length(numeric_props) >= 2) {
    save_plot("corr_heatmap", plotCorrHeatmap(data), w = 7, h = 6)
  }

  ## Cluster heatmap with dendrogram
  if (length(numeric_props) >= 2 && nrow(data) >= 3) {
    save_plot("cluster_heatmap", function() plotClusterHeatmap(data), w = 7, h = 7)
  }

  ## Violations summary
  viol_cols <- paste0(c("Lipinski", "Ghose", "Veber", "Egan", "Muegge"), " #violations")
  if (any(viol_cols %in% names(data))) {
    save_plot("violations", plotViolationsSummary(data), w = 7, h = 5)
  }

  ## Drug-likeness score distribution
  if (any(viol_cols %in% names(data))) {
    save_plot("score", plotDruglikenessScore(data), w = 7, h = 4)
  }

  paths
}

## ===================== Markdown report builder ===========================

#' Build the full Markdown report content
#'
#' @param datasets Named list of dataset info lists, each containing:
#'   \code{raw}, \code{filtered}, \code{filters}, \code{source_name}.
#' @param plot_paths Named list (per dataset) of plot file path lists.
#' @return A character scalar with the full Markdown document.
#' @keywords internal
buildReportMarkdown <- function(datasets, plot_paths) {

  lines <- character(0)
  today <- format(Sys.Date(), "%B %d, %Y")

  ## ---- Header ----
  lines <- c(lines,
    "---",
    "title: \"ADMETShiny Analysis Report\"",
    paste0("date: \"", today, "\""),
    "output: html_document",
    "---",
    "",
    ""
  )

  ## ---- Executive Summary ----
  n_datasets <- length(datasets)
  lines <- c(lines,
    "## 1. Executive Summary",
    "",
    paste0("This report was generated by **ADMETShiny v", ADMETSHINY_VERSION,
           " (", ADMETSHINY_CODENAME, ")** on ", today, ". ",
           "It summarizes the ADMET and drug-likeness analysis performed across ",
           n_datasets, " dataset(s)."),
    ""
  )

  summary_rows <- .safe_rbind(lapply(names(datasets), function(key) {
    ds <- datasets[[key]]
    stats <- computeDatasetStats(ds$raw, ds$filtered, ds$filters, ds$source_name)
    if (is.null(stats)) return(NULL)
    data.frame(
      Dataset = ds$source_name,
      `Molecules Uploaded` = stats$n_raw,
      `Molecules Filtered` = stats$n_filtered,
      `Filters Applied` = if (length(stats$filters_applied) > 0)
        paste(stats$filters_applied, collapse = ", ") else "None",
      `Mean Score` = stats$score$mean,
      check.names = FALSE, stringsAsFactors = FALSE, row.names = NULL
    )
  }))
  lines <- c(lines, .md_table(summary_rows), "")

  ## ---- Per-dataset sections ----
  for (i in seq_along(datasets)) {
    key <- names(datasets)[i]
    ds <- datasets[[key]]
    stats <- computeDatasetStats(ds$raw, ds$filtered, ds$filters, ds$source_name)
    if (is.null(stats)) next

    lines <- c(lines,
      "",
      paste0("## ", i + 1, ". ", ds$source_name),
      ""
    )

    ## Overview
    lines <- c(lines,
      "### Overview",
      "",
      paste0("- **Molecules uploaded:** ", stats$n_raw),
      paste0("- **Molecules after filtering:** ", stats$n_filtered,
             " (", stats$pct_filtered, "% of uploaded)"),
      paste0("- **Filters applied:** ",
             if (length(stats$filters_applied) > 0)
               paste(stats$filters_applied, collapse = ", ")
             else "None"),
      ""
    )

    ## Summary statistics
    lines <- c(lines,
      "### Summary Statistics of Key Properties",
      "",
      .md_table(stats$property_stats),
      ""
    )

    ## Drug-likeness filter results
    lines <- c(lines,
      "### Drug-likeness Filter Results",
      "",
      .md_table(stats$filter_stats),
      ""
    )

    ## ADMET classification
    lines <- c(lines,
      "### BOILED-Egg ADMET Classification",
      "",
      .md_table(stats$admet_stats),
      ""
    )

    ## Additional metrics
    if (!is.null(stats$additional)) {
      lines <- c(lines,
        "### Additional Drug-likeness Metrics",
        "",
        .md_table(stats$additional),
        ""
      )
    }

    ## Composite score
    lines <- c(lines,
      "### Composite Drug-likeness Score",
      "",
      "The composite score (0-100) represents the proportion of the five",
      "drug-likeness rules (Lipinski, Veber, Ghose, Egan, Muegge) that each",
      "compound satisfies, scaled to 0-100. A compound passing all five rules",
      "scores 100.",
      "",
      paste0("- **Mean score:** ", stats$score$mean, " / 100"),
      paste0("- **Median score:** ", stats$score$median, " / 100"),
      paste0("- **Standard deviation:** ", stats$score$sd),
      paste0("- **Excellent (score >= 80):** ", stats$score$excellent, "% of compounds"),
      paste0("- **Acceptable (score 60-79):** ", stats$score$acceptable, "% of compounds"),
      paste0("- **Poor (score < 60):** ", stats$score$poor, "% of compounds"),
      ""
    )

    ## Plots
    pp <- plot_paths[[key]]
    if (!is.null(pp) && length(pp) > 0) {
      lines <- c(lines, "### Visualizations", "")
      if (!is.null(pp$boiled_egg))
        lines <- c(lines, paste0("![BOILED-Egg (", ds$source_name, ")](", pp$boiled_egg, ")"), "")
      if (!is.null(pp$mw))
        lines <- c(lines, paste0("![Molecular Weight Distribution (", ds$source_name, ")](", pp$mw, ")"), "")
      if (!is.null(pp$tpsa))
        lines <- c(lines, paste0("![TPSA Distribution (", ds$source_name, ")](", pp$tpsa, ")"), "")
      if (!is.null(pp$logp))
        lines <- c(lines, paste0("![LogP Distribution (", ds$source_name, ")](", pp$logp, ")"), "")
      if (!is.null(pp$corr_heatmap))
        lines <- c(lines, paste0("![Correlation Heatmap (", ds$source_name, ")](", pp$corr_heatmap, ")"), "")
      if (!is.null(pp$cluster_heatmap))
        lines <- c(lines, paste0("![Cluster Heatmap with Dendrogram (", ds$source_name, ")](", pp$cluster_heatmap, ")"), "")
      if (!is.null(pp$violations))
        lines <- c(lines, paste0("![Violations Summary (", ds$source_name, ")](", pp$violations, ")"), "")
      if (!is.null(pp$score))
        lines <- c(lines, paste0("![Composite Score Distribution (", ds$source_name, ")](", pp$score, ")"), "")
    }
  }

  ## ---- Cross-dataset comparison ----
  if (n_datasets > 1) {
    lines <- c(lines,
      "",
      paste0("## ", n_datasets + 2, ". Cross-Dataset Comparison"),
      ""
    )

    ## Comparison table
    comp_rows <- .safe_rbind(lapply(names(datasets), function(key) {
      ds <- datasets[[key]]
      stats <- computeDatasetStats(ds$raw, ds$filtered, ds$filters, ds$source_name)
      if (is.null(stats)) return(NULL)
      data.frame(
        Dataset = ds$source_name,
        Uploaded = stats$n_raw,
        Filtered = stats$n_filtered,
        `% Pass` = stats$pct_filtered,
        `Filters` = length(stats$filters_applied),
        `Mean Score` = stats$score$mean,
        `Median Score` = stats$score$median,
        check.names = FALSE, stringsAsFactors = FALSE, row.names = NULL
      )
    }))
    lines <- c(lines, .md_table(comp_rows), "")

    ## Common filters
    all_filters <- unique(unlist(lapply(names(datasets), function(key) {
      datasets[[key]]$filters
    })))
    common_filters <- all_filters[sapply(all_filters, function(f) {
      all(sapply(names(datasets), function(key) f %in% datasets[[key]]$filters))
    })]
    if (length(common_filters) > 0) {
      lines <- c(lines,
        "**Filters applied to all datasets:**",
        "",
        paste(common_filters, collapse = ", "),
        ""
      )
    } else {
      lines <- c(lines,
        "**No common filters** were applied across all datasets.",
        ""
      )
    }
  }

  ## ---- References ----
  lines <- c(lines,
    "",
    paste0("## ", n_datasets + 3, ". References"),
    "",
    "1. Lipinski, C. A., Lombardo, F., Dominy, B. W., & Feeney, P. J. (1997).",
    "   *Advanced Drug Delivery Reviews*, 23(1-3), 3-25.",
    "2. Ghose, A. K., Viswanadhan, V. N., & Wendoloski, J. J. (1999).",
    "   *Journal of Combinatorial Chemistry*, 1(1), 55-68.",
    "3. Veber, D. F., et al. (2002). *Journal of Medicinal Chemistry*, 45(12), 2615-2623.",
    "4. Egan, W. J., Merz, K. M., & Baldwin, J. J. (2000).",
    "   *Journal of Medicinal Chemistry*, 43(21), 3867-3877.",
    "5. Muegge, I., Heald, S. L., & Brittelli, D. (2001).",
    "   *Journal of Medicinal Chemistry*, 44(12), 1841-1846.",
    "6. Daina, A., & Zoete, V. (2016). *ChemMedChem*, 11(11), 1117-1121.",
    "7. Seelig, A. (1998). *European Journal of Biochemistry*, 251(1-2), 252-261.",
    "8. Hughes, J. D., et al. (2008). *Bioorganic & Medicinal Chemistry Letters*,",
    "   18(17), 4872-4875.",
    "9. Gleeson, M. P., et al. (2011). *Journal of Medicinal Chemistry*,",
    "   54(13), 4459-4468.",
    "10. Teague, S. J., et al. (1999). *Angewandte Chemie International Edition*,",
    "    38(24), 3743-3748.",
    ""
  )

  paste(lines, collapse = "\n")
}

## ===================== Report renderer ===================================

#' Render the ADMETShiny report to HTML, PDF, or Word
#'
#' Generates a comprehensive report from one or more datasets and renders it
#' to the specified format using \pkg{rmarkdown}. Plots are generated as
#' temporary PNG files and embedded in the document.
#'
#' @param datasets Named list of dataset info lists, each containing:
#'   \code{raw} (data.frame or NULL), \code{filtered} (data.frame or NULL),
#'   \code{filters} (character vector), \code{source_name} (character).
#' @param format Character; one of \code{"html"}, \code{"pdf"}, \code{"doc"}.
#' @param output_file Character; path where the output file will be written.
#' @return Invisible \code{NULL}; called for the side-effect of writing the
#'   rendered report to \code{output_file}.
#' @keywords internal
renderReport <- function(datasets, format = "html", output_file) {

  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    stop("The 'rmarkdown' package is required to render reports. ",
         "Install it with: install.packages('rmarkdown')",
         call. = FALSE)
  }

  ## Filter out empty datasets
  datasets <- datasets[sapply(datasets, function(ds) !is.null(ds$raw))]

  if (length(datasets) == 0) {
    stop("No datasets loaded. Please upload data in at least one module.",
         call. = FALSE)
  }

  ## Use a persistent temp directory (not cleaned until after render) for
  ## both the Markdown source and the plot PNGs, so that rmarkdown can find
  ## the images during rendering.
  work_dir <- tempfile("admetshiny_report")
  dir.create(work_dir, recursive = TRUE)
  on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)

  plot_dir <- file.path(work_dir, "plots")
  dir.create(plot_dir, recursive = TRUE)

  ## Generate plots for each dataset
  plot_paths <- list()
  for (key in names(datasets)) {
    ds <- datasets[[key]]
    filtered <- if (!is.null(ds$filtered)) ds$filtered else ds$raw
    if (!is.null(filtered) && nrow(filtered) > 0) {
      plot_paths[[key]] <- generateReportPlots(filtered, key, plot_dir)
    }
  }

  ## Build Markdown content
  md_content <- buildReportMarkdown(datasets, plot_paths)

  ## Write to .md file in the work directory (same dir as the plots)
  md_file <- file.path(work_dir, "report.md")
  writeLines(md_content, md_file)

  ## Determine output format
  output_format <- switch(format,
    pdf  = rmarkdown::pdf_document(fig_caption = FALSE, latex_engine = "pdflatex"),
    doc  = rmarkdown::word_document(),
    html = rmarkdown::html_document(self_contained = TRUE),
    rmarkdown::html_document(self_contained = TRUE)
  )

  ## Render
  rmarkdown::render(
    input = md_file,
    output_format = output_format,
    output_file = basename(output_file),
    output_dir = dirname(output_file),
    quiet = TRUE,
    clean = TRUE
  )

  invisible(NULL)
}

## ===================== Console-friendly report ===========================

#' Generate an ADMETShiny report from the R console
#'
#' Generates a comprehensive ADMET and drug-likeness analysis report from one
#' or more datasets, directly from the R console (without launching the Shiny
#' app). The report includes per-dataset statistics, drug-likeness filter
#' results, BOILED-Egg ADMET classification, additional literature-supported
#' metrics, a composite drug-likeness score, cross-dataset comparison,
#' visualizations and references.
#'
#' @param data A data.frame, or a named list of data.frames. If a single
#'   data.frame is provided, it is treated as one dataset named
#'   \code{"Dataset 1"}. If a list is provided, each element is a separate
#'   dataset and the list names are used as dataset labels.
#' @param filters Character vector of drug-likeness filter names applied to
#'   the data (e.g. \code{c("Lipinski", "Veber")}). Use \code{character(0)}
#'   if no filters were applied. Default \code{character(0)}.
#' @param format Character; one of \code{"html"}, \code{"pdf"}, \code{"doc"}.
#'   Default \code{"html"}.
#' @param output_file Character; path where the output file will be written.
#'   If \code{NULL}, a default name is used in the current working directory.
#' @param source_name Character; human-readable name for the dataset. Only
#'   used when \code{data} is a single data.frame. Default \code{"Dataset"}.
#' @return Invisible \code{NULL}; called for the side-effect of writing the
#'   rendered report to \code{output_file}.
#' @export
#' @examples
#' \dontrun{
#' ## Generate a report from a single dataset
#' d <- read.csv("swissadme.csv", check.names = FALSE)
#' d <- fixSwissADME(d)
#' generateReport(d, filters = c("Lipinski", "Veber"), format = "html")
#'
#' ## Generate a report from multiple datasets
#' generateReport(
#'   data = list(SwissADME = d1, ADMETlab = d2),
#'   filters = list(SwissADME = c("Lipinski"), ADMETlab = c("Egan")),
#'   format = "pdf"
#' )
#' }
generateReport <- function(data, filters = character(0),
                           format = "html", output_file = NULL,
                           source_name = "Dataset") {

  ## Build the datasets list expected by renderReport
  if (is.data.frame(data)) {

    ## Single dataset
    datasets <- list(
      Dataset1 = list(
        raw         = data,
        filtered    = data,
        filters     = filters,
        source_name = source_name
      )
    )

  } else if (is.list(data) && !is.null(names(data))) {

    ## Named list of datasets
    datasets <- lapply(names(data), function(nm) {
      d <- data[[nm]]
      if (is.null(d) || !is.data.frame(d)) return(NULL)
      ## Get filters for this dataset if filters is a named list
      fltrs <- if (is.list(filters) && nm %in% names(filters)) {
        filters[[nm]]
      } else if (is.character(filters)) {
        filters
      } else {
        character(0)
      }
      list(
        raw         = d,
        filtered    = d,
        filters     = fltrs,
        source_name = nm
      )
    })
    names(datasets) <- names(data)
    datasets <- datasets[!sapply(datasets, is.null)]

  } else {
    stop("'data' must be a data.frame or a named list of data.frames.",
         call. = FALSE)
  }

  ## Default output file
  if (is.null(output_file)) {
    ext <- switch(format, pdf = "pdf", doc = "docx", "html")
    output_file <- file.path(getwd(),
                             paste0("ADMETShiny_Report_",
                                    format(Sys.Date(), "%Y%m%d"), ".", ext))
  }

  ## Render
  renderReport(datasets, format = format, output_file = output_file)

  message("Report generated: ", output_file)
  invisible(NULL)
}
