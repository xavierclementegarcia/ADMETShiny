# ---------------------------------------------------------------------------
# data_fixers.R
# Normalization of datasets exported from SwissADME, ADMETlab 3.0 and
# Deep-PK to the standard column schema used by the application.
#
# The standard schema uses a single generic `LogP` column for lipophilicity,
# regardless of the source. SwissADME-specific LogP variants (MLOGP, WLOGP,
# XLOGP3) are NOT created for non-SwissADME sources. When a source does not
# provide a required property (e.g. MR, #Heavy atoms, #Aromatic heavy atoms in
# ADMETlab 3.0), it is calculated from the SMILES using CDK.
# ---------------------------------------------------------------------------

#' Normalize a SwissADME dataset
#'
#' Cleans and standardizes a data.frame exported from SwissADME: coerces
#' the known numeric and character columns to the correct type, converts empty
#' strings to \code{NA}, and creates a generic \code{LogP} column from
#' \code{Consensus Log P} (or \code{WLOGP} as fallback) so that the
#' drug-likeness filters and plots work with a single lipophilicity column
#' across all data sources.
#'
#' @param data A data.frame read with \code{read.csv(..., check.names = FALSE)}.
#' @return A data.frame with numeric and categorical columns correctly typed,
#'   plus a \code{LogP} column.
#' @export
#' @seealso \code{\link{fixADMETlab}}, \code{\link{fixDeepPK}},
#'   \code{\link{applyFilters}}.
fixSwissADME <- function(data) {

  stopifnot(is.data.frame(data))

  numeric_columns <- c(
    "MW",
    "#Heavy atoms",
    "#Aromatic heavy atoms",
    "Fraction Csp3",
    "#Rotatable bonds",
    "#H-bond acceptors",
    "#H-bond donors",
    "MR",
    "TPSA",
    "iLOGP",
    "XLOGP3",
    "WLOGP",
    "MLOGP",
    "Silicos-IT Log P",
    "Consensus Log P",
    "ESOL Log S",
    "ESOL Solubility (mg/ml)",
    "ESOL Solubility (mol/l)",
    "Ali Log S",
    "Ali Solubility (mg/ml)",
    "Ali Solubility (mol/l)",
    "Silicos-IT LogSw",
    "Silicos-IT Solubility (mg/ml)",
    "Silicos-IT Solubility (mol/l)",
    "Log Kp (cm/s)",
    "Lipinski #violations",
    "Ghose #violations",
    "Veber #violations",
    "Egan #violations",
    "Muegge #violations",
    "Bioavailability Score",
    "PAINS #alerts",
    "Brenk #alerts",
    "Leadlikeness #violations",
    "Synthetic Accessibility"
  )

  character_columns <- c(
    "ESOL Class",
    "Ali Class",
    "Silicos-IT class",
    "GI absorption",
    "BBB permeant",
    "Pgp substrate",
    "CYP1A2 inhibitor",
    "CYP2C19 inhibitor",
    "CYP2C9 inhibitor",
    "CYP2D6 inhibitor",
    "CYP3A4 inhibitor"
  )

  numeric_columns <- intersect(numeric_columns, names(data))
  character_columns <- intersect(character_columns, names(data))

  data[numeric_columns] <- lapply(
    data[numeric_columns],
    function(x) {
      x <- trimws(x)
      x[x == ""] <- NA
      as.numeric(x)
    }
  )

  data[character_columns] <- lapply(data[character_columns], as.character)

  ## Create the generic LogP column from Consensus Log P (preferred) or
  ## WLOGP as fallback. This is the only LogP column the filters and plots
  ## require. The original SwissADME-specific MLOGP/WLOGP/XLOGP3 columns
  ## remain in the data.frame for reference but are not used by the
  ## drug-likeness filters.
  if ("Consensus Log P" %in% names(data)) {
    data$LogP <- data[["Consensus Log P"]]
  } else if ("WLOGP" %in% names(data)) {
    data$LogP <- data$WLOGP
  } else if ("MLOGP" %in% names(data)) {
    data$LogP <- data$MLOGP
  }

  ## Recalculate the #violations columns using the application's standard
  ## thresholds (LogP <= 5 for Lipinski, no aromatic heavy atoms for Muegge,
  ## etc.). SwissADME pre-computes these with its own thresholds (MLOGP <= 4.15
  ## for Lipinski, aromatic heavy atoms <= 7 for Muegge) which are
  ## inconsistent with the filters. Overwriting them ensures consistency
  ## between the #violations column and the individual property thresholds.
  data <- computeViolationColumns(data)

  required_base <- c("MW", "LogP", "#H-bond acceptors", "#H-bond donors")
  missing_cols <- setdiff(required_base, names(data))
  if (length(missing_cols) > 0) {
    warning("Missing expected columns: ",
            paste(missing_cols, collapse = ", "))
  }

  data
}

## ---------------------------------------------------------------------------

#' Normalize an ADMETlab 3.0 dataset
#'
#' Cleans and maps a data.frame exported from ADMETlab 3.0 to the application's
#' standard column schema: drops the SVG column, renames descriptor columns,
#' converts ADMET probabilities into categorical labels, and computes
#' drug-likeness violation columns.
#'
#' ADMETlab 3.0 does not natively provide Molar Refractivity (MR), the number
#' of heavy atoms, or the number of aromatic heavy atoms, which are required
#' by the Ghose and Muegge filters. These properties are calculated on-the-fly
#' from the SMILES using CDK (\code{\link{calcCDKDescriptors}}). Requires the
#' suggested package \pkg{rcdk} and a working Java JDK.
#'
#' The source's \code{logP} column is mapped to the generic \code{LogP}
#' column. No artificial \code{MLOGP}/\code{WLOGP}/\code{XLOGP3} columns are
#' created.
#'
#' @param data A data.frame read from an ADMETlab 3.0 CSV file.
#' @return A data.frame ready for filtering and plotting.
#' @export
#' @seealso \code{\link{fixSwissADME}}, \code{\link{fixDeepPK}},
#'   \code{\link{computeViolationColumns}}, \code{\link{calcCDKDescriptors}}.
fixADMETlab <- function(data) {

  ## 1. Remove the molstr (SVG) column if it exists
  if ("molstr" %in% names(data)) {
    data$molstr <- NULL
  }

  ## 2. Remove the first index column if it exists
  first_col_name <- names(data)[1]
  if (first_col_name == "0" || first_col_name == "") {
    data <- data[, -1, drop = FALSE]
  }

  ## 3. Rename columns from ADMETlab names to the application's standard names
  rename_map <- c(
    "smiles" = "CanonicalSMILES",
    "MW"     = "MW",
    "TPSA"   = "TPSA",
    "nHA"    = "#H-bond acceptors",
    "nHD"    = "#H-bond donors",
    "nRot"   = "#Rotatable bonds",
    "logP"   = "LogP"
  )

  for (old in names(rename_map)) {
    if (old %in% names(data)) {
      names(data)[names(data) == old] <- rename_map[[old]]
    }
  }

  ## 4. Convert key columns to numeric
  numeric_cols <- c("MW", "TPSA", "#H-bond acceptors", "#H-bond donors",
                    "#Rotatable bonds", "LogP", "logS", "logD",
                    "hia", "BBB", "pgp_sub")
  for (col in numeric_cols) {
    if (col %in% names(data)) {
      data[[col]] <- suppressWarnings(as.numeric(data[[col]]))
    }
  }

  ## 5. Calculate missing properties (MR, #Heavy atoms, #Aromatic heavy atoms)
  ##    from SMILES using CDK. ADMETlab 3.0 does not provide these natively,
  ##    but they are required by the Ghose and Muegge filters.
  missing_props <- setdiff(
    c("MR", "#Heavy atoms", "#Aromatic heavy atoms"),
    names(data)
  )

  if (length(missing_props) > 0 && "CanonicalSMILES" %in% names(data)) {
    if (!requireNamespace("rcdk", quietly = TRUE)) {
      stop("The 'rcdk' package is required to calculate MR, #Heavy atoms ",
           "and #Aromatic heavy atoms from SMILES for ADMETlab 3.0 data. ",
           "Install it with: install.packages('rcdk')",
           call. = FALSE)
    }

    smiles <- as.character(data$CanonicalSMILES)
    cdk_desc <- calcCDKDescriptors(smiles, which = c("mr", "heavy", "aroma"))

    ## Map CDK names to standard names
    cdk_rename <- c(
      AMR        = "MR",
      nAtom      = "#Heavy atoms",
      naAromAtom = "#Aromatic heavy atoms"
    )
    for (old in names(cdk_rename)) {
      if (old %in% names(cdk_desc)) {
        names(cdk_desc)[names(cdk_desc) == old] <- cdk_rename[[old]]
      }
    }

    ## Keep only the SMILES and the missing property columns
    keep <- intersect(c("SMILES", missing_props), names(cdk_desc))
    if (length(keep) > 1) {
      cdk_desc <- cdk_desc[, keep, drop = FALSE]
      data <- merge(data, cdk_desc,
                    by.x = "CanonicalSMILES", by.y = "SMILES",
                    all.x = TRUE)
    }
  }

  ## 6. Compute ADMET categorical properties from ADMETlab probabilities
  data <- computeADMETlabProperties(data)

  ## 7. Compute drug-likeness violation columns
  computeViolationColumns(data)
}

#' Compute ADMET categorical properties from ADMETlab probabilities
#'
#' Converts the numeric ADMETlab probabilities (0 to 1) into the categorical
#' labels used by the application (\code{"High"}/\code{"Low"} for GI
#' absorption, \code{"Yes"}/\code{"No"} for BBB permeability and P-gp
#' substrate) using a 0.5 threshold.
#'
#' @param data A data.frame with ADMETlab probability columns.
#' @return A data.frame with the added categorical columns.
#' @keywords internal
computeADMETlabProperties <- function(data) {

  if ("hia" %in% names(data)) {
    data[["GI absorption"]] <- ifelse(
      is.na(data$hia), NA,
      ifelse(data$hia >= 0.5, "High", "Low")
    )
  }

  if ("BBB" %in% names(data)) {
    data[["BBB permeant"]] <- ifelse(
      is.na(data$BBB), NA,
      ifelse(data$BBB >= 0.5, "Yes", "No")
    )
  }

  if ("pgp_sub" %in% names(data)) {
    data[["Pgp substrate"]] <- ifelse(
      is.na(data$pgp_sub), NA,
      ifelse(data$pgp_sub >= 0.5, "Yes", "No")
    )
  }

  data
}

## ---------------------------------------------------------------------------

#' Normalize a Deep-PK dataset
#'
#' Cleans and maps a data.frame exported from Deep-PK Learning to the
#' application's standard column schema. Because Deep-PK CSVs do not contain
#' physicochemical descriptors natively, all descriptors (MW, TPSA, LogP, MR,
#' HBA, HBD, rotatable bonds, heavy atoms, aromatic heavy atoms) are
#' calculated on-the-fly from the SMILES using CDK
#' (\code{\link{calcCDKDescriptors}}). If the Deep-PK CSV contains its own
#' LogP prediction, that value takes precedence over the CDK-calculated ALogP.
#' Requires the suggested package \pkg{rcdk} and a working Java JDK.
#'
#' @param data A data.frame read from a Deep-PK CSV file.
#' @return A data.frame ready for filtering and plotting.
#' @export
#' @seealso \code{\link{fixSwissADME}}, \code{\link{fixADMETlab}},
#'   \code{\link{calcCDKDescriptors}}, \code{\link{mapCDKDescriptors}}.
fixDeepPK <- function(data) {

  ## 1. Remove the first index column if it exists
  first_col_name <- names(data)[1]
  if (first_col_name == "0" || first_col_name == "") {
    data <- data[, -1, drop = FALSE]
  }

  ## 2. Rename SMILES -> CanonicalSMILES
  if ("SMILES" %in% names(data)) {
    names(data)[names(data) == "SMILES"] <- "CanonicalSMILES"
  }

  ## 3. Extract LogP from Deep-PK's own prediction column (if present)
  cols <- names(data)
  logp_col <- cols[grepl("Log\\(P\\)\\] Predictions", cols, ignore.case = TRUE)]
  has_deepkp_logp <- length(logp_col) > 0
  if (has_deepkp_logp) {
    data$LogP <- suppressWarnings(as.numeric(data[[logp_col[1]]]))
  }

  ## 4. Calculate ALL physicochemical descriptors from SMILES using CDK
  if (!requireNamespace("rcdk", quietly = TRUE)) {
    stop("The 'rcdk' package is required to process Deep-PK datasets.",
         call. = FALSE)
  }

  smiles <- as.character(data$CanonicalSMILES)
  cdk_desc <- calcCDKDescriptors(
    smiles,
    which = c("mw", "alogp", "tpsa", "hbd", "hba", "rotb", "heavy", "aroma", "mr")
  )
  mapped <- mapCDKDescriptors(cdk_desc)

  ## Prepare CDK results for merge
  mapped$CanonicalSMILES <- mapped$SMILES
  mapped$SMILES <- NULL

  ## If Deep-PK provided its own LogP, drop CDK's LogP from mapped (Deep-PK
  ## value takes precedence). Otherwise, keep CDK's LogP (ALogP).
  if (has_deepkp_logp) {
    mapped$LogP <- NULL
  }

  ## Keep only the standard descriptor columns
  keep_cols <- c("CanonicalSMILES", "MW", "LogP", "TPSA", "#H-bond donors",
                 "#H-bond acceptors", "#Rotatable bonds", "#Heavy atoms",
                 "#Aromatic heavy atoms", "MR")
  mapped <- mapped[, intersect(keep_cols, names(mapped)), drop = FALSE]

  data <- merge(data, mapped, by = "CanonicalSMILES", all.x = TRUE)

  ## 5. Map Deep-PK ADMET predictions to categorical labels
  clean_text <- function(x) gsub("<br/>|&nbsp;|<.*?>", " ", x)

  hia_col <- cols[grepl("Human Intestinal Absorption\\] Interpretation",
                        cols, ignore.case = TRUE)]
  if (length(hia_col) > 0) {
    hia <- clean_text(data[[hia_col[1]]])
    ## Deep-PK returns "Absorbed" or "Not absorbed".
    ## grepl("Absorbed", "Not absorbed") returns TRUE (substring match),
    ## so we must exclude "Not" explicitly.
    data[["GI absorption"]] <- ifelse(
      grepl("Absorbed", hia, ignore.case = TRUE) &
        !grepl("Not", hia, ignore.case = TRUE),
      "High", "Low"
    )
  }

  bbb_col <- cols[grepl("Blood-Brain Barrier\\] Interpretation",
                        cols, ignore.case = TRUE)]
  if (length(bbb_col) > 0) {
    bbb <- clean_text(data[[bbb_col[1]]])
    ## Deep-PK returns "Penetrable" or "Non-penetrable".
    ## grepl("Penetrable", "Non-penetrable") returns TRUE (substring match),
    ## so we must exclude "Non" explicitly.
    data[["BBB permeant"]] <- ifelse(
      grepl("Penetrable", bbb, ignore.case = TRUE) &
        !grepl("Non", bbb, ignore.case = TRUE),
      "Yes", "No"
    )
  }

  pgp_col <- cols[grepl("P-Glycoprotein Substrate\\] Interpretation",
                        cols, ignore.case = TRUE)]
  if (length(pgp_col) > 0) {
    pgp <- clean_text(data[[pgp_col[1]]])
    ## Deep-PK returns "Substrate" or "Non-substrate".
    ## grepl("Substrate", "Non-substrate") returns TRUE (substring match),
    ## so we must exclude "Non" explicitly.
    data[["Pgp substrate"]] <- ifelse(
      grepl("Substrate", pgp, ignore.case = TRUE) &
        !grepl("Non", pgp, ignore.case = TRUE),
      "Yes", "No"
    )
  }

  ## 6. Compute drug-likeness violation columns
  computeViolationColumns(data)
}
