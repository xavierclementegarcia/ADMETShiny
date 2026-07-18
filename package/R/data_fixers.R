# ---------------------------------------------------------------------------
# data_fixers.R
# Column mapping and dataset normalization for the ADMET Master Manager.
# The generic mapADMETColumns() function replaces the former platform-specific
# (replaces the former fixSwissADME/fixADMETlab/fixDeepPK functions).
# ---------------------------------------------------------------------------

## ---------------------------------------------------------------------------
## ADMET Master Manager helpers
## ---------------------------------------------------------------------------

#' Detect column types in a data.frame
#'
#' Returns a compact summary of every column in a data.frame, with its
#' inferred type (numeric / string), the number of unique non-NA values and
#' a small sample of the first few values. Used by the ADMET Master Manager
#' to help the user pick the right column mapping.
#'
#' @param data A data.frame.
#' @return A data.frame with columns: \code{column_name}, \code{detected_type},
#'   \code{n_unique}, \code{sample_values}.
#' @keywords internal
#' @seealso \code{\link{detectSMILESColumn}}, \code{\link{mapADMETColumns}}.
detectColumnTypes <- function(data) {

  stopifnot(is.data.frame(data))

  if (ncol(data) == 0) {
    return(data.frame(
      column_name    = character(0),
      detected_type  = character(0),
      n_unique       = integer(0),
      sample_values  = character(0),
      stringsAsFactors = FALSE
    ))
  }

  res <- data.frame(
    column_name   = names(data),
    detected_type = vapply(data, function(col) {
      ## Try to coerce to numeric: if more than half of the non-NA values
      ## parse successfully, treat the column as numeric.
      if (is.numeric(col)) return("numeric")
      x <- suppressWarnings(as.numeric(as.character(col)))
      non_na <- sum(!is.na(x))
      total  <- sum(!is.na(col))
      if (total > 0 && non_na / total >= 0.5) return("numeric")
      "string"
    }, character(1)),
    n_unique = vapply(data, function(col) {
      length(unique(as.character(col)[!is.na(as.character(col))]))
    }, integer(1)),
    sample_values = vapply(data, function(col) {
      vals <- as.character(col)[!is.na(as.character(col))]
      if (length(vals) == 0) return("")
      samp <- head(vals, 3)
      paste(samp, collapse = ", ")
    }, character(1)),
    stringsAsFactors = FALSE
  )

  res
}

## ---------------------------------------------------------------------------

#' Auto-detect SMILES column in a data.frame
#'
#' Tries to identify the column that contains SMILES strings. It first looks
#' for column names matching common conventions (\code{smiles},
#' \code{SMILES}, \code{CanonicalSMILES}, \code{canonical_smiles},
#' \code{IsomericSMILES}, etc.). If no name matches, it inspects the values
#' of every string column and picks the one whose values look most like
#' SMILES (contain carbon / aromatic / bracket / bond characters).
#'
#' @param data A data.frame.
#' @return Character name of the detected SMILES column, or \code{NULL} if
#'   none found.
#' @keywords internal
#' @seealso \code{\link{detectColumnTypes}}, \code{\link{mapADMETColumns}}.
detectSMILESColumn <- function(data) {

  stopifnot(is.data.frame(data))
  if (ncol(data) == 0) return(NULL)

  ## 1. Look for SMILES-like column names first.
  smiles_names <- c("smiles", "SMILES", "CanonicalSMILES", "canonical_smiles",
                    "canonical_smiles_can", "IsomericSMILES", "isomeric_smiles",
                    "MOL", "mol", "molstr", "structure", "Structure",
                    "canonical_smiles ")
  hit <- intersect(smiles_names, names(data))
  if (length(hit) > 0) return(hit[1])

  ## Fuzzy match: any column name containing "smiles" (case-insensitive).
  guess <- names(data)[grepl("smiles", names(data), ignore.case = TRUE)]
  if (length(guess) > 0) return(guess[1])

  ## 2. Inspect values: pick the string column whose values look most like
  ##    SMILES. A value "looks like SMILES" if it contains at least one
  ##    carbon-letter (C/c) and at least one bond / bracket / ring digit
  ##    character.
  score_col <- function(col) {
    vals <- as.character(col)[!is.na(as.character(col))]
    if (length(vals) == 0) return(0)
    vals <- head(vals, 100)
    looks_smiles <- vapply(vals, function(v) {
      has_carbon <- grepl("[Cc]", v)
      has_bond_or_bracket <-
        grepl("[\\[\\]\\(\\)=#\\\\/\\d]", v) ||
        grepl("c[1-9]", v) ||
        grepl("C[1-9]", v)
      has_carbon && has_bond_or_bracket
    }, logical(1))
    mean(looks_smiles, na.rm = TRUE)
  }

  scores <- vapply(data, score_col, numeric(1))
  best <- which.max(scores)
  if (length(best) == 0 || scores[best] < 0.3) return(NULL)
  names(data)[best]
}

## ---------------------------------------------------------------------------

#' Map user columns to standard ADMET schema
#'
#' Takes a raw data.frame and a user-specified column mapping, renames
#' columns to the application's standard schema, converts types, optionally
#' calculates missing descriptors with CDK, and computes violation columns
#' and ADMET properties.
#'
#' The \code{mapping} argument is a named character vector where each name
#' is the name of a column in \code{data} and each value is one of the
#' standard short field codes used by the ADMET Master Manager:
#' \code{"None"}, \code{"SMILES"}, \code{"Name"}, \code{"MW"},
#' \code{"LogP"}, \code{"WLOGP"}, \code{"TPSA"}, \code{"HBD"},
#' \code{"HBA"}, \code{"Rotatable Bonds"}, \code{"Molar Refractivity"},
#' \code{"Heavy Atoms"}, \code{"Aromatic Heavy Atoms"},
#' \code{"GI Absorption"}, \code{"GI Absorption_num"},
#' \code{"BBB Permeant"}, \code{"BBB Permeant_num"},
#' \code{"Pgp Substrate"}, \code{"Pgp Substrate_num"},
#' \code{"LogS"}, \code{"LogD"}.
#'
#' @param data A data.frame as uploaded by the user (CSV or Excel).
#' @param mapping A named character vector where names are user column names
#'   and values are standard field names. Use "None" for columns to skip.
#'   Example: c("mol_weight" = "MW", "logp" = "LogP", "smiles" = "SMILES")
#' @param calculate_cdk Logical. If TRUE and SMILES column is available,
#'   calculate missing descriptors (MW, LogP, TPSA, HBD, HBA, RB, MR,
#'   HeavyAtoms, AromAtoms) using CDK. Default TRUE.
#' @return A data.frame with standardized column names, violation columns,
#'   and ADMET properties.
#' @export
#' @seealso \code{\link{computeViolationColumns}}, \code{\link{computeADMETProperties}},
#'   \code{\link{calcCDKDescriptors}}.
mapADMETColumns <- function(data, mapping, calculate_cdk = TRUE) {

  stopifnot(is.data.frame(data))
  if (is.null(mapping) || length(mapping) == 0) {
    warning("No column mapping was provided; returning the data as-is.")
    return(data)
  }

  ## ----- 1. Rename columns based on mapping ----------------------------
  ## Translate the short field codes to the application's standard column
  ## names (matching the schema used by mapADMETColumns, mapCDKDescriptors 
  ## mapCDKDescriptors).
  code_to_standard <- c(
    "SMILES"                  = "SMILES",
    "Name"                    = "Name",
    "MW"                      = "MW",
    "LogP"                    = "LogP",
    "WLOGP"                   = "WLOGP",
    "TPSA"                    = "TPSA",
    "HBD"                     = "#H-bond donors",
    "HBA"                     = "#H-bond acceptors",
    "Rotatable Bonds"         = "#Rotatable bonds",
    "Molar Refractivity"      = "MR",
    "Heavy Atoms"             = "#Heavy atoms",
    "Aromatic Heavy Atoms"    = "#Aromatic heavy atoms",
    "GI Absorption"           = "GI absorption",
    "GI Absorption_num"       = "GI_absorption_num",
    "BBB Permeant"            = "BBB permeant",
    "BBB Permeant_num"        = "BBB_num",
    "Pgp Substrate"           = "Pgp substrate",
    "Pgp Substrate_num"       = "Pgp_num",
    "LogS"                    = "LogS",
    "LogD"                    = "LogD"
  )

  for (user_col in names(mapping)) {
    code <- mapping[[user_col]]
    if (is.na(code) || code == "None" || code == "") next
    new_name <- code_to_standard[[code]]
    if (is.na(new_name)) {
      warning("Unknown mapping code '", code, "' for column '", user_col,
              "'; skipping.")
      next
    }
    if (user_col %in% names(data)) {
      ## If the target name already exists and points to a different column,
      ## drop the old column first to avoid duplicate-name confusion.
      if (new_name %in% names(data) && new_name != user_col) {
        data[[new_name]] <- NULL
      }
      names(data)[names(data) == user_col] <- new_name
    }
  }

  ## ----- 2. Convert numeric fields ------------------------------------
  numeric_cols <- intersect(
    c("MW", "LogP", "WLOGP", "TPSA",
      "#H-bond donors", "#H-bond acceptors", "#Rotatable bonds",
      "MR", "#Heavy atoms", "#Aromatic heavy atoms",
      "LogS", "LogD",
      "GI_absorption_num", "BBB_num", "Pgp_num"),
    names(data)
  )
  for (col in numeric_cols) {
    data[[col]] <- suppressWarnings(as.numeric(data[[col]]))
  }

  ## ----- 3. Handle ADMET numeric -> categorical -----------------------
  if ("GI_absorption_num" %in% names(data)) {
    v <- data$GI_absorption_num
    data[["GI absorption"]] <- ifelse(
      is.na(v), NA,
      ifelse(v >= 0.5, "High", "Low")
    )
    data$GI_absorption_num <- NULL
  }
  if ("BBB_num" %in% names(data)) {
    v <- data$BBB_num
    data[["BBB permeant"]] <- ifelse(
      is.na(v), NA,
      ifelse(v >= 0.5, "Yes", "No")
    )
    data$BBB_num <- NULL
  }
  if ("Pgp_num" %in% names(data)) {
    v <- data$Pgp_num
    data[["Pgp substrate"]] <- ifelse(
      is.na(v), NA,
      ifelse(v >= 0.5, "Yes", "No")
    )
    data$Pgp_num <- NULL
  }

  ## ----- 4. Optionally calculate missing descriptors with CDK --------
  if (calculate_cdk && "SMILES" %in% names(data)) {
    missing_desc <- setdiff(
      c("MW", "LogP", "TPSA", "#H-bond donors", "#H-bond acceptors",
        "#Rotatable bonds", "MR", "#Heavy atoms", "#Aromatic heavy atoms"),
      names(data)
    )
    if (length(missing_desc) > 0) {
      if (!requireNamespace("rcdk", quietly = TRUE)) {
        warning("The 'rcdk' package is required to calculate missing ",
                "descriptors with CDK. Install it with: ",
                "install.packages('rcdk').")
      } else {
        tryCatch({
          smiles <- as.character(data$SMILES)
          ## Map missing standard names -> CDK "which" short codes.
          std_to_cdk <- c(
            "MW"                     = "mw",
            "LogP"                   = "alogp",
            "TPSA"                   = "tpsa",
            "#H-bond donors"         = "hbd",
            "#H-bond acceptors"      = "hba",
            "#Rotatable bonds"       = "rotb",
            "MR"                     = "mr",
            "#Heavy atoms"           = "heavy",
            "#Aromatic heavy atoms"  = "aroma"
          )
          cdk_which <- unique(std_to_cdk[missing_desc])
          cdk_desc <- calcCDKDescriptors(smiles, which = cdk_which)

          ## Rename CDK columns to standard names
          cdk_rename <- c(
            MW            = "MW",
            ALogP         = "LogP",
            AMR           = "MR",
            TopoPSA       = "TPSA",
            nHBDon        = "#H-bond donors",
            nHBAcc        = "#H-bond acceptors",
            nRotB         = "#Rotatable bonds",
            nAtom         = "#Heavy atoms",
            naAromAtom    = "#Aromatic heavy atoms"
          )
          for (old in names(cdk_rename)) {
            if (old %in% names(cdk_desc)) {
              names(cdk_desc)[names(cdk_desc) == old] <- cdk_rename[[old]]
            }
          }
          if ("ALogP2" %in% names(cdk_desc)) {
            cdk_desc$ALogP2 <- NULL
          }

          ## Use lookup to avoid row duplication
          cdk_lookup <- split(cdk_desc, cdk_desc$SMILES)
          for (col in missing_desc) {
            if (col %in% names(cdk_desc)) {
              data[[col]] <- sapply(smiles, function(s) {
                if (!is.na(s) && nchar(s) > 0 && s %in% names(cdk_lookup)) {
                  cdk_lookup[[s]][[col]][1]
                } else {
                  NA
                }
              })
            }
          }
        }, error = function(e) {
          warning("Could not calculate CDK descriptors: ", e$message)
        })
      }
    }
  }

  ## ----- 5. If WLOGP present and LogP not mapped, set LogP = WLOGP -----
  ## WLOGP is preferred for BOILED-Egg (official calibration), but only
  ## if the user didn't explicitly map a different column to LogP.
  if ("WLOGP" %in% names(data) && !"LogP" %in% names(data)) {
    data$LogP <- data$WLOGP
  }

  ## ----- 6. Compute drug-likeness violation columns -------------------
  data <- computeViolationColumns(data)

  ## ----- 7. Compute ADMET properties (only if not provided) -----------
  ## Don't overwrite user-provided GI / BBB / Pgp columns.
  need_admet <- (!"GI absorption" %in% names(data)) ||
                (!"BBB permeant"   %in% names(data)) ||
                (!"Pgp substrate"  %in% names(data))
  if (need_admet) {
    ## Temporarily stash any user-provided categorical columns so that
    ## computeADMETProperties doesn't overwrite them.
    stash_gi  <- if ("GI absorption" %in% names(data)) data[["GI absorption"]] else NULL
    stash_bbb <- if ("BBB permeant"  %in% names(data)) data[["BBB permeant"]]  else NULL
    stash_pgp <- if ("Pgp substrate" %in% names(data)) data[["Pgp substrate"]] else NULL

    data <- computeADMETProperties(data)

    if (!is.null(stash_gi))  data[["GI absorption"]] <- stash_gi
    if (!is.null(stash_bbb)) data[["BBB permeant"]]  <- stash_bbb
    if (!is.null(stash_pgp)) data[["Pgp substrate"]] <- stash_pgp
  }

  data
}
