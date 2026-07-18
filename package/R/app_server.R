# ---------------------------------------------------------------------------
# app_server.R
# Server function for the ADMETShiny application.
# ---------------------------------------------------------------------------

#' ADMETShiny application server
#'
#' The server function backing the ADMETShiny Shiny application. It is not
#' meant to be called directly by end users; use \code{\link{run_app}} to
#' launch the app.
#'
#' @param input,output,session Shiny input, output and session objects.
#' @keywords internal
app_server <- function(input, output, session) {

  ## Helper: safely extract id_col from Shiny input, handling NULL,
  ## character(0), NA, "None", and "" -- all of which mean "no label".
  .safe_id_col <- function(val) {
    if (is.null(val) || length(val) == 0 || !is.character(val) ||
        is.na(val[1]) || val[1] == "None" || val[1] == "") {
      return(NULL)
    }
    val[1]
  }

  ## --------------------------- Dark mode ---------------------------------
  is_dark_mode <- reactiveVal(FALSE)

  observeEvent(input$toggle_darkmode, {
    new_state <- !is_dark_mode()
    is_dark_mode(new_state)
    session$sendCustomMessage("toggle-dark-mode", new_state)
    updateActionButton(session, "toggle_darkmode",
                       icon = if (new_state) icon("sun") else icon("moon"))
  })

  ## ------------------------- Navbar navigation ---------------------------
  observeEvent(input$go_cdk,        updateNavbarPage(session, "main_nav", selected = "cdk"))
  observeEvent(input$go_master,     updateNavbarPage(session, "main_nav", selected = "master"))
  observeEvent(input$go_report,     updateNavbarPage(session, "main_nav", selected = "report"))
  observeEvent(input$back_home_cdk, updateNavbarPage(session, "main_nav", selected = "home"))

  ## ============================ CDK / webchem ============================
  smiles_table_rv <- reactiveVal(NULL)
  cdk_results_rv  <- reactiveVal(NULL)

  ## ---- Step 1: identifiers -> SMILES ----
  observeEvent(input$fetch_smiles, {
    req(input$cdk_ids)
    ids <- strsplit(input$cdk_ids, "\n")[[1]]

    withProgress(message = "Querying PubChem...", value = 0.3, {
      tryCatch({
        tabla <- getSmilesFromIdentifiers(ids, from = input$cdk_id_type)
        smiles_table_rv(tabla)
        showNotification(sprintf("%d SMILES obtained.", nrow(tabla)),
                         type = "message", duration = 5)
      }, error = function(e) {
        showNotification(paste("Error obtaining SMILES:", e$message),
                         type = "error", duration = 10)
        smiles_table_rv(NULL)
      })
    })
  })

  ## ---- Step 1 (alternative): manual SMILES entry ----
  observeEvent(input$use_manual_smiles, {
    req(input$cdk_manual_smiles)
    smiles_raw <- strsplit(input$cdk_manual_smiles, "\n")[[1]]
    smiles_raw <- trimws(smiles_raw)
    smiles_raw <- smiles_raw[smiles_raw != "" & !is.na(smiles_raw)]

    if (length(smiles_raw) == 0) {
      showNotification("Error: no SMILES entered.", type = "error", duration = 8)
      return(NULL)
    }

    ## Parse optional names
    names_vec <- NULL
    if (!is.null(input$cdk_manual_names) && nchar(trimws(input$cdk_manual_names)) > 0) {
      names_vec <- trimws(strsplit(input$cdk_manual_names, ",")[[1]])
      if (length(names_vec) != length(smiles_raw)) {
        showNotification(
          sprintf("Warning: %d SMILES but %d names. Names will be auto-generated.",
                  length(smiles_raw), length(names_vec)),
          type = "warning", duration = 8
        )
        names_vec <- NULL
      }
    }
    if (is.null(names_vec)) {
      names_vec <- paste0("Molecule_", seq_along(smiles_raw))
    }

    ## Build a data.frame consistent with the PubChem path
    tabla <- data.frame(
      query           = names_vec,
      cid             = NA_character_,
      CanonicalSMILES = smiles_raw,
      stringsAsFactors = FALSE
    )

    smiles_table_rv(tabla)
    showNotification(sprintf("%d SMILES loaded from manual entry.", nrow(tabla)),
                     type = "message", duration = 5)
  })

  ## ---- Step 1 (alternative): CSV upload with SMILES ----
  observeEvent(input$cdk_smiles_csv, {
    req(input$cdk_smiles_csv)
    tryCatch({
      d <- read.csv(input$cdk_smiles_csv$datapath,
                    check.names = FALSE,
                    stringsAsFactors = FALSE,
                    header = input$cdk_csv_has_header)

      ## Auto-detect columns that look like SMILES
      candidate_cols <- names(d)[grepl("smiles", names(d), ignore.case = TRUE)]
      if (length(candidate_cols) == 0) {
        ## Fall back to any character column
        char_cols <- names(d)[sapply(d, is.character)]
        candidate_cols <- char_cols
      }

      updateSelectInput(session, "cdk_csv_smiles_col",
                        choices = names(d),
                        selected = if (length(candidate_cols) > 0) candidate_cols[1] else names(d)[1])
    }, error = function(e) {
      showNotification(paste("Error reading CSV:", e$message),
                       type = "error", duration = 10)
    })
  })

  observeEvent(input$use_csv_smiles, {
    req(input$cdk_smiles_csv)
    req(input$cdk_csv_smiles_col)

    tryCatch({
      d <- read.csv(input$cdk_smiles_csv$datapath,
                    check.names = FALSE,
                    stringsAsFactors = FALSE,
                    header = input$cdk_csv_has_header)

      smiles_col <- input$cdk_csv_smiles_col
      if (!smiles_col %in% names(d)) {
        showNotification("Error: selected SMILES column not found in the file.",
                         type = "error", duration = 10)
        return(NULL)
      }

      smiles_raw <- trimws(as.character(d[[smiles_col]]))
      valid <- smiles_raw != "" & !is.na(smiles_raw)

      if (sum(valid) == 0) {
        showNotification("Error: no valid SMILES found in the selected column.",
                         type = "error", duration = 10)
        return(NULL)
      }

      ## Build a data.frame consistent with the PubChem path, preserving
      ## all original columns from the CSV
      d <- d[valid, , drop = FALSE]
      d$CanonicalSMILES <- smiles_raw[valid]

      ## Ensure a 'query' column exists for the molecule viewer
      if (!"query" %in% names(d)) {
        d$query <- if (smiles_col != "query") {
          as.character(d[[smiles_col]])
        } else {
          paste0("Molecule_", seq_len(nrow(d)))
        }
      }
      ## Ensure a 'cid' column exists (NA for non-PubChem sources)
      if (!"cid" %in% names(d)) {
        d$cid <- NA_character_
      }

      smiles_table_rv(d)

      if (sum(!valid) > 0) {
        showNotification(
          sprintf("%d SMILES loaded (%d empty/invalid rows skipped).",
                  nrow(d), sum(!valid)),
          type = "message", duration = 6
        )
      } else {
        showNotification(sprintf("%d SMILES loaded from CSV.", nrow(d)),
                         type = "message", duration = 5)
      }
    }, error = function(e) {
      showNotification(paste("Error loading SMILES from CSV:", e$message),
                       type = "error", duration = 10)
    })
  })

  ## ---- Load example dataset (drugs.csv) ----
  observeEvent(input$load_example_drugs, {
    drugs_file <- system.file("extdata", "drugs.csv", package = "admetshiny")
    if (drugs_file == "" || !file.exists(drugs_file)) {
      showNotification("Example dataset not found.", type = "error", duration = 8)
      return(NULL)
    }

    tryCatch({
      d <- read.csv(drugs_file, check.names = FALSE, stringsAsFactors = FALSE)

      ## The drugs.csv file has columns: name, smiles
      if (!"smiles" %in% tolower(names(d))) {
        showNotification("Error: no 'smiles' column found in example dataset.",
                         type = "error", duration = 10)
        return(NULL)
      }

      ## Find the SMILES column (case-insensitive)
      smiles_col <- names(d)[tolower(names(d)) == "smiles"][1]
      smiles_raw <- trimws(as.character(d[[smiles_col]]))
      valid <- smiles_raw != "" & !is.na(smiles_raw)

      if (sum(valid) == 0) {
        showNotification("Error: no valid SMILES found in example dataset.",
                         type = "error", duration = 10)
        return(NULL)
      }

      d <- d[valid, , drop = FALSE]
      d$CanonicalSMILES <- smiles_raw[valid]

      ## Ensure a 'query' column (use the 'name' column if present)
      if ("name" %in% tolower(names(d))) {
        name_col <- names(d)[tolower(names(d)) == "name"][1]
        d$query <- as.character(d[[name_col]])
      } else if (!"query" %in% names(d)) {
        d$query <- paste0("Molecule_", seq_len(nrow(d)))
      }
      ## Ensure a 'cid' column (NA for non-PubChem sources)
      if (!"cid" %in% names(d)) {
        d$cid <- NA_character_
      }

      smiles_table_rv(d)
      showNotification(sprintf("Example dataset loaded: %d drugs with SMILES.",
                               nrow(d)),
                       type = "message", duration = 5)
    }, error = function(e) {
      showNotification(paste("Error loading example dataset:", e$message),
                       type = "error", duration = 10)
    })
  })

  output$cdk_smiles_table <- renderDT({
    req(smiles_table_rv())
    smiles_table_rv()
  }, options = list(pageLength = 10, scrollX = TRUE))

  ## ---- SMILES download ----
  output$download_smiles <- downloadHandler(
    filename = function() {
      ext <- input$smiles_export_format
      paste0("PubChem_results.", ext)
    },
    content = function(file) {
      req(smiles_table_rv())
      if (input$smiles_export_format == "csv") {
        write.csv(smiles_table_rv(), file, row.names = FALSE)
      } else {
        if (!requireNamespace("openxlsx", quietly = TRUE)) {
          stop("The 'openxlsx' package is required for Excel export.",
               call. = FALSE)
        }
        openxlsx::write.xlsx(smiles_table_rv(), file, overwrite = TRUE)
      }
    }
  )

  ## ============== Molecule viewer ======================================
  observeEvent(smiles_table_rv(), {
    req(smiles_table_rv())
    df <- smiles_table_rv()

    ## Use row index as the selector value (works for all input methods,
    ## including manual SMILES and CSV upload where there is no PubChem CID)
    row_idx <- seq_len(nrow(df))

    display_names <- if ("Title" %in% names(df) && !all(is.na(df$Title))) {
      ifelse(!is.na(df$Title) & df$Title != "",
             paste0(df$Title, " (#", row_idx, ")"),
             paste0(df$query, " (#", row_idx, ")"))
    } else if ("cid" %in% names(df) && !all(is.na(df$cid))) {
      paste0(df$query, " (CID: ", df$cid, ")")
    } else {
      paste0(df$query, " (#", row_idx, ")")
    }

    choices <- setNames(as.character(row_idx), display_names)
    updateSelectInput(session, "cdk_mol_selector", choices = choices)
  })

  output$cdk_mol_info <- renderUI({
    req(input$cdk_mol_selector, smiles_table_rv())
    df <- smiles_table_rv()
    ## Select by row index (robust for all input methods)
    idx <- as.integer(input$cdk_mol_selector)
    if (is.na(idx) || idx < 1 || idx > nrow(df)) {
      return(tags$p("Select a molecule."))
    }
    row <- df[idx, , drop = FALSE]

    tagList(
      tags$p(tags$b("Name/Query: "), as.character(row$query)),
      if ("Title" %in% names(row) && !is.na(row$Title) && row$Title != "")
        tags$p(tags$b("Common name: "), as.character(row$Title)),
      if ("cid" %in% names(row) && !is.na(row$cid)) {
        tags$p(tags$b("CID: "), as.character(row$cid))
      },
      if ("MolecularFormula" %in% names(row) && !is.na(row$MolecularFormula))
        tags$p(tags$b("Formula: "), as.character(row$MolecularFormula)),
      if ("CanonicalSMILES" %in% names(row) && !is.na(row$CanonicalSMILES))
        tags$p(tags$b("SMILES: "), tags$br(),
               tags$code(style = "word-break: break-all; font-size: 11px;",
                         as.character(row$CanonicalSMILES))),
      if ("IUPACName" %in% names(row) && !is.na(row$IUPACName) &&
          row$IUPACName != "")
        tags$p(tags$b("IUPAC: "), as.character(row$IUPACName))
    )
  })

  output$cdk_molecule_image <- renderUI({
    req(input$cdk_mol_selector, smiles_table_rv())
    df <- smiles_table_rv()
    idx <- as.integer(input$cdk_mol_selector)
    if (is.na(idx) || idx < 1 || idx > nrow(df)) return(NULL)
    row <- df[idx, , drop = FALSE]

    ## If we have a PubChem CID, show the 2D image from PubChem
    if ("cid" %in% names(row) && !is.na(row$cid)) {
      cid <- row$cid
      url <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",
                    cid, "/PNG?image_size=large")
      tags$div(
        style = "width: 100%;",
        tags$img(src = url,
                 style = paste("max-width: 100%; max-height: 380px;",
                               "border: 1px solid #ddd; border-radius: 8px;",
                               "box-shadow: 0 2px 6px rgba(0,0,0,0.1);"),
                 alt = "Molecular structure"),
        tags$p(style = "color: #999; font-size: 11px; margin-top: 8px;",
               "Image from PubChem (PUG REST).")
      )
    } else {
      ## No CID available (manual SMILES or CSV upload): show SMILES text
      ## as a placeholder; the 2D structure will be rendered after CDK
      ## descriptor calculation in Step 2
      smiles_str <- if ("CanonicalSMILES" %in% names(row)) {
        as.character(row$CanonicalSMILES)
      } else {
        "(no SMILES available)"
      }
      tags$div(
        style = "width: 100%; text-align: center; padding: 40px 20px;",
        icon("image", style = "font-size: 48px; color: #ccc; margin-bottom: 15px;"),
        tags$p(style = "color: #999; font-size: 14px; margin-bottom: 10px;",
               "2D structure from PubChem is not available for manually",
               tags$br(),
               "entered SMILES or CSV-uploaded molecules."),
        tags$p(style = "color: #666; font-size: 12px; margin-bottom: 5px;",
               tags$b("SMILES:")),
        tags$code(style = "word-break: break-all; font-size: 12px; background: #f5f5f5; padding: 8px; border-radius: 4px; display: inline-block; max-width: 100%;",
                  smiles_str),
        tags$p(style = "color: #aaa; font-size: 11px; margin-top: 15px;",
               "Proceed to Step 2 to calculate CDK descriptors for this molecule.")
      )
    }
  })

  ## ---- Step 2: SMILES -> descriptors (rcdk / CDK) ----
  observeEvent(input$calc_cdk, {

    if (is.null(smiles_table_rv())) {
      showNotification("Error: first obtain the SMILES in Step 1.",
                       type = "error", duration = 10)
      return(NULL)
    }
    if (length(input$cdk_descriptors_core) == 0) {
      showNotification("Error: select at least one descriptor from the list.",
                       type = "error", duration = 10)
      return(NULL)
    }

    all_descriptors <- input$cdk_descriptors_core

    cols <- names(smiles_table_rv())
    ## Prefer CanonicalSMILES, then any column containing "smiles"
    if ("CanonicalSMILES" %in% cols) {
      smiles_col <- "CanonicalSMILES"
    } else {
      smiles_col <- cols[grepl("smiles", cols, ignore.case = TRUE)]
      if (length(smiles_col) == 0) {
        showNotification("Error: no SMILES column found in the loaded data.",
                         type = "error", duration = 10)
        return(NULL)
      }
      smiles_col <- smiles_col[1]
    }

    withProgress(message = "Calculating descriptors with CDK...", value = 0.5, {
      tryCatch({
        raw <- calcCDKDescriptors(smiles_table_rv()[[smiles_col]],
                                  which = all_descriptors)
        mapped <- mapCDKDescriptors(raw)

        mapped$TEMP_SMILES <- mapped$SMILES
        mapped$SMILES <- NULL

        ## Use lookup instead of merge to avoid row duplication
        orig_smiles <- as.character(smiles_table_rv()[[smiles_col]])
        cdk_lookup <- split(mapped, mapped$TEMP_SMILES)

        result_list <- lapply(seq_len(nrow(smiles_table_rv())), function(i) {
          smi <- orig_smiles[i]
          if (!is.na(smi) && smi %in% names(cdk_lookup)) {
            cdk_row <- cdk_lookup[[smi]][1, , drop = FALSE]
            cdk_row$TEMP_SMILES <- NULL
            ## Drop any original columns whose names collide with CDK-derived
            ## names to prevent duplicated column names in the cbind result.
            orig_row <- smiles_table_rv()[i, , drop = FALSE]
            dup_cols <- intersect(names(orig_row), names(cdk_row))
            if (length(dup_cols) > 0) {
              orig_row[[dup_cols[1]]] <- NULL
            }
            cbind(orig_row, cdk_row)
          } else {
            smiles_table_rv()[i, , drop = FALSE]
          }
        })

        final_df <- do.call(rbind, result_list)
        cdk_results_rv(final_df)

        showNotification(sprintf("Descriptors calculated for %d molecules.",
                                 nrow(final_df)),
                         type = "message", duration = 5)
      }, error = function(e) {
        showNotification(paste("Error calculating CDK descriptors:", e$message),
                         type = "error", duration = 15)
        cdk_results_rv(NULL)
      })
    })
  })

  output$cdk_results_table <- renderDT({
    req(cdk_results_rv())
    cdk_results_rv()
  }, options = list(pageLength = 10, scrollX = TRUE))

  ## ============== Step 3: CDK filtering =================================
  cdk_resultado <- eventReactive(input$cdk_run, {
    req(cdk_results_rv())
    tryCatch({
      out <- applyFilters(
        data = cdk_results_rv(), filters = input$cdk_filters,
        lipinski = list(mw = input$cdk_mw, logp = input$cdk_logp,
                        hba = input$cdk_hba, hbd = input$cdk_hbd,
                        violations = input$cdk_violations),
        veber = list(v_rb = input$cdk_v_rb, v_tpsa = input$cdk_v_tpsa,
                     v_hb_sum = input$cdk_v_hb_sum,
                     violations = input$cdk_veber_violations),
        ghose = list(g_mw_min = input$cdk_g_mw[1], g_mw_max = input$cdk_g_mw[2],
                     g_mr_min = input$cdk_g_mr[1], g_mr_max = input$cdk_g_mr[2],
                     g_logp_min = input$cdk_g_logp[1], g_logp_max = input$cdk_g_logp[2],
                     g_ha_min = input$cdk_g_ha[1], g_ha_max = input$cdk_g_ha[2],
                     violations = input$cdk_ghose_violations),
        egan = list(e_tpsa = input$cdk_e_tpsa, e_logp = input$cdk_e_logp,
                    violations = input$cdk_egan_violations),
        muegge = list(m_mw_min = input$cdk_m_mw[1], m_mw_max = input$cdk_m_mw[2],
                      m_logp_min = input$cdk_m_logp[1], m_logp_max = input$cdk_m_logp[2],
                      m_hba = input$cdk_m_hba, m_hbd = input$cdk_m_hbd,
                      m_rb = input$cdk_m_rb, m_tpsa = input$cdk_m_tpsa,
                      violations = input$cdk_muegge_violations)
      )
      if (nrow(out) == 0)
        showNotification("No compound meets the selected filters.",
                         type = "warning", duration = 6)
      out
    }, error = function(e) {
      showNotification(paste("Error applying filters:", e$message),
                       type = "error", duration = 8)
      cdk_results_rv()[0, ]
    })
  })

  output$cdk_tabla <- renderDT({
    req(cdk_resultado())
    cdk_resultado()
  }, options = list(pageLength = 10, scrollX = TRUE))

  output$cdk_download <- downloadHandler(
    filename = function() "CDK_Filtered_Dataset.csv",
    content  = function(file) write.csv(cdk_resultado(), file, row.names = FALSE)
  )

  ## ============== Dynamic plot controls for CDK =========================
  output$cdk_radar_id_selector <- renderUI({
    req(cdk_resultado())
    tagList(
      selectInput("cdk_radar_id_col", "Identifier column",
                  choices = names(cdk_resultado()),
                  selected = if ("query" %in% names(cdk_resultado())) "query" else names(cdk_resultado())[1]),
      selectizeInput("cdk_radar_ids", "Molecules to compare (max 5)",
                     choices = NULL, multiple = TRUE,
                     options = list(maxItems = 5, placeholder = "Select up to 5 molecules"))
    )
  })

  observeEvent(input$cdk_radar_id_col, {
    req(cdk_resultado(), input$cdk_radar_id_col)
    vals <- unique(as.character(cdk_resultado()[[input$cdk_radar_id_col]]))
    updateSelectizeInput(session, "cdk_radar_ids", choices = vals, server = TRUE)
  })

  output$cdk_smiles_col_selector <- renderUI({
    req(cdk_resultado())
    cols <- names(cdk_resultado())
    guess <- cols[grepl("smiles", cols, ignore.case = TRUE)]
    default <- if (length(guess) > 0) guess[1] else cols[1]
    selectInput("cdk_smiles_col", "SMILES column", choices = cols, selected = default)
  })

  output$cdk_label_col_selector <- renderUI({
    req(cdk_resultado())
    cols <- names(cdk_resultado())
    default <- if ("Title" %in% cols) "Title" else if ("query" %in% cols) "query" else cols[1]
    selectInput("cdk_label_col", "Label column", choices = cols, selected = default)
  })

  output$cdk_pca_controls <- renderUI({
    req(cdk_resultado())
    cols <- names(cdk_resultado())
    numeric_cols <- cols[sapply(cdk_resultado(), is.numeric)]
    tagList(
      selectInput("cdk_pca_variables", "Variables", choices = numeric_cols,
                  selected = intersect(c("MW","TPSA","LogP","#H-bond acceptors","#H-bond donors","#Rotatable bonds"), numeric_cols),
                  multiple = TRUE),
      selectInput("cdk_pca_color", "Colour by", choices = c("None", cols), selected = "None"),
      selectInput("cdk_pca_label", "Labels", choices = c("None", cols),
                  selected = if ("query" %in% cols) "query" else "None"),
      checkboxInput("cdk_pca_ellipse", "Show ellipses", TRUE),
      checkboxInput("cdk_pca_scale", "Scale variables", TRUE)
    )
  })

  output$cdk_tsne_controls <- renderUI({
    req(cdk_resultado())
    cols <- names(cdk_resultado())
    numeric_cols <- cols[sapply(cdk_resultado(), is.numeric)]
    tagList(
      selectInput("cdk_tsne_variables", "Variables", choices = numeric_cols,
                  selected = intersect(c("MW","TPSA","LogP","#H-bond acceptors","#H-bond donors","#Rotatable bonds"), numeric_cols),
                  multiple = TRUE),
      sliderInput("cdk_tsne_perplexity", "Perplexity", min = 5, max = 50, value = 30),
      sliderInput("cdk_tsne_iter", "Iterations", min = 500, max = 5000, value = 1000, step = 500),
      selectInput("cdk_tsne_color", "Colour by", choices = c("None", cols), selected = "None"),
      selectInput("cdk_tsne_label", "Labels", choices = c("None", cols),
                  selected = if ("query" %in% cols) "query" else "None"),
      checkboxInput("cdk_tsne_scale", "Scale variables", TRUE)
    )
  })

  output$cdk_umap_controls <- renderUI({
    req(cdk_resultado())
    cols <- names(cdk_resultado())
    numeric_cols <- cols[sapply(cdk_resultado(), is.numeric)]
    tagList(
      selectInput("cdk_umap_variables", "Variables", choices = numeric_cols,
                  selected = intersect(c("MW","TPSA","LogP","#H-bond acceptors","#H-bond donors","#Rotatable bonds"), numeric_cols),
                  multiple = TRUE),
      sliderInput("cdk_umap_n_neighbors", "n_neighbors", min = 2, max = 50, value = 15),
      sliderInput("cdk_umap_min_dist", "min_dist", min = 0, max = 1, value = 0.1, step = 0.05),
      selectInput("cdk_umap_color", "Colour by", choices = c("None", cols), selected = "None"),
      selectInput("cdk_umap_label", "Labels", choices = c("None", cols),
                  selected = if ("query" %in% cols) "query" else "None"),
      checkboxInput("cdk_umap_scale", "Scale variables", TRUE)
    )
  })

  output$cdk_parallel_controls <- renderUI({
    req(cdk_resultado())
    cols <- names(cdk_resultado())
    numeric_cols <- cols[sapply(cdk_resultado(), is.numeric)]
    tagList(
      selectInput("cdk_parallel_variables", "Variables", choices = numeric_cols,
                  selected = intersect(c("MW","TPSA","LogP","#H-bond acceptors","#H-bond donors","#Rotatable bonds"), numeric_cols),
                  multiple = TRUE),
      selectInput("cdk_parallel_color", "Colour by", choices = c("None", cols), selected = "None"),
      checkboxInput("cdk_parallel_scale", "Scale variables", TRUE)
    )
  })

  output$cdk_cluster_heatmap_controls <- renderUI({
    req(cdk_resultado())
    cols <- names(cdk_resultado())
    numeric_cols <- cols[sapply(cdk_resultado(), is.numeric)]
    ## Look for a suitable label column: prefer text columns that could be IDs
    char_cols <- cols[!sapply(cdk_resultado(), is.numeric)]
    default_id <- if (any(c("Title", "query", "Name", "name", "ID", "Compound", "Molecule") %in% cols)) {
      intersect(c("Title", "query", "Name", "name", "ID", "Compound", "Molecule"), cols)[1]
    } else if (length(char_cols) > 0) {
      char_cols[1]
    } else {
      "None"
    }
    tagList(
      selectInput("cdk_cluster_heatmap_variables", "Variables", choices = numeric_cols,
                  selected = c("MW","TPSA","LogP","#H-bond acceptors","#H-bond donors","#Rotatable bonds"),
                  multiple = TRUE),
      selectInput("cdk_cluster_heatmap_id_col", "Label column (optional)", choices = c("None", cols),
                  selected = default_id),
      selectInput("cdk_cluster_heatmap_method", "Clustering method",
                  choices = c("Ward D2" = "ward.D2", "Ward D" = "ward.D",
                              "Complete" = "complete", "Average (UPGMA)" = "average",
                              "Single" = "single", "McQuitty" = "mcquitty",
                              "Median" = "median", "Centroid" = "centroid"),
                  selected = "ward.D2"),
      checkboxInput("cdk_cluster_heatmap_scale", "Scale variables (z-score)", TRUE)
    )
  })

  output$cdk_violin_controls <- renderUI({
    req(cdk_resultado())
    cols <- names(cdk_resultado())
    numeric_cols <- cols[sapply(cdk_resultado(), is.numeric)]
    categorical_cols <- cols[!sapply(cdk_resultado(), is.numeric)]
    tagList(
      selectInput("cdk_violin_variable", "Variable", choices = numeric_cols, selected = "MW"),
      selectInput("cdk_violin_group", "Group by", choices = categorical_cols,
                  selected = if ("query" %in% categorical_cols) "query" else categorical_cols[1]),
      checkboxInput("cdk_violin_box", "Show boxplot", TRUE),
      checkboxInput("cdk_violin_points", "Show points", FALSE)
    )
  })

  ## ---- CDK: renderUI blocks for the new plot types ----
  output$cdk_histogram_custom_controls <- renderUI({
    req(cdk_resultado())
    cols <- names(cdk_resultado())
    numeric_cols <- cols[sapply(cdk_resultado(), is.numeric)]
    categorical_cols <- cols[!sapply(cdk_resultado(), is.numeric)]
    tagList(
      selectInput("cdk_hc_variable", "Variable", choices = numeric_cols, selected = "MW"),
      sliderInput("cdk_hc_bins", "Number of bins", min = 5, max = 100, value = 30),
      selectInput("cdk_hc_group", "Group by (overlay)", choices = c("None", categorical_cols), selected = "None"),
      checkboxInput("cdk_hc_density", "Show density curve", TRUE),
      checkboxInput("cdk_hc_rug", "Show rug plot", FALSE)
    )
  })

  ## ============== CDK plot rendering ====================================
  cdk_base_plot_types <- c("Radar plot (Chemical profile)",
                           "Tanimoto / AGNES (Structural similarity)")

  cdk_currentPlot <- reactive({
    req(cdk_resultado())
    validate(need(nrow(cdk_resultado()) > 0, "No data to plot."))
    pt <- input$cdk_plot_type

    if (pt %in% cdk_base_plot_types) {
      if (pt == "Radar plot (Chemical profile)") {
        req(input$cdk_radar_id_col)
        validate(need(length(input$cdk_radar_ids) >= 1,
                      "Select at least one molecule for the radar."))
        function() plotRadar(cdk_resultado(),
                             id_col = input$cdk_radar_id_col,
                             ids = input$cdk_radar_ids)
      } else {
        req(input$cdk_smiles_col)
        function() plotTanimoto(cdk_resultado(),
                                smiles_col = input$cdk_smiles_col,
                                label_col = input$cdk_label_col,
                                max_n = input$cdk_tanimoto_max_n,
                                method = input$cdk_agnes_method)
      }
    } else {
      switch(pt,
             "Boiled Egg" = plotBoiledEgg(cdk_resultado()),
             "Molecular Weight" = plotMW(cdk_resultado()),
             "TPSA" = plotTPSA(cdk_resultado()),
             "LogP" = plotLogP(cdk_resultado()),
             "Correlation Heatmap" = plotCorrHeatmap(cdk_resultado()),
             "Principal Component Analysys (PCA - Chemical space)" =
               apply_palette(plotPCA(data = cdk_resultado(),
                                     variables = input$cdk_pca_variables,
                                     color_by = input$cdk_pca_color,
                                     label_by = input$cdk_pca_label,
                                     scale_data = input$cdk_pca_scale,
                                     ellipse = input$cdk_pca_ellipse),
                             input$cdk_palette, cdk_resultado(), input$cdk_pca_color),
             "t-SNE (Chemical space)" =
               apply_palette(plotTSNE(data = cdk_resultado(),
                                      variables = input$cdk_tsne_variables,
                                      color_by = input$cdk_tsne_color,
                                      label_by = input$cdk_tsne_label,
                                      perplexity = input$cdk_tsne_perplexity,
                                      max_iter = input$cdk_tsne_iter,
                                      scale_data = input$cdk_tsne_scale),
                             input$cdk_palette, cdk_resultado(), input$cdk_tsne_color),
             "UMAP (Chemical space)" =
               apply_palette(plotUMAP(data = cdk_resultado(),
                                      variables = input$cdk_umap_variables,
                                      color_by = input$cdk_umap_color,
                                      label_by = input$cdk_umap_label,
                                      n_neighbors = input$cdk_umap_n_neighbors,
                                      min_dist = input$cdk_umap_min_dist,
                                      scale_data = input$cdk_umap_scale),
                             input$cdk_palette, cdk_resultado(), input$cdk_umap_color),
             "Parallel Coordinates" =
               apply_palette(plotParallel(data = cdk_resultado(),
                                          variables = input$cdk_parallel_variables,
                                          color_by = input$cdk_parallel_color,
                                          scale_data = input$cdk_parallel_scale),
                             input$cdk_palette, cdk_resultado(), input$cdk_parallel_color),
             "Cluster Heatmap (Dendrogram)" =
               function() plotClusterHeatmap(data = cdk_resultado(),
                                             variables = input$cdk_cluster_heatmap_variables,
                                             id_col = .safe_id_col(input$cdk_cluster_heatmap_id_col),
                                             method = input$cdk_cluster_heatmap_method,
                                             scale_data = input$cdk_cluster_heatmap_scale,
                                             palette = input$cdk_palette),
             "Violin Plot" =
               apply_palette(plotViolin(data = cdk_resultado(),
                                         variable = input$cdk_violin_variable,
                                         group_by = input$cdk_violin_group,
                                         show_box = input$cdk_violin_box,
                                         show_points = input$cdk_violin_points),
                             input$cdk_palette, cdk_resultado(), input$cdk_violin_group),
             "Custom Histogram" =
               plotHistogramCustom(data = cdk_resultado(),
                                   variable = input$cdk_hc_variable,
                                   bins = input$cdk_hc_bins,
                                   group_by = input$cdk_hc_group,
                                   show_density = input$cdk_hc_density,
                                   show_rug = input$cdk_hc_rug)
      )
    }
  })

  output$cdk_plot <- renderPlot({
    p <- cdk_currentPlot()
    if (is.function(p)) p() else print(p)
  })

  output$cdk_downloadPlot <- downloadHandler(
    filename = function() {
      paste0("cdk_", gsub("[^A-Za-z0-9]+", "_", input$cdk_plot_type), ".", input$cdk_format)
    },
    content = function(file) {
      p <- cdk_currentPlot()
      if (is.function(p)) {
        open_device <- switch(input$cdk_format,
          "png"  = function() grDevices::png(file, width = input$cdk_width, height = input$cdk_height, units = "in", res = as.numeric(input$cdk_dpi)),
          "jpeg" = function() grDevices::jpeg(file, width = input$cdk_width, height = input$cdk_height, units = "in", res = as.numeric(input$cdk_dpi)),
          "tiff" = function() grDevices::tiff(file, width = input$cdk_width, height = input$cdk_height, units = "in", res = as.numeric(input$cdk_dpi)),
          "pdf"  = function() grDevices::pdf(file, width = input$cdk_width, height = input$cdk_height),
          "svg"  = function() grDevices::svg(file, width = input$cdk_width, height = input$cdk_height),
          function() grDevices::png(file, width = input$cdk_width, height = input$cdk_height, units = "in", res = as.numeric(input$cdk_dpi))
        )
        open_device()
        p()
        grDevices::dev.off()
      } else {
        ggplot2::ggsave(filename = file, plot = p, device = input$cdk_format,
                        width = input$cdk_width, height = input$cdk_height,
                        units = "in", dpi = as.numeric(input$cdk_dpi))
      }
    }
  )

  ## ============================ ADMET Master ============================
  ## Four-step wizard module: upload any CSV/Excel dataset, map its columns
  ## to the application's standard schema, optionally calculate missing
  ## descriptors with CDK, apply the standard drug-likeness filters and
  ## produce the same plot catalogue as the other modules.
  master_raw_rv    <- reactiveVal(NULL)   ## raw uploaded data
  master_mapped_rv <- reactiveVal(NULL)   ## after mapping + standardization
  master_smiles_col_rv <- reactiveVal(NULL)  ## detected SMILES column (raw)

  ## Helper: convert a column name to a valid Shiny input id.
  .master_input_id <- function(idx) paste0("master_map_", idx)

  ## Helper: try to auto-map a user column name to a standard field code.
  .master_auto_select <- function(col_name) {
    lc <- tolower(trimws(col_name))
    if (lc %in% c("smiles", "canonicalsmiles", "isomericsmiles",
                  "mol", "molstr", "structure")) return("SMILES")
    if (lc %in% c("name", "id", "compound", "molecule", "title",
                  "molecule_name", "compound_name", "drug", "drug_name"))
      return("Name")
    if (lc == "mw" || lc == "molecular_weight" ||
        grepl("molecular.*weight", lc) || lc == "molwt") return("MW")
    if (lc == "wlogp") return("WLOGP")
    if (lc == "mlogp" || lc == "xlogp3" || lc == "ilogp") return("LogP")
    if (lc == "logp" || grepl("^log[\\s_.-]*p", lc) || lc == "clogp")
      return("LogP")
    if (lc == "tpsa" || lc == "topological_polar_surface_area")
      return("TPSA")
    if (lc %in% c("hbd", "nhd", "h_donors", "hdon", "nhd",
                  "h-bond donors", "h_bond_donors", "donors"))
      return("HBD")
    if (lc %in% c("hba", "nha", "h_acceptors", "hac", "nha",
                  "h-bond acceptors", "h_bond_acceptors", "acceptors"))
      return("HBA")
    if (grepl("rotat", lc)) return("Rotatable Bonds")
    if (lc %in% c("mr", "molar_refractivity", "molar refractivity",
                  "amr")) return("Molar Refractivity")
    if (grepl("arom", lc) && grepl("atom|heavy", lc))
      return("Aromatic Heavy Atoms")
    if (grepl("heavy", lc) && grepl("atom", lc)) return("Heavy Atoms")
    if (grepl("log[\\s_.-]*s", lc) || lc == "logs" || lc == "esol_log_s")
      return("LogS")
    if (grepl("log[\\s_.-]*d", lc) || lc == "logd") return("LogD")
    if (grepl("gi[\\s_.-]*abs", lc) || lc == "hia" ||
        lc == "absorption") return("GI Absorption")
    if (grepl("bbb", lc)) return("BBB Permeant")
    if (grepl("pgp", lc) || grepl("p-gp", lc) ||
        grepl("p.glycoprotein", lc)) return("Pgp Substrate")
    "None"
  }

  ## ---------------- Step 1: file upload ----------------
  observeEvent(input$master_file, {
    req(input$master_file)
    tryCatch({
      ext <- tolower(tools::file_ext(input$master_file$name))
      if (ext == "csv") {
        d <- read.csv(input$master_file$datapath,
                      check.names = FALSE, stringsAsFactors = FALSE)
      } else if (ext %in% c("xlsx", "xls")) {
        if (!requireNamespace("openxlsx", quietly = TRUE)) {
          stop("The 'openxlsx' package is required to read Excel files. ",
               "Install it with: install.packages('openxlsx').",
               call. = FALSE)
        }
        d <- openxlsx::read.xlsx(input$master_file$datapath, sheet = 1)
      } else {
        stop("Unsupported file format. Please upload a CSV or Excel file.")
      }
      ## Drop completely empty columns (often an artifact of CSV export)
      keep <- vapply(d, function(col) {
        !all(is.na(col)) && !all(as.character(col) == "")
      }, logical(1))
      d <- d[, keep, drop = FALSE]

      master_raw_rv(d)
      master_mapped_rv(NULL)
      master_smiles_col_rv(detectSMILESColumn(d))
      showNotification(sprintf("Loaded %d rows and %d columns.",
                               nrow(d), ncol(d)),
                       type = "message", duration = 5)
    }, error = function(e) {
      showNotification(paste("Error reading the file:", e$message),
                       type = "error", duration = 10)
      master_raw_rv(NULL)
      master_mapped_rv(NULL)
      master_smiles_col_rv(NULL)
    })
  })

  ## Step 1 outputs: preview, column info, SMILES detection banner
  output$master_preview <- renderDT({
    req(master_raw_rv())
    master_raw_rv()
  }, options = list(pageLength = 10, scrollX = TRUE))

  output$master_column_info <- renderDT({
    req(master_raw_rv())
    info <- detectColumnTypes(master_raw_rv())
    info
  }, options = list(pageLength = 25, scrollX = TRUE))

  output$master_smiles_detection <- renderUI({
    req(master_raw_rv())
    col <- master_smiles_col_rv()
    if (is.null(col)) {
      tags$div(
        class = "alert alert-warning",
        icon("triangle-exclamation"),
        tags$b(" No SMILES column detected automatically. "),
        "You can still map a column to SMILES manually in Step 2."
      )
    } else {
      tags$div(
        class = "alert alert-success",
        icon("check-circle"),
        tags$b(" SMILES column detected: "),
        tags$code(col),
        " (you can override this in Step 2)."
      )
    }
  })

  ## ---------------- Step 2: column mapping UI ----------------
  output$master_mapping_ui <- renderUI({
    req(master_raw_rv())
    cols <- names(master_raw_rv())
    detect_smiles <- master_smiles_col_rv()

    field_choices <- c(
      "None (skip)"                                = "None",
      "SMILES (string) - for CDK/Tanimoto"         = "SMILES",
      "Name/ID (string)"                           = "Name",
      "MW (numeric) - Lipinski/Ghose/Muegge"       = "MW",
      "LogP (numeric) - All filters, BOILED-Egg"   = "LogP",
      "WLOGP (numeric) - Official BOILED-Egg"      = "WLOGP",
      "TPSA (numeric) - Veber/Egan/BOILED-Egg"     = "TPSA",
      "HBD (numeric) - Lipinski/Veber/Muegge"      = "HBD",
      "HBA (numeric) - Lipinski/Veber/Muegge"      = "HBA",
      "Rotatable Bonds (numeric) - Veber/Muegge"   = "Rotatable Bonds",
      "Molar Refractivity (numeric) - Ghose"       = "Molar Refractivity",
      "Heavy Atoms (numeric) - Ghose"              = "Heavy Atoms",
      "Aromatic Heavy Atoms (numeric)"             = "Aromatic Heavy Atoms",
      "GI Absorption (categorical)"                = "GI Absorption",
      "GI Absorption (numeric 0-1)"                = "GI Absorption_num",
      "BBB Permeant (categorical)"                 = "BBB Permeant",
      "BBB Permeant (numeric 0-1)"                 = "BBB Permeant_num",
      "Pgp Substrate (categorical)"                = "Pgp Substrate",
      "Pgp Substrate (numeric 0-1)"                = "Pgp Substrate_num",
      "LogS (numeric)"                             = "LogS",
      "LogD (numeric)"                             = "LogD"
    )

    tagList(
      lapply(seq_along(cols), function(i) {
        col <- cols[i]
        default_val <- .master_auto_select(col)
        if ((is.null(default_val) || default_val == "None") &&
            !is.null(detect_smiles) && col == detect_smiles) {
          default_val <- "SMILES"
        }
        fluidRow(
          column(4, tags$div(style = "padding-top: 8px;",
                             tags$b(col))),
          column(8, selectInput(.master_input_id(i), NULL,
                                choices = field_choices,
                                selected = default_val))
        )
      })
    )
  })

  ## ---------------- Step 2: calculate button ----------------
  observeEvent(input$master_calculate, {
    req(master_raw_rv())
    raw <- master_raw_rv()
    cols <- names(raw)

    ## Build mapping vector from the dynamic selectInputs.
    mapping <- setNames(
      vapply(seq_along(cols), function(i) {
        val <- input[[.master_input_id(i)]]
        if (is.null(val)) "None" else val
      }, character(1)),
      cols
    )

    ## At least one column must be mapped to something other than None.
    if (all(mapping == "None")) {
      showNotification("Please map at least one column to a standard field.",
                       type = "error", duration = 6)
      return(NULL)
    }

    calc_cdk <- isTRUE(input$master_calculate_cdk)
    has_smiles <- "SMILES" %in% mapping
    needs_cdk <- calc_cdk && has_smiles

    tryCatch({
      if (needs_cdk) {
        withProgress(
          message = "Calculating missing descriptors with CDK...",
          detail  = "This may take a moment for large datasets.",
          value   = 0.5, {
            mapped <- mapADMETColumns(raw, mapping, calculate_cdk = TRUE)
          }
        )
      } else {
        mapped <- mapADMETColumns(raw, mapping, calculate_cdk = calc_cdk)
      }
      master_mapped_rv(mapped)
      showNotification(
        sprintf("Dataset standardized: %d rows, %d columns.",
                nrow(mapped), ncol(mapped)),
        type = "message", duration = 5
      )
    }, error = function(e) {
      showNotification(paste("Error during standardization:", e$message),
                       type = "error", duration = 10)
      master_mapped_rv(NULL)
    })
  })

  ## ---------------- Step 2: mapping summary table ----------------
  output$master_mapping_summary <- renderDT({
    req(master_mapped_rv())
    d <- master_mapped_rv()
    data.frame(
      column     = names(d),
      type       = vapply(d, function(col) {
                     if (is.numeric(col)) "numeric" else "string"
                   }, character(1)),
      n_non_NA   = vapply(d, function(col) {
                     sum(!is.na(col))
                   }, integer(1)),
      stringsAsFactors = FALSE
    )
  }, options = list(pageLength = 25, scrollX = TRUE))

  ## ---------------- Step 2: validation message ----------------
  output$master_validation <- renderUI({
    req(master_mapped_rv())
    d <- master_mapped_rv()
    cols <- names(d)

    checks <- list(
      list(name = "Lipinski filter",
           needs = c("MW", "LogP", "#H-bond acceptors", "#H-bond donors")),
      list(name = "Veber filter",
           needs = c("#Rotatable bonds", "TPSA",
                     "#H-bond acceptors", "#H-bond donors")),
      list(name = "Ghose filter",
           needs = c("MW", "MR", "LogP", "#Heavy atoms")),
      list(name = "Egan filter",    needs = c("TPSA", "LogP")),
      list(name = "Muegge filter",
           needs = c("MW", "LogP", "#H-bond acceptors", "#H-bond donors",
                     "#Rotatable bonds", "TPSA")),
      list(name = "BOILED-Egg",     needs = c("LogP", "TPSA")),
      list(name = "Tanimoto / AGNES", needs = c("SMILES")),
      list(name = "Radar plot",
           needs = c("MW", "LogP", "TPSA",
                     "#H-bond acceptors", "#H-bond donors"))
    )

    items <- vapply(checks, function(chk) {
      present <- chk$needs %in% cols
      if (all(present)) {
        paste0("<span style='color: #2e7d32;'>&#10004;</span> <b>",
               chk$name, "</b> &mdash; ready")
      } else {
        missing <- chk$needs[!present]
        paste0("<span style='color: #c62828;'>&#10006;</span> <b>",
               chk$name, "</b> &mdash; missing: ",
               paste(missing, collapse = ", "))
      }
    }, character(1))

    HTML(paste(items, collapse = "<br>"))
  })

  ## ---------------- Step 3: filter ----------------
  master_resultado <- eventReactive(input$master_run, {
    req(master_mapped_rv())
    tryCatch({
      out <- applyFilters(
        data = master_mapped_rv(), filters = input$master_filters,
        lipinski = list(mw = input$master_mw, logp = input$master_logp,
                        hba = input$master_hba, hbd = input$master_hbd,
                        violations = input$master_violations),
        veber = list(v_rb = input$master_v_rb, v_tpsa = input$master_v_tpsa,
                     v_hb_sum = input$master_v_hb_sum,
                     violations = input$master_veber_violations),
        ghose = list(g_mw_min = input$master_g_mw[1],
                     g_mw_max = input$master_g_mw[2],
                     g_mr_min = input$master_g_mr[1],
                     g_mr_max = input$master_g_mr[2],
                     g_logp_min = input$master_g_logp[1],
                     g_logp_max = input$master_g_logp[2],
                     g_ha_min = input$master_g_ha[1],
                     g_ha_max = input$master_g_ha[2],
                     violations = input$master_ghose_violations),
        egan = list(e_tpsa = input$master_e_tpsa,
                    e_logp = input$master_e_logp,
                    violations = input$master_egan_violations),
        muegge = list(m_mw_min = input$master_m_mw[1],
                      m_mw_max = input$master_m_mw[2],
                      m_logp_min = input$master_m_logp[1],
                      m_logp_max = input$master_m_logp[2],
                      m_hba = input$master_m_hba, m_hbd = input$master_m_hbd,
                      m_rb = input$master_m_rb, m_tpsa = input$master_m_tpsa,
                      violations = input$master_muegge_violations)
      )
      if (nrow(out) == 0) {
        showNotification("No compound meets the selected filters.",
                         type = "warning", duration = 6)
      }
      out
    }, error = function(e) {
      showNotification(paste("Error applying filters:", e$message),
                       type = "error", duration = 8)
      master_mapped_rv()[0, ]
    })
  })

  output$master_tabla <- renderDT({
    req(master_resultado())
    master_resultado()
  }, options = list(pageLength = 10, scrollX = TRUE))

  output$master_download <- downloadHandler(
    filename = function() "ADMET_Master_Filtered.csv",
    content  = function(file) write.csv(master_resultado(), file,
                                        row.names = FALSE)
  )

  ## ---------------- Step 4: dynamic plot controls ----------------
  output$master_radar_id_selector <- renderUI({
    req(master_resultado())
    cols <- names(master_resultado())
    default <- if ("Name" %in% cols) "Name" else cols[1]
    tagList(
      selectInput("master_radar_id_col", "Identifier column",
                  choices = cols, selected = default),
      selectizeInput("master_radar_ids", "Molecules to compare (max 5)",
                     choices = NULL, multiple = TRUE,
                     options = list(maxItems = 5,
                                    placeholder = "Select up to 5 molecules"))
    )
  })

  observeEvent(input$master_radar_id_col, {
    req(master_resultado(), input$master_radar_id_col)
    vals <- unique(as.character(
      master_resultado()[[input$master_radar_id_col]]))
    updateSelectizeInput(session, "master_radar_ids",
                         choices = vals, server = TRUE)
  })

  output$master_smiles_col_selector <- renderUI({
    req(master_resultado())
    cols <- names(master_resultado())
    guess <- cols[grepl("smiles", cols, ignore.case = TRUE)]
    default <- if (length(guess) > 0) guess[1] else cols[1]
    selectInput("master_smiles_col", "SMILES column",
                choices = cols, selected = default)
  })

  output$master_label_col_selector <- renderUI({
    req(master_resultado())
    cols <- names(master_resultado())
    default <- if ("Name" %in% cols) "Name" else cols[1]
    selectInput("master_label_col", "Label column",
                choices = cols, selected = default)
  })

  ## Dynamic BOILED-Egg LogP selector — only show LogP variants that exist
  output$master_boiled_egg_logp_ui <- renderUI({
    req(master_resultado())
    cols <- names(master_resultado())
    logp_cols <- cols[grepl("logp|LogP|WLOGP|MLOGP|XLOGP|iLOGP|Consensus",
                            cols, ignore.case = TRUE)]
    if (length(logp_cols) == 0) logp_cols <- "LogP"
    default <- if ("WLOGP" %in% logp_cols) "WLOGP" else logp_cols[1]
    selectInput("master_boiled_egg_logp", "LogP source for BOILED-Egg",
                choices = logp_cols, selected = default)
  })

  output$master_pca_controls <- renderUI({
    req(master_resultado())
    cols <- names(master_resultado())
    numeric_cols <- cols[sapply(master_resultado(), is.numeric)]
    tagList(
      selectInput("master_pca_variables", "Variables", choices = numeric_cols,
                  selected = intersect(c("MW","TPSA","LogP","#H-bond acceptors",
                                          "#H-bond donors","#Rotatable bonds"),
                                        numeric_cols),
                  multiple = TRUE),
      selectInput("master_pca_color", "Colour by",
                  choices = c("None", cols),
                  selected = if ("GI absorption" %in% cols) "GI absorption" else "None"),
      selectInput("master_pca_label", "Labels",
                  choices = c("None", cols),
                  selected = if ("Name" %in% cols) "Name" else "None"),
      checkboxInput("master_pca_ellipse", "Show ellipses", TRUE),
      checkboxInput("master_pca_scale", "Scale variables", TRUE)
    )
  })

  output$master_tsne_controls <- renderUI({
    req(master_resultado())
    cols <- names(master_resultado())
    numeric_cols <- cols[sapply(master_resultado(), is.numeric)]
    default_color <- if ("GI absorption" %in% cols) "GI absorption" else "None"
    default_label <- if ("Name" %in% cols) "Name" else "None"
    tagList(
      selectInput("master_tsne_variables", "Variables", choices = numeric_cols,
                  selected = intersect(c("MW","TPSA","LogP","#H-bond acceptors",
                                          "#H-bond donors","#Rotatable bonds"),
                                        numeric_cols),
                  multiple = TRUE),
      sliderInput("master_tsne_perplexity", "Perplexity",
                  min = 5, max = 50, value = 30),
      sliderInput("master_tsne_iter", "Iterations",
                  min = 500, max = 5000, value = 1000, step = 500),
      selectInput("master_tsne_color", "Colour by",
                  choices = c("None", cols), selected = default_color),
      selectInput("master_tsne_label", "Labels",
                  choices = c("None", cols), selected = default_label),
      checkboxInput("master_tsne_scale", "Scale variables", TRUE)
    )
  })

  output$master_umap_controls <- renderUI({
    req(master_resultado())
    cols <- names(master_resultado())
    numeric_cols <- cols[sapply(master_resultado(), is.numeric)]
    default_color <- if ("GI absorption" %in% cols) "GI absorption" else "None"
    default_label <- if ("Name" %in% cols) "Name" else "None"
    tagList(
      selectInput("master_umap_variables", "Variables", choices = numeric_cols,
                  selected = intersect(c("MW","TPSA","LogP","#H-bond acceptors",
                                          "#H-bond donors","#Rotatable bonds"),
                                        numeric_cols),
                  multiple = TRUE),
      sliderInput("master_umap_n_neighbors", "n_neighbors",
                  min = 2, max = 50, value = 15),
      sliderInput("master_umap_min_dist", "min_dist",
                  min = 0, max = 1, value = 0.1, step = 0.05),
      selectInput("master_umap_color", "Colour by",
                  choices = c("None", cols), selected = default_color),
      selectInput("master_umap_label", "Labels",
                  choices = c("None", cols), selected = default_label),
      checkboxInput("master_umap_scale", "Scale variables", TRUE)
    )
  })

  output$master_parallel_controls <- renderUI({
    req(master_resultado())
    cols <- names(master_resultado())
    numeric_cols <- cols[sapply(master_resultado(), is.numeric)]
    tagList(
      selectInput("master_parallel_variables", "Variables",
                  choices = numeric_cols,
                  selected = intersect(c("MW","TPSA","LogP","#H-bond acceptors",
                                          "#H-bond donors","#Rotatable bonds"),
                                        numeric_cols),
                  multiple = TRUE),
      selectInput("master_parallel_color", "Colour by",
                  choices = c("None", cols),
                  selected = if ("GI absorption" %in% cols) "GI absorption" else "None"),
      checkboxInput("master_parallel_scale", "Scale variables", TRUE)
    )
  })

  output$master_cluster_heatmap_controls <- renderUI({
    req(master_resultado())
    cols <- names(master_resultado())
    numeric_cols <- cols[sapply(master_resultado(), is.numeric)]
    char_cols <- cols[!sapply(master_resultado(), is.numeric)]
    default_id <- if (any(c("Name","name","ID","Compound","Molecule") %in% cols)) {
      intersect(c("Name","name","ID","Compound","Molecule"), cols)[1]
    } else if (length(char_cols) > 0) {
      char_cols[1]
    } else {
      "None"
    }
    tagList(
      selectInput("master_cluster_heatmap_variables", "Variables",
                  choices = numeric_cols,
                  selected = intersect(c("MW","TPSA","LogP","#H-bond acceptors",
                                          "#H-bond donors","#Rotatable bonds"),
                                        numeric_cols),
                  multiple = TRUE),
      selectInput("master_cluster_heatmap_id_col", "Label column (optional)",
                  choices = c("None", cols), selected = default_id),
      selectInput("master_cluster_heatmap_method", "Clustering method",
                  choices = c("Ward D2" = "ward.D2", "Ward D" = "ward.D",
                              "Complete" = "complete",
                              "Average (UPGMA)" = "average",
                              "Single" = "single", "McQuitty" = "mcquitty",
                              "Median" = "median", "Centroid" = "centroid"),
                  selected = "ward.D2"),
      checkboxInput("master_cluster_heatmap_scale",
                    "Scale variables (z-score)", TRUE)
    )
  })

  output$master_violin_controls <- renderUI({
    req(master_resultado())
    cols <- names(master_resultado())
    numeric_cols <- cols[sapply(master_resultado(), is.numeric)]
    categorical_cols <- cols[!sapply(master_resultado(), is.numeric)]
    if (length(numeric_cols) == 0) return(helpText("No numeric columns available."))
    default_group <- if ("GI absorption" %in% categorical_cols) {
      "GI absorption"
    } else if ("BBB permeant" %in% categorical_cols) {
      "BBB permeant"
    } else if (length(categorical_cols) > 0) {
      categorical_cols[1]
    } else {
      character(0)
    }
    tagList(
      selectInput("master_violin_variable", "Variable",
                  choices = numeric_cols,
                  selected = if ("MW" %in% numeric_cols) "MW" else numeric_cols[1]),
      selectInput("master_violin_group", "Group by",
                  choices = categorical_cols, selected = default_group),
      checkboxInput("master_violin_box", "Show boxplot", TRUE),
      checkboxInput("master_violin_points", "Show points", FALSE)
    )
  })

  output$master_histogram_custom_controls <- renderUI({
    req(master_resultado())
    cols <- names(master_resultado())
    numeric_cols <- cols[sapply(master_resultado(), is.numeric)]
    categorical_cols <- cols[!sapply(master_resultado(), is.numeric)]
    if (length(numeric_cols) == 0) return(helpText("No numeric columns available."))
    tagList(
      selectInput("master_hc_variable", "Variable",
                  choices = numeric_cols,
                  selected = if ("MW" %in% numeric_cols) "MW" else numeric_cols[1]),
      sliderInput("master_hc_bins", "Number of bins",
                  min = 5, max = 100, value = 30),
      selectInput("master_hc_group", "Group by (overlay)",
                  choices = c("None", categorical_cols), selected = "None"),
      checkboxInput("master_hc_density", "Show density curve", TRUE),
      checkboxInput("master_hc_rug", "Show rug plot", FALSE)
    )
  })

  ## ---------------- Step 4: plot rendering ----------------
  master_base_plot_types <- c("Radar plot (Chemical profile)",
                              "Tanimoto / AGNES (Structural similarity)")

  master_currentPlot <- reactive({
    req(master_resultado())
    validate(need(nrow(master_resultado()) > 0, "No data to plot."))
    pt <- input$master_plot_type

    if (pt %in% master_base_plot_types) {
      if (pt == "Radar plot (Chemical profile)") {
        req(input$master_radar_id_col)
        validate(need(length(input$master_radar_ids) >= 1,
                      "Select at least one molecule for the radar."))
        function() plotRadar(master_resultado(),
                             id_col = input$master_radar_id_col,
                             ids = input$master_radar_ids)
      } else {
        req(input$master_smiles_col)
        function() plotTanimoto(master_resultado(),
                                smiles_col = input$master_smiles_col,
                                label_col  = input$master_label_col,
                                max_n      = input$master_tanimoto_max_n,
                                method     = input$master_agnes_method)
      }
    } else {
      switch(pt,
        "Boiled Egg" = {
          logp_choice <- input$master_boiled_egg_logp
          if (is.null(logp_choice) || logp_choice == "")
            logp_choice <- "LogP"
          data_be <- master_resultado()
          if (logp_choice %in% names(data_be) && logp_choice != "LogP") {
            data_be$LogP <- data_be[[logp_choice]]
          }
          plotBoiledEgg(data_be, logp_source = logp_choice)
        },
        "Molecular Weight" = plotMW(master_resultado()),
        "TPSA"             = plotTPSA(master_resultado()),
        "LogP"             = plotLogP(master_resultado()),
        "Correlation Heatmap" = plotCorrHeatmap(master_resultado()),
        "Principal Component Analysys (PCA - Chemical space)" =
          apply_palette(plotPCA(data = master_resultado(),
                                variables  = input$master_pca_variables,
                                color_by   = input$master_pca_color,
                                label_by   = input$master_pca_label,
                                scale_data = input$master_pca_scale,
                                ellipse    = input$master_pca_ellipse),
                        input$master_palette, master_resultado(),
                        input$master_pca_color),
        "t-SNE (Chemical space)" =
          apply_palette(plotTSNE(data = master_resultado(),
                                 variables   = input$master_tsne_variables,
                                 color_by    = input$master_tsne_color,
                                 label_by    = input$master_tsne_label,
                                 perplexity  = input$master_tsne_perplexity,
                                 max_iter    = input$master_tsne_iter,
                                 scale_data  = input$master_tsne_scale),
                        input$master_palette, master_resultado(),
                        input$master_tsne_color),
        "UMAP (Chemical space)" =
          apply_palette(plotUMAP(data = master_resultado(),
                                 variables    = input$master_umap_variables,
                                 color_by     = input$master_umap_color,
                                 label_by     = input$master_umap_label,
                                 n_neighbors  = input$master_umap_n_neighbors,
                                 min_dist     = input$master_umap_min_dist,
                                 scale_data   = input$master_umap_scale),
                        input$master_palette, master_resultado(),
                        input$master_umap_color),
        "Parallel Coordinates" =
          apply_palette(plotParallel(data = master_resultado(),
                                     variables   = input$master_parallel_variables,
                                     color_by    = input$master_parallel_color,
                                     scale_data  = input$master_parallel_scale),
                        input$master_palette, master_resultado(),
                        input$master_parallel_color),
        "Cluster Heatmap (Dendrogram)" =
          function() plotClusterHeatmap(
            data       = master_resultado(),
            variables  = input$master_cluster_heatmap_variables,
            id_col     = .safe_id_col(input$master_cluster_heatmap_id_col),
            method     = input$master_cluster_heatmap_method,
            scale_data = input$master_cluster_heatmap_scale,
            palette    = input$master_palette
          ),
        "Violin Plot" =
          apply_palette(plotViolin(data = master_resultado(),
                                    variable    = input$master_violin_variable,
                                    group_by    = input$master_violin_group,
                                    show_box    = input$master_violin_box,
                                    show_points = input$master_violin_points),
                        input$master_palette, master_resultado(),
                        input$master_violin_group),
        "Custom Histogram" =
          plotHistogramCustom(data         = master_resultado(),
                              variable     = input$master_hc_variable,
                              bins         = input$master_hc_bins,
                              group_by     = input$master_hc_group,
                              show_density = input$master_hc_density,
                              show_rug     = input$master_hc_rug)
      )
    }
  })

  output$master_plot <- renderPlot({
    p <- master_currentPlot()
    if (is.function(p)) p() else print(p)
  })

  output$master_downloadPlot <- downloadHandler(
    filename = function() {
      paste0("admet_master_",
             gsub("[^A-Za-z0-9]+", "_", input$master_plot_type),
             ".", input$master_format)
    },
    content = function(file) {
      p <- master_currentPlot()
      if (is.function(p)) {
        open_device <- switch(input$master_format,
          "png"  = function() grDevices::png(file,
                              width  = input$master_width,
                              height = input$master_height,
                              units  = "in",
                              res    = as.numeric(input$master_dpi)),
          "jpeg" = function() grDevices::jpeg(file,
                               width  = input$master_width,
                               height = input$master_height,
                               units  = "in",
                               res    = as.numeric(input$master_dpi)),
          "tiff" = function() grDevices::tiff(file,
                               width  = input$master_width,
                               height = input$master_height,
                               units  = "in",
                               res    = as.numeric(input$master_dpi)),
          "pdf"  = function() grDevices::pdf(file,
                              width  = input$master_width,
                              height = input$master_height),
          "svg"  = function() grDevices::svg(file,
                              width  = input$master_width,
                              height = input$master_height),
          function() grDevices::png(file,
                            width  = input$master_width,
                            height = input$master_height,
                            units  = "in",
                            res    = as.numeric(input$master_dpi))
        )
        open_device()
        p()
        grDevices::dev.off()
      } else {
        ggplot2::ggsave(filename = file, plot = p,
                        device = input$master_format,
                        width  = input$master_width,
                        height = input$master_height,
                        units  = "in",
                        dpi    = as.numeric(input$master_dpi))
      }
    }
  )

  ## ============================================================
  ## Report Module
  ## ============================================================
  ##
  ## Collects the state of all four data-source modules (SwissADME,
  ## CDK & webchem, ADMETlab 3.0, Deep-PK) and generates a comprehensive
  ## report in HTML, PDF, or Word format. The report includes per-dataset
  ## statistics, drug-likeness filter results, BOILED-Egg ADMET
  ## classification, additional literature-supported metrics (Pfizer 3/75,
  ## GSK 4/400, lead-likeness), a composite drug-likeness score (0-100),
  ## cross-dataset comparison, visualizations and references.
  ##
  ## This section is purely additive: it reads from existing reactive
  ## values without modifying them.

  ## ---- Collect the state of all datasets ----
  report_state <- reactive({

    datasets <- list()

    ## CDK & webchem
    if (!is.null(cdk_results_rv())) {
      filt <- tryCatch(cdk_resultado(), error = function(e) NULL)
      datasets[["CDK"]] <- list(
        raw         = cdk_results_rv(),
        filtered    = if (!is.null(filt)) filt else cdk_results_rv(),
        filters     = if (!is.null(filt)) input$cdk_filters else character(0),
        source_name = "CDK & webchem"
      )
    }

    ## ADMET Master
    if (!is.null(master_mapped_rv())) {
      filt <- tryCatch(master_resultado(), error = function(e) NULL)
      datasets[["ADMET Master"]] <- list(
        raw         = master_mapped_rv(),
        filtered    = if (!is.null(filt)) filt else master_mapped_rv(),
        filters     = if (!is.null(filt)) input$master_filters else character(0),
        source_name = "ADMET Master Manager"
      )
    }

    datasets
  })

  ## ---- HTML preview (rendered on button click) ----
  report_html_rv <- reactiveVal(NULL)

  observeEvent(input$generate_report_btn, {

    st <- report_state()
    if (length(st) == 0) {
      showNotification("No datasets loaded. Please upload data in at least one module.",
                       type = "error", duration = 8)
      return(NULL)
    }

    withProgress(message = "Generating report preview...", value = 0.5, {
      tryCatch({
        tmpfile <- tempfile(fileext = ".html")
        renderReport(st, format = "html", output_file = tmpfile)
        report_html_rv(readLines(tmpfile, warn = FALSE))
        showNotification("Report preview generated successfully.",
                         type = "message", duration = 4)
      }, error = function(e) {
        showNotification(paste("Error generating report:", e$message),
                         type = "error", duration = 12)
      })
    })
  })

  ## ---- Preview container (placeholder or rendered HTML) ----
  output$report_preview_container <- renderUI({

    if (is.null(report_html_rv())) {
      return(tags$div(
        style = "text-align: center; padding: 80px; color: #999;",
        icon("file-lines", style = "font-size: 60px; margin-bottom: 20px;"),
        tags$h3("No report generated yet"),
        tags$p(style = "font-size: 16px;",
               "Click 'Generate Preview' to create a comprehensive report of your analyses."),
        br(),
        tags$p(style = "font-size: 13px; color: #aaa;",
               "The report includes statistics, drug-likeness scores, ADMET classification,",
               br(),
               "additional metrics, visualizations and cross-dataset comparison.")
      ))
    }

    tags$div(
      style = "width: 100%; height: 850px; overflow-y: auto; border: 1px solid #ddd; border-radius: 8px; padding: 20px; background: white;",
      HTML(paste(report_html_rv(), collapse = "\n"))
    )
  })

  ## ---- Download handlers ----

  output$download_report_pdf <- downloadHandler(
    filename = function() paste0("ADMETShiny_Report_", format(Sys.Date(), "%Y%m%d"), ".pdf"),
    content = function(file) {
      st <- report_state()
      if (length(st) == 0) {
        showNotification("No datasets loaded.", type = "error", duration = 8)
        return(NULL)
      }
      withProgress(message = "Generating PDF report...", value = 0.5, {
        tryCatch({
          renderReport(st, format = "pdf", output_file = file)
          showNotification("PDF report generated.", type = "message", duration = 4)
        }, error = function(e) {
          showNotification(paste("Error generating PDF:", e$message,
                                 "\n\nTip: Install LaTeX with tinytex::install_tinytex()"),
                           type = "error", duration = 15)
        })
      })
    }
  )

  output$download_report_doc <- downloadHandler(
    filename = function() paste0("ADMETShiny_Report_", format(Sys.Date(), "%Y%m%d"), ".docx"),
    content = function(file) {
      st <- report_state()
      if (length(st) == 0) {
        showNotification("No datasets loaded.", type = "error", duration = 8)
        return(NULL)
      }
      withProgress(message = "Generating Word report...", value = 0.5, {
        tryCatch({
          renderReport(st, format = "doc", output_file = file)
          showNotification("Word report generated.", type = "message", duration = 4)
        }, error = function(e) {
          showNotification(paste("Error generating Word document:", e$message),
                           type = "error", duration = 12)
        })
      })
    }
  )

  output$download_report_html <- downloadHandler(
    filename = function() paste0("ADMETShiny_Report_", format(Sys.Date(), "%Y%m%d"), ".html"),
    content = function(file) {
      st <- report_state()
      if (length(st) == 0) {
        showNotification("No datasets loaded.", type = "error", duration = 8)
        return(NULL)
      }
      withProgress(message = "Generating HTML report...", value = 0.5, {
        tryCatch({
          renderReport(st, format = "html", output_file = file)
          showNotification("HTML report generated.", type = "message", duration = 4)
        }, error = function(e) {
          showNotification(paste("Error generating HTML:", e$message),
                           type = "error", duration = 12)
        })
      })
    }
  )
}
