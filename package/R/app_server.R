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
  observeEvent(input$go_swissadme,  updateNavbarPage(session, "main_nav", selected = "swissadme"))
  observeEvent(input$go_cdk,        updateNavbarPage(session, "main_nav", selected = "cdk"))
  observeEvent(input$go_admetlab,   updateNavbarPage(session, "main_nav", selected = "admetlab"))
  observeEvent(input$go_deeppk,     updateNavbarPage(session, "main_nav", selected = "deeppk"))
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

  ## ---- Step 4: send CDK dataset to SwissADME module ----
  observeEvent(input$send_to_filter, {
    req(cdk_results_rv())
    datos_rv(cdk_results_rv())
    updateNavbarPage(session, "main_nav", selected = "swissadme")
    showNotification("CDK/webchem dataset loaded into the SwissADME Filtering module.",
                     type = "message", duration = 6)
  })

  ## ============================ SwissADME ================================
  datos_rv <- reactiveVal(NULL)

  observeEvent(input$archivo, {
    req(input$archivo)
    tryCatch({
      d <- read.csv(input$archivo$datapath, check.names = FALSE)
      datos_rv(fixSwissADME(d))
    }, error = function(e) {
      showNotification(paste("Error reading the file:", e$message),
                       type = "error", duration = 8)
    })
  })

  datos <- reactive({ datos_rv() })

  output$preview <- renderDT({ req(datos()); datos() })

  resultado <- eventReactive(input$run, {
    req(datos())
    tryCatch({
      out <- applyFilters(
        data = datos(), filters = input$filters,
        lipinski = list(mw = input$mw, logp = input$logp, hba = input$hba, hbd = input$hbd, violations = input$violations),
        veber = list(v_rb = input$v_rb, v_tpsa = input$v_tpsa, v_hb_sum = input$v_hb_sum, violations = input$veber_violations),
        ghose = list(g_mw_min = input$g_mw[1], g_mw_max = input$g_mw[2], g_mr_min = input$g_mr[1], g_mr_max = input$g_mr[2], g_logp_min = input$g_logp[1], g_logp_max = input$g_logp[2], g_ha_min = input$g_ha[1], g_ha_max = input$g_ha[2], violations = input$ghose_violations),
        egan = list(e_tpsa = input$e_tpsa, e_logp = input$e_logp, violations = input$egan_violations),
        muegge = list(m_mw_min = input$m_mw[1], m_mw_max = input$m_mw[2], m_logp_min = input$m_logp[1], m_logp_max = input$m_logp[2], m_hba = input$m_hba, m_hbd = input$m_hbd, m_rb = input$m_rb, m_tpsa = input$m_tpsa, violations = input$muegge_violations)
      )
      if (nrow(out) == 0)
        showNotification("No compound meets the selected filters.",
                         type = "warning", duration = 6)
      out
    }, error = function(e) {
      showNotification(paste("Error applying filters:", e$message),
                       type = "error", duration = 8)
      datos()[0, ]
    })
  })

  output$tabla <- renderDT({
    req(resultado())
    resultado()
  }, options = list(pageLength = 10, scrollX = TRUE))

  output$download <- downloadHandler(
    filename = function() "Filtered_Dataset.csv",
    content  = function(file) write.csv(resultado(), file, row.names = FALSE)
  )

  ## ------------- Dynamic controls: Radar & Tanimoto -------------------
  output$radar_id_selector <- renderUI({
    req(resultado())
    tagList(
      selectInput("radar_id_col", "Identifier column",
                  choices = names(resultado()),
                  selected = names(resultado())[1]),
      selectizeInput("radar_ids", "Molecules to compare (max 5)",
                     choices = NULL, multiple = TRUE,
                     options = list(maxItems = 5, placeholder = "Select up to 5 molecules"))
    )
  })

  observeEvent(input$radar_id_col, {
    req(resultado(), input$radar_id_col)
    vals <- unique(as.character(resultado()[[input$radar_id_col]]))
    updateSelectizeInput(session, "radar_ids", choices = vals, server = TRUE)
  })

  output$smiles_col_selector <- renderUI({
    req(resultado())
    cols <- names(resultado())
    guess <- cols[grepl("smiles", cols, ignore.case = TRUE)]
    default <- if (length(guess) > 0) guess[1] else cols[1]
    selectInput("smiles_col", "SMILES column", choices = cols, selected = default)
  })

  output$label_col_selector <- renderUI({
    req(resultado())
    cols <- names(resultado())
    default <- if ("Molecule" %in% cols) "Molecule" else cols[1]
    selectInput("label_col", "Label column", choices = cols, selected = default)
  })

  ## ------------- Dynamic controls: PCA / t-SNE / Parallel / Violin -----
  output$pca_controls <- renderUI({
    req(resultado())
    cols <- names(resultado())
    numeric_cols <- cols[sapply(resultado(), is.numeric)]
    tagList(
      selectInput("pca_variables", "Variables", choices = numeric_cols,
                  selected = c("MW","TPSA","LogP","#H-bond acceptors","#H-bond donors","#Rotatable bonds"),
                  multiple = TRUE),
      selectInput("pca_color", "Colour by", choices = c("None", cols), selected = "GI absorption"),
      selectInput("pca_label", "Labels", choices = c("None", cols), selected = "Molecule"),
      checkboxInput("pca_ellipse", "Show ellipses", TRUE),
      checkboxInput("pca_scale", "Scale variables", TRUE)
    )
  })

  output$tsne_controls <- renderUI({
    req(resultado())
    cols <- names(resultado())
    numeric_cols <- cols[sapply(resultado(), is.numeric)]
    default_color <- if ("GI absorption" %in% cols) "GI absorption" else "None"
    default_label <- if ("Molecule" %in% cols) "Molecule" else "None"
    tagList(
      selectInput("tsne_variables", "Variables", choices = numeric_cols,
                  selected = c("MW","TPSA","LogP","#H-bond acceptors","#H-bond donors","#Rotatable bonds"),
                  multiple = TRUE),
      sliderInput("tsne_perplexity", "Perplexity", min = 5, max = 50, value = 30),
      sliderInput("tsne_iter", "Iterations", min = 500, max = 5000, value = 1000, step = 500),
      selectInput("tsne_color", "Colour by", choices = c("None", cols), selected = default_color),
      selectInput("tsne_label", "Labels", choices = c("None", cols), selected = default_label),
      checkboxInput("tsne_scale", "Scale variables", TRUE)
    )
  })

  output$umap_controls <- renderUI({
    req(resultado())
    cols <- names(resultado())
    numeric_cols <- cols[sapply(resultado(), is.numeric)]
    default_color <- if ("GI absorption" %in% cols) "GI absorption" else "None"
    default_label <- if ("Molecule" %in% cols) "Molecule" else "None"
    tagList(
      selectInput("umap_variables", "Variables", choices = numeric_cols,
                  selected = c("MW","TPSA","LogP","#H-bond acceptors","#H-bond donors","#Rotatable bonds"),
                  multiple = TRUE),
      sliderInput("umap_n_neighbors", "n_neighbors", min = 2, max = 50, value = 15),
      sliderInput("umap_min_dist", "min_dist", min = 0, max = 1, value = 0.1, step = 0.05),
      selectInput("umap_color", "Colour by", choices = c("None", cols), selected = default_color),
      selectInput("umap_label", "Labels", choices = c("None", cols), selected = default_label),
      checkboxInput("umap_scale", "Scale variables", TRUE)
    )
  })

  output$parallel_controls <- renderUI({
    req(resultado())
    cols <- names(resultado())
    numeric_cols <- cols[sapply(resultado(), is.numeric)]
    tagList(
      selectInput("parallel_variables", "Variables", choices = numeric_cols,
                  selected = c("MW","TPSA","LogP","#H-bond acceptors","#H-bond donors","#Rotatable bonds"),
                  multiple = TRUE),
      selectInput("parallel_color", "Colour by", choices = c("None", cols), selected = "GI absorption"),
      checkboxInput("parallel_scale", "Scale variables", TRUE)
    )
  })

  output$cluster_heatmap_controls <- renderUI({
    req(resultado())
    cols <- names(resultado())
    numeric_cols <- cols[sapply(resultado(), is.numeric)]
    ## Look for a suitable label column: prefer text columns that could be IDs
    char_cols <- cols[!sapply(resultado(), is.numeric)]
    default_id <- if (any(c("Name", "name", "ID", "Compound", "Molecule") %in% cols)) {
      intersect(c("Name", "name", "ID", "Compound", "Molecule"), cols)[1]
    } else if (length(char_cols) > 0) {
      char_cols[1]
    } else {
      "None"
    }
    tagList(
      selectInput("cluster_heatmap_variables", "Variables", choices = numeric_cols,
                  selected = c("MW","TPSA","LogP","#H-bond acceptors","#H-bond donors","#Rotatable bonds"),
                  multiple = TRUE),
      selectInput("cluster_heatmap_id_col", "Label column (optional)", choices = c("None", cols),
                  selected = default_id),
      selectInput("cluster_heatmap_method", "Clustering method",
                  choices = c("Ward D2" = "ward.D2", "Ward D" = "ward.D",
                              "Complete" = "complete", "Average (UPGMA)" = "average",
                              "Single" = "single", "McQuitty" = "mcquitty",
                              "Median" = "median", "Centroid" = "centroid"),
                  selected = "ward.D2"),
      checkboxInput("cluster_heatmap_scale", "Scale variables (z-score)", TRUE)
    )
  })

  output$violin_controls <- renderUI({
    req(resultado())
    cols <- names(resultado())
    numeric_cols <- cols[sapply(resultado(), is.numeric)]
    categorical_cols <- cols[!sapply(resultado(), is.numeric)]
    tagList(
      selectInput("violin_variable", "Variable", choices = numeric_cols, selected = "MW"),
      selectInput("violin_group", "Group by", choices = categorical_cols, selected = "GI absorption"),
      checkboxInput("violin_box", "Show boxplot", TRUE),
      checkboxInput("violin_points", "Show points", FALSE)
    )
  })

  ## ---- SwissADME: renderUI blocks for the new plot types ----
  output$histogram_custom_controls <- renderUI({
    req(resultado())
    cols <- names(resultado())
    numeric_cols <- cols[sapply(resultado(), is.numeric)]
    categorical_cols <- cols[!sapply(resultado(), is.numeric)]
    tagList(
      selectInput("hc_variable", "Variable", choices = numeric_cols, selected = "MW"),
      sliderInput("hc_bins", "Number of bins", min = 5, max = 100, value = 30),
      selectInput("hc_group", "Group by (overlay)", choices = c("None", categorical_cols),
                  selected = "None"),
      checkboxInput("hc_density", "Show density curve", TRUE),
      checkboxInput("hc_rug", "Show rug plot", FALSE)
    )
  })

  ## ----------------------------- SwissADME plot ------------------------
  base_plot_types <- c("Radar plot (Chemical profile)",
                       "Tanimoto / AGNES (Structural similarity)")

  currentPlot <- reactive({
    req(resultado())
    validate(need(nrow(resultado()) > 0, "No data to plot."))
    pt <- input$plot_type

    if (pt %in% base_plot_types) {
      if (pt == "Radar plot (Chemical profile)") {
        req(input$radar_id_col)
        validate(need(length(input$radar_ids) >= 1, "Select at least one molecule."))
        function() plotRadar(resultado(),
                             id_col = input$radar_id_col, ids = input$radar_ids)
      } else {
        req(input$smiles_col)
        function() plotTanimoto(resultado(),
                                smiles_col = input$smiles_col,
                                label_col = input$label_col,
                                max_n = input$tanimoto_max_n,
                                method = input$agnes_method)
      }
    } else {
      switch(pt,
             "Boiled Egg" = {
               logp_choice <- input$boiled_egg_logp
               if (is.null(logp_choice) || logp_choice == "") logp_choice <- "LogP"
               data_be <- resultado()
               if (logp_choice %in% names(data_be) && logp_choice != "LogP") {
                 data_be$LogP <- data_be[[logp_choice]]
               }
               plotBoiledEgg(data_be)
             },
             "Molecular Weight" = plotMW(resultado()),
             "TPSA" = plotTPSA(resultado()),
             "LogP" = plotLogP(resultado()),
             "Correlation Heatmap" = plotCorrHeatmap(resultado()),
             "Principal Component Analysys (PCA - Chemical space)" =
               apply_palette(plotPCA(data = resultado(), variables = input$pca_variables,
                                     color_by = input$pca_color, label_by = input$pca_label,
                                     scale_data = input$pca_scale, ellipse = input$pca_ellipse),
                             input$swiss_palette, resultado(), input$pca_color),
             "t-SNE (Chemical space)" =
               apply_palette(plotTSNE(data = resultado(), variables = input$tsne_variables,
                                      color_by = input$tsne_color, label_by = input$tsne_label,
                                      perplexity = input$tsne_perplexity, max_iter = input$tsne_iter,
                                      scale_data = input$tsne_scale),
                             input$swiss_palette, resultado(), input$tsne_color),
             "UMAP (Chemical space)" =
               apply_palette(plotUMAP(data = resultado(), variables = input$umap_variables,
                                      color_by = input$umap_color, label_by = input$umap_label,
                                      n_neighbors = input$umap_n_neighbors,
                                      min_dist = input$umap_min_dist, scale_data = input$umap_scale),
                             input$swiss_palette, resultado(), input$umap_color),
             "Parallel Coordinates" =
               apply_palette(plotParallel(data = resultado(), variables = input$parallel_variables,
                                          color_by = input$parallel_color, scale_data = input$parallel_scale),
                             input$swiss_palette, resultado(), input$parallel_color),
             "Cluster Heatmap (Dendrogram)" =
               function() plotClusterHeatmap(data = resultado(),
                                             variables = input$cluster_heatmap_variables,
                                             id_col = .safe_id_col(input$cluster_heatmap_id_col),
                                             method = input$cluster_heatmap_method,
                                             scale_data = input$cluster_heatmap_scale,
                                             palette = input$swiss_palette),
             "Violin Plot" =
               apply_palette(plotViolin(data = resultado(), variable = input$violin_variable,
                                         group_by = input$violin_group, show_box = input$violin_box,
                                         show_points = input$violin_points),
                             input$swiss_palette, resultado(), input$violin_group),
             "Custom Histogram" =
               plotHistogramCustom(data = resultado(),
                                   variable = input$hc_variable,
                                   bins = input$hc_bins,
                                   group_by = input$hc_group,
                                   show_density = input$hc_density,
                                   show_rug = input$hc_rug)
      )
    }
  })

  output$plot <- renderPlot({
    p <- currentPlot()
    if (is.function(p)) p() else print(p)
  })

  output$downloadPlot <- downloadHandler(
    filename = function() paste0(gsub("[^A-Za-z0-9]+", "_", input$plot_type), ".", input$format),
    content = function(file) {
      p <- currentPlot()
      if (is.function(p)) {
        open_device <- switch(input$format,
          "png"  = function() grDevices::png(file, width = input$width, height = input$height, units = "in", res = as.numeric(input$dpi)),
          "jpeg" = function() grDevices::jpeg(file, width = input$width, height = input$height, units = "in", res = as.numeric(input$dpi)),
          "tiff" = function() grDevices::tiff(file, width = input$width, height = input$height, units = "in", res = as.numeric(input$dpi)),
          "pdf"  = function() grDevices::pdf(file, width = input$width, height = input$height),
          "svg"  = function() grDevices::svg(file, width = input$width, height = input$height),
          function() grDevices::png(file, width = input$width, height = input$height, units = "in", res = as.numeric(input$dpi))
        )
        open_device()
        p()
        grDevices::dev.off()
      } else {
        ggplot2::ggsave(filename = file, plot = p, device = input$format,
                        width = input$width, height = input$height,
                        units = "in", dpi = as.numeric(input$dpi))
      }
    }
  )

  ## ============================ ADMETlab 3.0 =============================
  admetlab_datos_rv <- reactiveVal(NULL)

  observeEvent(input$admetlab_archivo, {
    req(input$admetlab_archivo)
    tryCatch({
      d <- read.csv(input$admetlab_archivo$datapath, check.names = FALSE)
      withProgress(message = "Processing ADMETlab dataset (calculating CDK descriptors)...",
                   value = 0.5, {
        admetlab_datos_rv(fixADMETlab(d))
      })
      showNotification("ADMETlab dataset loaded and processed.",
                       type = "message", duration = 5)
    }, error = function(e) {
      showNotification(paste("Error reading ADMETlab:", e$message),
                       type = "error", duration = 8)
    })
  })

  admetlab_datos <- reactive({ admetlab_datos_rv() })

  output$admetlab_preview <- renderDT({
    req(admetlab_datos())
    admetlab_datos()
  }, options = list(pageLength = 10, scrollX = TRUE))

  admetlab_resultado <- eventReactive(input$admetlab_run, {
    req(admetlab_datos())
    tryCatch({
      out <- applyFilters(
        data = admetlab_datos(), filters = input$admetlab_filters,
        lipinski = list(mw = input$admetlab_mw, logp = input$admetlab_logp, hba = input$admetlab_hba, hbd = input$admetlab_hbd, violations = input$admetlab_violations),
        veber = list(v_rb = input$admetlab_v_rb, v_tpsa = input$admetlab_v_tpsa, v_hb_sum = input$admetlab_v_hb_sum, violations = input$admetlab_veber_violations),
        ghose = list(g_mw_min = input$admetlab_g_mw[1], g_mw_max = input$admetlab_g_mw[2], g_mr_min = input$admetlab_g_mr[1], g_mr_max = input$admetlab_g_mr[2], g_logp_min = input$admetlab_g_logp[1], g_logp_max = input$admetlab_g_logp[2], g_ha_min = input$admetlab_g_ha[1], g_ha_max = input$admetlab_g_ha[2], violations = input$admetlab_ghose_violations),
        egan = list(e_tpsa = input$admetlab_e_tpsa, e_logp = input$admetlab_e_logp, violations = input$admetlab_egan_violations),
        muegge = list(m_mw_min = input$admetlab_m_mw[1], m_mw_max = input$admetlab_m_mw[2], m_logp_min = input$admetlab_m_logp[1], m_logp_max = input$admetlab_m_logp[2], m_hba = input$admetlab_m_hba, m_hbd = input$admetlab_m_hbd, m_rb = input$admetlab_m_rb, m_tpsa = input$admetlab_m_tpsa, violations = input$admetlab_muegge_violations)
      )
      if (nrow(out) == 0) {
        showNotification("No compound meets the selected filters.",
                         type = "warning", duration = 6)
      }
      out
    }, error = function(e) {
      showNotification(paste("Error applying filters:", e$message),
                       type = "error", duration = 8)
      admetlab_datos()[0, ]
    })
  })

  output$admetlab_tabla <- renderDT({
    req(admetlab_resultado())
    admetlab_resultado()
  }, options = list(pageLength = 10, scrollX = TRUE))

  output$admetlab_download <- downloadHandler(
    filename = function() "ADMETlab_Filtered_Dataset.csv",
    content = function(file) write.csv(admetlab_resultado(), file, row.names = FALSE)
  )

  ## ---- ADMETlab dynamic controls ----
  output$admetlab_radar_id_selector <- renderUI({
    req(admetlab_resultado())
    tagList(
      selectInput("admetlab_radar_id_col", "Identifier column",
                  choices = names(admetlab_resultado()),
                  selected = names(admetlab_resultado())[1]),
      selectizeInput("admetlab_radar_ids", "Molecules to compare",
                     choices = NULL, multiple = TRUE, options = list(maxItems = 5))
    )
  })

  observeEvent(input$admetlab_radar_id_col, {
    req(admetlab_resultado(), input$admetlab_radar_id_col)
    updateSelectizeInput(session, "admetlab_radar_ids",
                         choices = unique(as.character(admetlab_resultado()[[input$admetlab_radar_id_col]])),
                         server = TRUE)
  })

  output$admetlab_smiles_col_selector <- renderUI({
    req(admetlab_resultado())
    cols <- names(admetlab_resultado())
    guess <- cols[grepl("smiles", cols, ignore.case = TRUE)]
    selectInput("admetlab_smiles_col", "SMILES column", choices = cols,
                selected = if (length(guess) > 0) guess[1] else cols[1])
  })

  output$admetlab_label_col_selector <- renderUI({
    req(admetlab_resultado())
    cols <- names(admetlab_resultado())
    selectInput("admetlab_label_col", "Label column", choices = cols, selected = cols[1])
  })

  output$admetlab_pca_controls <- renderUI({
    req(admetlab_resultado())
    cols <- names(admetlab_resultado())
    numeric_cols <- cols[sapply(admetlab_resultado(), is.numeric)]
    tagList(
      selectInput("admetlab_pca_variables", "Variables", choices = numeric_cols,
                  selected = intersect(c("MW","TPSA","LogP","#H-bond acceptors","#H-bond donors","#Rotatable bonds"), numeric_cols),
                  multiple = TRUE),
      selectInput("admetlab_pca_color", "Colour by", choices = c("None", cols), selected = "None"),
      selectInput("admetlab_pca_label", "Labels", choices = c("None", cols), selected = "None"),
      checkboxInput("admetlab_pca_ellipse", "Ellipses", TRUE),
      checkboxInput("admetlab_pca_scale", "Scale", TRUE)
    )
  })

  output$admetlab_tsne_controls <- renderUI({
    req(admetlab_resultado())
    cols <- names(admetlab_resultado())
    numeric_cols <- cols[sapply(admetlab_resultado(), is.numeric)]
    tagList(
      selectInput("admetlab_tsne_variables", "Variables", choices = numeric_cols,
                  selected = intersect(c("MW","TPSA","LogP","#H-bond acceptors","#H-bond donors","#Rotatable bonds"), numeric_cols),
                  multiple = TRUE),
      sliderInput("admetlab_tsne_perplexity", "Perplexity", min = 5, max = 50, value = 30),
      sliderInput("admetlab_tsne_iter", "Iterations", min = 500, max = 5000, value = 1000, step = 500),
      selectInput("admetlab_tsne_color", "Colour by", choices = c("None", cols), selected = "None"),
      selectInput("admetlab_tsne_label", "Labels", choices = c("None", cols), selected = "None"),
      checkboxInput("admetlab_tsne_scale", "Scale", TRUE)
    )
  })

  output$admetlab_umap_controls <- renderUI({
    req(admetlab_resultado())
    cols <- names(admetlab_resultado())
    numeric_cols <- cols[sapply(admetlab_resultado(), is.numeric)]
    tagList(
      selectInput("admetlab_umap_variables", "Variables", choices = numeric_cols,
                  selected = intersect(c("MW","TPSA","LogP","#H-bond acceptors","#H-bond donors","#Rotatable bonds"), numeric_cols),
                  multiple = TRUE),
      sliderInput("admetlab_umap_n_neighbors", "n_neighbors", min = 2, max = 50, value = 15),
      sliderInput("admetlab_umap_min_dist", "min_dist", min = 0, max = 1, value = 0.1, step = 0.05),
      selectInput("admetlab_umap_color", "Colour by", choices = c("None", cols), selected = "None"),
      selectInput("admetlab_umap_label", "Labels", choices = c("None", cols), selected = "None"),
      checkboxInput("admetlab_umap_scale", "Scale", TRUE)
    )
  })

  output$admetlab_parallel_controls <- renderUI({
    req(admetlab_resultado())
    cols <- names(admetlab_resultado())
    numeric_cols <- cols[sapply(admetlab_resultado(), is.numeric)]
    tagList(
      selectInput("admetlab_parallel_variables", "Variables", choices = numeric_cols,
                  selected = intersect(c("MW","TPSA","LogP","#H-bond acceptors","#H-bond donors","#Rotatable bonds"), numeric_cols),
                  multiple = TRUE),
      selectInput("admetlab_parallel_color", "Colour by", choices = c("None", cols), selected = "None"),
      checkboxInput("admetlab_parallel_scale", "Scale variables", TRUE)
    )
  })

  output$admetlab_cluster_heatmap_controls <- renderUI({
    req(admetlab_resultado())
    cols <- names(admetlab_resultado())
    numeric_cols <- cols[sapply(admetlab_resultado(), is.numeric)]
    ## Look for a suitable label column: prefer text columns that could be IDs
    char_cols <- cols[!sapply(admetlab_resultado(), is.numeric)]
    default_id <- if (any(c("Name", "name", "ID", "Compound", "Molecule") %in% cols)) {
      intersect(c("Name", "name", "ID", "Compound", "Molecule"), cols)[1]
    } else if (length(char_cols) > 0) {
      char_cols[1]
    } else {
      "None"
    }
    tagList(
      selectInput("admetlab_cluster_heatmap_variables", "Variables", choices = numeric_cols,
                  selected = c("MW","TPSA","LogP","#H-bond acceptors","#H-bond donors","#Rotatable bonds"),
                  multiple = TRUE),
      selectInput("admetlab_cluster_heatmap_id_col", "Label column (optional)", choices = c("None", cols),
                  selected = default_id),
      selectInput("admetlab_cluster_heatmap_method", "Clustering method",
                  choices = c("Ward D2" = "ward.D2", "Ward D" = "ward.D",
                              "Complete" = "complete", "Average (UPGMA)" = "average",
                              "Single" = "single", "McQuitty" = "mcquitty",
                              "Median" = "median", "Centroid" = "centroid"),
                  selected = "ward.D2"),
      checkboxInput("admetlab_cluster_heatmap_scale", "Scale variables (z-score)", TRUE)
    )
  })

  output$admetlab_violin_controls <- renderUI({
    req(admetlab_resultado())
    cols <- names(admetlab_resultado())
    numeric_cols <- cols[sapply(admetlab_resultado(), is.numeric)]
    categorical_cols <- cols[!sapply(admetlab_resultado(), is.numeric)]
    tagList(
      selectInput("admetlab_violin_variable", "Variable", choices = numeric_cols, selected = "MW"),
      selectInput("admetlab_violin_group", "Group by", choices = categorical_cols,
                  selected = if ("BBB permeant" %in% categorical_cols) "BBB permeant" else categorical_cols[1]),
      checkboxInput("admetlab_violin_box", "Boxplot", TRUE),
      checkboxInput("admetlab_violin_points", "Points", FALSE)
    )
  })

  ## ---- ADMETlab: renderUI blocks for the new plot types ----
  output$admetlab_histogram_custom_controls <- renderUI({
    req(admetlab_resultado())
    cols <- names(admetlab_resultado())
    numeric_cols <- cols[sapply(admetlab_resultado(), is.numeric)]
    categorical_cols <- cols[!sapply(admetlab_resultado(), is.numeric)]
    tagList(
      selectInput("admetlab_hc_variable", "Variable", choices = numeric_cols, selected = "MW"),
      sliderInput("admetlab_hc_bins", "Number of bins", min = 5, max = 100, value = 30),
      selectInput("admetlab_hc_group", "Group by (overlay)", choices = c("None", categorical_cols), selected = "None"),
      checkboxInput("admetlab_hc_density", "Show density curve", TRUE),
      checkboxInput("admetlab_hc_rug", "Show rug plot", FALSE)
    )
  })

  ## ---- ADMETlab plot rendering ----
  admetlab_base_plot_types <- c("Radar plot (Chemical profile)",
                                "Tanimoto / AGNES (Structural similarity)")

  admetlab_currentPlot <- reactive({
    req(admetlab_resultado())
    validate(need(nrow(admetlab_resultado()) > 0, "No data."))
    pt <- input$admetlab_plot_type

    if (pt %in% admetlab_base_plot_types) {
      if (pt == "Radar plot (Chemical profile)") {
        req(input$admetlab_radar_id_col)
        function() plotRadar(admetlab_resultado(),
                             id_col = input$admetlab_radar_id_col,
                             ids = input$admetlab_radar_ids)
      } else {
        req(input$admetlab_smiles_col)
        function() plotTanimoto(admetlab_resultado(),
                                smiles_col = input$admetlab_smiles_col,
                                label_col = input$admetlab_label_col,
                                max_n = input$admetlab_tanimoto_max_n,
                                method = input$admetlab_agnes_method)
      }
    } else {
      switch(pt,
             "Boiled Egg" = plotBoiledEgg(admetlab_resultado()),
             "TPSA" = plotTPSA(admetlab_resultado()),
             "LogP" = plotLogP(admetlab_resultado()),
             "Correlation Heatmap" = plotCorrHeatmap(admetlab_resultado()),
             "Principal Component Analysys (PCA - Chemical space)" =
               apply_palette(plotPCA(data = admetlab_resultado(), variables = input$admetlab_pca_variables,
                                     color_by = input$admetlab_pca_color, label_by = input$admetlab_pca_label,
                                     scale_data = input$admetlab_pca_scale, ellipse = input$admetlab_pca_ellipse),
                             input$admetlab_palette, admetlab_resultado(), input$admetlab_pca_color),
             "t-SNE (Chemical space)" =
               apply_palette(plotTSNE(data = admetlab_resultado(), variables = input$admetlab_tsne_variables,
                                      color_by = input$admetlab_tsne_color, label_by = input$admetlab_tsne_label,
                                      perplexity = input$admetlab_tsne_perplexity, max_iter = input$admetlab_tsne_iter,
                                      scale_data = input$admetlab_tsne_scale),
                             input$admetlab_palette, admetlab_resultado(), input$admetlab_tsne_color),
             "UMAP (Chemical space)" =
               apply_palette(plotUMAP(data = admetlab_resultado(), variables = input$admetlab_umap_variables,
                                      color_by = input$admetlab_umap_color, label_by = input$admetlab_umap_label,
                                      n_neighbors = input$admetlab_umap_n_neighbors,
                                      min_dist = input$admetlab_umap_min_dist,
                                      scale_data = input$admetlab_umap_scale),
                             input$admetlab_palette, admetlab_resultado(), input$admetlab_umap_color),
             "Parallel Coordinates" =
               apply_palette(plotParallel(data = admetlab_resultado(), variables = input$admetlab_parallel_variables,
                                          color_by = input$admetlab_parallel_color, scale_data = input$admetlab_parallel_scale),
                             input$admetlab_palette, admetlab_resultado(), input$admetlab_parallel_color),
             "Cluster Heatmap (Dendrogram)" =
               function() plotClusterHeatmap(data = admetlab_resultado(),
                                             variables = input$admetlab_cluster_heatmap_variables,
                                             id_col = .safe_id_col(input$admetlab_cluster_heatmap_id_col),
                                             method = input$admetlab_cluster_heatmap_method,
                                             scale_data = input$admetlab_cluster_heatmap_scale,
                                             palette = input$admetlab_palette),
             "Violin Plot" =
               apply_palette(plotViolin(data = admetlab_resultado(), variable = input$admetlab_violin_variable,
                                         group_by = input$admetlab_violin_group, show_box = input$admetlab_violin_box,
                                         show_points = input$admetlab_violin_points),
                             input$admetlab_palette, admetlab_resultado(), input$admetlab_violin_group),
             "Custom Histogram" =
               plotHistogramCustom(data = admetlab_resultado(),
                                   variable = input$admetlab_hc_variable,
                                   bins = input$admetlab_hc_bins,
                                   group_by = input$admetlab_hc_group,
                                   show_density = input$admetlab_hc_density,
                                   show_rug = input$admetlab_hc_rug)
      )
    }
  })

  output$admetlab_plot <- renderPlot({
    p <- admetlab_currentPlot()
    if (is.function(p)) p() else print(p)
  })

  output$admetlab_downloadPlot <- downloadHandler(
    filename = function() {
      paste0("admetlab_", gsub("[^A-Za-z0-9]+", "_", input$admetlab_plot_type), ".", input$admetlab_format)
    },
    content = function(file) {
      p <- admetlab_currentPlot()
      if (is.function(p)) {
        open_device <- switch(input$admetlab_format,
          "png" = function() grDevices::png(file, width = input$admetlab_width, height = input$admetlab_height, units = "in", res = as.numeric(input$admetlab_dpi)),
          "pdf" = function() grDevices::pdf(file, width = input$admetlab_width, height = input$admetlab_height),
          function() grDevices::png(file, width = input$admetlab_width, height = input$admetlab_height, units = "in", res = as.numeric(input$admetlab_dpi))
        )
        open_device()
        p()
        grDevices::dev.off()
      } else {
        ggplot2::ggsave(filename = file, plot = p, device = input$admetlab_format,
                        width = input$admetlab_width, height = input$admetlab_height,
                        units = "in", dpi = as.numeric(input$admetlab_dpi))
      }
    }
  )

  ## ============================ Deep-PK ==================================
  deeppk_datos_rv <- reactiveVal(NULL)

  observeEvent(input$deeppk_archivo, {
    req(input$deeppk_archivo)
    tryCatch({
      d <- read.csv(input$deeppk_archivo$datapath, check.names = FALSE)
      withProgress(message = "Processing Deep-PK dataset (calculating CDK descriptors)...",
                   value = 0.5, {
        deeppk_datos_rv(fixDeepPK(d))
      })
      showNotification("Deep-PK dataset loaded and processed.",
                       type = "message", duration = 5)
    }, error = function(e) {
      showNotification(paste("Error reading Deep-PK:", e$message),
                       type = "error", duration = 8)
    })
  })

  deeppk_datos <- reactive({ deeppk_datos_rv() })

  output$deeppk_preview <- renderDT({
    req(deeppk_datos())
    deeppk_datos()
  }, options = list(pageLength = 10, scrollX = TRUE))

  deeppk_resultado <- eventReactive(input$deeppk_run, {
    req(deeppk_datos())
    tryCatch({
      out <- applyFilters(
        data = deeppk_datos(), filters = input$deeppk_filters,
        lipinski = list(mw = input$deeppk_mw, logp = input$deeppk_logp, hba = input$deeppk_hba, hbd = input$deeppk_hbd, violations = input$deeppk_violations),
        veber = list(v_rb = input$deeppk_v_rb, v_tpsa = input$deeppk_v_tpsa, v_hb_sum = input$deeppk_v_hb_sum, violations = input$deeppk_veber_violations),
        ghose = list(g_mw_min = input$deeppk_g_mw[1], g_mw_max = input$deeppk_g_mw[2], g_mr_min = input$deeppk_g_mr[1], g_mr_max = input$deeppk_g_mr[2], g_logp_min = input$deeppk_g_logp[1], g_logp_max = input$deeppk_g_logp[2], g_ha_min = input$deeppk_g_ha[1], g_ha_max = input$deeppk_g_ha[2], violations = input$deeppk_ghose_violations),
        egan = list(e_tpsa = input$deeppk_e_tpsa, e_logp = input$deeppk_e_logp, violations = input$deeppk_egan_violations),
        muegge = list(m_mw_min = input$deeppk_m_mw[1], m_mw_max = input$deeppk_m_mw[2], m_logp_min = input$deeppk_m_logp[1], m_logp_max = input$deeppk_m_logp[2], m_hba = input$deeppk_m_hba, m_hbd = input$deeppk_m_hbd, m_rb = input$deeppk_m_rb, m_tpsa = input$deeppk_m_tpsa, violations = input$deeppk_muegge_violations)
      )
      if (nrow(out) == 0) {
        showNotification("No compound meets the selected filters.",
                         type = "warning", duration = 6)
      }
      out
    }, error = function(e) {
      showNotification(paste("Error applying filters:", e$message),
                       type = "error", duration = 8)
      deeppk_datos()[0, ]
    })
  })

  output$deeppk_tabla <- renderDT({
    req(deeppk_resultado())
    deeppk_resultado()
  }, options = list(pageLength = 10, scrollX = TRUE))

  output$deeppk_download <- downloadHandler(
    filename = function() "DeepPK_Filtered.csv",
    content = function(file) write.csv(deeppk_resultado(), file, row.names = FALSE)
  )

  ## ---- Deep-PK dynamic controls ----
  output$deeppk_radar_id_selector <- renderUI({
    req(deeppk_resultado())
    tagList(
      selectInput("deeppk_radar_id_col", "ID", choices = names(deeppk_resultado()),
                  selected = names(deeppk_resultado())[1]),
      selectizeInput("deeppk_radar_ids", "Molecules", choices = NULL, multiple = TRUE,
                     options = list(maxItems = 5))
    )
  })

  observeEvent(input$deeppk_radar_id_col, {
    req(deeppk_resultado(), input$deeppk_radar_id_col)
    updateSelectizeInput(session, "deeppk_radar_ids",
                         choices = unique(as.character(deeppk_resultado()[[input$deeppk_radar_id_col]])),
                         server = TRUE)
  })

  output$deeppk_smiles_col_selector <- renderUI({
    req(deeppk_resultado())
    cols <- names(deeppk_resultado())
    guess <- cols[grepl("smiles", cols, ignore.case = TRUE)]
    selectInput("deeppk_smiles_col", "SMILES", choices = cols,
                selected = if (length(guess) > 0) guess[1] else cols[1])
  })

  output$deeppk_label_col_selector <- renderUI({
    req(deeppk_resultado())
    cols <- names(deeppk_resultado())
    selectInput("deeppk_label_col", "Label", choices = cols, selected = cols[1])
  })

  output$deeppk_pca_controls <- renderUI({
    req(deeppk_resultado())
    cols <- names(deeppk_resultado())
    num_cols <- cols[sapply(deeppk_resultado(), is.numeric)]
    tagList(
      selectInput("deeppk_pca_variables", "Variables", choices = num_cols,
                  selected = intersect(c("MW","TPSA","LogP","#H-bond acceptors","#H-bond donors","#Rotatable bonds"), num_cols),
                  multiple = TRUE),
      selectInput("deeppk_pca_color", "Colour", choices = c("None", cols), selected = "None"),
      selectInput("deeppk_pca_label", "Labels", choices = c("None", cols), selected = "None"),
      checkboxInput("deeppk_pca_ellipse", "Ellipses", TRUE),
      checkboxInput("deeppk_pca_scale", "Scale", TRUE)
    )
  })

  output$deeppk_tsne_controls <- renderUI({
    req(deeppk_resultado())
    cols <- names(deeppk_resultado())
    num_cols <- cols[sapply(deeppk_resultado(), is.numeric)]
    tagList(
      selectInput("deeppk_tsne_variables", "Variables", choices = num_cols,
                  selected = intersect(c("MW","TPSA","LogP","#H-bond acceptors","#H-bond donors","#Rotatable bonds"), num_cols),
                  multiple = TRUE),
      sliderInput("deeppk_tsne_perplexity", "Perplexity", min = 5, max = 50, value = 30),
      sliderInput("deeppk_tsne_iter", "Iterations", min = 500, max = 5000, value = 1000, step = 500),
      selectInput("deeppk_tsne_color", "Colour", choices = c("None", cols), selected = "None"),
      selectInput("deeppk_tsne_label", "Labels", choices = c("None", cols), selected = "None"),
      checkboxInput("deeppk_tsne_scale", "Scale", TRUE)
    )
  })

  output$deeppk_umap_controls <- renderUI({
    req(deeppk_resultado())
    cols <- names(deeppk_resultado())
    num_cols <- cols[sapply(deeppk_resultado(), is.numeric)]
    tagList(
      selectInput("deeppk_umap_variables", "Variables", choices = num_cols,
                  selected = intersect(c("MW","TPSA","LogP","#H-bond acceptors","#H-bond donors","#Rotatable bonds"), num_cols),
                  multiple = TRUE),
      sliderInput("deeppk_umap_n_neighbors", "n_neighbors", min = 2, max = 50, value = 15),
      sliderInput("deeppk_umap_min_dist", "min_dist", min = 0, max = 1, value = 0.1, step = 0.05),
      selectInput("deeppk_umap_color", "Colour", choices = c("None", cols), selected = "None"),
      selectInput("deeppk_umap_label", "Labels", choices = c("None", cols), selected = "None"),
      checkboxInput("deeppk_umap_scale", "Scale", TRUE)
    )
  })

  output$deeppk_parallel_controls <- renderUI({
    req(deeppk_resultado())
    cols <- names(deeppk_resultado())
    num_cols <- cols[sapply(deeppk_resultado(), is.numeric)]
    tagList(
      selectInput("deeppk_parallel_variables", "Variables", choices = num_cols,
                  selected = intersect(c("MW","TPSA","LogP","#H-bond acceptors","#H-bond donors","#Rotatable bonds"), num_cols),
                  multiple = TRUE),
      selectInput("deeppk_parallel_color", "Colour", choices = c("None", cols), selected = "None"),
      checkboxInput("deeppk_parallel_scale", "Scale", TRUE)
    )
  })

  output$deeppk_cluster_heatmap_controls <- renderUI({
    req(deeppk_resultado())
    cols <- names(deeppk_resultado())
    numeric_cols <- cols[sapply(deeppk_resultado(), is.numeric)]
    ## Look for a suitable label column: prefer text columns that could be IDs
    char_cols <- cols[!sapply(deeppk_resultado(), is.numeric)]
    default_id <- if (any(c("Name", "name", "ID", "Compound", "Molecule") %in% cols)) {
      intersect(c("Name", "name", "ID", "Compound", "Molecule"), cols)[1]
    } else if (length(char_cols) > 0) {
      char_cols[1]
    } else {
      "None"
    }
    tagList(
      selectInput("deeppk_cluster_heatmap_variables", "Variables", choices = numeric_cols,
                  selected = c("MW","TPSA","LogP","#H-bond acceptors","#H-bond donors","#Rotatable bonds"),
                  multiple = TRUE),
      selectInput("deeppk_cluster_heatmap_id_col", "Label column (optional)", choices = c("None", cols),
                  selected = default_id),
      selectInput("deeppk_cluster_heatmap_method", "Clustering method",
                  choices = c("Ward D2" = "ward.D2", "Ward D" = "ward.D",
                              "Complete" = "complete", "Average (UPGMA)" = "average",
                              "Single" = "single", "McQuitty" = "mcquitty",
                              "Median" = "median", "Centroid" = "centroid"),
                  selected = "ward.D2"),
      checkboxInput("deeppk_cluster_heatmap_scale", "Scale variables (z-score)", TRUE)
    )
  })

  output$deeppk_violin_controls <- renderUI({
    req(deeppk_resultado())
    cols <- names(deeppk_resultado())
    num_cols <- cols[sapply(deeppk_resultado(), is.numeric)]
    cat_cols <- cols[!sapply(deeppk_resultado(), is.numeric)]
    tagList(
      selectInput("deeppk_violin_variable", "Variable", choices = num_cols, selected = "MW"),
      selectInput("deeppk_violin_group", "Group by", choices = cat_cols,
                  selected = if ("BBB permeant" %in% cat_cols) "BBB permeant" else cat_cols[1]),
      checkboxInput("deeppk_violin_box", "Boxplot", TRUE),
      checkboxInput("deeppk_violin_points", "Points", FALSE)
    )
  })

  ## ---- Deep-PK: renderUI blocks for the new plot types ----
  output$deeppk_histogram_custom_controls <- renderUI({
    req(deeppk_resultado())
    cols <- names(deeppk_resultado())
    numeric_cols <- cols[sapply(deeppk_resultado(), is.numeric)]
    categorical_cols <- cols[!sapply(deeppk_resultado(), is.numeric)]
    tagList(
      selectInput("deeppk_hc_variable", "Variable", choices = numeric_cols, selected = "MW"),
      sliderInput("deeppk_hc_bins", "Number of bins", min = 5, max = 100, value = 30),
      selectInput("deeppk_hc_group", "Group by (overlay)", choices = c("None", categorical_cols), selected = "None"),
      checkboxInput("deeppk_hc_density", "Show density curve", TRUE),
      checkboxInput("deeppk_hc_rug", "Show rug plot", FALSE)
    )
  })

  ## ---- Deep-PK plot rendering ----
  deeppk_base_plot_types <- c("Radar plot (Chemical profile)",
                              "Tanimoto / AGNES (Structural similarity)")

  deeppk_currentPlot <- reactive({
    req(deeppk_resultado())
    validate(need(nrow(deeppk_resultado()) > 0, "No data."))
    pt <- input$deeppk_plot_type

    if (pt %in% deeppk_base_plot_types) {
      if (pt == "Radar plot (Chemical profile)") {
        req(input$deeppk_radar_id_col)
        function() plotRadar(deeppk_resultado(),
                             id_col = input$deeppk_radar_id_col,
                             ids = input$deeppk_radar_ids)
      } else {
        req(input$deeppk_smiles_col)
        function() plotTanimoto(deeppk_resultado(),
                                smiles_col = input$deeppk_smiles_col,
                                label_col = input$deeppk_label_col,
                                max_n = input$deeppk_tanimoto_max_n,
                                method = input$deeppk_agnes_method)
      }
    } else {
      switch(pt,
             "Boiled Egg" = plotBoiledEgg(deeppk_resultado()),
             "TPSA" = plotTPSA(deeppk_resultado()),
             "LogP" = plotLogP(deeppk_resultado()),
             "Correlation Heatmap" = plotCorrHeatmap(deeppk_resultado()),
             "Principal Component Analysys (PCA - Chemical space)" =
               apply_palette(plotPCA(data = deeppk_resultado(), variables = input$deeppk_pca_variables,
                                     color_by = input$deeppk_pca_color, label_by = input$deeppk_pca_label,
                                     scale_data = input$deeppk_pca_scale, ellipse = input$deeppk_pca_ellipse),
                             input$deeppk_palette, deeppk_resultado(), input$deeppk_pca_color),
             "t-SNE (Chemical space)" =
               apply_palette(plotTSNE(data = deeppk_resultado(), variables = input$deeppk_tsne_variables,
                                      color_by = input$deeppk_tsne_color, label_by = input$deeppk_tsne_label,
                                      perplexity = input$deeppk_tsne_perplexity, max_iter = input$deeppk_tsne_iter,
                                      scale_data = input$deeppk_tsne_scale),
                             input$deeppk_palette, deeppk_resultado(), input$deeppk_tsne_color),
             "UMAP (Chemical space)" =
               apply_palette(plotUMAP(data = deeppk_resultado(), variables = input$deeppk_umap_variables,
                                      color_by = input$deeppk_umap_color, label_by = input$deeppk_umap_label,
                                      n_neighbors = input$deeppk_umap_n_neighbors,
                                      min_dist = input$deeppk_umap_min_dist,
                                      scale_data = input$deeppk_umap_scale),
                             input$deeppk_palette, deeppk_resultado(), input$deeppk_umap_color),
             "Parallel Coordinates" =
               apply_palette(plotParallel(data = deeppk_resultado(), variables = input$deeppk_parallel_variables,
                                          color_by = input$deeppk_parallel_color, scale_data = input$deeppk_parallel_scale),
                             input$deeppk_palette, deeppk_resultado(), input$deeppk_parallel_color),
             "Cluster Heatmap (Dendrogram)" =
               function() plotClusterHeatmap(data = deeppk_resultado(),
                                             variables = input$deeppk_cluster_heatmap_variables,
                                             id_col = .safe_id_col(input$deeppk_cluster_heatmap_id_col),
                                             method = input$deeppk_cluster_heatmap_method,
                                             scale_data = input$deeppk_cluster_heatmap_scale,
                                             palette = input$deeppk_palette),
             "Violin Plot" =
               apply_palette(plotViolin(data = deeppk_resultado(), variable = input$deeppk_violin_variable,
                                         group_by = input$deeppk_violin_group, show_box = input$deeppk_violin_box,
                                         show_points = input$deeppk_violin_points),
                             input$deeppk_palette, deeppk_resultado(), input$deeppk_violin_group),
             "Custom Histogram" =
               plotHistogramCustom(data = deeppk_resultado(),
                                   variable = input$deeppk_hc_variable,
                                   bins = input$deeppk_hc_bins,
                                   group_by = input$deeppk_hc_group,
                                   show_density = input$deeppk_hc_density,
                                   show_rug = input$deeppk_hc_rug)
      )
    }
  })

  output$deeppk_plot <- renderPlot({
    p <- deeppk_currentPlot()
    if (is.function(p)) p() else print(p)
  })

  output$deeppk_downloadPlot <- downloadHandler(
    filename = function() {
      paste0("deeppk_", gsub("[^A-Za-z0-9]+", "_", input$deeppk_plot_type), ".", input$deeppk_format)
    },
    content = function(file) {
      p <- deeppk_currentPlot()
      if (is.function(p)) {
        open_device <- switch(input$deeppk_format,
          "png" = function() grDevices::png(file, width = input$deeppk_width, height = input$deeppk_height, units = "in", res = as.numeric(input$deeppk_dpi)),
          "pdf" = function() grDevices::pdf(file, width = input$deeppk_width, height = input$deeppk_height),
          function() grDevices::png(file, width = input$deeppk_width, height = input$deeppk_height, units = "in", res = as.numeric(input$deeppk_dpi))
        )
        open_device()
        p()
        grDevices::dev.off()
      } else {
        ggplot2::ggsave(filename = file, plot = p, device = input$deeppk_format,
                        width = input$deeppk_width, height = input$deeppk_height,
                        units = "in", dpi = as.numeric(input$deeppk_dpi))
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

    ## SwissADME
    if (!is.null(datos_rv())) {
      filt <- tryCatch(resultado(), error = function(e) NULL)
      datasets[["SwissADME"]] <- list(
        raw         = datos_rv(),
        filtered    = if (!is.null(filt)) filt else datos_rv(),
        filters     = if (!is.null(filt)) input$filters else character(0),
        source_name = "SwissADME Manager"
      )
    }

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

    ## ADMETlab 3.0
    if (!is.null(admetlab_datos_rv())) {
      filt <- tryCatch(admetlab_resultado(), error = function(e) NULL)
      datasets[["ADMETlab"]] <- list(
        raw         = admetlab_datos_rv(),
        filtered    = if (!is.null(filt)) filt else admetlab_datos_rv(),
        filters     = if (!is.null(filt)) input$admetlab_filters else character(0),
        source_name = "ADMETlab 3.0 Manager"
      )
    }

    ## Deep-PK
    if (!is.null(deeppk_datos_rv())) {
      filt <- tryCatch(deeppk_resultado(), error = function(e) NULL)
      datasets[["DeepPK"]] <- list(
        raw         = deeppk_datos_rv(),
        filtered    = if (!is.null(filt)) filt else deeppk_datos_rv(),
        filters     = if (!is.null(filt)) input$deeppk_filters else character(0),
        source_name = "Deep-PK Manager"
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
