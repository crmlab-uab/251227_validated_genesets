#!/usr/bin/env Rscript
# build_tf_list.R
# Author: C. Ryan Miller, MD, PhD (rmiller@uab.edu)
# Created: 2025-12-31
# Purpose: Generate curated transcription factor list from HGNC gene groups (primary source)
#          Validates against existing TF database and BioMart GO annotations
#          All columns track provenance via suffix convention
#
# Usage: Rscript build_tf_list.R [--skip-mouse]
#
# Outputs:
#   - tf_human_curated.csv: Main TF list with all annotations
#   - tf_human_curated.provenance.yaml: Column provenance metadata

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})

# Load centralized config
.script_dir <- (function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("--file=", "", file_arg[1]))))
  }
  return(normalizePath("."))
})()
source(file.path(.script_dir, "lib", "config_loader.R"))

# Parse command-line arguments
option_list <- list(
  make_option(c("--skip-mouse"), action = "store_true", default = FALSE,
              help = "Skip mouse ortholog mapping (faster, no BioMart queries)")
)
opt <- parse_args(OptionParser(option_list = option_list))

cat("=== Transcription Factor List Generator ===\n")
cat("Primary source: HGNC Gene Groups\n")
cat("Validation sources: Lambert et al. 2018 TF Census, JASPAR 2020, BioMart GO:0003700\n\n")

# Initialize column provenance tracker
column_provenance <- list()
add_provenance <- function(col, source, description) {
  column_provenance[[col]] <<- list(
    source = source,
    description = description,
    added = format(Sys.time(), "%Y-%m-%d")
  )
}

# ============================================================================
# STEP 1: Load HGNC Gene Groups and Filter for TFs (Primary Source)
# ============================================================================
cat("Step 1: Loading HGNC gene groups file...\n")

# Check for local HGNC gene groups file
hgnc_groups_file <- file.path(repo_root, "curated/kinases/inputs/hgnc_gene_groups.csv")
if (!file.exists(hgnc_groups_file)) {
  # Try TF inputs dir
  hgnc_groups_file <- input_path("hgnc_gene_groups.csv")
}
if (!file.exists(hgnc_groups_file)) {
  cat("  Downloading HGNC gene groups from genenames.org...\n")
  download.file(
    "https://www.genenames.org/cgi-bin/genegroup/download-all/csv",
    hgnc_groups_file,
    quiet = TRUE
  )
}

check_file_nonempty(hgnc_groups_file, "HGNC gene groups file not found")

# Load HGNC gene groups
hgnc_all <- fread(hgnc_groups_file, header = TRUE, sep = "\t")
cat("  Loaded", nrow(hgnc_all), "total gene-group entries\n")

# Define TF-related group patterns
# Include: transcription factor families, DNA-binding domains, homeoboxes, zinc fingers (specific types)
tf_patterns <- c(
  "transcription factor",
  "^basic helix.loop.helix",
  "^basic leucine zipper",
  "^forkhead box",
  "^ETS ",
  "^homeobox",
  "^nuclear receptor",
  "^GATA zinc finger",
  "^SOX ",
  "^PAX ",
  "^POU class",
  "^LIM class",
  "^KLF ",
  "^SP transcription",
  "^IRF ",
  "^STAT ",
  "^SMAD ",
  "^CREB ",
  "^ATF ",
  "^NF.kappa",
  "^NFAT ",
  "^E2F ",
  "^RFX ",
  "^TEA domain",
  "^HIF ",
  "^MYC ",
  "^MAX ",
  "^MAD ",
  "^REL ",
  "^JUN ",
  "^FOS ",
  "^MAF ",
  "^AP.1 ",
  "^RUNX ",
  "^EBF ",
  "^TCF.LEF",
  "^GLI ",
  "^TFAP",
  "^TBX ",
  "^MEF2",
  "^MYOD",
  "^NEUROD",
  "^NEUROG",
  "^ASCL",
  "^HAND",
  "^TWIST",
  "^SNAI",
  "^ZEB ",
  "^SIX ",
  "^EYA ",
  "^DACH",
  "^DLX ",
  "^MSX ",
  "^OTX ",
  "^EMX ",
  "^NKX ",
  "^PITX",
  "^LHX ",
  "^ISL ",
  "^PROX",
  "^PBX ",
  "^MEIS",
  "^TEAD",
  "^YAP ",
  "^TAZ ",
  "^HIPPO",
  "^NOTCH",
  "^HES ",
  "^HEY ",
  "^ID ",
  "^SREBP",
  "^PPAR",
  "^LXR ",
  "^FXR ",
  "^RXR ",
  "^RAR ",
  "^VDR ",
  "^THR ",
  "^GR ",
  "^ER ",
  "^AR ",
  "^PR ",
  "^MR "
)

# Build regex pattern
tf_regex <- paste0("(", paste(tf_patterns, collapse = "|"), ")", collapse = "")

# Filter for TF-related groups
tf_groups <- hgnc_all[grepl(tf_regex, `Group name`, ignore.case = TRUE)]
cat("  Found", nrow(tf_groups), "entries in TF-related groups\n")

# Show group distribution
n_groups <- length(unique(tf_groups$`Group name`))
cat("  Spanning", n_groups, "HGNC TF groups\n")

# Select and rename columns
tfs <- tf_groups[, .(
  HGNC_symbol = `Approved symbol`,
  HGNC_ID = `HGNC ID`,
  Name_hgnc = `Approved name`,
  Group_hgnc = `Group name`,
  Group_ID_hgnc = `Group ID`,
  Status_hgnc = Status,
  Locus_type_hgnc = `Locus type`
)]

# Remove duplicates (genes can be in multiple groups - keep first occurrence)
tfs <- tfs[!duplicated(HGNC_symbol)]
cat("  Unique TFs from HGNC groups:", nrow(tfs), "\n")

add_provenance("HGNC_symbol", "hgnc", "Authoritative human gene symbol from HGNC gene groups")
add_provenance("HGNC_ID", "hgnc", "HGNC database identifier")
add_provenance("Name_hgnc", "hgnc", "Full gene name from HGNC")
add_provenance("Group_hgnc", "hgnc", "HGNC gene group name (TF-related)")
add_provenance("Group_ID_hgnc", "hgnc", "HGNC gene group numeric ID")
add_provenance("Status_hgnc", "hgnc", "HGNC gene approval status")
add_provenance("Locus_type_hgnc", "hgnc", "HGNC locus type classification")

# ============================================================================
# STEP 2: Validate against Lambert et al. 2018 TF Census
# ============================================================================
cat("\nStep 2: Validating against Lambert et al. 2018 TF Census (PMID:29425488)...\n")

# The 221003_TranscriptionFactors_v_1.01.csv is from humantfs.ccbr.utoronto.ca
# This is the official TF list from Lambert et al. 2018 "The Human Transcription Factors" (Cell)
existing_tf_file <- input_path("221003_TranscriptionFactors_v_1.01.csv")
if (file.exists(existing_tf_file)) {
  existing_tfs <- fread(existing_tf_file)
  cat("  Loaded Lambert TF Census with", nrow(existing_tfs), "entries\n")

  # Merge with existing TF data
  tf_cols <- existing_tfs[, .(
    HGNC_symbol,
    Ensembl_ID_lambert = Ensembl_ID,
    DBD_lambert = DBD,
    TF_assessment_lambert = `TF assessment`
  )]

  tfs <- merge(tfs, tf_cols, by = "HGNC_symbol", all.x = TRUE)

  # Add InLambert flag
  tfs[, InLambert := ifelse(!is.na(DBD_lambert), "Y", "N")]

  n_in_lambert <- sum(tfs$InLambert == "Y")
  cat("  Matched in Lambert TF Census:", n_in_lambert, "/", nrow(tfs), "\n")

  # Also find TFs in Lambert but NOT in HGNC groups
  lambert_only <- existing_tfs[!HGNC_symbol %in% tfs$HGNC_symbol]
  if (nrow(lambert_only) > 0) {
    cat("  TFs in Lambert but not in HGNC groups:", nrow(lambert_only), "\n")
    cat("    (Will add these to final list)\n")

    # Add these TFs to our list
    lambert_additions <- lambert_only[, .(
      HGNC_symbol,
      Ensembl_ID_lambert = Ensembl_ID,
      DBD_lambert = DBD,
      TF_assessment_lambert = `TF assessment`,
      InLambert = "Y"
    )]
    lambert_additions[, `:=`(
      HGNC_ID = NA_character_,
      Name_hgnc = NA_character_,
      Group_hgnc = "Lambert TF Census (not in HGNC groups)",
      Group_ID_hgnc = NA_integer_,
      Status_hgnc = NA_character_,
      Locus_type_hgnc = NA_character_
    )]

    tfs <- rbind(tfs, lambert_additions, fill = TRUE)
    cat("  Total TFs after adding Lambert entries:", nrow(tfs), "\n")
  }

  add_provenance("Ensembl_ID_lambert", "lambert", "Ensembl gene ID from Lambert et al. 2018 TF Census")
  add_provenance("DBD_lambert", "lambert", "DNA-binding domain classification from Lambert TF Census")
  add_provenance("TF_assessment_lambert", "lambert", "TF assessment/confidence from Lambert TF Census")
  add_provenance("InLambert", "validation", "Present in Lambert et al. 2018 TF Census (Y/N)")
} else {
  cat("  Warning: Lambert TF Census file not found, skipping validation\n")
  tfs[, `:=`(Ensembl_ID_lambert = NA_character_, DBD_lambert = NA_character_,
             TF_assessment_lambert = NA_character_, InLambert = "N")]
}

# ============================================================================
# STEP 3: Validate with BioMart GO:0003700
# ============================================================================
cat("\nStep 3: Validating with BioMart GO:0003700 (DNA-binding TF activity)...\n")

if (requireNamespace("biomaRt", quietly = TRUE)) {
  library(biomaRt)

  tryCatch({
    human_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

    # Query for genes with GO:0003700 (DNA-binding transcription factor activity)
    # Also include GO:0000981 (RNA pol II specific)
    go_tfs <- getBM(
      attributes = c("external_gene_name", "ensembl_gene_id", "go_id"),
      filters = "go",
      values = c("GO:0003700", "GO:0000981"),
      mart = human_mart
    )
    setDT(go_tfs)

    # Get unique genes
    go_tfs_unique <- go_tfs[, .(
      Ensembl_ID_biomart = ensembl_gene_id[1],
      GO_TF_biomart = paste(unique(go_id[go_id != ""]), collapse = ";")
    ), by = external_gene_name]
    setnames(go_tfs_unique, "external_gene_name", "HGNC_symbol")

    cat("  BioMart GO:0003700/GO:0000981 genes:", nrow(go_tfs_unique), "\n")

    # Merge with our TF list
    tfs <- merge(tfs, go_tfs_unique, by = "HGNC_symbol", all.x = TRUE)

    # Add InGO_TF flag
    tfs[, InGO_TF := ifelse(!is.na(GO_TF_biomart), "Y", "N")]

    n_in_go <- sum(tfs$InGO_TF == "Y")
    cat("  Matched with GO TF terms:", n_in_go, "/", nrow(tfs), "\n")

    # Find GO TFs not in our list
    go_only <- go_tfs_unique[!HGNC_symbol %in% tfs$HGNC_symbol]
    go_only <- go_only[!is.na(HGNC_symbol) & HGNC_symbol != ""]
    if (nrow(go_only) > 0) {
      cat("  GO TFs not in HGNC/TFDB:", nrow(go_only), "\n")
      cat("    (Will add these to final list)\n")

      # Add these to our list
      go_additions <- go_only[, .(
        HGNC_symbol,
        Ensembl_ID_biomart,
        GO_TF_biomart,
        InGO_TF = "Y"
      )]
      go_additions[, `:=`(
        HGNC_ID = NA_character_,
        Name_hgnc = NA_character_,
        Group_hgnc = "BioMart GO:0003700 (not in HGNC groups)",
        Group_ID_hgnc = NA_integer_,
        Status_hgnc = NA_character_,
        Locus_type_hgnc = NA_character_,
        Ensembl_ID_tfdb = NA_character_,
        DBD_tfdb = NA_character_,
        TF_assessment_tfdb = NA_character_,
        InTFDB = "N"
      )]

      tfs <- rbind(tfs, go_additions, fill = TRUE)
      cat("  Total TFs after adding GO entries:", nrow(tfs), "\n")
    }

    add_provenance("Ensembl_ID_biomart", "biomart", "Ensembl gene ID from BioMart")
    add_provenance("GO_TF_biomart", "biomart", "GO terms for TF activity from BioMart")
    add_provenance("InGO_TF", "validation", "Has GO:0003700 or GO:0000981 annotation (Y/N)")

  }, error = function(e) {
    cat("  Warning: BioMart query failed:", e$message, "\n")
    tfs[, `:=`(Ensembl_ID_biomart = NA_character_, GO_TF_biomart = NA_character_,
               InGO_TF = "N")]
  })
} else {
  cat("  Warning: biomaRt package not available, skipping GO validation\n")
  tfs[, `:=`(Ensembl_ID_biomart = NA_character_, GO_TF_biomart = NA_character_,
             InGO_TF = "N")]
}

# ============================================================================
# STEP 4: Validate against JASPAR 2020 (DNA-binding motif database)
# ============================================================================
cat("\nStep 4: Validating against JASPAR 2020 database...\n")

# Query JASPAR API for human/vertebrate TFs with binding motifs
jaspar_tfs <- NULL
tryCatch({
  if (!requireNamespace("httr", quietly = TRUE)) {
    stop("Package 'httr' is required for JASPAR API. Install with: install.packages('httr')")
  }
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required for JASPAR API. Install with: install.packages('jsonlite')")
  }

  # Fetch all CORE vertebrate matrices from JASPAR
  # Using pagination to get all entries
  jaspar_url <- "https://jaspar.elixir.no/api/v1/matrix/?page_size=1000&collection=CORE&tax_group=vertebrates&format=json"

  cat("  Querying JASPAR API for vertebrate TFs...\n")
  response <- httr::GET(jaspar_url, httr::timeout(60))

  if (httr::status_code(response) == 200) {
    jaspar_data <- jsonlite::fromJSON(httr::content(response, as = "text", encoding = "UTF-8"))

    # Extract TF names from matrices
    if (!is.null(jaspar_data$results) && nrow(jaspar_data$results) > 0) {
      jaspar_tfs <- data.table(
        HGNC_symbol = jaspar_data$results$name,
        JASPAR_matrix_id = jaspar_data$results$matrix_id
      )

      # Clean up - some JASPAR names have variants like "::alt" or special chars
      jaspar_tfs[, HGNC_symbol := gsub("::.*$", "", HGNC_symbol)]
      jaspar_tfs[, HGNC_symbol := gsub("\\(var\\..*\\)$", "", HGNC_symbol)]
      jaspar_tfs[, HGNC_symbol := trimws(HGNC_symbol)]

      # Aggregate multiple matrices per gene
      jaspar_tfs <- jaspar_tfs[, .(
        JASPAR_matrix_ids = paste(unique(JASPAR_matrix_id), collapse = ";")
      ), by = HGNC_symbol]

      cat("  JASPAR TFs with binding motifs:", nrow(jaspar_tfs), "\n")

      # Merge with our TF list
      tfs <- merge(tfs, jaspar_tfs, by = "HGNC_symbol", all.x = TRUE)

      # Add InJASPAR flag
      tfs[, InJASPAR := ifelse(!is.na(JASPAR_matrix_ids), "Y", "N")]

      n_in_jaspar <- sum(tfs$InJASPAR == "Y")
      cat("  Matched in JASPAR:", n_in_jaspar, "/", nrow(tfs), "\n")

      # Find JASPAR TFs not in our list
      jaspar_only <- jaspar_tfs[!HGNC_symbol %in% tfs$HGNC_symbol]
      jaspar_only <- jaspar_only[!is.na(HGNC_symbol) & HGNC_symbol != ""]
      if (nrow(jaspar_only) > 0) {
        cat("  JASPAR TFs not in HGNC/Lambert/GO:", nrow(jaspar_only), "\n")
        cat("    (Will add these to final list)\n")

        # Add these to our list
        jaspar_additions <- jaspar_only[, .(
          HGNC_symbol,
          JASPAR_matrix_ids,
          InJASPAR = "Y"
        )]
        jaspar_additions[, `:=`(
          HGNC_ID = NA_character_,
          Name_hgnc = NA_character_,
          Group_hgnc = "JASPAR 2020 (not in other sources)",
          Group_ID_hgnc = NA_integer_,
          Status_hgnc = NA_character_,
          Locus_type_hgnc = NA_character_,
          Ensembl_ID_lambert = NA_character_,
          DBD_lambert = NA_character_,
          TF_assessment_lambert = NA_character_,
          InLambert = "N",
          Ensembl_ID_biomart = NA_character_,
          GO_TF_biomart = NA_character_,
          InGO_TF = "N"
        )]

        tfs <- rbind(tfs, jaspar_additions, fill = TRUE)
        cat("  Total TFs after adding JASPAR entries:", nrow(tfs), "\n")
      }

      add_provenance("JASPAR_matrix_ids", "jaspar", "JASPAR matrix ID(s) for known binding motifs")
      add_provenance("InJASPAR", "validation", "Has binding motif in JASPAR 2020 database (Y/N)")
    }
  } else {
    cat("  Warning: JASPAR API returned status", httr::status_code(response), "\n")
  }
}, error = function(e) {
  cat("  Warning: JASPAR API query failed:", e$message, "\n")
})

if (is.null(jaspar_tfs)) {
  cat("  JASPAR validation skipped\n")
  tfs[, `:=`(JASPAR_matrix_ids = NA_character_, InJASPAR = "N")]
}

# ============================================================================
# STEP 5: Create Consensus Ensembl ID with Status
# ============================================================================
cat("\nStep 5: Creating consensus Ensembl ID with source status...\n")

# Determine Ensembl ID and its source status
# Priority: Lambert (curated) > BioMart (automated)
# Status values:
#   - "confirmed": Both sources agree
#   - "lambert": Only in Lambert TF Census (curated, high confidence)
#   - "biomart": Only in BioMart GO annotation
#   - NA: No Ensembl ID available

tfs[, `:=`(
  # Check which sources have Ensembl IDs
  has_lambert = !is.na(Ensembl_ID_lambert) & Ensembl_ID_lambert != "",
  has_biomart = !is.na(Ensembl_ID_biomart) & Ensembl_ID_biomart != ""
)]

# Create consensus Ensembl ID (prefer Lambert)
tfs[, Ensembl_ID := ifelse(has_lambert, Ensembl_ID_lambert, Ensembl_ID_biomart)]

# Create status column
tfs[, Ensembl_ID_status := ifelse(
  has_lambert & has_biomart & Ensembl_ID_lambert == Ensembl_ID_biomart,
  "confirmed",
  ifelse(has_lambert & has_biomart & Ensembl_ID_lambert != Ensembl_ID_biomart,
         "lambert",  # Prefer curated source when they disagree
         ifelse(has_lambert, "lambert",
                ifelse(has_biomart, "biomart", NA_character_)))
)]

# Clean up temporary columns
tfs[, `:=`(has_lambert = NULL, has_biomart = NULL)]

# Report stats
n_confirmed <- sum(tfs$Ensembl_ID_status == "confirmed", na.rm = TRUE)
n_lambert_only <- sum(tfs$Ensembl_ID_status == "lambert", na.rm = TRUE)
n_biomart_only <- sum(tfs$Ensembl_ID_status == "biomart", na.rm = TRUE)
n_none <- sum(is.na(tfs$Ensembl_ID_status))

cat("  Ensembl ID status:\n")
cat("    confirmed (both sources agree):", n_confirmed, "\n")
cat("    lambert (curated source only):", n_lambert_only, "\n")
cat("    biomart (GO annotation only):", n_biomart_only, "\n")
cat("    none:", n_none, "\n")

add_provenance("Ensembl_ID", "derived", "Consensus Ensembl ID (Lambert > BioMart > HGNC priority)")
add_provenance("Ensembl_ID_status", "derived", "Source status: confirmed, lambert, biomart, hgnc")

# ============================================================================
# STEP 5b: Query HGNC for TFs still missing Ensembl IDs
# ============================================================================
missing_ensembl <- tfs[is.na(Ensembl_ID) | Ensembl_ID == ""]
if (nrow(missing_ensembl) > 0) {
  cat("\nStep 5b: Querying HGNC for", nrow(missing_ensembl), "TFs missing Ensembl IDs...\n")

  tryCatch({
    if (!requireNamespace("httr", quietly = TRUE)) {
      stop("Package 'httr' is required for HGNC API")
    }
    if (!requireNamespace("jsonlite", quietly = TRUE)) {
      stop("Package 'jsonlite' is required for HGNC API")
    }

    # Query HGNC for each missing symbol (case-insensitive search)
    hgnc_results <- list()
    symbols_to_query <- unique(missing_ensembl$HGNC_symbol)

    # Batch symbols into chunks to avoid too many API calls
    cat("  Querying HGNC REST API...\n")
    for (sym in symbols_to_query) {
      # Try exact match first, then case-insensitive
      url <- paste0("https://rest.genenames.org/fetch/symbol/", toupper(sym))

      tryCatch({
        response <- httr::GET(url,
                              httr::add_headers(Accept = "application/json"),
                              httr::timeout(10))

        if (httr::status_code(response) == 200) {
          data <- jsonlite::fromJSON(httr::content(response, as = "text", encoding = "UTF-8"))

          if (!is.null(data$response$docs) && length(data$response$docs) > 0) {
            doc <- data$response$docs[1, ]
            if (!is.null(doc$ensembl_gene_id) && !is.na(doc$ensembl_gene_id) && doc$ensembl_gene_id != "") {
              hgnc_results[[sym]] <- list(
                HGNC_symbol = sym,
                Ensembl_ID_hgnc = doc$ensembl_gene_id,
                HGNC_symbol_official = doc$symbol
              )
            }
          }
        }
        Sys.sleep(0.1)  # Rate limiting
      }, error = function(e) {
        # Silently skip failures
      })
    }

    if (length(hgnc_results) > 0) {
      hgnc_dt <- rbindlist(hgnc_results)
      cat("  Found Ensembl IDs for", nrow(hgnc_dt), "TFs via HGNC\n")

      # Update TFs with HGNC results
      for (i in seq_len(nrow(hgnc_dt))) {
        sym <- hgnc_dt$HGNC_symbol[i]
        ensembl <- hgnc_dt$Ensembl_ID_hgnc[i]
        tfs[HGNC_symbol == sym & (is.na(Ensembl_ID) | Ensembl_ID == ""),
            `:=`(Ensembl_ID = ensembl, Ensembl_ID_status = "hgnc")]
      }

      # Update stats
      n_hgnc <- sum(tfs$Ensembl_ID_status == "hgnc", na.rm = TRUE)
      n_still_missing <- sum(is.na(tfs$Ensembl_ID) | tfs$Ensembl_ID == "")
      cat("  Ensembl IDs from HGNC:", n_hgnc, "\n")
      cat("  Still missing:", n_still_missing, "\n")
    } else {
      cat("  No additional Ensembl IDs found via HGNC\n")
    }
  }, error = function(e) {
    cat("  Warning: HGNC query failed:", e$message, "\n")
  })
}

# ============================================================================
# STEP 6: Add Mouse Ortholog Mapping
# ============================================================================
if (!opt$`skip-mouse`) {
  cat("\nStep 6: Mapping mouse orthologs via BioMart...\n")

  if (requireNamespace("biomaRt", quietly = TRUE)) {
    tryCatch({
      human_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
      hgnc_symbols <- unique(tfs$HGNC_symbol[!is.na(tfs$HGNC_symbol) & tfs$HGNC_symbol != ""])

      orthologs <- getBM(
        attributes = c("external_gene_name", "mmusculus_homolog_associated_gene_name",
                       "mmusculus_homolog_ensembl_gene"),
        filters = "external_gene_name",
        values = hgnc_symbols,
        mart = human_mart
      )
      setDT(orthologs)

      orthologs <- orthologs[mmusculus_homolog_associated_gene_name != ""]
      orthologs <- orthologs[!duplicated(external_gene_name)]

      setnames(orthologs, c("external_gene_name", "mmusculus_homolog_associated_gene_name",
                            "mmusculus_homolog_ensembl_gene"),
               c("HGNC_symbol", "Mouse_Symbol", "Ensembl_Mouse"))

      tfs <- merge(tfs, orthologs, by = "HGNC_symbol", all.x = TRUE)

      n_mouse <- sum(!is.na(tfs$Mouse_Symbol))
      cat("  Mouse orthologs found:", n_mouse, "/", nrow(tfs), "\n")

      add_provenance("Mouse_Symbol", "biomart", "Mouse ortholog gene symbol from Ensembl BioMart")
      add_provenance("Ensembl_Mouse", "biomart", "Mouse Ensembl gene ID from BioMart")

    }, error = function(e) {
      cat("  Warning: BioMart query failed:", e$message, "\n")
      tfs[, `:=`(Mouse_Symbol = NA_character_, Ensembl_Mouse = NA_character_)]
    })
  } else {
    cat("  Warning: biomaRt package not available, skipping mouse mapping\n")
    tfs[, `:=`(Mouse_Symbol = NA_character_, Ensembl_Mouse = NA_character_)]
  }
} else {
  cat("\nStep 6: Skipped mouse ortholog mapping (--skip-mouse)\n")
  tfs[, `:=`(Mouse_Symbol = NA_character_, Ensembl_Mouse = NA_character_)]
}

# ============================================================================
# STEP 7: Reorder columns and write output
# ============================================================================
cat("\nStep 7: Writing output...\n")

# Define column order
col_order <- c(
  # Core identifiers
  "HGNC_symbol", "HGNC_ID", "Ensembl_ID", "Ensembl_ID_status", "Name_hgnc",
  # HGNC group info
  "Group_hgnc", "Group_ID_hgnc", "Status_hgnc", "Locus_type_hgnc",
  # Lambert TF Census validation (DBD and assessment, not Ensembl_ID)
  "DBD_lambert", "TF_assessment_lambert", "InLambert",
  # BioMart GO validation (GO terms, not Ensembl_ID)
  "GO_TF_biomart", "InGO_TF",
  # JASPAR validation
  "JASPAR_matrix_ids", "InJASPAR",
  # Mouse orthologs
  "Mouse_Symbol", "Ensembl_Mouse"
)

# Select columns that exist
existing_cols <- intersect(col_order, names(tfs))
tfs <- tfs[, ..existing_cols]

# Create output directory
ensure_dir(paths$output_dir)

# Write main output
output_file <- output_path("tf_human_curated.csv")
fwrite(tfs, output_file)
cat("  ✓ TF list:", output_file, "\n")
cat("    (", nrow(tfs), "transcription factors,", ncol(tfs), "columns)\n")

# Write provenance YAML
if (requireNamespace("yaml", quietly = TRUE)) {
  provenance_file <- output_path("tf_human_curated.provenance.yaml")

  provenance_yaml <- list(
    file = basename(output_file),
    generated = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    primary_source = "HGNC Gene Groups (TF-related)",
    validation_sources = c("Lambert et al. 2018 TF Census", "JASPAR 2020", "BioMart GO:0003700"),
    total_tfs = nrow(tfs),
    sources = list(
      hgnc = list(
        name = "HGNC Gene Names",
        url = "https://www.genenames.org/",
        role = "Primary source"
      ),
      lambert = list(
        name = "Lambert et al. 2018 - The Human Transcription Factors",
        file = "221003_TranscriptionFactors_v_1.01.csv",
        url = "https://humantfs.ccbr.utoronto.ca/",
        pubmed = "PMID:29425488",
        doi = "10.1016/j.cell.2018.01.029",
        journal = "Cell 172(4):650-665 (2018)",
        role = "Validation source - curated TF census with DBD classification"
      ),
      jaspar = list(
        name = "JASPAR 2020",
        url = "https://jaspar.elixir.no/",
        api = "https://jaspar.elixir.no/api/v1/",
        pubmed = "PMID:31701148",
        role = "Validation source - TFs with experimentally-determined binding motifs"
      ),
      biomart = list(
        name = "Ensembl BioMart",
        url = "https://www.ensembl.org/biomart",
        go_terms = c("GO:0003700", "GO:0000981"),
        role = "Validation source and mouse orthologs"
      )
    ),
    columns = column_provenance
  )

  yaml::write_yaml(provenance_yaml, provenance_file)
  cat("  ✓ Provenance:", provenance_file, "\n")
}

# Summary statistics
cat("\n=== Summary ===\n")
cat("Total transcription factors:", nrow(tfs), "\n")
cat("In Lambert TF Census:", sum(tfs$InLambert == "Y", na.rm = TRUE), "\n")
cat("In JASPAR 2020:", sum(tfs$InJASPAR == "Y", na.rm = TRUE), "\n")
cat("With GO:0003700/GO:0000981:", sum(tfs$InGO_TF == "Y", na.rm = TRUE), "\n")
cat("With Ensembl ID:", sum(!is.na(tfs$Ensembl_ID) & tfs$Ensembl_ID != ""), "\n")
cat("With mouse ortholog:", sum(!is.na(tfs$Mouse_Symbol)), "\n")

# DBD distribution (top 10)
if ("DBD_lambert" %in% names(tfs)) {
  cat("\nDNA-binding domain distribution (top 10):\n")
  dbd_dist <- tfs[!is.na(DBD_lambert) & DBD_lambert != "", .N, by = DBD_lambert][order(-N)][1:10]
  print(dbd_dist)
}

cat("\n✓ Transcription factor list generation complete.\n")
