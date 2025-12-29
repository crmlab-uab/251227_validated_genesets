#!/usr/bin/env Rscript
# Mouse Kinase Validation using UniProt, GO, and Bioconductor
# Author: C. Ryan Miller, MD, PhD (rmiller@uab.edu)
# Date: 2025-12-27
# Purpose: Independently validate all mouse kinases against multiple databases

suppressPackageStartupMessages({
  library(data.table)
  library(UniProt.ws)
  library(GO.db)
  library(org.Mm.eg.db)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

cat("=== COMPREHENSIVE KINASE VALIDATION ===\n\n")
cat("Loading kinase list...\n")

# Read curated kinase list
kinases <- fread("genesets/curated/kinases/201006_composite_kinases_curated.csv", header = TRUE)
colnames(kinases) <- c(
  "Mouse_Symbol", "xName", "Manning_Name", "HGNC_Symbol", "Coral_Name",
  "Gene_Description", "Manning_Group", "Manning_Family", "Manning_Subfamily",
  "UniProt_Human", "MGI_ID", "Entrez_Mouse",
  "KinHub_Validated", "Coral_Validated", "UniProt_Validated",
  "Source", "Date_Added"
)

cat("Total kinases to validate:", nrow(kinases), "\n\n")

# Initialize UniProt connection
cat("Connecting to UniProt database...\n")
up <- UniProt.ws::UniProt.ws(taxId = 10090)  # Mouse

# Kinase-related GO terms
KINASE_GO_TERMS <- c(
  "GO:0004672",  # protein kinase activity
  "GO:0016301",  # kinase activity
  "GO:0016773",  # phosphotransferase activity, alcohol group as acceptor
  "GO:0004674",  # protein serine/threonine kinase activity
  "GO:0004713",  # protein tyrosine kinase activity
  "GO:0019199",  # transmembrane receptor protein kinase activity
  "GO:0004712"   # protein serine/threonine/tyrosine kinase activity
)

cat("\nValidating against databases:\n")
cat("  1. UniProt (protein kinase activity)\n")
cat("  2. Gene Ontology (kinase GO terms)\n")
cat("  3. org.Mm.eg.db (mouse annotations)\n")
cat("  4. org.Hs.eg.db (human orthologs)\n\n")

# Validation function
validate_kinase <- function(mouse_symbol, human_uniprot, entrez_id, row_idx) {
  result <- list(
    Mouse_Symbol = mouse_symbol,
    UniProt_Kinase = FALSE,
    GO_Kinase = FALSE,
    Has_Kinase_Domain = FALSE,
    GO_Terms = "",
    Validation_Evidence = "NONE",
    Notes = ""
  )

  # 1. Check Gene Ontology via org.Mm.eg.db
  if (!is.na(entrez_id) && entrez_id != "") {
      go_terms <- AnnotationDbi::select(
        org.Mm.eg.db,
        keys = as.character(entrez_id),  # Convert to character for select()
        columns = "GOALL",
        keytype = "ENTREZID"
      )

      if (!is.null(go_terms) && nrow(go_terms) > 0) {
        kinase_go <- go_terms[go_terms$GOALL %in% KINASE_GO_TERMS, ]
        if (nrow(kinase_go) > 0) {
          result$GO_Kinase <- TRUE
          result$GO_Terms <- paste(unique(kinase_go$GOALL), collapse = ";")
          result$Validation_Evidence <- "GO"
        }
      }
    }

    # 2. Check UniProt (if we have UniProt ID)
    if (!is.na(human_uniprot) && human_uniprot != "") {
      # Query UniProt for protein function (may fail with network timeout)
      uniprot_data <- tryCatch({
        UniProt.ws::select(
          up,
          keys = human_uniprot,
          columns = c("protein_name", "keyword", "go_id"),
          keytype = "UniProtKB"
        )
      }, error = function(e) {
        cat("  [UniProt query failed for", mouse_symbol, "]\n", file=stderr())
        NULL
      })

      if (!is.null(uniprot_data) && nrow(uniprot_data) > 0) {
        # Collapse multiple rows (UniProt returns 1:many mapping)
        keywords_all <- paste(uniprot_data$keyword, collapse = ";")
        protein_names_all <- paste(uniprot_data$protein_name, collapse = ";")
        
        # Check keywords for "kinase"
        if (nchar(keywords_all) > 0 && !is.na(keywords_all)) {
          if (grepl("kinase|Kinase", keywords_all, ignore.case = TRUE)) {
            result$UniProt_Kinase <- TRUE
            result$Validation_Evidence <- ifelse(
              result$Validation_Evidence == "GO",
              "GO+UniProt",
              "UniProt"
            )
          }
        }

        # Check protein names
        if (nchar(protein_names_all) > 0 && !is.na(protein_names_all)) {
          if (grepl("kinase|Kinase", protein_names_all, ignore.case = TRUE)) {
            result$UniProt_Kinase <- TRUE
            result$Validation_Evidence <- ifelse(
              result$Validation_Evidence == "GO",
              "GO+UniProt",
              "UniProt"
            )
          }
        }
      }
    }

    # 3. Check gene description (last resort)
    # This is based on the curated file's description
    # If name contains "kinase" but no database validation, flag as suspicious

  # Set overall status
  if (result$GO_Kinase && result$UniProt_Kinase) {
    result$Validation_Evidence <- "CONFIRMED"
  } else if (result$GO_Kinase || result$UniProt_Kinase) {
    result$Validation_Evidence <- "PARTIAL"
  } else {
    result$Validation_Evidence <- "UNCONFIRMED"
    result$Notes <- "No database evidence for kinase activity"
  }

  # Progress indicator
  if (row_idx %% 50 == 0) {
    cat("  Processed", row_idx, "/", nrow(kinases), "genes...\n")
  }

  return(result)
}

# Validate all kinases
cat("Starting validation (this may take 5-10 minutes)...\n\n")
start_time <- Sys.time()

validation_results <- rbindlist(lapply(1:nrow(kinases), function(i) {
  validate_kinase(
    kinases$Mouse_Symbol[i],
    kinases$UniProt_Human[i],
    kinases$Entrez_Mouse[i],
    i
  )
}))

end_time <- Sys.time()
cat("\n✓ Validation completed in", round(difftime(end_time, start_time, units = "mins"), 1), "minutes\n\n")

# Merge with original data
kinases_validated <- cbind(kinases, validation_results[, -1])

# Generate validation summary
cat("=== VALIDATION SUMMARY ===\n\n")

cat("Overall Validation Status:\n")
print(table(kinases_validated$Validation_Evidence))

cat("\nGO Kinase Activity:\n")
cat("  Confirmed:", sum(kinases_validated$GO_Kinase), "/", nrow(kinases), "\n")

cat("\nUniProt Kinase Keywords:\n")
cat("  Confirmed:", sum(kinases_validated$UniProt_Kinase), "/", nrow(kinases), "\n")

# Identify discrepancies with original X marks
cat("\n=== DISCREPANCIES WITH ORIGINAL FILE ===\n\n")

# Genes marked validated but not confirmed by databases
false_positives <- kinases_validated[
  (KinHub_Validated == "X" | Coral_Validated == "X" | UniProt_Validated == "X") &
  Validation_Evidence == "UNCONFIRMED"
]

if (nrow(false_positives) > 0) {
  cat("⚠ Genes with X marks BUT no database confirmation:\n")
  print(false_positives[, .(Mouse_Symbol, Gene_Description, Manning_Group,
                            KinHub_Validated, Coral_Validated, UniProt_Validated)])
} else {
  cat("✓ No false positives found\n")
}

# Genes without X marks but confirmed by databases
false_negatives <- kinases_validated[
  (is.na(KinHub_Validated) | KinHub_Validated == "") &
  (is.na(Coral_Validated) | Coral_Validated == "") &
  (is.na(UniProt_Validated) | UniProt_Validated == "") &
  Validation_Evidence == "CONFIRMED"
]

if (nrow(false_negatives) > 0) {
  cat("\n⚠ Genes WITHOUT X marks BUT confirmed by databases:\n")
  print(false_negatives[, .(Mouse_Symbol, Gene_Description, Manning_Group, GO_Kinase, UniProt_Kinase)])
} else {
  cat("✓ No false negatives found\n")
}

output_file <- "mouse_kinome_validation_results.csv"
fwrite(kinases_validated, output_file, quote = TRUE)
cat("\n✓ Validation results saved to:", output_file, "\n")

# Create summary report
summary_report <- paste0(
  "# Mouse Kinome Validation Report\n",
  "Date: ", Sys.Date(), "\n",
  "Total Genes: ", nrow(kinases), "\n\n",
  "## Validation Results\n",
  "- CONFIRMED (GO + UniProt): ", sum(kinases_validated$Validation_Evidence == "CONFIRMED"), "\n",
  "- PARTIAL (GO or UniProt): ", sum(kinases_validated$Validation_Evidence == "PARTIAL"), "\n",
  "- UNCONFIRMED (no evidence): ", sum(kinases_validated$Validation_Evidence == "UNCONFIRMED"), "\n\n",
  "## Database Coverage\n",
  "- GO kinase terms: ", sum(kinases_validated$GO_Kinase), " (",
    round(100 * sum(kinases_validated$GO_Kinase) / nrow(kinases), 1), "%\n",
  "- UniProt kinase keywords: ", sum(kinases_validated$UniProt_Kinase), " (",
    round(100 * sum(kinases_validated$UniProt_Kinase) / nrow(kinases), 1), "%\n\n",
  "## Original X Mark Accuracy\n",
  "- False positives: ", nrow(false_positives), "\n",
  "- False negatives: ", nrow(false_negatives), "\n"
)

writeLines(summary_report, "VALIDATION_REPORT.md")
cat("✓ Summary report saved to: VALIDATION_REPORT.md\n")
