#!/usr/bin/env Rscript
# scan_invalid_entrez.R
# Author: C. Ryan Miller, MD, PhD (rmiller@uab.edu)
# Updated: 2025-12-31 (refactored to use config_loader.R)
# Purpose: Scan for invalid Entrez IDs in kinase list
# Usage: Rscript lib/scan_invalid_entrez.R

library(data.table)
library(org.Mm.eg.db)

# Load centralized config
.script_dir <- (function() {
 cmd_args <- commandArgs(trailingOnly = FALSE)
 file_arg <- grep("--file=", cmd_args, value = TRUE)
 if (length(file_arg) > 0) {
   return(dirname(normalizePath(sub("--file=", "", file_arg[1]))))
 }
 return(normalizePath("."))
})()
if (file.exists(file.path(.script_dir, "config_loader.R"))) {
  source(file.path(.script_dir, "config_loader.R"))
} else if (file.exists(file.path(.script_dir, "lib", "config_loader.R"))) {
  source(file.path(.script_dir, "lib", "config_loader.R"))
}

cat("=== Scanning for Invalid Entrez IDs ===\n\n")

# Find kinase input file
candidates <- list.files(paths$input_dir, pattern='201006_composite_kinases_curated.*\\.csv$', full.names=TRUE, ignore.case=TRUE)
if (length(candidates) == 0) stop('Missing input snapshot: please place 201006_composite_kinases_curated__YYMMDD.csv in ', paths$input_dir)
kinases_file <- sort(candidates, decreasing=TRUE)[1]

cat("Loading kinases from:", kinases_file, "\n")
kinases <- fread(kinases_file, header=TRUE)
cat("  Found", nrow(kinases), "kinases\n")

# Filter out excluded genes from config
exclude_genes <- get_exclude_genes()
symbol_col_temp <- intersect(c("Mouse_symbol", "Mouse_Symbol", "external_gene_name", "symbol"), names(kinases))[1]
if (length(exclude_genes) > 0 && !is.na(symbol_col_temp)) {
  n_before <- nrow(kinases)
  excluded_found <- kinases[[symbol_col_temp]][kinases[[symbol_col_temp]] %in% exclude_genes]
  kinases <- kinases[!kinases[[symbol_col_temp]] %in% exclude_genes]
  n_excluded <- n_before - nrow(kinases)
  if (n_excluded > 0) {
    cat("  Excluded", n_excluded, "genes per config:", paste(excluded_found, collapse=", "), "\n")
  }
}
cat("  Processing", nrow(kinases), "kinases\n")

cat("Checking against org.Mm.eg.db...\n\n")
all_entrez <- keys(org.Mm.eg.db, keytype='ENTREZID')

# Find the Entrez ID column (handles different naming conventions)
entrez_col <- intersect(c("Entrez_ID_mouse", "Entrez_Mouse", "entrez_id", "entrezgene_id"), names(kinases))[1]
symbol_col <- intersect(c("Mouse_symbol", "Mouse_Symbol", "external_gene_name", "symbol"), names(kinases))[1]

if (is.na(entrez_col)) {
  cat("No Entrez ID column found in input file. Skipping validation.\n")
  cat("Available columns:", paste(names(kinases), collapse=", "), "\n")
} else {
  cat("Using column:", entrez_col, "for Entrez IDs\n")
  bad <- kinases[!is.na(get(entrez_col)) & get(entrez_col) != '' & !(as.character(get(entrez_col)) %in% all_entrez)]

  if (nrow(bad) > 0) {
    cat("Found", nrow(bad), "invalid Entrez IDs:\n")
    if (!is.na(symbol_col)) {
      print(bad[, .SD, .SDcols = c(symbol_col, entrez_col)])
    } else {
      print(bad[[entrez_col]])
    }
  } else {
    cat("All Entrez IDs are valid.\n")
  }
}
