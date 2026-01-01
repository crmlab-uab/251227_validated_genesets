#!/usr/bin/env Rscript
# 06_validate.R
# Author: C. Ryan Miller, MD, PhD (rmiller@uab.edu)
# Created: 2025-12-31
# Purpose: Run post-generation validation checks on kinase gene sets
# Usage: Rscript 06_validate.R [--skip-genome] [--skip-comprehensive]

suppressPackageStartupMessages(library(data.table))

# Load centralized config (provides repo_root, paths, cfg, helper functions)
.script_dir <- (function() {
 cmd_args <- commandArgs(trailingOnly = FALSE)
 file_arg <- grep("--file=", cmd_args, value = TRUE)
 if (length(file_arg) > 0) {
   return(dirname(normalizePath(sub("--file=", "", file_arg[1]))))
 }
 return(normalizePath("."))
})()
source(file.path(.script_dir, "lib", "config_loader.R"))

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
skip_genome <- "--skip-genome" %in% args
skip_comprehensive <- "--skip-comprehensive" %in% args

cat("=== Kinase Validation Pipeline ===\n\n")
cat("Script directory:", .script_dir, "\n")
cat("Output directory:", paths$output_dir, "\n\n")

validation_passed <- TRUE
validation_results <- list()

# --- Step 1: Summarize kinase groups ---
cat("[1/4] Summarizing kinase groups...\n")
tryCatch({
  source(file.path(.script_dir, "lib", "summarize_kinase_groups.R"))
  validation_results$summarize <- "PASSED"
  cat("     Done.\n\n")
}, error = function(e) {
  cat("     ERROR:", conditionMessage(e), "\n\n")
  validation_results$summarize <<- paste("FAILED:", conditionMessage(e))
  validation_passed <<- FALSE
})

# --- Step 2: Scan for invalid Entrez IDs ---
cat("[2/4] Scanning for invalid Entrez IDs...\n")
tryCatch({
  source(file.path(.script_dir, "lib", "scan_invalid_entrez.R"))
  validation_results$entrez <- "PASSED"
  cat("     Done.\n\n")
}, error = function(e) {
  cat("     ERROR:", conditionMessage(e), "\n\n")
  validation_results$entrez <<- paste("FAILED:", conditionMessage(e))
  validation_passed <<- FALSE
})

# --- Step 3: Validate against genome (optional, requires GTF) ---
if (!skip_genome) {
  cat("[3/4] Validating against genome annotations...\n")
  gtf_file <- "/data/bRNA3F/data/genome/mouse/gencode.vM37.primary_assembly.annotation.gtf.gz"
  if (file.exists(gtf_file)) {
    tryCatch({
      source(file.path(.script_dir, "lib", "validate_against_genome.R"))
      validation_results$genome <- "PASSED"
      cat("     Done.\n\n")
    }, error = function(e) {
      cat("     ERROR:", conditionMessage(e), "\n\n")
      validation_results$genome <<- paste("FAILED:", conditionMessage(e))
      validation_passed <<- FALSE
    })
  } else {
    cat("     SKIPPED: GTF file not found at", gtf_file, "\n\n")
    validation_results$genome <- "SKIPPED (GTF not found)"
  }
} else {
  cat("[3/4] Genome validation: SKIPPED (--skip-genome)\n\n")
  validation_results$genome <- "SKIPPED (user request)"
}

# --- Step 4: Comprehensive database validation (optional, slow) ---
if (!skip_comprehensive) {
  cat("[4/4] Running comprehensive database validation...\n")
  cat("     (This may take 5-10 minutes due to UniProt queries)\n")
  tryCatch({
    source(file.path(.script_dir, "lib", "comprehensive_kinase_validation.R"))
    validation_results$comprehensive <- "PASSED"
    cat("     Done.\n\n")
  }, error = function(e) {
    cat("     ERROR:", conditionMessage(e), "\n\n")
    validation_results$comprehensive <<- paste("FAILED:", conditionMessage(e))
    validation_passed <<- FALSE
  })
} else {
  cat("[4/4] Comprehensive validation: SKIPPED (--skip-comprehensive)\n\n")
  validation_results$comprehensive <- "SKIPPED (user request)"
}

# --- Summary ---
cat("=== Validation Summary ===\n\n")
for (name in names(validation_results)) {
  status <- validation_results[[name]][1]  # Ensure scalar
  icon <- if (grepl("^PASSED", status)) "OK" else if (grepl("^SKIPPED", status)) "--" else "!!"
  cat(sprintf("  [%s] %s: %s\n", icon, name, status))
}

cat("\n")
if (validation_passed) {
  cat("All validations passed.\n")
} else {
  cat("Some validations failed. Check output above for details.\n")
}

# List generated validation outputs
cat("\nValidation outputs:\n")
val_files <- list.files(paths$output_dir, pattern = "(validation|summary)", full.names = TRUE)
if (length(val_files) > 0) {
  for (f in val_files) {
    cat("  ", f, "\n")
  }
} else {
  cat("  (none)\n")
}
