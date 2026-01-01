#!/usr/bin/env Rscript
# export_curated_gmt.R
# Export curated gene sets (kinases, phosphatases, TFs) as GMT files
# Author: Claude Code
# Created: 2026-01-01
#
# Usage: Rscript export_curated_gmt.R
#        Run from the repository root directory

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

cat("=== Exporting Curated Gene Sets to GMT ===\n\n")

# Output directory for GMT files
gmt_dir <- "curated/gmt"
if (!dir.exists(gmt_dir)) {
  dir.create(gmt_dir, recursive = TRUE)
}

# Helper function to write GMT line
write_gmt_line <- function(name, description, genes, file, append = FALSE) {
  genes <- genes[!is.na(genes) & genes != ""]
  line <- paste(c(name, description, genes), collapse = "\t")
  cat(line, "\n", file = file, append = append, sep = "")
}

# ============================================================================
# KINASES
# ============================================================================
kinase_file <- "curated/kinases/outputs/kinases_human_curated.csv"
if (file.exists(kinase_file)) {
  cat("Processing kinases...\n")
  kinases <- read_csv(kinase_file, show_col_types = FALSE)

  # Human kinases - all
  human_gmt <- file.path(gmt_dir, "kinases_human.gmt")
  write_gmt_line(
    "KINASES_HUMAN_ALL",
    "All human kinases from HGNC (n=546)",
    kinases$HGNC_symbol,
    human_gmt
  )

  # Human kinases - protein substrate only
  protein_kinases <- kinases %>% filter(Substrate_protein == "Y")
  write_gmt_line(
    "KINASES_HUMAN_PROTEIN",
    "Human protein kinases (n=388)",
    protein_kinases$HGNC_symbol,
    human_gmt, append = TRUE
  )

  # Mouse kinases
  mouse_gmt <- file.path(gmt_dir, "kinases_mouse.gmt")
  mouse_kinases <- kinases %>% filter(!is.na(Mouse_Symbol))
  write_gmt_line(
    "KINASES_MOUSE_ALL",
    paste0("All mouse kinases (orthologs, n=", nrow(mouse_kinases), ")"),
    mouse_kinases$Mouse_Symbol,
    mouse_gmt
  )

  cat("  Human:", nrow(kinases), "genes ->", human_gmt, "\n")
  cat("  Mouse:", nrow(mouse_kinases), "genes ->", mouse_gmt, "\n")
} else {
  cat("  Kinase file not found, skipping\n")
}

# ============================================================================
# PHOSPHATASES
# ============================================================================
phos_file <- "curated/phosphatases/outputs/phosphatases_human_curated.csv"
if (file.exists(phos_file)) {
  cat("Processing phosphatases...\n")
  phosphatases <- read_csv(phos_file, show_col_types = FALSE)

  # Human phosphatases - all
  human_gmt <- file.path(gmt_dir, "phosphatases_human.gmt")
  write_gmt_line(
    "PHOSPHATASES_HUMAN_ALL",
    paste0("All human phosphatases from HGNC (n=", nrow(phosphatases), ")"),
    phosphatases$HGNC_symbol,
    human_gmt
  )

  # Human phosphatases - protein substrate only
  protein_phos <- phosphatases %>% filter(Substrate_protein == "Y")
  write_gmt_line(
    "PHOSPHATASES_HUMAN_PROTEIN",
    paste0("Human protein phosphatases (n=", nrow(protein_phos), ")"),
    protein_phos$HGNC_symbol,
    human_gmt, append = TRUE
  )

  # Human phosphatases - catalytic only
  catalytic_phos <- phosphatases %>% filter(Is_catalytic == "Y")
  write_gmt_line(
    "PHOSPHATASES_HUMAN_CATALYTIC",
    paste0("Human catalytic phosphatases (n=", nrow(catalytic_phos), ")"),
    catalytic_phos$HGNC_symbol,
    human_gmt, append = TRUE
  )

  # Mouse phosphatases
  mouse_gmt <- file.path(gmt_dir, "phosphatases_mouse.gmt")
  mouse_phos <- phosphatases %>% filter(!is.na(Mouse_Symbol))
  write_gmt_line(
    "PHOSPHATASES_MOUSE_ALL",
    paste0("All mouse phosphatases (orthologs, n=", nrow(mouse_phos), ")"),
    mouse_phos$Mouse_Symbol,
    mouse_gmt
  )

  cat("  Human:", nrow(phosphatases), "genes ->", human_gmt, "\n")
  cat("  Mouse:", nrow(mouse_phos), "genes ->", mouse_gmt, "\n")
} else {
  cat("  Phosphatase file not found, skipping\n")
}

# ============================================================================
# TRANSCRIPTION FACTORS
# ============================================================================
tf_file <- "curated/tf/outputs/tf_human_curated.csv"
if (file.exists(tf_file)) {
  cat("Processing transcription factors...\n")
  tfs <- read_csv(tf_file, show_col_types = FALSE)

  # Human TFs - all
  human_gmt <- file.path(gmt_dir, "tf_human.gmt")
  write_gmt_line(
    "TF_HUMAN_ALL",
    paste0("All human transcription factors (n=", nrow(tfs), ")"),
    tfs$HGNC_symbol,
    human_gmt
  )

  # Mouse TFs
  mouse_gmt <- file.path(gmt_dir, "tf_mouse.gmt")
  mouse_tfs <- tfs %>% filter(!is.na(Mouse_Symbol))
  write_gmt_line(
    "TF_MOUSE_ALL",
    paste0("All mouse transcription factors (orthologs, n=", nrow(mouse_tfs), ")"),
    mouse_tfs$Mouse_Symbol,
    mouse_gmt
  )

  cat("  Human:", nrow(tfs), "genes ->", human_gmt, "\n")
  cat("  Mouse:", nrow(mouse_tfs), "genes ->", mouse_gmt, "\n")
} else {
  cat("  TF file not found, skipping\n")
}

# ============================================================================
# COMBINED GMT (all sets in one file)
# ============================================================================
cat("\nCreating combined GMT files...\n")

# Human combined
human_combined <- file.path(gmt_dir, "curated_genesets_human.gmt")
system(paste("cat", file.path(gmt_dir, "*_human.gmt"), ">", human_combined))

# Mouse combined
mouse_combined <- file.path(gmt_dir, "curated_genesets_mouse.gmt")
system(paste("cat", file.path(gmt_dir, "*_mouse.gmt"), ">", mouse_combined))

cat("  Combined human:", human_combined, "\n")
cat("  Combined mouse:", mouse_combined, "\n")

# Summary
cat("\n=== GMT Export Complete ===\n")
cat("Output directory:", gmt_dir, "\n")
system(paste("ls -la", gmt_dir))
