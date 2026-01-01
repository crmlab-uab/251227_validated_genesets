#!/usr/bin/env Rscript
# 05_export_gmt.R
# Author: C. Ryan Miller
# Updated: 2025-12-31 (refactored to use config_loader.R)
# Purpose: Export validated kinase gene sets as GMT files for GSEA
# Usage: Rscript 05_export_gmt.R

library(data.table)

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

# Ensure output directory exists
ensure_dir(paths$output_dir)

# --- Export Human Kinases GMT ---
human_file <- output_path("kinases_human_annotated.csv")
if (file.exists(human_file)) {
  human_kinases <- fread(human_file)
  human_sym_col <- intersect(c("external_gene_name", "hgnc_symbol", "symbol"), names(human_kinases))[1]
  if (!is.na(human_sym_col)) {
    human_gmt <- output_path("kinases_human_allsources.gmt")
    cat("Exporting human kinome GMT...\n", file = stderr())
    cat("HUMAN_KINOME\tAll curated human kinases\t",
        paste(human_kinases[[human_sym_col]], collapse = "\t"),
        "\n",
        file = human_gmt)
    check_file_nonempty(human_gmt, paste0("Human GMT output missing or empty: ", human_gmt))
    if (requireNamespace("tools", quietly = TRUE)) {
      md5 <- tools::md5sum(human_gmt)
      cat(sprintf("%s  %s\n", unname(md5), basename(human_gmt)), file = paste0(human_gmt, ".md5"))
    }
    cat(sprintf("\u2713 Human GMT: %s (%d genes)\n", human_gmt, nrow(human_kinases)), file = stderr())
  }
} else {
  cat("Human annotated file not found, skipping human GMT export\n", file = stderr())
}

# --- Export Mouse Kinases GMT ---
mouse_file <- output_path("kinases_mouse_orthologs.csv")
if (!file.exists(mouse_file)) {
  # Fall back to old naming pattern
  candidates <- list.files(paths$output_dir, pattern = "mouse_kinome.*\\.csv$", full.names = TRUE)
  if (length(candidates) > 0) {
    mouse_file <- sort(candidates, decreasing = TRUE)[1]
  }
}

if (file.exists(mouse_file)) {
  mouse_kinases <- fread(mouse_file)
  mouse_sym_col <- intersect(c("Mouse_symbol", "external_gene_name", "symbol"), names(mouse_kinases))[1]
  if (!is.na(mouse_sym_col)) {
    mouse_gmt <- output_path("kinases_mouse_allsources.gmt")
    cat("Exporting mouse kinome GMT...\n", file = stderr())
    cat("MOUSE_KINOME\tAll curated mouse kinases\t",
        paste(mouse_kinases[[mouse_sym_col]], collapse = "\t"),
        "\n",
        file = mouse_gmt)
    check_file_nonempty(mouse_gmt, paste0("Mouse GMT output missing or empty: ", mouse_gmt))
    if (requireNamespace("tools", quietly = TRUE)) {
      md5 <- tools::md5sum(mouse_gmt)
      cat(sprintf("%s  %s\n", unname(md5), basename(mouse_gmt)), file = paste0(mouse_gmt, ".md5"))
    }
    cat(sprintf("\u2713 Mouse GMT: %s (%d genes)\n", mouse_gmt, nrow(mouse_kinases)), file = stderr())
  }
} else {
  cat("Mouse orthologs file not found, skipping mouse GMT export\n", file = stderr())
}

cat("Done.\n", file = stderr())
