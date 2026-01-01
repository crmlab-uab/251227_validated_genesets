# export_kinase_gmt.R
# Export mouse kinases as GMT for GSEA, with clarified human symbol columns
# Author: C. Ryan Miller, MD, PhD (rmiller@uab.edu)
# Updated: 2025-12-31 (refactored to use config_loader.R)
#
# Column mapping based on curated file:
# Human_Symbol: Manning (canonical)
# Human_Symbol2: KinHub
# Human_Symbol3: Coral
# Human_Symbol4: (spare/legacy, often matches canonical)
#
# This file is designed to be sourced by 05_export_gmt.R or used standalone

library(data.table)

# Load config if not already loaded
if (!exists("repo_root") || !exists("paths")) {
  config_loader_path <- file.path(dirname(sys.frame(1)$ofile %||% "."), "config_loader.R")
  if (file.exists(config_loader_path)) {
    source(config_loader_path)
  } else {
    cur <- normalizePath(".", mustWork = FALSE)
    for (i in 1:10) {
      cand <- file.path(cur, "genesets_config.yaml")
      if (file.exists(cand)) {
        repo_root <- cur
        break
      }
      parent <- dirname(cur)
      if (parent == cur) break
      cur <- parent
    }
    if (!exists("repo_root")) {
      repo_root <- normalizePath(".", mustWork = FALSE)
    }
    paths <- list(
      output_dir = file.path(repo_root, "curated/kinases/outputs")
    )
    output_path <- function(f) file.path(paths$output_dir, f)
    ensure_dir <- function(d) {
      if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
      invisible(d)
    }
  }
}

ensure_dir(paths$output_dir)

if (!exists("input_file") || !file.exists(input_file)) {
  candidates <- list.files(paths$output_dir, pattern = "mouse_kinome_definitive.*\\.csv$",
                           full.names = TRUE, ignore.case = TRUE)
  if (length(candidates) == 0 && exists("paths") && !is.null(paths$input_dir)) {
    candidates <- list.files(paths$input_dir, pattern = "mouse_kinome_definitive.*\\.csv$",
                             full.names = TRUE, ignore.case = TRUE)
  }
  if (length(candidates) > 0) {
    input_file <- sort(candidates, decreasing = TRUE)[1]
  } else {
    stop("No mouse_kinome_definitive*.csv file found. Run 04_map_human_to_mouse.R first.")
  }
}

if (!file.exists(input_file) || file.info(input_file)$size == 0) {
  stop(paste0("Required input file missing or empty: ", input_file,
              "\nCheck that previous pipeline steps completed successfully."))
}

gmt_mouse_file <- output_path("mouse_kinome_all_sources.gmt")
kinases <- fread(input_file)

symbol_col <- intersect(c("Mouse_symbol", "external_gene_name", "symbol"), names(kinases))[1]
if (is.na(symbol_col)) {
  stop("Cannot find symbol column in input file")
}

cat("Exporting mouse kinome GMT...\n")
cat("MOUSE_KINOME\tAll curated mouse kinases\t",
    paste(kinases[[symbol_col]], collapse = "\t"),
    "\n",
    file = gmt_mouse_file)

cat("\u2713 GMT export complete:", gmt_mouse_file, "\n")
