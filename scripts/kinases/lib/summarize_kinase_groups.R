#!/usr/bin/env Rscript
# summarize_kinase_groups.R
# Author: C. Ryan Miller
# Updated: 2025-12-31 (refactored to use config_loader.R)
# Purpose: Summarize kinases by Manning group classification
# Usage: Rscript lib/summarize_kinase_groups.R

library(data.table)

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

# Load annotated kinases (prefer human, fall back to mouse)
human_file <- output_path("kinases_human_annotated.csv")
mouse_file <- output_path("kinases_mouse_orthologs.csv")

if (file.exists(human_file)) {
  kinases <- fread(human_file)
  species <- "human"
} else if (file.exists(mouse_file)) {
  kinases <- fread(mouse_file)
  species <- "mouse"
} else {
  stop("No annotated kinase file found. Run the pipeline first.")
}

cat(sprintf("=== Kinase Group Summary (%s) ===\n\n", species))

# Summarize by group (if Group column exists)
if ("Group name" %in% names(kinases)) {
  group_col <- "Group name"
} else if ("Manning_Group" %in% names(kinases)) {
  group_col <- "Manning_Group"
} else if ("Group" %in% names(kinases)) {
  group_col <- "Group"
} else {
  cat("No group column found. Available columns:\n")
  print(names(kinases))
  quit(status = 1)
}

group_summary <- kinases[, .N, by = group_col][order(-N)]
setnames(group_summary, group_col, "Group")

cat("Kinases by group:\n")
print(group_summary)

cat(sprintf("\nTotal: %d kinases in %d groups\n", nrow(kinases), nrow(group_summary)))

# Write summary
output_file <- output_path("kinases_group_summary.csv")
fwrite(group_summary, output_file)
cat(sprintf("\n✓ Group summary written to: %s\n", output_file))
