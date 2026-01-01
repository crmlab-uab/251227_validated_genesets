#!/usr/bin/env Rscript
# 02_reconcile_sources.R
# Author: C. Ryan Miller, MD, PhD (rmiller@uab.edu)
# Created: 2025-12-31
# Purpose: Reconcile Manning tableS1 and KinHub gene lists, flagging discrepancies
# Usage: Rscript 02_reconcile_sources.R

suppressPackageStartupMessages(library(data.table))

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

cat("=== Manning vs KinHub Reconciliation ===\n\n")

# Find input files
manning_file <- cfg_get("sources.manning.path", NULL)
if (is.null(manning_file)) {
  manning_candidates <- list.files(paths$input_dir, pattern = "Manning.*\\.csv$", full.names = TRUE)
  if (length(manning_candidates) > 0) {
    manning_file <- sort(manning_candidates, decreasing = TRUE)[1]
  }
} else {
  manning_file <- file.path(repo_root, manning_file)
}

kinhub_candidates <- list.files(paths$input_dir, pattern = "KinHub.*\\.csv$", full.names = TRUE)
if (length(kinhub_candidates) == 0) {
  cat("No KinHub CSV found in", paths$input_dir, "\n")
  cat("Skipping reconciliation.\n")
  quit(status = 0)
}
kinhub_file <- sort(kinhub_candidates, decreasing = TRUE)[1]

if (is.null(manning_file) || !file.exists(manning_file)) {
  cat("No Manning file found. Skipping reconciliation.\n")
  quit(status = 0)
}

cat("Manning file:", manning_file, "\n")
cat("KinHub file:", kinhub_file, "\n\n")

# Load files
manning <- fread(manning_file, encoding = "UTF-8")
kinhub <- fread(kinhub_file, encoding = "UTF-8")

cat("Manning entries:", nrow(manning), "\n")
cat("KinHub entries:", nrow(kinhub), "\n\n")

# Standardize column names
setnames(manning, 1, "Manning_Name")
setnames(kinhub, c(1, 3), c("HGNC_Name", "Manning_Name_KH"))

# Get unique gene lists
manning_names <- unique(manning$Manning_Name)
kinhub_manning <- unique(kinhub$Manning_Name_KH)

# Build reconciliation table
reconciliation <- data.table(
  Manning_Name = union(manning_names, kinhub_manning)
)

# Add source flags
reconciliation[, In_Manning := Manning_Name %in% manning_names]
reconciliation[, In_KinHub := Manning_Name %in% kinhub_manning]

# Merge Manning metadata
manning_meta <- manning[, .(
  Manning_Name,
  Manning_Group = Group,
  Manning_Family = Family,
  Manning_Pseudogene = `Pseudogene?`,
  Manning_Novelty = Novelty
)]
reconciliation <- merge(reconciliation, manning_meta, by = "Manning_Name", all.x = TRUE)

# Merge KinHub metadata
kinhub_meta <- kinhub[, .(
  Manning_Name_KH,
  KinHub_HGNC = kinhub[[1]],  # HGNC_Name
  KinHub_Group = Group,
  KinHub_Family = Family,
  KinHub_UniProt = UniprotID
)]
setnames(kinhub_meta, "Manning_Name_KH", "Manning_Name")
reconciliation <- merge(reconciliation, kinhub_meta, by = "Manning_Name", all.x = TRUE)

# Determine status and exclusion reason
reconciliation[, Status := fcase(
  In_Manning & In_KinHub, "BOTH",
  In_Manning & !In_KinHub & grepl("Y|R", Manning_Pseudogene), "MANNING_ONLY_PSEUDO",
  In_Manning & !In_KinHub & grepl("Pseudogene", Manning_Novelty, ignore.case = TRUE), "MANNING_ONLY_PSEUDO",
  In_Manning & !In_KinHub, "MANNING_ONLY",
  !In_Manning & In_KinHub, "KINHUB_ONLY",
  default = "UNKNOWN"
)]

# Add exclusion flag and reason
reconciliation[, Exclude := Status == "MANNING_ONLY_PSEUDO"]

reconciliation[, Exclusion_Reason := fcase(
  Status == "MANNING_ONLY_PSEUDO", "Pseudogene (not in KinHub)",
  Status == "KINHUB_ONLY", "KinHub addition (not in Manning 2002)",
  Status == "MANNING_ONLY", "Manning-only (review needed)",
  default = ""
)]

# Sort
setorder(reconciliation, -In_Manning, -In_KinHub, Manning_Name)

# Print summary
cat("=== Reconciliation Summary ===\n\n")
status_counts <- table(reconciliation$Status)
for (s in names(status_counts)) {
  cat(sprintf("  %-20s: %d\n", s, status_counts[s]))
}
cat(sprintf("\n  Total genes: %d\n", nrow(reconciliation)))

# Show genes needing review
manning_only_review <- reconciliation[Status == "MANNING_ONLY"]
kinhub_only <- reconciliation[Status == "KINHUB_ONLY"]

if (nrow(manning_only_review) > 0) {
  cat("\n=== Manning-only (non-pseudogene, needs review) ===\n")
  print(manning_only_review[, .(Manning_Name, Manning_Group, Manning_Family, Manning_Novelty)])
}

if (nrow(kinhub_only) > 0) {
  cat("\n=== KinHub-only (not in Manning 2002) ===\n")
  print(kinhub_only[, .(Manning_Name, KinHub_HGNC, KinHub_Group, KinHub_Family)])
}

# Save reconciliation file
output_file <- file.path(paths$input_dir, "manning_kinhub_reconciliation.csv")
fwrite(reconciliation, output_file)
cat("\n✓ Reconciliation saved to:", output_file, "\n")

# Also save a simple exclusion list for downstream use
exclusion_list <- reconciliation[Exclude == TRUE, .(
  Manning_Name,
  Exclusion_Reason,
  Manning_Group,
  Manning_Family
)]
exclusion_file <- file.path(paths$input_dir, "pseudogene_exclusions.csv")
fwrite(exclusion_list, exclusion_file)
cat("✓ Pseudogene exclusion list saved to:", exclusion_file, "\n")
cat("  (", nrow(exclusion_list), "pseudogenes flagged for exclusion)\n")

cat("\nReconciliation complete.\n")
