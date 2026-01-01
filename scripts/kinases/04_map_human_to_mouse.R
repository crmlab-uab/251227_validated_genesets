#!/usr/bin/env Rscript
# 04_map_human_to_mouse.R
# Author: C. Ryan Miller
# Created: 2025-12-28 02:17 CST
# Updated: 2025-12-31 (refactored to use config_loader.R)
# Purpose: Map human kinases to mouse orthologs via BioMart
# Usage: Rscript 04_map_human_to_mouse.R
#
# NOTE: This is a thin runner. All mapping logic is consolidated in scripts/kinases/lib/
# Do not duplicate mapping logic here; update lib/ scripts for mapping changes.
# Output: mouse_kinome_definitive.csv with Mouse_Symbol, Ensembl_Gene_ID, Entrez_ID, UniProt_ID, MGI_ID, Description

library(data.table)
library(biomaRt)
library(httr)

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

# Get curated kinase list from inputs
candidates <- list.files(paths$input_dir, pattern = "201006_composite_kinases_curated.*\\.csv$",
                         full.names = TRUE, ignore.case = TRUE)
if (length(candidates) == 0) {
  stop("Missing input snapshot: place 201006_composite_kinases_curated__YYMMDD.csv in ", paths$input_dir)
}
curated_file <- sort(candidates, decreasing = TRUE)[1]
check_file_nonempty(curated_file, paste0("Curated kinases input missing or empty: ", curated_file))
cat("Using curated kinase input:", curated_file, "\n", file = stderr())

curated <- fread(curated_file)
mouse_symbols <- unique(curated$Mouse_symbol)

# Query Ensembl BioMart for annotation IDs
cat("Querying BioMart for mouse kinase annotations...\n", file = stderr())
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
bm <- getBM(
  attributes = c(
    "external_gene_name", "ensembl_gene_id", "entrezgene_id",
    "uniprotswissprot", "mgi_id", "description"
  ),
  filters = "external_gene_name",
  values = mouse_symbols,
  mart = ensembl
)
bm <- as.data.table(bm)
setnames(bm, c("external_gene_name", "ensembl_gene_id", "entrezgene_id", "uniprotswissprot", "mgi_id", "description"),
         c("Mouse_symbol", "Ensembl_Gene_ID", "Entrez_ID", "UniProt_ID", "MGI_ID", "Description"))

# Merge with curated to ensure only true kinases
bm <- bm[Mouse_symbol %in% mouse_symbols]

# Remove duplicates, keep one row per Mouse_symbol (prefer with all IDs present)
setorder(bm, Mouse_symbol, -Ensembl_Gene_ID, -UniProt_ID, -MGI_ID, -Entrez_ID)
bm <- bm[!duplicated(Mouse_symbol)]

# Write output to canonical outputs directory
ensure_dir(paths$output_dir)
outfile <- output_path("kinases_mouse_orthologs.csv")

fwrite(bm, outfile)
check_file_nonempty(outfile, paste0("Output file missing or empty: ", outfile))

# Write MD5 checksum
if (requireNamespace("tools", quietly = TRUE)) {
  md5 <- tools::md5sum(outfile)
  md5file <- paste0(outfile, ".md5")
  cat(sprintf("%s  %s\n", unname(md5), basename(outfile)), file = md5file)
}

cat(sprintf("\u2713 Definitive mouse kinome table saved to: %s (%d kinases)\n", outfile, nrow(bm)), file = stderr())
