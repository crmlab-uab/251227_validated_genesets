# 01_fetch_geneset_BioMart.R
# Author: bRNA3F AI Agent
# Created: 2025-12-27 02:40 CST
# Updated: 2025-12-31 (refactored to use config_loader.R)
# Purpose: Fetch and annotate the definitive kinome for mouse or human from BioMart
# Usage: Rscript 01_fetch_geneset_BioMart.R --species mouse|human

library(optparse)
library(data.table)
library(biomaRt)

# Load centralized config (provides repo_root, paths, cfg, helper functions)
# Determine script directory for sourcing config_loader.R
.script_dir <- (function() {
 cmd_args <- commandArgs(trailingOnly = FALSE)
 file_arg <- grep("--file=", cmd_args, value = TRUE)
 if (length(file_arg) > 0) {
   return(dirname(normalizePath(sub("--file=", "", file_arg[1]))))
 }
 return(normalizePath("."))
})()
source(file.path(.script_dir, "lib", "config_loader.R"))

# Parse command-line arguments (override config species if provided)
option_list <- list(
  make_option(c("-s", "--species"), type = "character", default = NULL,
              help = "Species: mouse or human", metavar = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))

# Determine species: command-line > config > error
species <- if (!is.null(opt$species)) tolower(opt$species) else cfg_get("species", NULL)
if (is.null(species) || !(species %in% c("mouse", "human"))) {
  stop("--species must be 'mouse' or 'human' (via argument or config)")
}

cat(sprintf("Fetching kinome for: %s\n", species), file = stderr())

# Set BioMart dataset and output paths
date_tag <- format(Sys.time(), "%y%m%d")
ensure_dir(paths$input_dir)

if (species == "mouse") {
  dataset <- "mmusculus_gene_ensembl"
  outfile_canonical <- input_path("kinases_mouse_biomart.csv")
  outfile_ts <- input_path(paste0("kinases_mouse_biomart__", date_tag, ".csv"))
  id_cols <- c("ensembl_gene_id", "external_gene_name", "description", "go_id")
} else {
  dataset <- "hsapiens_gene_ensembl"
  outfile_canonical <- input_path("kinases_human_biomart.csv")
  outfile_ts <- input_path(paste0("kinases_human_biomart__", date_tag, ".csv"))
  id_cols <- c("ensembl_gene_id", "external_gene_name", "description", "go_id", "hgnc_symbol", "hgnc_id")
}

# Fetch all protein-coding genes with kinase activity
kinase_go_terms <- c("GO:0004672", "GO:0004674", "GO:0016301", "GO:0016773")
cat("Querying BioMart for kinase genes...\n", file = stderr())

mart <- useMart("ensembl", dataset = dataset)
all_kinases <- getBM(
  attributes = id_cols,
  filters = "go",
  values = kinase_go_terms,
  mart = mart
)
setDT(all_kinases)
all_kinases <- unique(all_kinases[, ..id_cols])

# Post-process: collapse GO IDs per gene
if (species == "mouse") {
  all_kinases[, go_id := as.character(go_id)]
  collapsed <- all_kinases[, .(
    go_id = paste(sort(unique(go_id[!is.na(go_id) & go_id != ""])), collapse = ";"),
    description = unique(description),
    external_gene_name = unique(external_gene_name)
  ), by = ensembl_gene_id]

  # Filter out blank gene names, Rik, Gm genes
  collapsed <- collapsed[
    external_gene_name != "" &
    !grepl("^(Gm[0-9]+|.*Rik)$", external_gene_name, ignore.case = TRUE)
  ]

  # Parse MGI ID from description
  collapsed[, mgi_id := sub('.*\\[Source:MGI Symbol;Acc:(MGI:[0-9]+)\\].*', '\\1', description)]
  collapsed[, description := sub(' ?\\[Source:MGI Symbol;Acc:MGI:[0-9]+\\]$', '', description)]
  collapsed[, mgi_id := sub('^MGI:', '', mgi_id)]

  setorder(collapsed, external_gene_name)
  setcolorder(collapsed, c("ensembl_gene_id", "external_gene_name", "description", "mgi_id", "go_id"))
} else {
  # Human
  all_kinases[, go_id := as.character(go_id)]
  collapsed <- all_kinases[, .(
    go_id = paste(sort(unique(go_id[!is.na(go_id) & go_id != ""])), collapse = ";"),
    description = unique(description),
    external_gene_name = unique(external_gene_name),
    hgnc_symbol = unique(hgnc_symbol),
    hgnc_id = unique(hgnc_id)
  ), by = ensembl_gene_id]

  # Filter out blank gene names
  collapsed <- collapsed[
    external_gene_name != "" &
    !grepl("^(Gm[0-9]+|.*Rik)$", external_gene_name, ignore.case = TRUE)
  ]

  # Parse HGNC ID from description
  collapsed[, hgnc_id := sub('.*\\[Source:HGNC Symbol;Acc:HGNC:([0-9]+)\\].*', '\\1', description)]
  collapsed[, description := sub(' ?\\[Source:HGNC Symbol;Acc:HGNC:[0-9]+\\]$', '', description)]

  setorder(collapsed, external_gene_name)
  setcolorder(collapsed, c("ensembl_gene_id", "external_gene_name", "description", "hgnc_symbol", "hgnc_id", "go_id"))
}

# Write outputs to canonical inputs directory
fwrite(collapsed, outfile_canonical)
fwrite(collapsed, outfile_ts)

# Validate outputs
for (f in c(outfile_canonical, outfile_ts)) {
  check_file_nonempty(f, paste0("Output file missing or empty: ", f, "\nCheck BioMart connection."))
}

# Write MD5 checksum
if (requireNamespace("tools", quietly = TRUE)) {
  md5 <- tools::md5sum(outfile_canonical)
  md5file <- paste0(outfile_canonical, ".md5")
  cat(sprintf("%s  %s\n", unname(md5), basename(outfile_canonical)), file = md5file)
}

cat(sprintf("Done. Outputs:\n  %s\n  %s\n  (%d unique kinases)\n",
            outfile_canonical, outfile_ts, nrow(collapsed)), file = stderr())
