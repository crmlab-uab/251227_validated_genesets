# fetch_kinome.R
# Author: bRNA3F AI Agent
# Created: 2025-12-27 02:40 CST
# Purpose: Fetch and annotate the definitive kinome for mouse or human (protein, metabolic, lipid kinases) from authoritative sources.
# Usage: Rscript fetch_kinome.R --species mouse|human

library(optparse)
library(data.table)
library(biomaRt)

# Helper for error reporting
check_file_exists <- function(f, msg=NULL) {
  if (!file.exists(f)) stop(ifelse(is.null(msg), paste0('Missing required file: ', f), msg), call.=FALSE)
}
check_file_nonempty <- function(f, msg=NULL) {
  if (!file.exists(f) || file.info(f)$size == 0) stop(ifelse(is.null(msg), paste0('File missing or empty: ', f), msg), call.=FALSE)
}



# Robust repo root detection
cmdArgs <- commandArgs(trailingOnly = FALSE)
fileArgIdx <- grep("--file=", cmdArgs)
fileArg <- if (length(fileArgIdx) > 0) sub("--file=", "", cmdArgs[fileArgIdx[1]]) else NA_character_
script_dir <- if (!is.na(fileArg) && nzchar(fileArg)) dirname(normalizePath(fileArg)) else normalizePath(getwd())
find_repo_root <- function(start_dir) {
  cur <- normalizePath(start_dir)
  for (i in 1:10) {
    if (file.exists(file.path(cur, "genesets_config.yaml"))) return(cur)
    next_cur <- dirname(cur)
    if (next_cur == cur) break
    cur <- next_cur
  }
  stop("Could not find repository root from ", start_dir)
}
repo_root <- find_repo_root(script_dir)

# Always read species from YAML first
suppressMessages(library(yaml))
config_file <- Sys.getenv('KINASES_CONFIG', unset = file.path(repo_root, 'genesets_config.yaml'))
cat(sprintf("[DEBUG] Checking config file: %s\n", config_file))
species <- NULL
if (file.exists(config_file)) {
  cfg <- yaml::read_yaml(config_file)
  cat(sprintf("[DEBUG] cfg$species: %s\n", as.character(cfg$species)))
  if (!is.null(cfg$species)) species <- tolower(cfg$species)
} else {
  cat(sprintf("[DEBUG] Config file not found: %s\n", config_file))
}
flush.console()

# Then override with command-line argument if provided
library(optparse)
option_list <- list(
  make_option(c("-s", "--species"), type="character", default=NULL, help="Species: mouse or human", metavar="character")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (!is.null(opt$species)) {
  cat(sprintf("[DEBUG] Overriding species from argument: %s\n", opt$species))
  species <- tolower(opt$species)
}
cat(sprintf("[DEBUG] Final species value: %s\n", as.character(species)))
flush.console()
if (is.null(species) || !(species %in% c("mouse", "human"))) {
  stop("--species must be 'mouse' or 'human' (argument or YAML config)")
}

cat(sprintf("Fetching kinome for: %s\n", species), file=stderr())

# Set BioMart dataset and output file (write snapshots to canonical inputs folder)
if (species == "mouse") {
  dataset <- "mmusculus_gene_ensembl"
  date_tag <- format(Sys.time(), "%y%m%d")
  outdir <- file.path(repo_root, "curated","kinases","inputs")
  dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
  if (!dir.exists(outdir)) stop(paste0('Failed to create output directory: ', outdir))
  # canonical stable filename (used by downstream steps) and a timestamped snapshot
  outfile_canonical <- file.path(outdir, "kinases_mouse.csv")
  outfile_ts <- file.path(outdir, paste0("mouse_kinome__", date_tag, ".csv"))
  id_cols <- c(
    "ensembl_gene_id",        # Ensembl Gene ID
    "external_gene_name",     # Gene Symbol
    "description",            # Gene Description
    "go_id"                   # GO term (for reference)
  )
} else {
  dataset <- "hsapiens_gene_ensembl"
  date_tag <- format(Sys.time(), "%y%m%d")
  outdir <- file.path(repo_root, "curated","kinases","inputs")
  dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
  if (!dir.exists(outdir)) stop(paste0('Failed to create output directory: ', outdir))
  # Always write canonical output for downstream scripts
  outfile_canonical <- file.path(outdir, "kinases_human.csv")
  outfile_ts <- file.path(outdir, paste0("kinases_human__", date_tag, ".csv"))
  id_cols <- c(
    "ensembl_gene_id",        # Ensembl Gene ID
    "external_gene_name",     # Gene Symbol
    "description",            # Gene Description
    "go_id",                  # GO term (for reference)
    "hgnc_symbol",            # HGNC Symbol
    "hgnc_id"                 # HGNC ID
  )
}


# 1. Fetch all protein-coding genes with kinase activity (GO:0004672, GO:0004674, GO:0016301, GO:0016773)
mart <- useMart("ensembl", dataset=dataset)
kinase_go_terms <- c("GO:0004672", "GO:0004674", "GO:0016301", "GO:0016773")
cat("Querying BioMart for kinase genes...\n", file=stderr())
all_kinases <- getBM(attributes = c(id_cols, "go_id"),
                     filters = "go",
                     values = kinase_go_terms,
                     mart = mart)
setDT(all_kinases)
all_kinases <- unique(all_kinases[, ..id_cols])

# 2. Optionally merge with KinHub, Coral, Manning, UniProt lists (future extension)
# (For now, BioMart GO-based is the most systematic and up-to-date)

# 3. Post-process: collapse GO IDs per gene, parse MGI ID for mouse
if (species == "mouse") {
  # Collapse GO IDs per gene
  all_kinases[, go_id := as.character(go_id)]
  collapsed <- all_kinases[, .(
    go_id = paste(sort(unique(go_id[!is.na(go_id) & go_id != ""])), collapse = ";"),
    description = unique(description),
    external_gene_name = unique(external_gene_name)
  ), by = ensembl_gene_id]
  # Filter: remove blank gene names, Rik, Gm genes
  collapsed <- collapsed[external_gene_name != "" &
                         !grepl("^(Gm[0-9]+|.*Rik)$", external_gene_name, ignore.case=TRUE)]
  # Parse MGI ID from description (format: [Source:MGI Symbol;Acc:MGI:number])
  collapsed[, mgi_id := sub('.*\\[Source:MGI Symbol;Acc:(MGI:[0-9]+)\\].*', '\\1', description)]
  # Remove the [Source:MGI Symbol;Acc:MGI:number] suffix from description
  collapsed[, description := sub(' ?\\[Source:MGI Symbol;Acc:MGI:[0-9]+\\]$', '', description)]
  # Remove 'MGI:' prefix from mgi_id
  collapsed[, mgi_id := sub('^MGI:', '', mgi_id)]
  # Sort by external_gene_name
  setorder(collapsed, external_gene_name)
  # Reorder columns
  setcolorder(collapsed, c("ensembl_gene_id", "external_gene_name", "description", "mgi_id", "go_id"))
  # write canonical file for downstream use and also keep a dated snapshot
  # Always write canonical output for downstream scripts
  fwrite(collapsed, outfile_canonical)
  fwrite(collapsed, outfile_ts)
  for (f in c(outfile_canonical, outfile_ts)) {
    check_file_nonempty(f, paste0('Output file missing or empty: ', f, '\nCheck BioMart connection and output path.'))
  }
  # write md5 checksum for canonical file
  if (requireNamespace("tools", quietly=TRUE)) {
    md5 <- tools::md5sum(outfile_canonical)
    md5file <- paste0(outfile_canonical, ".md5")
    cat(sprintf("%s  %s\n", unname(md5), basename(outfile_canonical)), file=md5file)
  }
  cat(sprintf("Done. Outputs: %s and %s (%d unique kinases)\n", outfile_canonical, outfile_ts, nrow(collapsed)), file=stderr())
  # Also write a numbered output copy for downstream steps (01_ prefix)
  outdir_out <- file.path(repo_root, "genesets","curated","kinases","outputs")
  dir.create(outdir_out, recursive=TRUE, showWarnings=FALSE)
  numbered <- file.path(outdir_out, paste0("01_", basename(outfile_canonical)))
  file.copy(outfile_canonical, numbered, overwrite=TRUE)
  if (requireNamespace("tools", quietly=TRUE)) {
    md5 <- tools::md5sum(numbered)
    cat(sprintf("%s  %s\n", unname(md5), basename(numbered)), file=paste0(numbered, ".md5"))
  }
} else {
  # For human: collapse GO IDs per gene
  all_kinases[, go_id := as.character(go_id)]
  collapsed <- all_kinases[, .(
    go_id = paste(sort(unique(go_id[!is.na(go_id) & go_id != ""])), collapse = ";"),
    description = unique(description),
    external_gene_name = unique(external_gene_name),
    hgnc_symbol = unique(hgnc_symbol),
    hgnc_id = unique(hgnc_id)
  ), by = ensembl_gene_id]
  # Filter: remove blank gene names, Rik, Gm genes
  collapsed <- collapsed[external_gene_name != "" &
                         !grepl("^(Gm[0-9]+|.*Rik)$", external_gene_name, ignore.case=TRUE)]
  # Parse HGNC ID from description (format: [Source:HGNC Symbol;Acc:HGNC:number])
  collapsed[, hgnc_id := sub('.*\\[Source:HGNC Symbol;Acc:HGNC:([0-9]+)\\].*', '\\1', description)]
  # Remove the [Source:HGNC Symbol;Acc:HGNC:number] suffix from description
  collapsed[, description := sub(' ?\\[Source:HGNC Symbol;Acc:HGNC:[0-9]+\\]$', '', description)]
  # Sort by external_gene_name
  setorder(collapsed, external_gene_name)
  setcolorder(collapsed, c("ensembl_gene_id", "external_gene_name", "description", "hgnc_symbol", "hgnc_id", "go_id"))
  # for human, write canonical file and a dated snapshot (both with kinases_human* root)
  fwrite(collapsed, outfile_canonical)
  fwrite(collapsed, outfile_ts)
  for (f in c(outfile_canonical, outfile_ts)) {
    check_file_nonempty(f, paste0('Output file missing or empty: ', f, '\nCheck BioMart connection and output path.'))
  }
  # write md5 checksum for canonical file
  if (requireNamespace("tools", quietly=TRUE)) {
    md5 <- tools::md5sum(outfile_canonical)
    md5file <- paste0(outfile_canonical, ".md5")
    cat(sprintf("%s  %s\n", unname(md5), basename(outfile_canonical)), file=md5file)
  }
  cat(sprintf("Done. Outputs: %s and %s (%d unique kinases)\n", outfile_canonical, outfile_ts, nrow(collapsed)), file=stderr())
}
