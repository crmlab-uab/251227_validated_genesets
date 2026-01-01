#!/usr/bin/env Rscript
# validate_against_genome.R
# Author: C. Ryan Miller, MD, PhD (rmiller@uab.edu)
# Updated: 2025-12-31 (refactored to use config_loader.R)
# Purpose: Validate kinase list against downloaded mouse genome annotations
# Usage: Rscript lib/validate_against_genome.R

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

cat("=== Kinase Genome Validation ===\n\n")

# Find kinase input file
candidates <- list.files(paths$input_dir, pattern='201006_composite_kinases_curated.*\\.csv$', full.names=TRUE, ignore.case=TRUE)
if (length(candidates) == 0) stop('Missing input snapshot: place 201006_composite_kinases_curated__YYMMDD.csv in ', paths$input_dir)
kinase_file <- sort(candidates, decreasing=TRUE)[1]
gtf_file <- "/data/bRNA3F/data/genome/mouse/gencode.vM37.primary_assembly.annotation.gtf.gz"

# Load kinase list
cat("Loading kinase list...\n")
kinases <- fread(kinase_file, header = TRUE)
cat("  Found", nrow(kinases), "kinases\n")

# Filter out excluded genes from config
exclude_genes <- get_exclude_genes()
if (length(exclude_genes) > 0) {
  # Find the symbol column first
  sym_col <- intersect(c("Mouse_symbol", "Mouse_Symbol", "external_gene_name", "symbol"), names(kinases))[1]
  if (!is.na(sym_col)) {
    n_before <- nrow(kinases)
    excluded_found <- kinases[[sym_col]][kinases[[sym_col]] %in% exclude_genes]
    kinases <- kinases[!kinases[[sym_col]] %in% exclude_genes]
    n_excluded <- n_before - nrow(kinases)
    if (n_excluded > 0) {
      cat("  Excluded", n_excluded, "genes per config:", paste(excluded_found, collapse=", "), "\n")
    }
  }
}
cat("  Processing", nrow(kinases), "kinases\n\n")

# Parse GTF to extract gene information
cat("Parsing genome GTF (this may take 30-60 seconds)...\n")
gtf_genes <- fread(
  cmd = paste0("zcat ", gtf_file, " | grep -v '^#' | awk '$3==\"gene\"'"),
  sep = "\t",
  header = FALSE,
  col.names = c("chr", "source", "type", "start", "end", "score", "strand", "frame", "attributes")
)

cat("  Found", nrow(gtf_genes), "genes in genome\n\n")

# Extract gene attributes
cat("Extracting gene names and MGI IDs...\n")
extract_attr <- function(attr_string, key) {
  pattern <- paste0(key, ' "([^\"]+)"')
  match <- regmatches(attr_string, regexpr(pattern, attr_string, perl = TRUE))
  if (length(match) == 0) return(NA)
  gsub(paste0(key, ' "([^\"]+)"'), "\\1", match, perl = TRUE)
}

gtf_genes[, gene_symbol := sapply(attributes, function(x) extract_attr(x, "gene_name"))]
gtf_genes[, gene_id := sapply(attributes, function(x) extract_attr(x, "gene_id"))]
gtf_genes[, gene_type := sapply(attributes, function(x) extract_attr(x, "gene_type"))]
gtf_genes[, mgi_id := sapply(attributes, function(x) extract_attr(x, "mgi_id"))]

# Remove version numbers from Ensembl IDs
gtf_genes[, gene_id_clean := gsub("\\..*$", "", gene_id)]

cat("  Extracted", sum(!is.na(gtf_genes$gene_symbol)), "gene symbols\n")
cat("  Extracted", sum(!is.na(gtf_genes$mgi_id)), "MGI IDs\n\n")

# Create lookup table (one row per gene symbol)
genome_lookup <- gtf_genes[!is.na(gene_symbol), 
  .(gene_id = gene_id_clean[1], 
    gene_type = gene_type[1],
    mgi_id = mgi_id[1],
    chr = chr[1]),
  by = gene_symbol
]

cat("=== Validation Results ===\n\n")

# Find symbol column (handles different naming conventions)
symbol_col <- intersect(c("Mouse_symbol", "Mouse_Symbol", "external_gene_name", "symbol"), names(kinases))[1]
mgi_col <- intersect(c("MGI_ID", "MGI_Mouse", "mgi_id"), names(kinases))[1]

if (is.na(symbol_col)) {
  stop("No mouse symbol column found in input file. Available columns: ", paste(names(kinases), collapse=", "))
}
cat("Using column:", symbol_col, "for gene symbols\n")

# Validate each kinase
symbols <- kinases[[symbol_col]]
validation_results <- data.table(
  Mouse_Symbol = symbols,
  In_Genome = symbols %in% genome_lookup$gene_symbol,
  Ensembl_ID = NA_character_,
  Gene_Type = NA_character_,
  MGI_ID_Genome = NA_character_,
  MGI_ID_Match = NA,
  Chromosome = NA_character_,
  Status = NA_character_
)

for (i in 1:nrow(kinases)) {
  symbol <- symbols[i]
  mgi_curated <- if (!is.na(mgi_col)) kinases[[mgi_col]][i] else NA
  
  if (length(symbol) > 0 && !is.na(symbol) && symbol != "" && symbol %in% genome_lookup$gene_symbol) {
    genome_data <- genome_lookup[gene_symbol == symbol]
    validation_results[i, Ensembl_ID := genome_data$gene_id[1]]
    validation_results[i, Gene_Type := genome_data$gene_type[1]]
    validation_results[i, MGI_ID_Genome := genome_data$mgi_id[1]]
    validation_results[i, Chromosome := genome_data$chr[1]]
    
    # Check MGI ID match
    if (!is.na(mgi_curated) && mgi_curated != "" && !is.na(genome_data$mgi_id[1])) {
      validation_results[i, MGI_ID_Match := (mgi_curated == genome_data$mgi_id[1])]
    }
    
    # Determine status
    if (genome_data$gene_type[1] == "protein_coding") {
      validation_results[i, Status := "VALIDATED"]
    } else {
      validation_results[i, Status := paste0("WARNING: ", genome_data$gene_type[1])]
    }
  } else {
    validation_results[i, Status := "NOT_FOUND"]
  }
  
  # Progress
  if (i %% 50 == 0) cat("  Processed", i, "/", nrow(kinases), "kinases...\n")
}

cat("\n=== Summary ===\n\n")
cat("In genome:", sum(validation_results$In_Genome), "/", nrow(validation_results), "\n")
cat("Protein coding:", sum(validation_results$Gene_Type == "protein_coding", na.rm = TRUE), "\n")
cat("Not found:", sum(validation_results$Status == "NOT_FOUND"), "\n")
cat("Non-protein-coding:", sum(grepl("WARNING", validation_results$Status)), "\n\n")

# Show genes not found in genome
not_found <- validation_results[Status == "NOT_FOUND"]
if (nrow(not_found) > 0) {
  cat("=== Genes NOT FOUND in genome ===\n")
  # Use dynamic column selection
  not_found_kinases <- kinases[kinases[[symbol_col]] %in% not_found$Mouse_Symbol]
  print(not_found_kinases[, .SD, .SDcols = intersect(names(not_found_kinases),
        c(symbol_col, "HGNC", "gene_name", "Group", mgi_col, "Entrez_ID_mouse"))])
  cat("\n")
}

# Show non-protein-coding genes
non_coding <- validation_results[grepl("WARNING", Status)]
if (nrow(non_coding) > 0) {
  cat("=== Non-protein-coding genes ===\n")
  print(non_coding[, .(Mouse_Symbol, Gene_Type, Status)])
  cat("\n")
}

# Show MGI ID mismatches
mgi_mismatch <- validation_results[!is.na(MGI_ID_Match) & MGI_ID_Match == FALSE]
if (nrow(mgi_mismatch) > 0 && !is.na(mgi_col)) {
  cat("=== MGI ID mismatches ===\n")
  mismatch_kinases <- kinases[kinases[[symbol_col]] %in% mgi_mismatch$Mouse_Symbol]
  mismatch_data <- data.table(
    Mouse_Symbol = mismatch_kinases[[symbol_col]],
    MGI_Curated = mismatch_kinases[[mgi_col]]
  )
  mismatch_data <- merge(mismatch_data, validation_results[, .(Mouse_Symbol, MGI_ID_Genome)], by = "Mouse_Symbol")
  print(mismatch_data)
  cat("\n")
}

# Save full validation results
output_file <- output_path("genome_validation_results.csv")
full_results <- cbind(kinases, validation_results[, -1])
fwrite(full_results, output_file)
cat("✓ Full results saved to:", output_file, "\n\n")

# Special check: Find Cilk1/ICK missing data
cilk1_data <- genome_lookup[gene_symbol == "Cilk1"]
if (nrow(cilk1_data) > 0) {
  cat("=== Cilk1 (ICK) genome data ===\n")
  print(cilk1_data)
  cat("\n")
}

cat("✓ Genome validation complete!\n")
