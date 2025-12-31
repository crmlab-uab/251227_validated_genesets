#!/usr/bin/env Rscript
# Validate kinase list against downloaded mouse genome annotations
# Author: C. Ryan Miller, MD, PhD (rmiller@uab.edu)
# Created: 2025-12-27 03:10 CST

library(data.table)

cat("=== Kinase Genome Validation ===\n\n")

# Paths
inputs_dir <- 'curated/kinases/inputs'
candidates <- list.files(inputs_dir, pattern='201006_composite_kinases_curated.*\\.csv$', full.names=TRUE, ignore.case=TRUE)
if (length(candidates) == 0) stop('Missing input snapshot: please place 201006_composite_kinases_curated__YYMMDD.csv in ', inputs_dir)
kinase_file <- sort(candidates, decreasing=TRUE)[1]
gtf_file <- "/data/bRNA3F/data/genome/mouse/gencode.vM37.primary_assembly.annotation.gtf.gz"

# Load kinase list
cat("Loading kinase list...\n")
kinases <- fread(kinase_file, header = TRUE)
cat("  Found", nrow(kinases), "kinases\n\n")

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

# Validate each kinase
validation_results <- data.table(
  Mouse_Symbol = kinases$Mouse_Symbol,
  In_Genome = kinases$Mouse_Symbol %in% genome_lookup$gene_symbol,
  Ensembl_ID = NA_character_,
  Gene_Type = NA_character_,
  MGI_ID_Genome = NA_character_,
  MGI_ID_Match = NA,
  Chromosome = NA_character_,
  Status = NA_character_
)

for (i in 1:nrow(kinases)) {
  symbol <- kinases$Mouse_Symbol[i]
  mgi_curated <- kinases$MGI_Mouse[i]
  
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
  print(kinases[Mouse_Symbol %in% not_found$Mouse_Symbol, 
                .(Mouse_Symbol, Human_Symbol, Description, Group, MGI_Mouse, Entrez_Mouse)])
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
if (nrow(mgi_mismatch) > 0) {
  cat("=== MGI ID mismatches ===\n")
  mismatch_data <- kinases[Mouse_Symbol %in% mgi_mismatch$Mouse_Symbol, .(Mouse_Symbol, MGI_Mouse)]
  mismatch_data <- merge(mismatch_data, validation_results[, .(Mouse_Symbol, MGI_ID_Genome)], by = "Mouse_Symbol")
  print(mismatch_data)
  cat("\n")
}

# Save full validation results
output_file <- "/data/251227_validated_genesets/kinases/genome_validation_results.csv"
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
