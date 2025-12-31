# Author: C. Ryan Miller
# Created: 2025-12-28 02:17 CST
# Commit: 26ec675324c74c530b55664519f79329a3d427d8

library(data.table)
library(KEGGREST)
library(dplyr)
library(org.Hs.eg.db)
library(AnnotationDbi)


# Load config and set input/output dirs from YAML if available
library(yaml)
config_file <- Sys.getenv('KINASES_CONFIG', unset = 'genesets_config.yaml')
if (file.exists(config_file)) {
  cfg <- yaml::read_yaml(config_file)
  input_dir <- if (!is.null(cfg$input_dir)) cfg$input_dir else 'curated/kinases/inputs'
  output_dir <- if (!is.null(cfg$output_dir)) cfg$output_dir else 'curated/kinases/outputs'
} else {
  input_dir <- 'curated/kinases/inputs'
  output_dir <- 'curated/kinases/outputs'
}
kinome_file <- file.path(input_dir, 'human_kinome_with_group_family_GOlast_metabolic_lipid_flags.csv')
dt <- fread(kinome_file)


# Get all human lipid metabolism pathway IDs from KEGG
lipid_pathways <- keggList("pathway", "hsa")
lipid_ids <- names(lipid_pathways)[grepl("lipid", lipid_pathways, ignore.case = TRUE)]

# Get all human genes in lipid metabolism pathways (Entrez IDs)
lipid_genes <- unique(unlist(lapply(lipid_ids, function(pid) {
  kegg_genes <- keggGet(pid)[[1]]$GENE
  if (is.null(kegg_genes)) return(NULL)
  gsub(" .*", "", kegg_genes)
})))

# Map KEGG Entrez IDs to HGNC symbols using org.Hs.eg.db
entrez_to_hgnc <- AnnotationDbi::select(org.Hs.eg.db, keys=lipid_genes, keytype="ENTREZID", columns=c("SYMBOL"))
kegg_lipid_hgnc <- unique(na.omit(entrez_to_hgnc$SYMBOL))

# Add/overwrite Lipid column (Y/N) using robust mapping
dt$Lipid <- ifelse(
  dt$hgnc_symbol %in% kegg_lipid_hgnc |
  dt$external_gene_name %in% kegg_lipid_hgnc,
  "Y", "N"
)


# Move go_id to end
if ("go_id" %in% names(dt)) {
  go_col <- dt$go_id
  dt$go_id <- NULL
  dt$go_id <- go_col
}
output_file <- file.path(output_dir, 'human_kinome_with_group_family_GOlast_metabolic_lipid_flags.csv')
fwrite(dt, output_file)
cat("Annotated lipid kinases using KEGG. Output file:", output_file, "\n")
