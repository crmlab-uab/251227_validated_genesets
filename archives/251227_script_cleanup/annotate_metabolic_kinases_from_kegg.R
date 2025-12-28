# Author: C. Ryan Miller
# Created: 2025-12-28 02:17 CST
# Commit: 26ec675324c74c530b55664519f79329a3d427d8

library(data.table)
library(KEGGREST)
library(dplyr)
library(org.Hs.eg.db)
library(AnnotationDbi)


# Load kinome table

kinome_file <- "human_kinome_with_group_family_GOlast_custom2.csv"
dt <- fread(kinome_file)

# Rename 'Custom' column to 'Protein' if present
if ("Custom" %in% names(dt)) {
  names(dt)[names(dt) == "Custom"] <- "Protein"
}
# Set Protein column to Y for 'Protein kinase', N otherwise
if ("Protein" %in% names(dt)) {
  dt$Protein <- ifelse(dt$Protein == "Protein kinase", "Y", "N")
}

# Get all human metabolic pathway IDs from KEGG
metabolic_pathways <- keggList("pathway", "hsa")
metabolic_ids <- names(metabolic_pathways)[grepl("Metabolic pathways|metabolism", metabolic_pathways, ignore.case = TRUE)]

# Get all human genes in metabolic pathways (Entrez IDs)
metabolic_genes <- unique(unlist(lapply(metabolic_ids, function(pid) {
  kegg_genes <- keggGet(pid)[[1]]$GENE
  if (is.null(kegg_genes)) return(NULL)
  # KEGG GENE is a vector: c("1234 description", ...)
  gsub(" .*", "", kegg_genes)
})))

# Map KEGG Entrez IDs to HGNC symbols using org.Hs.eg.db
entrez_to_hgnc <- AnnotationDbi::select(org.Hs.eg.db, keys=metabolic_genes, keytype="ENTREZID", columns=c("SYMBOL"))
kegg_metabolic_hgnc <- unique(na.omit(entrez_to_hgnc$SYMBOL))

# Add Metabolic column (Y/N) using robust mapping
dt$Metabolic <- ifelse(
  dt$hgnc_symbol %in% kegg_metabolic_hgnc |
  dt$external_gene_name %in% kegg_metabolic_hgnc,
  "Y", "N"
)

# Prepare for Lipid annotation (optional, for downstream script)
if (!"Lipid" %in% names(dt)) {
  dt$Lipid <- "N"
}


# Move go_id to end
if ("go_id" %in% names(dt)) {
  go_col <- dt$go_id
  dt$go_id <- NULL
  dt$go_id <- go_col
}
fwrite(dt, "human_kinome_with_group_family_GOlast_metabolic_lipid_flags.csv")

cat("Annotated metabolic kinases using KEGG. Output file: human_kinome_with_group_family_GOlast_metabolic_lipid_flags.csv\n")
