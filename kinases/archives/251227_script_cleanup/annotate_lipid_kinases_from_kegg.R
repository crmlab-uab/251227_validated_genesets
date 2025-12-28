library(data.table)
library(KEGGREST)
library(dplyr)
library(org.Hs.eg.db)
library(AnnotationDbi)

# Load kinome table (with Metabolic column already present)
kinome_file <- "human_kinome_with_group_family_GOlast_metabolic_lipid_flags.csv"
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
fwrite(dt, "human_kinome_with_group_family_GOlast_metabolic_lipid_flags.csv")

cat("Annotated lipid kinases using KEGG. Output file: human_kinome_with_group_family_GOlast_metabolic_lipid_flags.csv\n")
