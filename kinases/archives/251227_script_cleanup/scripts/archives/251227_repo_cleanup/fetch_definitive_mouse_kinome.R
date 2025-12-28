# Author: C. Ryan Miller
# Created: 2025-12-28 02:17 CST
# Commit: 26ec675324c74c530b55664519f79329a3d427d8

# Fetch definitive mouse kinome from authoritative sources
# Output: mouse_kinome_definitive.csv with Mouse_Symbol, Ensembl_Gene_ID, Entrez_ID, UniProt_ID, MGI_ID, Description

library(data.table)
library(biomaRt)
library(httr)

# 1. Get kinase list from KinHub (curated, gold standard)
kinhub_url <- "http://www.kinhub.org/kinases.html"
# For automation, use the KinHub CSV (if available) or curated list from literature
# Here, use the 201006_composite_kinases_curated.csv as KinHub gold standard
curated <- fread("201006_composite_kinases_curated.csv")
mouse_symbols <- unique(curated$Mouse_Symbol)

# 2. Query Ensembl BioMart for all annotation IDs for these symbols
ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
bm <- getBM(
  attributes = c(
    "external_gene_name", "ensembl_gene_id", "entrezgene_id", "uniprotswissprot", "mgi_id", "description"
  ),
  filters = "external_gene_name",
  values = mouse_symbols,
  mart = ensembl
)
bm <- as.data.table(bm)
setnames(bm, c("external_gene_name", "ensembl_gene_id", "entrezgene_id", "uniprotswissprot", "mgi_id", "description"),
         c("Mouse_Symbol", "Ensembl_Gene_ID", "Entrez_ID", "UniProt_ID", "MGI_ID", "Description"))

# 3. Merge with KinHub curated to ensure only true kinases
bm <- bm[Mouse_Symbol %in% mouse_symbols]

# 4. Remove duplicates, keep one row per Mouse_Symbol (prefer with all IDs present)
setorder(bm, Mouse_Symbol, -(!is.na(Ensembl_Gene_ID)), -(!is.na(UniProt_ID)), -(!is.na(MGI_ID)), -(!is.na(Entrez_ID)))
bm <- bm[!duplicated(Mouse_Symbol)]

# 5. Output
fwrite(bm, "mouse_kinome_definitive.csv")
cat("✓ Definitive mouse kinome table saved to mouse_kinome_definitive.csv\n")