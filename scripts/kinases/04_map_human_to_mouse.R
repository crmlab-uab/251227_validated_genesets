# Author: C. Ryan Miller
# Created: 2025-12-28 02:17 CST
# Commit: 26ec675324c74c530b55664519f79329a3d427d8

# Fetch definitive mouse kinome from authoritative sources
# Author: C. Ryan Miller
# Created: 2025-12-28 02:17 CST
# Commit: 26ec675324c74c530b55664519f79329a3d427d8

# Fetch definitive mouse kinome from authoritative sources
#
# NOTE: This is a thin runner. All mapping logic is consolidated in scripts/kinases/lib/ (see map_human_to_mouse_uniprot.R, merge_kinase_uniprot_validation.R).
# Do not duplicate mapping logic here; update lib/ scripts for mapping changes.
# Output: mouse_kinome_definitive.csv with Mouse_Symbol, Ensembl_Gene_ID, Entrez_ID, UniProt_ID, MGI_ID, Description

library(data.table)
library(biomaRt)
library(httr)

# Helper for error reporting
check_file_exists <- function(f, msg=NULL) {
  if (!file.exists(f)) stop(ifelse(is.null(msg), paste0('Missing required file: ', f), msg), call.=FALSE)
}
check_file_nonempty <- function(f, msg=NULL) {
  if (!file.exists(f) || file.info(f)$size == 0) stop(ifelse(is.null(msg), paste0('File missing or empty: ', f), msg), call.=FALSE)
}

# 1. Get kinase list from KinHub (curated, gold standard)
kinhub_url <- "http://www.kinhub.org/kinases.html"
# For automation, use the KinHub CSV (if available) or curated list from literature
# Here, use the 201006_composite_kinases_curated CSV snapshot from inputs as KinHub gold standard
inputs_dir <- '../../../curated/kinases/inputs'
candidates <- list.files(inputs_dir, pattern='201006_composite_kinases_curated.*\\.csv$', full.names=TRUE, ignore.case=TRUE)
if (length(candidates) == 0) stop('Missing input snapshot: place 201006_composite_kinases_curated__YYMMDD.csv in ', inputs_dir)
curated_file <- sort(candidates, decreasing=TRUE)[1]
check_file_nonempty(curated_file, paste0('Curated kinases input missing or empty: ', curated_file))
curated <- fread(curated_file)

mouse_symbols <- unique(curated$Mouse_symbol)

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
         c("Mouse_symbol", "Ensembl_Gene_ID", "Entrez_ID", "UniProt_ID", "MGI_ID", "Description"))

# 3. Merge with KinHub curated to ensure only true kinases
bm <- bm[Mouse_symbol %in% mouse_symbols]


# 4. Remove duplicates, keep one row per Mouse_symbol (prefer with all IDs present)
setorder(bm, Mouse_symbol, -Ensembl_Gene_ID, -UniProt_ID, -MGI_ID, -Entrez_ID)
bm <- bm[!duplicated(Mouse_symbol)]

# 5. Output
outfile <- "mouse_kinome_definitive.csv"
fwrite(bm, outfile)
check_file_nonempty(outfile, paste0('Output file missing or empty: ', outfile))
cat("\u2713 Definitive mouse kinome table saved to mouse_kinome_definitive.csv\n")
