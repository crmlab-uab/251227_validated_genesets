# Clean and enrich comprehensive mouse kinome table
# 1. Parse MGI ID from Description
# 2. Fill in UniProt_ID if missing (using biomaRt)
# 3. Optionally remove GO columns

library(data.table)
library(biomaRt)

# Load data
kinases <- fread('comprehensive_mouse_kinome_biomart.csv')

# 1. Fetch MGI_IDs for all Ensembl_Gene_IDs using biomaRt
library(biomaRt)
# Debug: print a sample of Ensembl_Gene_IDs and BioMart result
# Standardize Ensembl_Gene_ID and Mouse_Symbol in both tables before merging
ensembl <- useMart('ensembl', dataset='mmusculus_gene_ensembl')
kinases[, Ensembl_Gene_ID := trimws(as.character(Ensembl_Gene_ID))]
kinases[, Mouse_Symbol := trimws(as.character(Mouse_Symbol))]
# 1. Try mapping by Ensembl_Gene_ID
query_ids <- unique(kinases$Ensembl_Gene_ID)
bm1 <- getBM(attributes = c('ensembl_gene_id', 'mgi_id'),
            filters = 'ensembl_gene_id',
            values = query_ids,
            mart = ensembl)
bm1 <- as.data.table(bm1)
bm1[, ensembl_gene_id := trimws(as.character(ensembl_gene_id))]
setnames(bm1, c('ensembl_gene_id', 'mgi_id'), c('Ensembl_Gene_ID', 'MGI_ID_ensembl'))
# 2. Try mapping by Mouse_Symbol (external_gene_name)
query_symbols <- unique(kinases$Mouse_Symbol)
bm2 <- getBM(attributes = c('external_gene_name', 'mgi_id'),
            filters = 'external_gene_name',
            values = query_symbols,
            mart = ensembl)
bm2 <- as.data.table(bm2)
bm2[, external_gene_name := trimws(as.character(external_gene_name))]
setnames(bm2, c('external_gene_name', 'mgi_id'), c('Mouse_Symbol', 'MGI_ID_symbol'))
# Merge both mappings into kinases
kinases <- merge(kinases, bm1, by = 'Ensembl_Gene_ID', all.x = TRUE)
kinases <- merge(kinases, bm2, by = 'Mouse_Symbol', all.x = TRUE)
# Fill MGI_ID from Ensembl mapping, then symbol mapping
kinases[, MGI_ID := fifelse(!is.na(MGI_ID_ensembl) & MGI_ID_ensembl != '', MGI_ID_ensembl,
                     fifelse(!is.na(MGI_ID_symbol) & MGI_ID_symbol != '', MGI_ID_symbol, NA_character_))]
kinases[, c('MGI_ID_ensembl', 'MGI_ID_symbol') := NULL]
# Remove [ ... ] annotation from Description
kinases[, Description := gsub(' \\[.*\\]', '', Description)]
kinases[, MGI_ID := sub('.*Acc:(MGI:[0-9]+).*', '\1', Description)]
kinases[!grepl('^MGI:', MGI_ID), MGI_ID := NA]

# 2. Fill in missing UniProt_IDs using biomaRt
missing_uniprot <- kinases[is.na(UniProt_ID) | UniProt_ID == '', unique(Ensembl_Gene_ID)]
if (length(missing_uniprot) > 0) {
  ensembl <- useMart('ensembl', dataset='mmusculus_gene_ensembl')
  bm <- getBM(attributes = c('ensembl_gene_id', 'uniprotswissprot'),
              filters = 'ensembl_gene_id',
              values = missing_uniprot,
              mart = ensembl)
  bm <- as.data.table(bm)
  setnames(bm, 'uniprotswissprot', 'UniProt_ID_fill')
  bm <- bm[!is.na(UniProt_ID_fill) & UniProt_ID_fill != '', .(UniProt_ID_fill = first(UniProt_ID_fill)), by=ensembl_gene_id]
  kinases <- merge(kinases, bm, by.x='Ensembl_Gene_ID', by.y='ensembl_gene_id', all.x=TRUE)
  kinases[UniProt_ID == '' | is.na(UniProt_ID), UniProt_ID := UniProt_ID_fill]
  kinases[, UniProt_ID_fill := NULL]
}

# 3. Remove GO columns if not needed
kinases[, c('GO_ID','GO_Term') := NULL]

fwrite(kinases, 'comprehensive_mouse_kinome_cleaned.csv')
cat('✓ Cleaned comprehensive mouse kinome table saved to comprehensive_mouse_kinome_cleaned.csv\n')

# Deduplicate to unique gene list (one row per Ensembl_Gene_ID, Mouse_Symbol, Entrez_ID, UniProt_ID, Description, MGI_ID)
kinases_unique <- unique(kinases, by=c('Ensembl_Gene_ID','Mouse_Symbol','Entrez_ID','UniProt_ID','Description','MGI_ID'))
# Reorder columns: MGI_ID as 5th, Description last
setcolorder(kinases_unique, c('Ensembl_Gene_ID','Mouse_Symbol','Entrez_ID','UniProt_ID','MGI_ID','Description'))
fwrite(kinases_unique, 'comprehensive_mouse_kinome_unique.csv')
cat('✓ Unique gene list saved to comprehensive_mouse_kinome_unique.csv\\n')
