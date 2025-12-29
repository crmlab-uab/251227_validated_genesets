# Map mouse gene symbols to UniProt IDs using biomaRt and a local mapping table
# Author: C. Ryan Miller, MD, PhD (rmiller@uab.edu)

library(data.table)
library(biomaRt)

# Load kinase list
kinases <- fread('genesets/curated/kinases/201006_composite_kinases_curated.csv', header=TRUE)

# Query biomaRt for mouse gene symbol to UniProt mapping
cat('Querying biomaRt for mouse gene symbol to UniProt mapping...\n')
ensembl <- useMart('ensembl', dataset='mmusculus_gene_ensembl')
biomart_map <- getBM(
  attributes = c('mgi_symbol', 'uniprotswissprot', 'entrezgene_id'),
  filters = 'mgi_symbol',
  values = unique(kinases$Mouse_Symbol),
  mart = ensembl
)
setDT(biomart_map)
setnames(biomart_map, c('mgi_symbol', 'uniprotswissprot', 'entrezgene_id'), c('Mouse_Symbol', 'UniProt_Mouse_biomart', 'Entrez_Mouse'))

# Optionally, load a local mapping table (e.g., downloaded from UniProt)
local_map_file <- '../uniprot_mouse_idmapping_selected.tab'
if (file.exists(local_map_file)) {
  cat('Loading local UniProt mapping table...\n')
  local_map <- fread(local_map_file, header=TRUE)
  # Expect columns: Gene_Name, UniProtKB, etc.
  setnames(local_map, c('Gene_Name', 'UniProtKB'), c('Mouse_Symbol', 'UniProt_Mouse_local'))
  # Merge with kinases
  kinases <- merge(kinases, local_map[, .(Mouse_Symbol, UniProt_Mouse_local)], by='Mouse_Symbol', all.x=TRUE)
}

# Merge biomaRt results
kinases <- merge(kinases, biomart_map[, .(Mouse_Symbol, UniProt_Mouse_biomart)], by='Mouse_Symbol', all.x=TRUE)

# Save output
output_file <- '../251227_curated_kinases_biomart.csv'
fwrite(kinases, output_file)
cat(paste0('✓ Mouse UniProt mapping (biomaRt/local) complete. Output: ', output_file, '\n'))
