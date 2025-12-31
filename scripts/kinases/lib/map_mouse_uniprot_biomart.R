# Map mouse gene symbols to UniProt IDs using biomaRt and a local mapping table
# Author: C. Ryan Miller, MD, PhD (rmiller@uab.edu)

library(data.table)
library(biomaRt)


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
candidates <- list.files(input_dir, pattern='201006_composite_kinases_curated.*\\.csv$', full.names=TRUE, ignore.case=TRUE)
if (length(candidates) == 0) stop('Missing input snapshot: please place 201006_composite_kinases_curated__YYMMDD.csv in ', input_dir)
kinases_file <- sort(candidates, decreasing=TRUE)[1]
kinases <- fread(kinases_file, header=TRUE)

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
local_map_file <- file.path(input_dir, 'uniprot_mouse_idmapping_selected.tab')
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
output_file <- file.path(output_dir, '251227_curated_kinases_biomart.csv')
fwrite(kinases, output_file)
cat(paste0('✓ Mouse UniProt mapping (biomaRt/local) complete. Output: ', output_file, '\n'))
