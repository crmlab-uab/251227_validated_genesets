# Fetch UniProt mouse gene name to UniProtKB mapping for kinases only
# Author: C. Ryan Miller, MD, PhD (rmiller@uab.edu)

library(httr)
library(data.table)


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
kinase_file <- sort(candidates, decreasing=TRUE)[1]
kinases <- fread(kinase_file, header=TRUE)
gene_symbols <- unique(na.omit(kinases$Mouse_Symbol))

# Prepare output file
output_file <- file.path(output_dir, 'uniprot_mouse_kinase_idmapping.tab')

cat('Fetching UniProt mouse gene name to UniProtKB mapping for kinases only...\n')

# Query UniProt for each gene symbol
results <- list()
for (symbol in gene_symbols) {
  query <- paste0('gene_exact:', symbol, ' AND organism_id:10090')
  url <- paste0('https://rest.uniprot.org/uniprotkb/search?query=', URLencode(query), '&fields=accession,gene_names&format=tsv&size=500')
  resp <- GET(url)
  if (status_code(resp) != 200) next
  dat <- fread(text=content(resp, as='text', encoding='UTF-8'), sep='\t', header=TRUE)
  if (nrow(dat) > 0) {
    dat[, Mouse_Symbol := symbol]
    results[[length(results)+1]] <- dat
  }
  cat('Fetched', symbol, '\n')
}

if (length(results) > 0) {
  mapping <- rbindlist(results, fill=TRUE)
  fwrite(mapping, output_file, sep='\t')
  cat('✓ Kinase mapping file saved to', output_file, '\n')
} else {
  cat('No kinase mapping data retrieved.\n')
}
