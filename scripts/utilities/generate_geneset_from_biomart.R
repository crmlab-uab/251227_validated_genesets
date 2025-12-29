#!/usr/bin/env Rscript
# Generic generator: query BioMart by filter (go/interpro/custom), apply optional augmentation, write CSV+md5
suppressPackageStartupMessages({
  library(data.table)
  library(biomaRt)
  library(readxl)
})

args <- commandArgs(trailingOnly = TRUE)
parse_args <- function(a){
  out <- list()
  for (x in a) {
    kv <- strsplit(x, "=")[[1]]
    if (length(kv)==2) out[[kv[1]]] <- kv[2]
  }
  return(out)
}
opts <- parse_args(args)

filter_type <- ifelse(!is.null(opts$filter_type), opts$filter_type, 'go')
filter_value <- ifelse(!is.null(opts$filter_value), opts$filter_value, stop('filter_value required'))
dataset <- ifelse(!is.null(opts$dataset), opts$dataset, 'hsapiens_gene_ensembl')
outdir <- ifelse(!is.null(opts$outdir), opts$outdir, 'genesets/curated/custom')
outfile <- file.path(outdir, ifelse(!is.null(opts$outfile), opts$outfile, paste0('geneset_', filter_type, '_', gsub('[^A-Za-z0-9]','_',filter_value), '.csv')))
apply_manning_flag <- ifelse(!is.null(opts$apply_manning), opts$apply_manning == 'TRUE', FALSE)
manning_path <- ifelse(!is.null(opts$manning_path), opts$manning_path, NULL)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

source('scripts/utilities/functions_genesets.R', local=FALSE)

cat('Connecting to Ensembl BioMart (dataset=', dataset, ')\n')
mart <- tryCatch({ useEnsembl(biomart='genes', dataset=dataset) }, error=function(e){ message('useEnsembl failed, using useMart fallback'); useMart('ensembl', dataset) })

attributes <- c('ensembl_gene_id','hgnc_symbol','external_gene_name','uniprotswissprot','entrezgene_id','description')

cat('Query type:', filter_type, 'value:', filter_value, '\n')
vals <- strsplit(filter_value, ',')[[1]]
vals <- trimws(vals)
results <- list()
if (tolower(filter_type) == 'go') {
  for (v in vals) {
    cat(' Querying GO term:', v, '\n')
    tmp <- getBM(attributes=attributes, filters='go', values=v, mart=mart)
    results[[v]] <- as.data.table(tmp)
  }
} else if (tolower(filter_type) == 'interpro') {
  for (v in vals) {
    cat(' Querying InterPro:', v, '\n')
    tmp <- getBM(attributes=attributes, filters='interpro', values=v, mart=mart)
    results[[v]] <- as.data.table(tmp)
  }
} else if (tolower(filter_type) == 'custom') {
  # expected format: filter_name:filter_value (only single entry supported)
  parts <- strsplit(filter_value, ':')[[1]]
  if (length(parts) != 2) stop('custom filter_value must be filter:value')
  tmp <- getBM(attributes=attributes, filters=parts[1], values=parts[2], mart=mart)
  results[[parts[2]]] <- as.data.table(tmp)
} else stop('Unsupported filter_type')

# combine and deduplicate by ensembl_gene_id
if (length(results) == 0) {
  res <- data.table()
} else {
  res <- unique(rbindlist(results, use.names=TRUE, fill=TRUE), by=c('ensembl_gene_id','hgnc_symbol','external_gene_name'))
}

dt <- as.data.table(res)
setnames(dt, old=c('hgnc_symbol','external_gene_name','uniprotswissprot','entrezgene_id'), new=c('HGNC_gene_name','external_gene_name','UniProtID_Human','Entrez_ID'), skip_absent=TRUE)
if (!'external_gene_name' %in% names(dt) && 'HGNC_gene_name' %in% names(dt)) dt[, external_gene_name := HGNC_gene_name]

if (apply_manning_flag) {
  if (is.null(manning_path)) stop('apply_manning requested but no manning_path provided')
  cat('Reading Manning file:', manning_path, '\n')
  manning_dt <- tryCatch(fread(manning_path), error=function(e) as.data.table(readxl::read_excel(manning_path)))
  dt <- apply_manning(dt, manning_dt)
}

# optional dedupe
if (exists('dedupe_genes')) dt <- tryCatch(dedupe_genes(dt), error=function(e){ message('dedupe failed: ', e$message); dt })

dir.create(dirname(outfile), recursive=TRUE, showWarnings=FALSE)
fwrite(dt, outfile)
md5 <- tools::md5sum(outfile)
md5file <- paste0(outfile, '.md5')
write(sprintf('%s  %s', md5, basename(outfile)), md5file)
cat('Wrote', outfile, 'rows=', nrow(dt), 'md5=', md5, '\n')
