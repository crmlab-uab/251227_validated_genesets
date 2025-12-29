#!/usr/bin/env Rscript
# Kinases-specific wrapper: query BioMart by GO:0004672 and augment with Manning/KEGG sources
suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
})

source('scripts/utilities/generate_geneset_from_biomart.R', local=FALSE)

## Curated lists to capture kinases: multiple GO terms + InterPro kinase domain IDs
go_terms <- c('GO:0004672', 'GO:0004713', 'GO:0004674')
interpro_ids <- c('IPR000719')

# find Manning file
manning_path <- list.files('genesets', pattern='manning.*(csv|xls|xlsx)$', recursive=TRUE, full.names=TRUE, ignore.case=TRUE)[1]
if (is.na(manning_path) || is.null(manning_path)) stop('Manning file not found under genesets')

outdir <- 'genesets/curated/kinases'
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

# 1) Domain-based (InterPro)
cat('Generating domain-based kinases (InterPro)\n')
cmd1 <- c(
  'filter_type=interpro',
  paste0('filter_value=', paste(interpro_ids, collapse=',')),
  'apply_manning=TRUE',
  paste0('manning_path=', shQuote(manning_path)),
  'outfile=kinases_human_domain.csv',
  paste0('outdir=', outdir)
)
system2('Rscript', args=c('scripts/utilities/generate_geneset_from_biomart.R', cmd1), stdout=TRUE, stderr=TRUE)

# 2) GO-based (functional)
cat('Generating GO-based kinases\n')
cmd2 <- c(
  'filter_type=go',
  paste0('filter_value=', paste(go_terms, collapse=',')),
  'apply_manning=TRUE',
  paste0('manning_path=', shQuote(manning_path)),
  'outfile=kinases_human_go.csv',
  paste0('outdir=', outdir)
)
system2('Rscript', args=c('scripts/utilities/generate_geneset_from_biomart.R', cmd2), stdout=TRUE, stderr=TRUE)

# 3) Union of both (de-duplicated by Ensembl gene id)
cat('Combining domain + GO results into union file\n')
library(data.table)
f1 <- file.path(outdir, 'kinases_human_domain.csv')
f2 <- file.path(outdir, 'kinases_human_go.csv')
dt1 <- if (file.exists(f1)) fread(f1) else data.table()
dt2 <- if (file.exists(f2)) fread(f2) else data.table()
union_dt <- unique(rbindlist(list(dt1, dt2), use.names=TRUE, fill=TRUE), by='ensembl_gene_id')
union_out <- file.path(outdir, 'kinases_human_union.csv')
fwrite(union_dt, union_out)
# md5s
write(sprintf('%s  %s', tools::md5sum(f1), basename(f1)), file.path(outdir,'checksums','kinases_human_domain.csv.md5'))
write(sprintf('%s  %s', tools::md5sum(f2), basename(f2)), file.path(outdir,'checksums','kinases_human_go.csv.md5'))
write(sprintf('%s  %s', tools::md5sum(union_out), basename(union_out)), file.path(outdir,'checksums','kinases_human_union.csv.md5'))

# Commit outputs
system(sprintf('git add -f %s %s %s %s', shQuote(f1), shQuote(f2), shQuote(union_out), shQuote(file.path(outdir,'checksums'))), ignore.stdout=TRUE)
system('git commit -m "Generate kinases_human domain/go/union outputs from BioMart, apply Manning, add checksums" || true')
