#!/usr/bin/env Rscript
# Kinases-specific wrapper: query BioMart by GO:0004672 and augment with Manning/KEGG sources
suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
})

## Note: do NOT source the generic generator here; it's invoked as a separate Rscript process

## Curated lists to capture kinases: multiple GO terms + InterPro kinase domain IDs
go_terms <- c('GO:0004672', 'GO:0004713', 'GO:0004674')
interpro_ids <- c('IPR000719')

# find Manning file
manning_path <- list.files('genesets', pattern='manning.*(csv|xls|xlsx)$', recursive=TRUE, full.names=TRUE, ignore.case=TRUE)[1]
if (is.na(manning_path) || is.null(manning_path)) stop('Manning file not found under genesets')

outdir <- 'genesets/curated/kinases'
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
# temp dir for intermediate files
tmpdir <- file.path(outdir, 'temp')
dir.create(tmpdir, recursive=TRUE, showWarnings=FALSE)

# 1) Domain-based (InterPro) -> write to temp
cat('Generating domain-based kinases (InterPro) into temp\n')
cmd1 <- c(
  'filter_type=interpro',
  paste0('filter_value=', paste(interpro_ids, collapse=',')),
  'apply_manning=TRUE',
  paste0('manning_path=', shQuote(manning_path)),
  'outfile=kinases_human_domain.csv',
  paste0('outdir=', tmpdir)
)
system2('Rscript', args=c('scripts/utilities/generate_geneset_from_biomart.R', cmd1), stdout=TRUE, stderr=TRUE)

# 2) GO-based (functional) -> write to temp
cat('Generating GO-based kinases into temp\n')
cmd2 <- c(
  'filter_type=go',
  paste0('filter_value=', paste(go_terms, collapse=',')),
  'apply_manning=TRUE',
  paste0('manning_path=', shQuote(manning_path)),
  'outfile=kinases_human_go.csv',
  paste0('outdir=', tmpdir)
)
system2('Rscript', args=c('scripts/utilities/generate_geneset_from_biomart.R', cmd2), stdout=TRUE, stderr=TRUE)

# 3) Union of both (de-duplicated by Ensembl gene id) -> final
cat('Combining domain + GO temp results into union file\n')
library(data.table)
f1 <- file.path(tmpdir, 'kinases_human_domain.csv')
f2 <- file.path(tmpdir, 'kinases_human_go.csv')
dt1 <- if (file.exists(f1)) fread(f1) else data.table()
dt2 <- if (file.exists(f2)) fread(f2) else data.table()
union_dt <- unique(rbindlist(list(dt1, dt2), use.names=TRUE, fill=TRUE), by='ensembl_gene_id')
union_out <- file.path(outdir, 'kinases_human_union.csv')
fwrite(union_dt, union_out)

# md5s: temp checksums for intermediates, final checksum for union
cs_temp <- file.path(tmpdir,'checksums'); dir.create(cs_temp, recursive=TRUE, showWarnings=FALSE)
cs_out <- file.path(outdir,'checksums'); dir.create(cs_out, recursive=TRUE, showWarnings=FALSE)
if (file.exists(f1)) write(sprintf('%s  %s', tools::md5sum(f1), basename(f1)), file.path(cs_temp,'kinases_human_domain.csv.md5'))
if (file.exists(f2)) write(sprintf('%s  %s', tools::md5sum(f2), basename(f2)), file.path(cs_temp,'kinases_human_go.csv.md5'))
if (file.exists(union_out)) write(sprintf('%s  %s', tools::md5sum(union_out), basename(union_out)), file.path(cs_out,'kinases_human_union.csv.md5'))

# Commit outputs (include temp intermediates if desired)
system(sprintf('git add -f %s %s %s %s %s', shQuote(f1), shQuote(f2), shQuote(union_out), shQuote(cs_temp), shQuote(cs_out)), ignore.stdout=TRUE)
system('git commit -m "Generate kinases_human domain/go (temp) and union outputs from BioMart; add checksums" || true')
