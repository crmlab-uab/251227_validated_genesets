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

# find Manning file under curated inputs
## compute repo root from script location so paths are absolute to repo root
cmdArgs <- commandArgs(trailingOnly = FALSE)
fileArgIdx <- grep("--file=", cmdArgs)
fileArg <- if (length(fileArgIdx) > 0) sub("--file=", "", cmdArgs[fileArgIdx[1]]) else NA_character_
script_dir <- if (!is.na(fileArg) && nzchar(fileArg)) dirname(normalizePath(fileArg)) else normalizePath(".")
repo_root <- normalizePath(file.path(script_dir, "..", "..", ".."))

inputs_dir <- file.path(repo_root, 'genesets','curated','kinases','inputs')
if (!dir.exists(inputs_dir)) stop('Inputs directory not found: ', inputs_dir)
manning_path <- list.files(inputs_dir, pattern='manning.*(csv|xls|xlsx)$', recursive=FALSE, full.names=TRUE, ignore.case=TRUE)[1]
if (is.na(manning_path) || is.null(manning_path)) stop('Manning file not found under ', inputs_dir)

# canonical outputs dir
outdir <- file.path(repo_root, 'genesets','curated','kinases','outputs')
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
# temp dir for intermediate files (project-wide)
tmpdir <- file.path(repo_root, 'output','temp')
dir.create(tmpdir, recursive=TRUE, showWarnings=FALSE)

# Create session log directory and redirect output to a session log
sess_dir <- file.path(repo_root, 'sessions', format(Sys.time(), '%y%m%d_%H%M%S'))
dir.create(sess_dir, recursive=TRUE, showWarnings=FALSE)

logf <- file.path(sess_dir, '04_generate_from_biomart.log')
cat(sprintf('Writing log to %s\n', logf))
log_con <- file(logf, open="wt")
sink(log_con, append=TRUE, split=TRUE)
sink(log_con, type='message', append=TRUE)
on.exit({
  sink(type='message')
  sink()
  close(log_con)
}, add=TRUE)

# 1) Domain-based (InterPro) -> write to temp
cat('Generating domain-based kinases (InterPro) into output/temp\n')
  cmd1 <- c(
  'filter_type=interpro',
  paste0('filter_value=', paste(interpro_ids, collapse=',')),
  'apply_manning=TRUE',
  paste0('manning_path=', shQuote(manning_path)),
    'outfile=kinases_human_domain_interpro.csv',
    paste0('outdir=', tmpdir)
)
system2('Rscript', args=c('scripts/utilities/generate_geneset_from_biomart.R', cmd1), stdout=TRUE, stderr=TRUE)

# 2) GO-based (functional) -> write to temp
cat('Generating GO-based kinases into output/temp\n')
  cmd2 <- c(
  'filter_type=go',
  paste0('filter_value=', paste(go_terms, collapse=',')),
  'apply_manning=TRUE',
  paste0('manning_path=', shQuote(manning_path)),
    'outfile=kinases_human_go_filtered.csv',
    paste0('outdir=', tmpdir)
)
system2('Rscript', args=c('scripts/utilities/generate_geneset_from_biomart.R', cmd2), stdout=TRUE, stderr=TRUE)


# 3) Union of both (de-duplicated by Ensembl gene id) -> final
cat('Combining domain + GO temp results into union file\n')
library(data.table)
# intermediates: add numbered prefix 04_
f1 <- file.path(tmpdir, 'kinases_human_domain_interpro.csv')
f2 <- file.path(tmpdir, 'kinases_human_go_filtered.csv')
dt1 <- if (file.exists(f1)) fread(f1) else data.table()
dt2 <- if (file.exists(f2)) fread(f2) else data.table()
union_dt <- unique(rbindlist(list(dt1, dt2), use.names=TRUE, fill=TRUE), by='ensembl_gene_id')

# Clean description: remove [ ... ]
if ("description" %in% names(union_dt)) {
  union_dt[, description := gsub("\\s*\\[[^]]*]", "", description)]
}

# Add metabolic flag from Mammalian Metabolic Enzyme Database
metabolic_csv <- file.path(repo_root, "genesets","curated","kinases","val_sources","Mammalian_Metabolic_Final.csv")
metabolic_verified_syms <- character(0)
if (file.exists(metabolic_csv)) {
  metabolic_dt <- tryCatch(fread(metabolic_csv), error=function(e) NULL)
  if (!is.null(metabolic_dt) && "Gene Symbol" %in% names(metabolic_dt)) {
    metabolic_verified_syms <- toupper(trimws(metabolic_dt[["Gene Symbol"]]))
    cat(sprintf("Loaded %d metabolic enzyme gene symbols from Mammalian Metabolic Enzyme Database.\n", length(metabolic_verified_syms)))
  }
}
if ("external_gene_name" %in% names(union_dt)) {
  union_dt[, Metabolic := ifelse(toupper(external_gene_name) %in% metabolic_verified_syms, "Verified", "")]
}

date_tag <- format(Sys.time(), "%y%m%d")
# write numbered intermediates
int1 <- file.path(tmpdir, paste0('04_', basename(f1)))
int2 <- file.path(tmpdir, paste0('04_', basename(f2)))
if (nrow(dt1)>0) fwrite(dt1, int1)
if (nrow(dt2)>0) fwrite(dt2, int2)
if (file.exists(int1)) write(sprintf('%s  %s', tools::md5sum(int1), basename(int1)), file.path(dirname(int1), paste0(basename(int1), '.md5')))
if (file.exists(int2)) write(sprintf('%s  %s', tools::md5sum(int2), basename(int2)), file.path(dirname(int2), paste0(basename(int2), '.md5')))
# union output with 04_ prefix
union_out <- file.path(outdir, paste0('04_kinases_human_union__', date_tag, '.csv'))
fwrite(union_dt, union_out)

# md5s: temp checksums for intermediates, final checksum for union
if (file.exists(f1)) write(sprintf('%s  %s', tools::md5sum(f1), basename(f1)), file.path(dirname(f1), paste0(basename(f1), '.md5')))
if (file.exists(f2)) write(sprintf('%s  %s', tools::md5sum(f2), basename(f2)), file.path(dirname(f2), paste0(basename(f2), '.md5')))
if (file.exists(union_out)) write(sprintf('%s  %s', tools::md5sum(union_out), basename(union_out)), file.path(dirname(union_out), paste0(basename(union_out), '.md5')))

# Commit outputs (include temp intermediates if desired)
system(sprintf('git add -f %s %s %s', shQuote(f1), shQuote(f2), shQuote(union_out)), ignore.stdout=TRUE)
system('git commit -m "Generate kinases: domain/go intermediates and union output; add checksums" || true')
