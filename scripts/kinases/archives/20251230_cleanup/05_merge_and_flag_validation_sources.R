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

# --- Provenance annotation setup (must come after repo_root is defined) ---
# Load KinHub and HGNC symbol sets for provenance tracking
kinhub_csv <- file.path(repo_root, 'genesets','curated','kinases','val_sources','kinhub_mapping.csv')
kinhub_syms <- character(0)
if (file.exists(kinhub_csv)) {
  kinhub_dt <- tryCatch(fread(kinhub_csv), error=function(e) NULL)
  if (!is.null(kinhub_dt) && 'HGNC' %in% names(kinhub_dt)) {
    kinhub_syms <- unique(na.omit(trimws(kinhub_dt[['HGNC']])))
  }
}
hgnc_csv <- file.path(repo_root, 'genesets','curated','kinases','val_sources','hgnc_kinase_groups.csv')
hgnc_syms <- character(0)
if (file.exists(hgnc_csv)) {
  hgnc_dt <- tryCatch(fread(hgnc_csv), error=function(e) NULL)
  if (!is.null(hgnc_dt) && 'symbol' %in% names(hgnc_dt)) {
    hgnc_syms <- unique(na.omit(trimws(hgnc_dt[['symbol']])))
  }
}

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
f1 <- file.path(tmpdir, 'kinases_human_domain_interpro.csv')
f2 <- file.path(tmpdir, 'kinases_human_go_filtered.csv')
dt1 <- if (file.exists(f1)) fread(f1) else data.table()
dt2 <- if (file.exists(f2)) fread(f2) else data.table()
union_dt <- unique(rbindlist(list(dt1, dt2), use.names=TRUE, fill=TRUE), by='ensembl_gene_id')

# Clean description: remove [ ... ]
if ("description" %in% names(union_dt)) {
  union_dt[, description := gsub("\\s*\\[[^]]*]", "", description)]
}

# Remove duplicate Group/Family/Subfamily columns (keep only canonical, drop Manning_*)
manning_cols_to_remove <- intersect(c("Manning_Group","Manning_Family","Manning_Subfamily"), names(union_dt))
if (length(manning_cols_to_remove) > 0) union_dt[, (manning_cols_to_remove) := NULL]


# Add metabolic flag from Mammalian Metabolic Enzyme Database with mouse-to-human ortholog mapping
metabolic_csv <- file.path(repo_root, "genesets","curated","kinases","val_sources","Mammalian_Metabolic_Final.csv")
ortholog_csv <- file.path(repo_root, "archives", "251227_repo_cleanup2", "human_mouse_orthologs.csv")
metabolic_verified_syms <- character(0)
if (file.exists(metabolic_csv)) {
  metabolic_dt <- tryCatch(fread(metabolic_csv), error=function(e) NULL)
  if (!is.null(metabolic_dt) && "Gene Symbol" %in% names(metabolic_dt)) {
    mouse_syms <- unique(na.omit(trimws(metabolic_dt[["Gene Symbol"]])))
    mouse_syms_upper <- toupper(mouse_syms)
    # Map mouse symbols to human using ortholog table
    human_syms <- character(0)
    if (file.exists(ortholog_csv)) {
      ortho_dt <- tryCatch(fread(ortholog_csv), error=function(e) NULL)
      if (!is.null(ortho_dt) && "Mouse_Symbol" %in% names(ortho_dt) && "HGNC_Symbol" %in% names(ortho_dt)) {
        ortho_dt[, Mouse_Symbol_upper := toupper(trimws(Mouse_Symbol))]
        ortho_dt[, HGNC_Symbol_upper := toupper(trimws(HGNC_Symbol))]
        # Map mouse symbols to human symbols
        matched_human <- ortho_dt[Mouse_Symbol_upper %in% mouse_syms_upper & HGNC_Symbol_upper != "", HGNC_Symbol_upper]
        human_syms <- unique(na.omit(matched_human))
      }
    }
    # Combine mouse symbols (for direct match) and mapped human symbols
    metabolic_verified_syms <- unique(c(mouse_syms_upper, human_syms))
    cat(sprintf("Loaded %d metabolic enzyme gene symbols from Mammalian Metabolic Enzyme Database (mouse), mapped to %d human symbols.\n", length(mouse_syms_upper), length(human_syms)))
  }
}
# Flexible matching: check external_gene_name, hgnc_symbol, symbol columns
symbol_cols <- intersect(c("external_gene_name","hgnc_symbol","symbol","gene_symbol"), names(union_dt))
if (!"Metabolic" %in% names(union_dt)) union_dt[, Metabolic := ""]
if (!"Provenance_KinHub" %in% names(union_dt)) union_dt[, Provenance_KinHub := FALSE]
if (!"Provenance_HGNC" %in% names(union_dt)) union_dt[, Provenance_HGNC := FALSE]

# Add Provenance_Manning column (default FALSE)
if (!"Provenance_Metabolic" %in% names(union_dt)) union_dt[, Provenance_Metabolic := FALSE]
if (!"Provenance_Manning" %in% names(union_dt)) union_dt[, Provenance_Manning := FALSE]

# Load Manning 'Name' column for provenance mapping
manning_names <- character(0)
if (file.exists(manning_path)) {
  manning_dt <- tryCatch(fread(manning_path), error=function(e) NULL)
  if (!is.null(manning_dt) && "Name" %in% names(manning_dt)) {
    manning_names <- unique(na.omit(trimws(manning_dt[["Name"]])))
    manning_names <- toupper(manning_names)
  }
}

# For each row, set provenance columns based on symbol matches
for (i in seq_len(nrow(union_dt))) {
  syms <- unique(na.omit(toupper(trimws(as.character(union_dt[i, ..symbol_cols])))))
  # KinHub provenance
  if (length(intersect(syms, toupper(kinhub_syms))) > 0) union_dt[i, Provenance_KinHub := TRUE]
  # HGNC provenance
  if (length(intersect(syms, toupper(hgnc_syms))) > 0) union_dt[i, Provenance_HGNC := TRUE]
  # Metabolic provenance and flag
  if (length(intersect(syms, metabolic_verified_syms)) > 0) {
    union_dt[i, Metabolic := "Verified"]
    union_dt[i, Provenance_Metabolic := TRUE]
  }
  # Manning provenance (Name column, mapped to official gene names)
  if (length(intersect(syms, manning_names)) > 0) {
    union_dt[i, Provenance_Manning := TRUE]
  }
}

# Write 04_kinases_human_union__DATE.csv as before
date_tag <- format(Sys.time(), "%y%m%d")
int1 <- file.path(tmpdir, paste0('04_', basename(f1)))
int2 <- file.path(tmpdir, paste0('04_', basename(f2)))
if (nrow(dt1)>0) fwrite(dt1, int1)
if (nrow(dt2)>0) fwrite(dt2, int2)
if (file.exists(int1)) write(sprintf('%s  %s', tools::md5sum(int1), basename(int1)), file.path(dirname(int1), paste0(basename(int1), '.md5')))
if (file.exists(int2)) write(sprintf('%s  %s', tools::md5sum(int2), basename(int2)), file.path(dirname(int2), paste0(basename(int2), '.md5')))
union_out <- file.path(outdir, paste0('04_kinases_human_union__', date_tag, '.csv'))
fwrite(union_dt, union_out)
if (file.exists(f1)) write(sprintf('%s  %s', tools::md5sum(f1), basename(f1)), file.path(dirname(f1), paste0(basename(f1), '.md5')))
if (file.exists(f2)) write(sprintf('%s  %s', tools::md5sum(f2), basename(f2)), file.path(dirname(f2), paste0(basename(f2), '.md5')))
if (file.exists(union_out)) write(sprintf('%s  %s', tools::md5sum(union_out), basename(union_out)), file.path(dirname(union_out), paste0(basename(union_out), '.md5')))

# NEW: Write 05_kinases_human.with_val_sources__DATE.csv as a fully annotated version (all columns from union_dt)
out_05 <- file.path(outdir, paste0('05_kinases_human.with_val_sources__', date_tag, '.csv'))
fwrite(union_dt, out_05)
if (file.exists(out_05)) write(sprintf('%s  %s', tools::md5sum(out_05), basename(out_05)), file.path(dirname(out_05), paste0(basename(out_05), '.md5')))

# Commit outputs (include temp intermediates if desired)
system(sprintf('git add -f %s %s %s %s', shQuote(f1), shQuote(f2), shQuote(union_out), shQuote(out_05)), ignore.stdout=TRUE)
system('git commit -m "Generate kinases: domain/go intermediates, union, and fully annotated 05 output; add checksums" || true')
