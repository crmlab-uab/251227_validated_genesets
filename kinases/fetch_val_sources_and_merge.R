#!/usr/bin/env Rscript
# Generic import/merge for validation sources (CSV, GMT, webpage-derived mappings)
# Scans `val_sources/` directory for files and merges any CSV mappings into the main kinases table.

suppressWarnings(suppressMessages({
  library(data.table)
  library(stringr)
}))

wd <- normalizePath(".", mustWork=TRUE)
setwd(wd)
cat("Working dir:", wd, "\n")

val_dir <- "val_sources"
if (!dir.exists(val_dir)) {
  cat("No val_sources directory found; skipping val_sources merge\n")
  quit(status = 0)
}

# find files
files <- list.files(val_dir, full.names = TRUE)
if (length(files) == 0) {
  cat("No validation source files found in", val_dir, "\n")
  quit(status = 0)
}

# read main kinases file (default)
kin_file <- "kinases_human.csv"
if (!file.exists(kin_file)) stop("kinases_human.csv not found in working dir")
kin <- fread(kin_file, na.strings=c("","NA"))

merged <- copy(kin)

for (f in files) {
  ext <- tolower(tools::file_ext(f))
  cat("Processing validation source:", f, "(ext=", ext, ")\n")
  if (ext == "csv" || ext == "tsv" || ext == "txt") {
    dt <- tryCatch(fread(f, na.strings=c("","NA")), error=function(e) NULL)
    if (is.null(dt)) { cat("Failed to read", f, "\n"); next }
    # try to detect common key columns
    colnames(dt) <- trimws(colnames(dt))
    key_by_ensembl <- intersect(tolower(colnames(dt)), c("ensembl","ensembl_gene_id","ensemblid"))
    key_by_name <- intersect(tolower(colnames(dt)), c("name","symbol","gene","gene_symbol","external_gene_name"))
    if (length(key_by_ensembl) > 0 && "ensembl_gene_id" %in% names(merged)) {
      kcol <- names(dt)[which(tolower(names(dt))==key_by_ensembl[1])]
      dt[, (kcol) := toupper(trimws(get(kcol)))]
      merged[, ensembl_gene_id := toupper(trimws(ensembl_gene_id))]
      merged <- merge(merged, dt, by.x = "ensembl_gene_id", by.y = kcol, all.x = TRUE, sort=FALSE, suffixes=c("",".vs"))
      cat("Merged by Ensembl using", kcol, "from", f, "\n")
    } else if (length(key_by_name) > 0 && "external_gene_name" %in% names(merged)) {
      kcol <- names(dt)[which(tolower(names(dt))==key_by_name[1])]
      dt[, (kcol) := toupper(trimws(get(kcol)))]
      merged[, external_gene_name := toupper(trimws(external_gene_name))]
      merged <- merge(merged, dt, by.x = "external_gene_name", by.y = kcol, all.x = TRUE, sort=FALSE, suffixes=c("",".vs"))
      cat("Merged by gene symbol using", kcol, "from", f, "\n")
    } else {
      cat("No suitable key found in", f, "— skipping merge\n")
    }
  } else if (ext == "gmt") {
    cat("GMT import detected for", f, ": GMT import support will be added; skipping for now\n")
    # TODO: parse GMT and create mapping table for merging
  } else if (ext %in% c("html","htm")) {
    cat("HTML page detected:", f, "— if this is a KinHub page use the dedicated script to parse it. Skipping here.\n")
  } else {
    cat("Unsupported file type:", f, "— skipping\n")
  }
}

outf <- "kinases/kinases_human.with_val_sources.csv"
fwrite(merged, outf)
cat("Wrote merged output:", outf, "\n")
