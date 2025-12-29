#!/usr/bin/env Rscript
# Generic import/merge for validation sources (CSV, GMT, webpage-derived mappings)
# Scans `val_sources/` directory for files and merges any CSV mappings into the main sources table (baseline: kinases_human.csv).

suppressWarnings(suppressMessages({
  library(data.table)
  library(stringr)
  library(xml2)
  library(rvest)
  library(readxl)
}))

wd <- normalizePath(".", mustWork=TRUE)
setwd(wd)
cat("Working dir:", wd, "\n")

possible_dirs <- c("genesets/curated", "val_sources", file.path("kinases","val_sources"))
found <- possible_dirs[file.exists(possible_dirs)]
val_dir <- if (length(found) > 0) found[1] else NA
if (is.na(val_dir) || !dir.exists(val_dir)) {
  cat("No val_sources directory found; skipping val_sources merge\n")
  quit(status = 0)
}

# find files
files <- list.files(val_dir, full.names = TRUE)
if (length(files) == 0) {
  cat("No validation source files found in", val_dir, "\n")
  quit(status = 0)
}

# read main kinases file (default location inside kinases/)
kin_file <- "kinases_human.csv"
if (!file.exists(kin_file)) {
  if (file.exists(file.path("kinases", "kinases_human.csv"))) kin_file <- file.path("kinases", "kinases_human.csv")
}
if (!file.exists(kin_file)) stop("Sources baseline kinases_human.csv not found (expected kinases/kinases_human.csv)\nPlease run the builder to regenerate kinases_human.csv and retry.")
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
  } else if (ext == "xls" || ext == "xlsx") {
    # read Excel spreadsheet
    dt <- tryCatch(as.data.table(readxl::read_excel(f)), error=function(e) NULL)
    if (is.null(dt)) { cat("Failed to read excel", f, "\n"); next }
    colnames(dt) <- trimws(colnames(dt))
    key_by_ensembl <- intersect(tolower(colnames(dt)), c("ensembl","ensembl_gene_id","ensemblid"))
    key_by_name <- intersect(tolower(colnames(dt)), c("name","symbol","gene","gene_symbol","external_gene_name"))
    if (length(key_by_ensembl) > 0 && "ensembl_gene_id" %in% names(merged)) {
      kcol <- names(dt)[which(tolower(names(dt))==key_by_ensembl[1])]
      dt[, (kcol) := toupper(trimws(get(kcol)))]
      merged[, ensembl_gene_id := toupper(trimws(ensembl_gene_id))]
      merged <- merge(merged, dt, by.x = "ensembl_gene_id", by.y = kcol, all.x = TRUE, sort=FALSE, suffixes=c("",".vs"))
      cat("Merged Excel by Ensembl using", kcol, "from", f, "\n")
    } else if (length(key_by_name) > 0 && "external_gene_name" %in% names(merged)) {
      kcol <- names(dt)[which(tolower(names(dt))==key_by_name[1])]
      dt[, (kcol) := toupper(trimws(get(kcol)))]
      merged[, external_gene_name := toupper(trimws(external_gene_name))]
      merged <- merge(merged, dt, by.x = "external_gene_name", by.y = kcol, all.x = TRUE, sort=FALSE, suffixes=c("",".vs"))
      cat("Merged Excel by gene symbol using", kcol, "from", f, "\n")
    } else {
      cat("No suitable key found in", f, "— skipping merge\n")
    }
  } else if (ext == "gmt") {
    # parse GMT: <set_name> <description> <gene1> <gene2> ...
    lines <- readLines(f, warn = FALSE)
    gmt_dt <- data.table()
    for (ln in lines) {
      parts <- strsplit(ln, "\t")[[1]]
      if (length(parts) < 3) next
      set_name <- parts[1]
      genes <- parts[-c(1,2)]
      genes <- toupper(trimws(genes))
      tmp <- data.table(external_gene_name = genes, val_source = set_name, val_type = "GMT")
      gmt_dt <- rbind(gmt_dt, tmp)
    }
    if (nrow(gmt_dt) > 0) {
      merged[, external_gene_name := toupper(trimws(external_gene_name))]
      # collapse val_source per gene
      gmt_agg <- gmt_dt[, .(val_sources = paste(unique(val_source), collapse=",")), by = external_gene_name]
      merged <- merge(merged, gmt_agg, by = "external_gene_name", all.x = TRUE, sort = FALSE)
      cat("GMT merged:", f, "(added val_sources column)\n")
    }
  } else if (ext %in% c("html","htm")) {
    # Parse KinHub HTML specifically if filename indicates kinhub
    base <- tolower(basename(f))
    doc <- tryCatch(read_html(f), error=function(e) NULL)
    if (is.null(doc)) { cat("Failed to read HTML file:", f, "\n"); next }
    # If KinHub, build kinhub_mapping_raw.tsv expected by KinHub merger
    if (grepl("kinhub", base)) {
      cat("Parsing KinHub HTML and creating kinases/kinhub_mapping_raw.tsv\n")
      tables <- html_table(doc, fill=TRUE)
      # try to find a table with headers including 'Kinase' or 'HGNC'
      sel <- NULL
      for (i in seq_along(tables)) {
        h <- tolower(colnames(tables[[i]]))
        if (any(grepl("kinase|hgnc|symbol", h))) { sel <- i; break }
      }
      if (is.null(sel)) sel <- 1
      kt <- as.data.table(tables[[sel]])
      # Heuristics: find symbol column
      cn <- tolower(colnames(kt))
      sym_col <- names(kt)[which(grepl("symbol|kinase|name|gene", cn))[1]]
      group_col <- names(kt)[which(grepl("group|family|family name|group name", cn))[1]]
      # Build mapping: HGNC (symbol), Group, Family, SubFamily (best-effort)
      kt_up <- copy(kt)
      if (!is.null(sym_col)) kt_up[, (sym_col) := toupper(trimws(get(sym_col)))]
      out_km <- data.table(HGNC = if (!is.null(sym_col)) kt_up[[sym_col]] else NA_character_,
                           Group = if (!is.null(group_col)) as.character(kt_up[[group_col]]) else NA_character_,
                           Family = NA_character_,
                           SubFamily = NA_character_)
      fwrite(out_km, file = file.path("kinases","kinhub_mapping_raw.tsv"), sep = "\t", na = "")
      cat("Wrote kinases/kinhub_mapping_raw.tsv (rows:", nrow(out_km), ")\n")
      # Now call the existing KinHub merge script inside kinases/
      owd <- getwd()
      setwd("kinases")
      tryCatch({ source("fetch_kinhub_and_merge.R") }, error=function(e) cat("KinHub script error:", e$message, "\n"))
      setwd(owd)
    } else {
      # Generic HTML table parsing
      cat("Attempting to parse HTML table from", f, "\n")
      tabs <- html_table(doc, fill=TRUE)
      if (length(tabs) >= 1) {
        dt <- as.data.table(tabs[[1]])
        colnames(dt) <- trimws(colnames(dt))
        key_by_ensembl <- intersect(tolower(colnames(dt)), c("ensembl","ensembl_gene_id","ensemblid"))
        key_by_name <- intersect(tolower(colnames(dt)), c("name","symbol","gene","gene_symbol","external_gene_name"))
        if (length(key_by_ensembl) > 0 && "ensembl_gene_id" %in% names(merged)) {
          kcol <- names(dt)[which(tolower(names(dt))==key_by_ensembl[1])]
          dt[, (kcol) := toupper(trimws(get(kcol)))]
          merged[, ensembl_gene_id := toupper(trimws(ensembl_gene_id))]
          merged <- merge(merged, dt, by.x = "ensembl_gene_id", by.y = kcol, all.x = TRUE, sort=FALSE, suffixes=c("",".vs"))
          cat("Merged HTML table by Ensembl using", kcol, "from", f, "\n")
        } else if (length(key_by_name) > 0 && "external_gene_name" %in% names(merged)) {
          kcol <- names(dt)[which(tolower(names(dt))==key_by_name[1])]
          dt[, (kcol) := toupper(trimws(get(kcol)))]
          merged[, external_gene_name := toupper(trimws(external_gene_name))]
          merged <- merge(merged, dt, by.x = "external_gene_name", by.y = kcol, all.x = TRUE, sort=FALSE, suffixes=c("",".vs"))
          cat("Merged HTML table by gene symbol using", kcol, "from", f, "\n")
        } else {
          cat("No suitable key found in HTML table; skipping\n")
        }
      } else {
        cat("No HTML tables found in", f, "\n")
      }
    }
  } else {
    cat("Unsupported file type:", f, "— skipping\n")
  }
}

outf <- "kinases/kinases_human.with_val_sources.csv"
fwrite(merged, outf)
cat("Wrote merged output:", outf, "\n")
