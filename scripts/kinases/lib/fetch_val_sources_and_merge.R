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

# Prefer kinases-specific validation and inputs directories; fall back to generic locations
possible_dirs <- c(
  file.path("genesets","curated","kinases","val_sources"),
  file.path("genesets","curated","kinases","inputs"),
  file.path("genesets","curated","kinases"),
  file.path("genesets","curated"),
  "val_sources",
  file.path("kinases","val_sources")
)
found <- possible_dirs[file.exists(possible_dirs)]
val_dir <- if (length(found) > 0) found[1] else NA
if (is.na(val_dir) || !dir.exists(val_dir)) {
  cat("No val_sources or kinases inputs directory found; skipping val_sources merge\n")
  quit(status = 0)
}

# find files (ignore subdirectories)
files_all <- list.files(val_dir, full.names = TRUE)
files <- files_all[!file.info(files_all)$isdir]
if (length(files) == 0) {
  cat("No validation source files found in", val_dir, "\n")
  quit(status = 0)
}

# locate main kinases baseline (try multiple canonical locations)
kin <- NULL
# Prefer the canonical inputs location created by the builder, then repo-root fallbacks, then any outputs snapshots
canonical_input <- file.path("genesets","curated","kinases","inputs","kinases_human.csv")
candidate_baselines <- c(canonical_input, "kinases_human.csv", file.path("kinases","kinases_human.csv"))
out_dir <- file.path("genesets","curated","kinases","outputs")
out_candidates <- character(0)
if (dir.exists(out_dir)) {
  out_candidates <- list.files(out_dir, pattern = "^kinases_human.*\\.csv$", full.names = TRUE)
  if (length(out_candidates) > 0) {
    fi <- file.info(out_candidates)
    out_candidates <- out_candidates[order(fi$mtime, decreasing = TRUE)]
  }
}
candidate_baselines <- c(candidate_baselines, out_candidates)
for (kb in candidate_baselines) {
  if (!is.null(kb) && nzchar(kb) && file.exists(kb)) {
    kin_try <- tryCatch(fread(kb, na.strings = c("", "NA")), error = function(e) NULL)
    if (!is.null(kin_try)) { kin <- kin_try; cat("Using baseline:", kb, "\n"); break }
  }
}
if (is.null(kin)) stop("Sources baseline kinases_human.csv not found in expected locations. Please run the builder to regenerate kinases_human.csv and retry.")

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
      # write kinHub raw mapping into canonical val_sources dir
      val_dir <- file.path("genesets","curated","kinases","val_sources")
      dir.create(val_dir, recursive=TRUE, showWarnings=FALSE)
      km_path <- file.path(val_dir, "kinhub_mapping_raw.tsv")
      fwrite(out_km, file = km_path, sep = "\t", na = "")
      cat("Wrote", km_path, "(rows:", nrow(out_km), ")\n")
      # Now invoke the canonical KinHub merge runner (02) from repo root; it will read from val_sources and write to outputs
      merge_runner <- file.path("scripts","kinases","bin","02_fetch_validation_sources.R")
      if (file.exists(merge_runner)) {
        tryCatch({ source(merge_runner) }, error=function(e) cat("KinHub merge runner error:", e$message, "\n"))
      } else {
        cat("No KinHub merge runner found at", merge_runner, "— skipping KinHub merge step\n")
      }
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

out_dir <- file.path("genesets","curated","kinases","outputs")
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
stamp <- format(Sys.time(), "%y%m%d")
outf <- file.path(out_dir, paste0("05_kinases_human.with_val_sources__", stamp, ".csv"))
fwrite(merged, outf)
# write md5
md5f <- paste0(outf, ".md5")
write(sprintf("%s  %s", tools::md5sum(outf), basename(outf)), md5f)
cat("Wrote merged output:", outf, "and checksum", md5f, "\n")
