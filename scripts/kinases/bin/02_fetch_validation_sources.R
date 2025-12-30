#!/usr/bin/env Rscript
# Unified kinase validation source fetch/merge script
# Usage: Rscript fetch_validation_sources.R --source=kinhub|hgnc [--config=genesets_config.yaml]

suppressWarnings(suppressMessages({
  library(data.table)
  library(httr)
  library(jsonlite)
  if (requireNamespace("yaml", quietly=TRUE)) library(yaml)
}))

# --- Argument parsing ---
args <- commandArgs(trailingOnly=TRUE)
get_arg <- function(key, default=NULL) {
  val <- grep(paste0("^--",key,"="), args, value=TRUE)
  if (length(val)) sub(paste0("^--",key,"="), "", val[1]) else default
}
source_type <- get_arg("source", "kinhub")
cfg_file <- get_arg("config", NA)

# --- Config and repo root ---
cmdArgs <- commandArgs(trailingOnly = FALSE)
fileArgIdx <- grep("--file=", cmdArgs)
fileArg <- if (length(fileArgIdx) > 0) sub("--file=", "", cmdArgs[fileArgIdx[1]]) else NA_character_
script_dir <- if (!is.na(fileArg) && nzchar(fileArg)) dirname(normalizePath(fileArg)) else normalizePath(".")
repo_root <- normalizePath(file.path(script_dir, "..", "..", ".."))

# Load config if available
cfg <- list()
if (!is.na(cfg_file) && file.exists(cfg_file)) {
  cfg <- yaml::read_yaml(cfg_file)
} else {
  # try to find config upward
  search_up <- function(start, filename = "genesets_config.yaml", max_up = 6) {
    cur <- normalizePath(start)
    for (i in 0:max_up) {
      cand <- file.path(cur, filename)
      if (file.exists(cand)) return(cand)
      parent <- dirname(cur)
      if (parent == cur) break
      cur <- parent
    }
    return(NA_character_)
  }
  cfgf <- search_up(script_dir, "genesets_config.yaml", 8)
  if (!is.na(cfgf) && nzchar(cfgf)) cfg <- yaml::read_yaml(cfgf)
}
cfg_get <- function(path, default) {
  parts <- strsplit(path, "\\.")[[1]]
  cur <- cfg
  for (p in parts) { if (is.null(cur[[p]])) return(default); cur <- cur[[p]] }
  cur
}

# --- Input file selection ---
input_dir <- file.path(repo_root, "genesets","curated","kinases","inputs")
outdir <- file.path(repo_root, "genesets","curated","kinases","val_sources")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

canonical <- file.path(input_dir, "kinases_human.csv")
if (file.exists(canonical)) {
  kin_file <- canonical
} else {
  candidates <- list.files(input_dir, pattern = "^(kinases_human|human_kinome)__.*\\.csv$", full.names = TRUE)
  if (length(candidates) == 0) stop("No kinases_human*.csv or human_kinome*.csv file found in ", input_dir)
  candidate_info <- file.info(candidates)
  kin_file <- rownames(candidate_info)[which.max(candidate_info$mtime)]
}
cat("Using baseline kinome input:", kin_file, "\n", file=stderr())
kin <- tryCatch(fread(kin_file), error = function(e) stop("Failed to read baseline kinome CSV: ", e$message))

# --- Symbol column detection ---
symbol_cols <- c("external_gene_name","hgnc_symbol","symbol","gene","Gene","gene_symbol")
possible <- setdiff(names(kin), c("ensembl_gene_id","hgnc_id","description","go_id","mgi_id"))
sym_col <- intersect(symbol_cols, names(kin))[1]
if (is.null(sym_col) || is.na(sym_col)) {
  possible <- setdiff(names(kin), c("ensembl_gene_id","hgnc_id","description","go_id","mgi_id"))
  if (length(possible) == 0) stop("Cannot detect gene symbol column in baseline kinome CSV")
  sym_col <- possible[1]
}

# --- Source-specific logic ---
if (tolower(source_type) == "kinhub") {
  # --- KinHub logic ---
  km_file <- file.path(outdir, "kinhub_mapping_raw.tsv")
  # Fetch KinHub mapping if missing
  if (!file.exists(km_file)) {
    cat("kinhub_mapping_raw.tsv not found, fetching KinHub HTML...\n", file=stderr())
    fetch_script <- file.path(repo_root, "scripts", "kinases", "bin", "fetch_kinhub_html_and_parse.R")
    if (!file.exists(fetch_script)) stop("fetch_kinhub_html_and_parse.R not found in scripts/kinases/bin")
    cmd <- sprintf('Rscript "%s" "%s"', fetch_script, km_file)
    status <- system(cmd)
    if (status != 0 || !file.exists(km_file)) stop("Failed to fetch/parse KinHub HTML to TSV")
  }
  km <- fread(km_file, header=FALSE, sep="\t", na.strings=c("","NA"))
  setnames(km, c("V1","V2","V3","V4"), c("HGNC","Group","Family","SubFamily"))
  km[, HGNC := toupper(trimws(HGNC))]
  # Optionally fetch Ensembl IDs for mapping
  unique_syms <- unique(na.omit(km$HGNC))
  fetch_by_symbol <- function(sym){
    url <- paste0("https://rest.genenames.org/fetch/symbol/", URLencode(sym))
    res <- tryCatch(httr::GET(url, httr::add_headers(Accept="application/json"), httr::timeout(10)), error=function(e) NULL)
    if(is.null(res) || httr::status_code(res)!=200) return(NA_character_)
    jd <- tryCatch(content(res, as="text", encoding="UTF-8"), error=function(e) NULL)
    if(is.null(jd)) return(NA_character_)
    jd2 <- tryCatch(jsonlite::fromJSON(jd, simplifyVector=FALSE), error=function(e) NULL)
    if(is.null(jd2) || is.null(jd2$response$docs) || length(jd2$response$docs)<1) return(NA_character_)
    doc <- jd2$response$docs[[1]]
    if(!is.null(doc$ensembl_gene_id)) return(as.character(doc$ensembl_gene_id))
    if(!is.null(doc$ensembl_gene_ids)) return(paste(unlist(doc$ensembl_gene_ids), collapse=","))
    return(NA_character_)
  }
  sym_map <- data.table(HGNC=character(), ensembl=character())
  for(sym in unique_syms){
    Sys.sleep(0.05)
    en <- fetch_by_symbol(sym)
    sym_map <- rbind(sym_map, data.table(HGNC=sym, ensembl=en))
  }
  km2 <- merge(km, sym_map, by="HGNC", all.x=TRUE)
  fwrite(km2, file.path(outdir, "kinhub_mapping.csv"))
  cat(sprintf("KinHub mapping written: %s\n", file.path(outdir, "kinhub_mapping.csv")), file=stderr())
} else if (tolower(source_type) == "hgnc") {
  # --- HGNC logic ---
  headers <- add_headers(Accept = "application/json", `User-Agent` = "bRNA3F/1.0 (+https://github.com/bRNA3F)")
  base_url <- "https://rest.genenames.org/fetch/symbol/"
  symbols <- unique(na.omit(trimws(kin[[sym_col]])))
  cat(sprintf("Fetching HGNC metadata for %d unique symbols...\n", length(symbols)), file=stderr())
  results <- vector("list", length(symbols))
  for (i in seq_along(symbols)) {
    sym <- symbols[i]
    url <- paste0(base_url, URLencode(sym, reserved=TRUE))
    res <- tryCatch(GET(url, headers, timeout(30)), error = function(e) e)
    if (!inherits(res, "error") && inherits(res, "response") && status_code(res) == 200) {
      dat <- tryCatch(content(res, as="text", encoding="UTF-8"), error=function(e) NULL)
      if (!is.null(dat)) {
        json <- tryCatch(fromJSON(dat), error=function(e) NULL)
        if (!is.null(json) && !is.null(json$response$docs) && length(json$response$docs) > 0) {
          doc <- json$response$docs[[1]]
          doc$queried_symbol <- sym
          results[[i]] <- as.data.table(doc)
        } else {
          results[[i]] <- data.table(queried_symbol=sym, hgnc_found=FALSE)
        }
      } else {
        results[[i]] <- data.table(queried_symbol=sym, hgnc_found=FALSE)
      }
    } else {
      results[[i]] <- data.table(queried_symbol=sym, hgnc_found=FALSE)
    }
    if (i %% 25 == 0) cat(sprintf("%d/%d symbols processed...\n", i, length(symbols)), file=stderr())
  }
  out_dt <- rbindlist(results, use.names=TRUE, fill=TRUE)
  fwrite(out_dt, file.path(outdir, "hgnc_kinase_groups.csv"))
  cat(sprintf("HGNC metadata written: %s\n", file.path(outdir, "hgnc_kinase_groups.csv")), file=stderr())
} else {
  stop("Unknown --source argument: must be 'kinhub' or 'hgnc'")
}
cat("Done.\n", file=stderr())
