#!/usr/bin/env Rscript
# 02_fetch_validation_sources.R
# Author: bRNA3F AI Agent
# Updated: 2025-12-31 (refactored to use config_loader.R)
# Purpose: Fetch external validation sources (KinHub, HGNC) for kinase validation
# Usage: Rscript 02_fetch_validation_sources.R --source=kinhub|hgnc

suppressWarnings(suppressMessages({
  library(data.table)
  library(httr)
  library(jsonlite)
}))

# Load centralized config
.script_dir <- (function() {
 cmd_args <- commandArgs(trailingOnly = FALSE)
 file_arg <- grep("--file=", cmd_args, value = TRUE)
 if (length(file_arg) > 0) {
   return(dirname(normalizePath(sub("--file=", "", file_arg[1]))))
 }
 return(normalizePath("."))
})()
source(file.path(.script_dir, "lib", "config_loader.R"))

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(key, default = NULL) {
  val <- grep(paste0("^--", key, "="), args, value = TRUE)
  if (length(val)) sub(paste0("^--", key, "="), "", val[1]) else default
}
source_type <- get_arg("source", "kinhub")

# Ensure directories exist
ensure_dir(paths$input_dir)
ensure_dir(paths$output_dir)

# Find baseline kinome input file
canonical <- input_path("kinases_human_biomart.csv")
if (file.exists(canonical)) {
  kin_file <- canonical
} else {
  candidates <- list.files(paths$input_dir, pattern = "^kinases_human.*\\.csv$", full.names = TRUE)
  if (length(candidates) == 0) {
    stop("No kinases_human*.csv file found in ", paths$input_dir)
  }
  candidate_info <- file.info(candidates)
  kin_file <- rownames(candidate_info)[which.max(candidate_info$mtime)]
}
check_file_nonempty(kin_file, paste0("Baseline kinome input missing or empty: ", kin_file))
cat("Using baseline kinome input:", kin_file, "\n", file = stderr())
kin <- tryCatch(fread(kin_file), error = function(e) stop("Failed to read baseline kinome CSV: ", e$message))

# Symbol column detection
symbol_cols <- c("external_gene_name", "hgnc_symbol", "symbol", "gene", "Gene", "gene_symbol")
sym_col <- intersect(symbol_cols, names(kin))[1]
if (is.null(sym_col) || is.na(sym_col)) {
  possible <- setdiff(names(kin), c("ensembl_gene_id", "hgnc_id", "description", "go_id", "mgi_id"))
  if (length(possible) == 0) stop("Cannot detect gene symbol column in baseline kinome CSV")
  sym_col <- possible[1]
}

# Source-specific logic
if (tolower(source_type) == "kinhub") {
  # KinHub logic
  km_file <- output_path("kinhub_mapping_raw.tsv")

  # Fetch KinHub mapping if missing
  if (!file.exists(km_file)) {
    cat("kinhub_mapping_raw.tsv not found, fetching KinHub HTML...\n", file = stderr())
    fetch_script <- file.path(repo_root, "scripts", "kinases", "lib", "fetch_kinhub_html_and_parse.R")
    if (!file.exists(fetch_script)) stop("fetch_kinhub_html_and_parse.R not found in scripts/kinases/lib")
    cmd <- sprintf('Rscript "%s" "%s"', fetch_script, km_file)
    status <- system(cmd)
    if (status != 0 || !file.exists(km_file)) stop("Failed to fetch/parse KinHub HTML to TSV")
  }
  check_file_nonempty(km_file, paste0("KinHub mapping file missing or empty: ", km_file))

  km <- fread(km_file, header = FALSE, sep = "\t", na.strings = c("", "NA"))
  setnames(km, c("V1", "V2", "V3", "V4"), c("HGNC", "Group", "Family", "SubFamily"))
  km[, HGNC := toupper(trimws(HGNC))]

  # Fetch Ensembl IDs for mapping
  unique_syms <- unique(na.omit(km$HGNC))
  fetch_by_symbol <- function(sym) {
    url <- paste0("https://rest.genenames.org/fetch/symbol/", URLencode(sym))
    res <- tryCatch(httr::GET(url, httr::add_headers(Accept = "application/json"), httr::timeout(10)), error = function(e) NULL)
    if (is.null(res) || httr::status_code(res) != 200) return(NA_character_)
    jd <- tryCatch(content(res, as = "text", encoding = "UTF-8"), error = function(e) NULL)
    if (is.null(jd)) return(NA_character_)
    jd2 <- tryCatch(jsonlite::fromJSON(jd, simplifyVector = FALSE), error = function(e) NULL)
    if (is.null(jd2) || is.null(jd2$response$docs) || length(jd2$response$docs) < 1) return(NA_character_)
    doc <- jd2$response$docs[[1]]
    if (!is.null(doc$ensembl_gene_id)) return(as.character(doc$ensembl_gene_id))
    if (!is.null(doc$ensembl_gene_ids)) return(paste(unlist(doc$ensembl_gene_ids), collapse = ","))
    return(NA_character_)
  }

  sym_map <- data.table(HGNC = character(), ensembl = character())
  for (sym in unique_syms) {
    Sys.sleep(0.05)
    en <- fetch_by_symbol(sym)
    sym_map <- rbind(sym_map, data.table(HGNC = sym, ensembl = en))
  }

  km2 <- merge(km, sym_map, by = "HGNC", all.x = TRUE)
  kinmap_out <- output_path("kinases_human_kinhub.csv")
  fwrite(km2, kinmap_out)
  check_file_nonempty(kinmap_out, paste0("KinHub mapping output missing or empty: ", kinmap_out))
  cat(sprintf("KinHub mapping written: %s\n", kinmap_out), file = stderr())

} else if (tolower(source_type) == "hgnc") {
  # HGNC logic
  headers <- add_headers(Accept = "application/json", `User-Agent` = "bRNA3F/1.0 (+https://github.com/bRNA3F)")
  base_url <- "https://rest.genenames.org/fetch/symbol/"
  symbols <- unique(na.omit(trimws(kin[[sym_col]])))

  cat(sprintf("Fetching HGNC metadata for %d unique symbols...\n", length(symbols)), file = stderr())
  results <- vector("list", length(symbols))

  for (i in seq_along(symbols)) {
    sym <- symbols[i]
    url <- paste0(base_url, URLencode(sym, reserved = TRUE))
    res <- tryCatch(GET(url, headers, timeout(30)), error = function(e) e)

    if (!inherits(res, "error") && inherits(res, "response") && status_code(res) == 200) {
      dat <- tryCatch(content(res, as = "text", encoding = "UTF-8"), error = function(e) NULL)
      if (!is.null(dat)) {
        json <- tryCatch(fromJSON(dat), error = function(e) NULL)
        if (!is.null(json) && !is.null(json$response$docs) && length(json$response$docs) > 0) {
          doc <- json$response$docs[[1]]
          doc$queried_symbol <- sym
          results[[i]] <- as.data.table(doc)
        } else {
          results[[i]] <- data.table(queried_symbol = sym, hgnc_found = FALSE)
        }
      } else {
        results[[i]] <- data.table(queried_symbol = sym, hgnc_found = FALSE)
      }
    } else {
      results[[i]] <- data.table(queried_symbol = sym, hgnc_found = FALSE)
    }
    if (i %% 25 == 0) cat(sprintf("%d/%d symbols processed...\n", i, length(symbols)), file = stderr())
  }

  out_dt <- rbindlist(results, use.names = TRUE, fill = TRUE)
  hgnc_out <- output_path("kinases_human_hgnc.csv")
  fwrite(out_dt, hgnc_out)
  check_file_nonempty(hgnc_out, paste0("HGNC metadata output missing or empty: ", hgnc_out))
  cat(sprintf("HGNC metadata written: %s\n", hgnc_out), file = stderr())

} else {
  stop("Unknown --source argument: must be 'kinhub' or 'hgnc'")
}

cat("Done.\n", file = stderr())
