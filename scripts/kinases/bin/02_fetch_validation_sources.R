#!/usr/bin/env Rscript
# read genesets config if present
if (requireNamespace("yaml", quietly=TRUE)) {
  # robust config discovery: search upward from the script location and cwd
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  fileArgIdx <- grep("--file=", cmdArgs)
  fileArg <- if (length(fileArgIdx) > 0) sub("--file=", "", cmdArgs[fileArgIdx[1]]) else NA_character_
  script_dir <- if (!is.na(fileArg) && nzchar(fileArg)) dirname(normalizePath(fileArg)) else getwd()
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
  if (is.na(cfgf)) cfgf <- search_up(getwd(), "genesets_config.yaml", 8)
  if (!is.na(cfgf) && nzchar(cfgf)) cfg <- yaml::read_yaml(cfgf) else cfg <- list()
} else cfg <- list()
# repository root (scripts/kinases/bin -> ../../../)
repo_root <- normalizePath(file.path(script_dir, "..", "..", ".."))
cfg_get <- function(path, default) {
  parts <- strsplit(path, "\\.")[[1]]
  cur <- cfg
  for (p in parts) { if (is.null(cur[[p]])) return(default); cur <- cur[[p]] }
  cur
}
# override input filenames from config
kin_f <- cfg_get("input_files.base_gene_list", "kinases_human.csv")

# Step guard: skip if both merge_kinhub and merge_val_sources are disabled
if (!cfg_get("steps.merge_kinhub", FALSE) && !cfg_get("steps.merge_val_sources", FALSE)) {
  cat("Skipping KinHub/validation merge step per genesets_config.yaml\n")
  quit(status=0)
}

suppressWarnings(suppressMessages({
  library(data.table)
  library(httr)
  library(jsonlite)
}))
# Prefer kinases val_sources location for KinHub mapping file
val_dir <- file.path(repo_root, "genesets","curated","kinases","val_sources")
if (!dir.exists(val_dir)) dir.create(val_dir, recursive=TRUE, showWarnings=FALSE)
km_file <- file.path(val_dir, "kinhub_mapping_raw.tsv")
# If KinHub mapping TSV is missing and enabled in config, fetch and parse HTML
kinhub_enabled <- cfg_get("validation.kinhub_enabled", TRUE)
if (!file.exists(km_file)) {
  if (kinhub_enabled) {
    cat("kinhub_mapping_raw.tsv not found, fetching KinHub HTML...\n", file=stderr())
    fetch_script <- file.path(repo_root, "scripts", "kinases", "bin", "fetch_kinhub_html_and_parse.R")
    if (!file.exists(fetch_script)) stop("fetch_kinhub_html_and_parse.R not found in scripts/kinases/bin")
    cmd <- sprintf('Rscript "%s" "%s"', fetch_script, km_file)
    status <- system(cmd)
    if (status != 0 || !file.exists(km_file)) stop("Failed to fetch/parse KinHub HTML to TSV")
  } else {
    stop("Expected kinhub_mapping_raw.tsv in genesets/curated/kinases/val_sources — KinHub fetch/parse is disabled in config (validation.kinhub_enabled: false)")
  }
}
km <- fread(km_file, header=FALSE, sep="\t", na.strings=c("","NA"))
setnames(km, c("V1","V2","V3","V4"), c("HGNC","Group","Family","SubFamily"))
km[, HGNC := toupper(trimws(HGNC))]
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
  # try fields
  if(!is.null(doc$ensembl_gene_id)) return(as.character(doc$ensembl_gene_id))
  if(!is.null(doc$ensembl_gene_ids)) return(paste(unlist(doc$ensembl_gene_ids), collapse=","))
  # sometimes ensembl is in 'ensembl_gene_id' nested
  return(NA_character_)
}
# iterate
sym_map <- data.table(HGNC=character(), ensembl=character())
for(sym in unique_syms){
  Sys.sleep(0.05)
  en <- fetch_by_symbol(sym)
  sym_map <- rbind(sym_map, data.table(HGNC=sym, ensembl=en))
}
# merge back
km2 <- merge(km, sym_map, by="HGNC", all.x=TRUE)
val_dir_out <- file.path(repo_root, "genesets","curated","kinases","val_sources")
dir.create(val_dir_out, recursive=TRUE, showWarnings=FALSE)
fwrite(km2, file.path(val_dir_out, "kinhub_mapping.csv"))
# Determine baseline kinases file: prefer canonical inputs location, then configured base_gene_list, then outputs snapshots
kin <- NULL
canonical_input <- file.path(repo_root, "genesets","curated","kinases","inputs","kinases_human.csv")
if (file.exists(canonical_input)) {
  kin <- fread(canonical_input, na.strings=c("","NA"))
  cat("Using canonical input baseline:", canonical_input, "\n")
} else {
  # check configured kin_f in cwd or repo_root
  kin_f_repo <- file.path(repo_root, kin_f)
  if (!is.null(kin_f) && file.exists(kin_f)) {
    kin <- fread(kin_f, na.strings=c("","NA"))
    cat("Using configured baseline (cwd):", kin_f, "\n")
  } else if (!is.null(kin_f) && file.exists(kin_f_repo)) {
    kin <- fread(kin_f_repo, na.strings=c("","NA"))
    cat("Using configured baseline (repo_root):", kin_f_repo, "\n")
  } else {
    # search for canonical outputs
    possible <- list.files(file.path(repo_root, "genesets","curated","kinases","outputs"), pattern = "kinases_human.*\\.csv$", full.names=TRUE)
    if (length(possible) > 0) {
      kin <- fread(possible[1], na.strings=c("","NA"))
      cat("Using baseline from outputs:", possible[1], "\n")
    }
  }
}
if (is.null(kin)) stop("Sources baseline kinases_human.csv not found — expected genesets/curated/kinases/inputs/kinases_human.csv or configured kin_f or genesets/curated/kinases/outputs snapshots")
# standardize ensembl ids
kin[, ensembl_gene_id := toupper(trimws(ensembl_gene_id))]
km2[, ensembl := toupper(trimws(ensembl))]
merged <- merge(kin, km2, by.x="ensembl_gene_id", by.y="ensembl", all.x=TRUE, sort=FALSE)
# report mismatches where both have Manning/Group (guard if columns missing)
if ("Manning_Group" %in% names(merged)) {
  merged[, kin_group := ifelse(is.na(Manning_Group), NA_character_, Manning_Group)]
} else {
  merged[, kin_group := NA_character_]
}
if ("Group" %in% names(merged)) {
  merged[, kh_group := ifelse(is.na(Group), NA_character_, Group)]
} else {
  merged[, kh_group := NA_character_]
}
# counts
total <- nrow(merged)
matched_khub <- sum(!is.na(merged$kh_group))
# Report status (sources merged)
cat(sprintf("Total rows: %d\nRows matched to kinhub by Ensembl: %d\n", total, matched_khub))
# list conflicts where both present and differ
conflicts <- merged[!is.na(kin_group) & !is.na(kh_group) & kin_group!=kh_group, .(external_gene_name, ensembl_gene_id, kin_group, kh_group)]
cat(sprintf("Conflicts (group) count: %d\n", nrow(conflicts)))
if(nrow(conflicts)>0) print(head(conflicts,20))
out_dir <- file.path("genesets","curated","kinases","outputs")
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
date_tag <- format(Sys.time(), "%y%m%d")
kinhub_out <- file.path(out_dir, paste0("02_kinases_human.with_kinhub__", date_tag, ".csv"))
km_out <- file.path(out_dir, paste0("02_kinhub_mapping__", date_tag, ".csv"))
fwrite(merged, file = kinhub_out)
fwrite(km2, file = km_out)
# write md5s
write(sprintf('%s  %s', tools::md5sum(kinhub_out), basename(kinhub_out)), paste0(kinhub_out, '.md5'))
write(sprintf('%s  %s', tools::md5sum(km_out), basename(km_out)), paste0(km_out, '.md5'))
cat(sprintf("Wrote %s and %s (KinHub fields merged into baseline)\n", kinhub_out, km_out))
