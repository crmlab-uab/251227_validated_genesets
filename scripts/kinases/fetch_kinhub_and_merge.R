# Author: C. Ryan Miller
# Created: 2025-12-28 02:17 CST
# Commit: 26ec675324c74c530b55664519f79329a3d427d8

#!/usr/bin/env Rscript
# read genesets config if present
if (requireNamespace("yaml", quietly=TRUE)) {
  cfg_candidates <- c("../genesets_config.yaml","../config/genesets_config.yaml","genesets_config.yaml")
  cfgf <- cfg_candidates[file.exists(cfg_candidates)][1]
  if (!is.na(cfgf) && nzchar(cfgf)) cfg <- yaml::read_yaml(cfgf) else cfg <- list()
} else cfg <- list()
cfg_get <- function(path, default) {
  parts <- strsplit(path, "\\.")[[1]]
  cur <- cfg
  for (p in parts) { if (is.null(cur[[p]])) return(default); cur <- cur[[p]] }
  cur
}
# override input filenames from config
kin_f <- cfg_get("input_files.base_gene_list", "kinases_human.csv")

# Step guard: skip if merge_kinhub is disabled
if (!cfg_get("steps.merge_kinhub", FALSE)) {
  cat("Skipping KinHub merge step per genesets_config.yaml\n")
  quit(status=0)
}

suppressWarnings(suppressMessages({
  library(data.table)
  library(httr)
  library(jsonlite)
}))
setwd("/data/251227_validated_genesets/kinases")
km <- fread("kinhub_mapping_raw.tsv", header=FALSE, sep="\t", na.strings=c("","NA"))
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
fwrite(km2, "kinhub_mapping.csv")
# Now merge with sources baseline (kinases_human.csv) by ensembl_gene_id
kin <- fread(kin_f, na.strings=c("","NA"))
# standardize ensembl ids
kin[, ensembl_gene_id := toupper(trimws(ensembl_gene_id))]
km2[, ensembl := toupper(trimws(ensembl))]
merged <- merge(kin, km2, by.x="ensembl_gene_id", by.y="ensembl", all.x=TRUE, sort=FALSE)
# report mismatches where both have Manning/Group
merged[, kin_group := ifelse(is.na(Manning_Group), NA_character_, Manning_Group)]
merged[, kh_group := ifelse(is.na(Group), NA_character_, Group)]
# counts
total <- nrow(merged)
matched_khub <- sum(!is.na(merged$kh_group))
# Report status (sources merged)
cat(sprintf("Total rows: %d\nRows matched to kinhub by Ensembl: %d\n", total, matched_khub))
# list conflicts where both present and differ
conflicts <- merged[!is.na(kin_group) & !is.na(kh_group) & kin_group!=kh_group, .(external_gene_name, ensembl_gene_id, kin_group, kh_group)]
cat(sprintf("Conflicts (group) count: %d\n", nrow(conflicts)))
if(nrow(conflicts)>0) print(head(conflicts,20))
# write merged file with kinhub fields prefixed
write.csv(merged, file="kinases_human.with_kinhub.csv", row.names=FALSE)
cat("Wrote kinases_human.with_kinhub.csv and kinhub_mapping.csv (KinHub fields merged into sources baseline)\n")
