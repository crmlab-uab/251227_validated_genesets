# Author: C. Ryan Miller
# Created: 2025-12-28 02:17 CST
# Commit: 26ec675324c74c530b55664519f79329a3d427d8

#!/usr/bin/env Rscript
library(data.table)
library(httr)
library(jsonlite)

cache_file <- "hgnc_lookup_cache.rds"
load_cache <- function() {
  if (file.exists(cache_file)) return(readRDS(cache_file))
  return(list())
}
save_cache <- function(cache) saveRDS(cache, cache_file)

hgnc_query <- function(sym_or_id) {
  # sym_or_id may be a symbol or HGNC:ID
  cache <- load_cache()
  key <- toupper(trimws(as.character(sym_or_id)))
  if (key=="NA" || key=="" || is.na(key)) return(list(symbol=NA,hgnc_id=NA,aliases=character(),prev_symbols=character(),ensembl=NA))
  if (!is.null(cache[[key]])) return(cache[[key]])

  base <- "https://rest.genenames.org"
  out <- list(symbol=NA, hgnc_id=NA, aliases=character(), prev_symbols=character(), ensembl=NA)
  # try search by symbol
  url <- paste0(base, "/search/symbol/", URLencode(sym_or_id, reserved=TRUE))
  resp <- try(GET(url, accept_json()), silent=TRUE)
  if (!inherits(resp, "try-error") && status_code(resp)==200) {
    j <- fromJSON(rawToChar(resp$content), simplifyVector = TRUE)
    if (!is.null(j$response$docs) && length(j$response$docs)>=1) {
      doc <- j$response$docs[[1]]
      out$symbol <- ifelse(!is.null(doc$symbol), doc$symbol, NA)
      out$hgnc_id <- ifelse(!is.null(doc$hgnc_id), doc$hgnc_id, NA)
      out$aliases <- ifelse(!is.null(doc$alias_symbol), doc$alias_symbol, character())
      out$prev_symbols <- ifelse(!is.null(doc$prev_symbol), doc$prev_symbol, character())
      out$ensembl <- ifelse(!is.null(doc$ensembl_gene_id), doc$ensembl_gene_id, NA)
    }
  }

  # if still empty, try fetch by HGNC id
  if ((is.na(out$hgnc_id) || out$hgnc_id=="") && grepl("^HGNC:", key)) {
    url2 <- paste0(base, "/fetch/hgnc_id/", URLencode(key, reserved=TRUE))
    resp2 <- try(GET(url2, accept_json()), silent=TRUE)
    if (!inherits(resp2, "try-error") && status_code(resp2)==200) {
      j2 <- fromJSON(rawToChar(resp2$content), simplifyVector = TRUE)
      if (!is.null(j2$response$docs) && length(j2$response$docs)>=1) {
        doc <- j2$response$docs[[1]]
        out$symbol <- ifelse(!is.null(doc$symbol), doc$symbol, out$symbol)
        out$hgnc_id <- ifelse(!is.null(doc$hgnc_id), doc$hgnc_id, out$hgnc_id)
        out$aliases <- unique(c(out$aliases, ifelse(!is.null(doc$alias_symbol), doc$alias_symbol, character())))
        out$prev_symbols <- unique(c(out$prev_symbols, ifelse(!is.null(doc$prev_symbol), doc$prev_symbol, character())))
        out$ensembl <- ifelse(!is.null(doc$ensembl_gene_id), doc$ensembl_gene_id, out$ensembl)
      }
    }
  }

  cache[[key]] <- out
  save_cache(cache)
  return(out)
}

main <- function() {
  # work in the kinases folder (script expected to be placed in that folder)
  kin_dir <- normalizePath(".", mustWork = TRUE)
  setwd(kin_dir)

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

  kin_f <- cfg_get("input_files.base_gene_list", "kinases_human.csv")
  # Manning CSV moved into kinases/data/; default relative to kinases working dir
  manning_f <- cfg_get("input_files.manning", "data/manning_2002_TableS1.csv")
  # Step guard: skip if augment_aliases is disabled
  if (!cfg_get("steps.augment_aliases", TRUE)) {
    cat("Skipping alias/Ensembl augmentation per genesets_config.yaml\n")
    quit(status=0)
  }
  if (!file.exists(kin_f)) stop("kinases_human.csv not found in working dir")
  if (!file.exists(manning_f)) stop("manning_2002_TableS1.csv not found in working dir")

  kin <- fread(kin_f, na.strings=c("","NA"))
  man <- fread(manning_f, na.strings=c("","NA"))

  # normalize manning symbol column name guesses
  man_col_candidates <- c("Symbol","Gene","gene_symbol","gene","GENE")
  man_sym_col <- intersect(names(man), man_col_candidates)
  if (length(man_sym_col)==0) {
    man_sym_col <- names(man)[1]
  } else {
    man_sym_col <- man_sym_col[1]
  }
  setnames(man, man_sym_col, "manning_symbol")

  kin[, HGNC_key := toupper(trimws(HGNC_gene_name))]
  kin[, HGNC_ID_val := ifelse(is.na(HGNC_ID), NA, as.character(HGNC_ID))]

  man[, man_symbol_up := toupper(trimws(manning_symbol))]
  # resolve HGNC info for Manning symbols (cached)
  man[, hgnc_res := lapply(manning_symbol, function(s) {
    tryCatch(hgnc_query(s), error=function(e) list(symbol=NA,hgnc_id=NA,aliases=character(),prev_symbols=character(),ensembl=NA))
  })]

  man[, hgnc_id := sapply(hgnc_res, function(x) ifelse(is.null(x$hgnc_id)||x$hgnc_id=="", NA, x$hgnc_id))]
  man[, hgnc_symbol := sapply(hgnc_res, function(x) ifelse(is.null(x$symbol)||x$symbol=="", NA, x$symbol))]
  man[, ensembl_m := sapply(hgnc_res, function(x) ifelse(is.null(x$ensembl)||x$ensembl=="", NA, x$ensembl))]
  man[, aliases := lapply(hgnc_res, function(x) unique(c(ifelse(is.null(x$alias_symbol), character(), x$alias_symbol), ifelse(is.null(x$prev_symbol), character(), x$prev_symbol), ifelse(is.null(x$symbol), character(), x$symbol)))) ]

  kin_matches <- copy(kin)
  kin_matches[, Manning_Group_new := NA_character_]
  kin_matches[, Manning_Family_new := NA_character_]
  kin_matches[, Manning_Subfamily_new := NA_character_]

  matched_count <- 0L
  for (i in seq_len(nrow(man))) {
    hgnc_id <- man[i, hgnc_id]
    ens_m <- toupper(as.character(man[i, ensembl_m]))
    alias_vec <- toupper(unlist(man[i, aliases]))

    matched_idx <- integer(0)
    if (!is.na(hgnc_id) && hgnc_id!="NA") {
      matched_idx <- which(!is.na(kin$HGNC_ID) & kin$HGNC_ID==hgnc_id)
    }
    if (length(matched_idx)==0 && !is.na(ens_m) && ens_m!="NA") {
      matched_idx <- which(toupper(kin$ensembl_gene_id)==ens_m)
    }
    if (length(matched_idx)==0 && length(alias_vec)>0) {
      for (a in alias_vec) {
        idx1 <- which(toupper(kin$HGNC_gene_name)==a)
        idx2 <- which(toupper(kin$external_gene_name)==a)
        if (length(idx1)>0) matched_idx <- idx1
        if (length(idx2)>0) matched_idx <- unique(c(matched_idx, idx2))
        if (length(matched_idx)>0) break
      }
    }

    if (length(matched_idx)>0) {
      pick <- matched_idx[1]
      annotated_idxs <- matched_idx[!is.na(kin$Manning_Group[matched_idx])]
      if (length(annotated_idxs)>0) pick <- annotated_idxs[1]
      kin_matches[pick, Manning_Group_new := man[i, ifelse("Group"%in%names(man), get("Group"), NA_character_) ]]
      kin_matches[pick, Manning_Family_new := man[i, ifelse("Family"%in%names(man), get("Family"), NA_character_) ]]
      kin_matches[pick, Manning_Subfamily_new := man[i, ifelse("SubFamily"%in%names(man), get("SubFamily"), NA_character_) ]]
      matched_count <- matched_count + 1L
    }
  }

  kin_out <- copy(kin)
  # prefer existing Manning columns; fill from new when empty
  for (col in c("Manning_Group","Manning_Family","Manning_Subfamily")) {
    newcol_src <- paste0(col, "_new")
    if (!(col %in% names(kin_out))) kin_out[, (col) := NA_character_]
    if (newcol_src %in% names(kin_matches)) {
      kin_out[, (col) := ifelse(is.na(get(col)) | get(col)=="", kin_matches[[newcol_src]], get(col))]
    }
  }

  outf <- "kinases_human.with_aliases_matches.csv"
  fwrite(kin_out, outf)
  cat(sprintf("Wrote %s\n", outf))
  cat(sprintf("Manning alias/Ensembl matches added: %d\n", matched_count))
}

main()
