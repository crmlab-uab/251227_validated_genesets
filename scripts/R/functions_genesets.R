# Author: C. Ryan Miller
# Created: 2025-12-28 02:25 CST
# Commit: 85c49284e5027d5940219124b28324bb8eba3205

## Composable functions for genesets pipeline
library(data.table)
library(httr)
library(jsonlite)

read_config <- function(cfgf="genesets_config.yaml") {
  if (file.exists(cfgf) && requireNamespace("yaml", quietly=TRUE)) {
    return(yaml::read_yaml(cfgf))
  }
  return(list())
}

hgnc_query_cached <- function(query, cache_file="kinases/hgnc_lookup_cache.rds") {
  cache <- if (file.exists(cache_file)) readRDS(cache_file) else list()
  key <- toupper(trimws(as.character(query)))
  if (key=="" || is.na(key)) return(list(symbol=NA, hgnc_id=NA, aliases=character(), prev_symbols=character(), ensembl=NA))
  if (!is.null(cache[[key]])) return(cache[[key]])
  base <- "https://rest.genenames.org"
  out <- list(symbol=NA, hgnc_id=NA, aliases=character(), prev_symbols=character(), ensembl=NA)
  # try symbol search
  url <- paste0(base, "/search/symbol/", URLencode(query, reserved=TRUE))
  resp <- try(GET(url, accept_json()), silent=TRUE)
  if (!inherits(resp, "try-error") && status_code(resp)==200) {
    j <- fromJSON(rawToChar(resp$content), simplifyVector=TRUE)
    if (!is.null(j$response$docs) && length(j$response$docs)>=1) {
      doc <- j$response$docs[[1]]
      if (!is.list(doc)) doc <- as.list(doc)
      out$symbol <- ifelse(!is.null(doc$symbol), doc$symbol, NA)
      out$hgnc_id <- ifelse(!is.null(doc$hgnc_id), doc$hgnc_id, NA)
      out$aliases <- ifelse(!is.null(doc$alias_symbol), doc$alias_symbol, character())
      out$prev_symbols <- ifelse(!is.null(doc$prev_symbol), doc$prev_symbol, character())
      out$ensembl <- ifelse(!is.null(doc$ensembl_gene_id), doc$ensembl_gene_id, NA)
    }
  }
  if ((is.na(out$hgnc_id) || out$hgnc_id=="") && grepl("^HGNC:", key)) {
    url2 <- paste0(base, "/fetch/hgnc_id/", URLencode(key, reserved=TRUE))
    resp2 <- try(GET(url2, accept_json()), silent=TRUE)
    if (!inherits(resp2, "try-error") && status_code(resp2)==200) {
      j2 <- fromJSON(rawToChar(resp2$content), simplifyVector=TRUE)
      if (!is.null(j2$response$docs) && length(j2$response$docs)>=1) {
        doc <- j2$response$docs[[1]]
        if (!is.list(doc)) doc <- as.list(doc)
        out$symbol <- ifelse(!is.null(doc$symbol), doc$symbol, out$symbol)
        out$hgnc_id <- ifelse(!is.null(doc$hgnc_id), doc$hgnc_id, out$hgnc_id)
        out$aliases <- unique(c(out$aliases, ifelse(!is.null(doc$alias_symbol), doc$alias_symbol, character())))
        out$prev_symbols <- unique(c(out$prev_symbols, ifelse(!is.null(doc$prev_symbol), doc$prev_symbol, character())))
        out$ensembl <- ifelse(!is.null(doc$ensembl_gene_id), doc$ensembl_gene_id, out$ensembl)
      }
    }
  }
  cache[[key]] <- out
  dir.create(dirname(cache_file), showWarnings=FALSE, recursive=TRUE)
  saveRDS(cache, cache_file)
  return(out)
}

fetch_genes <- function(species="human", dataset=NULL) {
  library(biomaRt)
  library(KEGGREST)
  if (is.null(dataset)) dataset <- ifelse(tolower(species)=="mouse","mmusculus_gene_ensembl","hsapiens_gene_ensembl")
  mart <- useMart("ensembl", dataset=dataset)
  kinase_go_terms <- c("GO:0004672", "GO:0004674", "GO:0016301", "GO:0016773")
  all_kinases <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description", "go_id"), filters = "go", values = kinase_go_terms, mart = mart)
  setDT(all_kinases)
  all_kinases <- unique(all_kinases[, .(ensembl_gene_id, external_gene_name, description)])
  return(all_kinases)
}

map_hgnc <- function(genes_dt, cache_file="kinases/hgnc_lookup_cache.rds") {
  genes <- copy(genes_dt)
  genes[, HGNC_gene_name := NA_character_]
  genes[, HGNC_ID := NA_character_]
  uniq <- unique(na.omit(toupper(trimws(genes$external_gene_name))))
  for (s in uniq) {
    res <- hgnc_query_cached(s, cache_file=cache_file)
    if (!is.null(res$symbol) && !is.na(res$symbol)) {
      genes[toupper(trimws(external_gene_name))==s, HGNC_gene_name := res$symbol]
    }
    if (!is.null(res$hgnc_id) && !is.na(res$hgnc_id)) {
      genes[toupper(trimws(external_gene_name))==s, HGNC_ID := res$hgnc_id]
    }
  }
  return(genes)
}

augment_aliases <- function(genes_dt, manning_dt, cache_file="kinases/hgnc_lookup_cache.rds") {
  # Try to expand matches using HGNC aliases and Ensembl crosswalk
  kin <- copy(genes_dt)
  kin[, HGNC_gene_name_up := toupper(trimws(HGNC_gene_name))]
  man <- copy(manning_dt)
  # normalize manning symbol col
  man_name_col <- intersect(names(man), c("Name","Symbol","Gene"))[1]
  if (is.na(man_name_col)) man_name_col <- names(man)[1]
  setnames(man, man_name_col, "manning_symbol")
  man[, man_up := toupper(trimws(manning_symbol))]
  # build alias map from HGNC for manning symbols
  alias_map <- list()
  for (s in unique(na.omit(man$man_up))) {
    res <- hgnc_query_cached(s, cache_file=cache_file)
    aliases <- unique(c(res$aliases, res$prev_symbols, res$symbol))
    for (a in aliases) alias_map[[toupper(a)]] <- res$symbol
  }
  # apply: if external gene name matches alias, set HGNC_gene_name to approved
  kin[, ext_up := toupper(trimws(external_gene_name))]
  kin[, HGNC_gene_name := ifelse(ext_up %in% names(alias_map) & (is.na(HGNC_gene_name) | HGNC_gene_name==""), alias_map[ext_up], HGNC_gene_name)]
  kin[, HGNC_ID := HGNC_ID] # no-op to preserve column
  kin[, c("HGNC_gene_name_up","ext_up") := NULL]
  return(kin)
}

apply_manning <- function(genes_dt, manning_dt) {
  kin <- copy(genes_dt)
  # normalize manning names
  man_name_col <- intersect(names(manning_dt), c("Name","Symbol","Gene"))[1]
  if (is.na(man_name_col)) man_name_col <- names(manning_dt)[1]
  setnames(manning_dt, man_name_col, "manning_symbol")
  man_small <- unique(manning_dt[, .(manning_symbol, Group, Family, Subfamily)])
  man_small[, man_up := toupper(trimws(manning_symbol))]
  kin[, join_sym := toupper(trimws(ifelse(!is.na(HGNC_gene_name) & HGNC_gene_name!="", HGNC_gene_name, external_gene_name)))]
  man_small[, man_up := toupper(trimws(manning_symbol))]
  setnames(man_small, "man_up", "Name_up")
  setnames(man_small, "manning_symbol", "Name")
  man_small[, Name_up := toupper(trimws(Name))]
  merged <- merge(kin, man_small, by.x="join_sym", by.y="Name_up", all.x=TRUE, sort=FALSE)
  # copy Manning columns to standard names
  if ("Group" %in% names(merged)) merged[, Manning_Group := Group]
  if ("Family" %in% names(merged)) merged[, Manning_Family := Family]
  if ("Subfamily" %in% names(merged)) merged[, Manning_Subfamily := Subfamily]
  merged[, join_sym := NULL]
  return(merged)
}

dedupe_genes <- function(genes_dt) {
  dt <- copy(genes_dt)
  dt[, ext_up := toupper(trimws(external_gene_name))]
  dt[, ensembl_num := suppressWarnings(as.numeric(gsub("^ENSG0*", "", ensembl_gene_id)))]
  keep_idx <- dt[, {
    rows <- .I
    has_manning <- which(!is.na(Manning_Group) | !is.na(Manning_Family) | !is.na(Manning_Subfamily))
    if (length(has_manning) > 0) {
      vals <- ensembl_num[has_manning]
      if (all(is.na(vals))) .(keep = rows[has_manning[1]]) else .(keep = rows[has_manning[which.min(ifelse(is.na(vals), Inf, vals))]])
    } else {
      vals_all <- ensembl_num
      if (all(is.na(vals_all))) .(keep = rows[1]) else .(keep = rows[which.min(ifelse(is.na(vals_all), Inf, vals_all))])
    }
  }, by=ext_up]$keep
  res <- dt[keep_idx]
  res[, c("ext_up","ensembl_num") := NULL]
  return(res)
}

merge_kinhub_dt <- function(genes_dt, kinhub_dt) {
  # kinhub_dt expected to have HGNC or ensembl column
  kin <- copy(genes_dt)
  kh <- copy(kinhub_dt)
  if (!"ensembl" %in% names(kh) && "HGNC" %in% names(kh)) {
    # attempt to map HGNC to ensembl via HGNC REST
    kh[, ensembl := NA_character_]
    for (i in seq_len(nrow(kh))) {
      sym <- kh$HGNC[i]
      res <- hgnc_query_cached(sym)
      kh$ensembl[i] <- res$ensembl
    }
  }
  kin[, ensembl := toupper(trimws(ensembl_gene_id))]
  kh[, ensembl := toupper(trimws(ensembl))]
  merged <- merge(kin, kh, by.x="ensembl", by.y="ensembl", all.x=TRUE, sort=FALSE)
  return(merged)
}

export_gmt <- function(genes_dt, outfile) {
  # Minimal GMT exporter: one set 'KINASES'
  genes <- unique(na.omit(genes_dt$external_gene_name))
  txt <- c(paste0("KINASES\tKinase gene set\t", paste(genes, collapse="\t")))
  writeLines(txt, outfile)
  return(outfile)
}
