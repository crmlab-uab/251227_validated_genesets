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
# Manning CSV moved into kinases/data/; default is relative to this script's working dir (kinases)
manning_f <- cfg_get("input_files.manning", "data/manning_2002_TableS1.csv")

# Step guard: skip if annotate_manning is disabled
if (!cfg_get("steps.annotate_manning", TRUE)) {
  cat("Skipping Manning annotation step per genesets_config.yaml\n")
  quit(status=0)
}

# Merge Manning classification into kinases_human.csv
suppressWarnings(suppressMessages({
  library(data.table)
  library(lubridate)
}))
wd <- dirname("/data/251227_validated_genesets/kinases/add_manning_annotation.R")
setwd(wd)
cat("Working dir:", getwd(), "\n")
# Files
manning_file <- manning_f
kin_file <- kin_f
if (!file.exists(manning_file)) stop("Manning file not found: ", manning_file)
if (!file.exists(kin_file)) stop("Kinases file not found: ", kin_file)
# Read
manning <- fread(manning_file, na.strings=c("", "NA"), showProgress=FALSE)
kin <- fread(kin_file, na.strings=c("", "NA"), showProgress=FALSE)

# Resolve HGNC IDs to approved HUGO symbols (use HGNC REST API) for better matching
suppressWarnings(suppressMessages({
  if (!requireNamespace("httr", quietly=TRUE)) install.packages("httr", repos = "https://cloud.r-project.org")
  if (!requireNamespace("jsonlite", quietly=TRUE)) install.packages("jsonlite", repos = "https://cloud.r-project.org")
  library(httr)
  library(jsonlite)
}))

# Extract HGNC IDs from `HGNC ID` column (format HGNC:###)
if ("HGNC ID" %in% names(kin)) {
  kin[, HGNC_ID_clean := trimws(`HGNC ID`)]
} else if ("HGNC_ID" %in% names(kin)) {
  kin[, HGNC_ID_clean := trimws(HGNC_ID)]
} else {
  kin[, HGNC_ID_clean := NA_character_]
}

# Function to query HGNC REST for one id
fetch_hgnc <- function(hgnc_id) {
  # Return a list with approved symbol, alias_symbol (vector), prev_symbol (vector)
  if (is.na(hgnc_id) || hgnc_id == "") return(list(approved=NA_character_, aliases=character(0), prevs=character(0)))
  url <- paste0("https://rest.genenames.org/fetch/hgnc_id/", URLencode(hgnc_id))
  res <- tryCatch(httr::GET(url, httr::add_headers(Accept = "application/json"), httr::timeout(10)), error = function(e) NULL)
  if (is.null(res) || httr::status_code(res) != 200) return(list(approved=NA_character_, aliases=character(0), prevs=character(0)))
  jtxt <- tryCatch(httr::content(res, as = "text", encoding = "UTF-8"), error=function(e) NULL)
  if (is.null(jtxt)) return(list(approved=NA_character_, aliases=character(0), prevs=character(0)))
  jd <- tryCatch(jsonlite::fromJSON(jtxt, simplifyVector = FALSE), error=function(e) NULL)
  if (is.null(jd)) return(list(approved=NA_character_, aliases=character(0), prevs=character(0)))
  docs <- jd$response$docs
  if (is.null(docs) || length(docs) < 1) return(list(approved=NA_character_, aliases=character(0), prevs=character(0)))
  doc1 <- if (is.data.frame(docs)) docs[1,] else docs[[1]]
  approved <- NA_character_
  aliases <- character(0)
  prevs <- character(0)
  if (is.list(doc1) && !is.null(doc1$symbol)) approved <- as.character(doc1$symbol)
  if (is.list(doc1) && !is.null(doc1$alias_symbol)) aliases <- unlist(doc1$alias_symbol)
  if (is.list(doc1) && !is.null(doc1$prev_symbol)) prevs <- unlist(doc1$prev_symbol)
  # Ensure character vectors
  aliases <- as.character(aliases)
  prevs <- as.character(prevs)
  return(list(approved=approved, aliases=aliases, prevs=prevs))
}

# Build unique HGNC ID list and query; also build alias -> approved mapping
unique_ids <- unique(na.omit(kin$HGNC_ID_clean))
id_info <- list()
alias_map <- list()
if (length(unique_ids) > 0) {
  for (id in unique_ids) {
    Sys.sleep(0.05)
    info <- fetch_hgnc(id)
    id_info[[id]] <- info
    appr <- ifelse(is.null(info$approved) || is.na(info$approved), NA_character_, info$approved)
    # Map aliases and previous symbols to approved
    all_aliases <- unique(c(info$aliases, info$prevs))
    if (length(all_aliases) > 0 && !is.na(appr)) {
      for (a in all_aliases) {
        if (!is.na(a) && nzchar(a)) alias_map[[toupper(a)]] <- appr
      }
    }
  }
}

# Create simple id -> approved named vector
id_map <- sapply(id_info, function(x) ifelse(is.null(x$approved) || is.na(x$approved), NA_character_, x$approved))

# Apply mapping back to kin table
kin[, HGNC_approved := ifelse(!is.na(HGNC_ID_clean) & HGNC_ID_clean %in% names(id_map), id_map[HGNC_ID_clean], NA_character_)]
kin[, HGNC_approved := unlist(HGNC_approved)]

# If external name matches an alias/previous symbol, map to approved
kin[, ext_upper := toupper(trimws(external_gene_name))]
kin[, mapped_from_alias := ifelse(ext_upper %in% names(alias_map), alias_map[ext_upper], NA_character_)]

# Final join_sym: prefer HGNC_approved, else mapped_from_alias, else external_gene_name
kin[, join_sym := ifelse(!is.na(HGNC_approved) & HGNC_approved != "", HGNC_approved,
                         ifelse(!is.na(mapped_from_alias) & mapped_from_alias != "", mapped_from_alias, external_gene_name))]
# Normalize symbol columns
setnames(manning, old=names(manning), new=names(manning))
# Manning symbols are in column 'Name'
if (!"Name" %in% names(manning)) stop("Manning file missing 'Name' column")
if (!"external_gene_name" %in% names(kin)) {
  # try other common names
  possible <- c("external_gene_name","gene","Symbol","external_gene_id")
  found <- intersect(possible, names(kin))
  if (length(found)==0) stop("kinases_human.csv missing expected gene symbol column")
}
# Prepare merge columns
manning_small <- unique(manning[, .(Name, Group, Family, Subfamily)])
# Remove any existing Manning_* columns from kin to avoid duplicate-suffixed columns from prior runs
old_mcols <- grep("^Manning_", names(kin), value = TRUE)
if (length(old_mcols) > 0) kin[, (old_mcols) := NULL]
# Rename to avoid collision
setnames(manning_small, c("Name","Group","Family","Subfamily"), c("Name","Manning_Group","Manning_Family","Manning_Subfamily"))
# Trim whitespace and uppercase (join_sym already set using HGNC_approved when available)
kin[, join_sym := toupper(trimws(join_sym))]
manning_small[, Name := toupper(trimws(Name))]
# Merge
merged <- merge(kin, manning_small, by.x="join_sym", by.y="Name", all.x=TRUE, sort=FALSE)
# Drop helper after initial merge only; we'll reuse for fuzzy matching

# Fuzzy matching for unmatched symbols
merged[, join_sym_orig := external_gene_name]
unmatched_idx <- which(is.na(merged$Manning_Group) & is.na(merged$Manning_Family) & is.na(merged$Manning_Subfamily))
if (length(unmatched_idx) > 0) {
  manning_names_list <- toupper(trimws(manning$Name))
  # Prepare a lookup table from original manning table
  # We'll attempt strategies: 1) remove trailing digits from symbol 2) exact prefix matches
  extra_matches <- 0
  for (i in unmatched_idx) {
    sym <- toupper(trimws(merged$join_sym[i]))
    if (is.na(sym) || sym == "") next
    # 1) strip trailing digits (e.g., ABL1 -> ABL)
    sym_nod <- sub("\\d+$", "", sym)
    found <- NA_character_
    if (sym_nod != sym && sym_nod %in% manning_names_list) found <- sym_nod
    # 2) try exact match to any Manning name
    if (is.na(found) && sym %in% manning_names_list) found <- sym
    # 3) try prefix: Manning name is prefix of symbol or vice versa
    if (is.na(found)) {
      hits <- manning_names_list[startsWith(sym, manning_names_list) | startsWith(manning_names_list, sym)]
      if (length(hits) == 1) found <- hits
      else if (length(hits) > 1) found <- hits[which.max(nchar(hits))]
    }
    if (!is.na(found)) {
      # get Manning fields from original manning table (case-insensitive)
      mrow <- manning[toupper(trimws(Name)) == found]
      if (nrow(mrow) >= 1) {
        merged$Manning_Group[i] <- as.character(mrow$Group[1])
        merged$Manning_Family[i] <- as.character(mrow$Family[1])
        merged$Manning_Subfamily[i] <- as.character(mrow$Subfamily[1])
        extra_matches <- extra_matches + 1
      }
    }
  }
  cat(sprintf("Fuzzy matching added %d extra matches\n", extra_matches))
}

# Drop helper
merged[, join_sym := NULL]
merged[, join_sym_orig := NULL]

# Deduplicate rows by external_gene_name: prefer Manning-annotated rows; otherwise pick row
# with the smallest numeric Ensembl ID (explicit). Record removed rows for audit.
dedupe_report <- list()
merged[, ext_up := toupper(trimws(external_gene_name))]
if (any(duplicated(merged$ext_up))) {
  # compute numeric ENSG value (NA when unparsable)
  merged[, ensembl_num := suppressWarnings(as.numeric(gsub("^ENSG0*", "", ensembl_gene_id)))]
  # identify keep indices per group
  keep_idx <- merged[, {
    rows <- .I
    has_manning <- which(!is.na(Manning_Group) | !is.na(Manning_Family) | !is.na(Manning_Subfamily))
    if (length(has_manning) > 0) {
      # among Manning-annotated pick smallest ensembl_num if available
      vals <- ensembl_num[has_manning]
      if (all(is.na(vals))) {
        .(keep = rows[has_manning[1]])
      } else {
        rel <- which.min(ifelse(is.na(vals), Inf, vals))
        .(keep = rows[has_manning[rel]])
      }
    } else {
      # no Manning rows: pick smallest ensembl_num among all rows
      vals_all <- ensembl_num
      if (all(is.na(vals_all))) {
        .(keep = rows[1])
      } else {
        rel <- which.min(ifelse(is.na(vals_all), Inf, vals_all))
        .(keep = rows[rel])
      }
    }
  }, by=ext_up]$keep

  # Build removed_rows report
  removed_rows <- list()
  all_idx <- seq_len(nrow(merged))
  to_remove <- setdiff(all_idx, keep_idx)
  if (length(to_remove) > 0) {
    # group removed by symbol
    for (sym in unique(merged$ext_up[to_remove])) {
      rows <- which(merged$ext_up == sym)
      rm_rows <- setdiff(rows, keep_idx)
      if (length(rm_rows) > 0) removed_rows[[sym]] <- merged[rm_rows, .(external_gene_name, ensembl_gene_id, description)]
    }
    # log removals
    for (nm in names(removed_rows)) {
      cat(sprintf("Dedup group %s: removed %s\n", nm, paste(removed_rows[[nm]]$ensembl_gene_id, collapse=",")))
    }
    merged <- merged[keep_idx]
  }
  dedupe_report$removed <- removed_rows
}
merged[, c("ext_up","ensembl_num") := NULL]

  # === Lookup HGNC IDs by Ensembl gene ID for kept rows (prefer authoritative mapping)
  # Function to fetch HGNC info by Ensembl gene id
  fetch_hgnc_by_ensembl <- function(ensembl_id) {
    if (is.na(ensembl_id) || ensembl_id == "") return(list(hgnc_id=NA_character_, approved=NA_character_))
    url <- paste0("https://rest.genenames.org/fetch/ensembl_gene_id/", URLencode(ensembl_id))
    res <- tryCatch(httr::GET(url, httr::add_headers(Accept = "application/json"), httr::timeout(10)), error = function(e) NULL)
    if (is.null(res) || httr::status_code(res) != 200) return(list(hgnc_id=NA_character_, approved=NA_character_))
    jtxt <- tryCatch(httr::content(res, as = "text", encoding = "UTF-8"), error=function(e) NULL)
    if (is.null(jtxt)) return(list(hgnc_id=NA_character_, approved=NA_character_))
    jd <- tryCatch(jsonlite::fromJSON(jtxt, simplifyVector = FALSE), error=function(e) NULL)
    if (is.null(jd)) return(list(hgnc_id=NA_character_, approved=NA_character_))
    docs <- jd$response$docs
    if (is.null(docs) || length(docs) < 1) return(list(hgnc_id=NA_character_, approved=NA_character_))
    doc1 <- if (is.data.frame(docs)) docs[1,] else docs[[1]]
    hgnc_id <- NA_character_
    approved <- NA_character_
    if (is.list(doc1) && !is.null(doc1$hgnc_id)) hgnc_id <- as.character(doc1$hgnc_id)
    if (is.list(doc1) && !is.null(doc1$symbol)) approved <- as.character(doc1$symbol)
    return(list(hgnc_id=hgnc_id, approved=approved))
  }

  # Query HGNC for all kept Ensembl IDs
  ensembl_ids <- unique(na.omit(merged$ensembl_gene_id))
  ensembl_map_hgnc <- list()
  if (length(ensembl_ids) > 0) {
    for (eid in ensembl_ids) {
      Sys.sleep(0.05)
      info <- fetch_hgnc_by_ensembl(eid)
      ensembl_map_hgnc[[eid]] <- info
    }
  }

  # Apply mappings: overwrite/create `HGNC_ID_clean` and `HGNC_approved` from authoritative lookup
  merged[, HGNC_ID_clean := ifelse(!is.na(ensembl_gene_id) & ensembl_gene_id %in% names(ensembl_map_hgnc),
                                   sapply(ensembl_map_hgnc[ensembl_gene_id], function(x) ifelse(is.null(x$hgnc_id), NA_character_, x$hgnc_id)),
                                   NA_character_)]
  merged[, HGNC_approved := ifelse(!is.na(ensembl_gene_id) & ensembl_gene_id %in% names(ensembl_map_hgnc),
                                   sapply(ensembl_map_hgnc[ensembl_gene_id], function(x) ifelse(is.null(x$approved), NA_character_, x$approved)),
                                   NA_character_)]

  # Remove any original parsed HGNC ID columns (from description or prior parsing)
  old_hgnc_cols <- intersect(c("HGNC ID","HGNC_ID"), names(merged))
  if (length(old_hgnc_cols) > 0) merged[, (old_hgnc_cols) := NULL]

  # Consolidate duplicate .x/.y columns produced from earlier merges: keep base name, coalesce values
  dup_suffs <- grep("\\.(x|y)$", names(merged), value=TRUE)
  if (length(dup_suffs) > 0) {
    bases <- unique(sub("\\.(x|y)$", "", dup_suffs))
    for (b in bases) {
      xcol <- paste0(b, ".x")
      ycol <- paste0(b, ".y")
      # Ensure base column exists
      if (!b %in% names(merged)) merged[, (b) := NA_character_]
      # Fill base with .x then .y where NA
      if (xcol %in% names(merged)) merged[is.na(get(b)) & !is.na(get(xcol)), (b) := get(xcol)]
      if (ycol %in% names(merged)) merged[is.na(get(b)) & !is.na(get(ycol)), (b) := get(ycol)]
      # Remove suffix columns
      to_drop <- intersect(c(xcol, ycol), names(merged))
      if (length(to_drop) > 0) merged[, (to_drop) := NULL]
    }
  }

  # Rename HGNC columns to requested names and drop helper alias columns
  if ("HGNC_ID_clean" %in% names(merged)) setnames(merged, "HGNC_ID_clean", "HGNC_ID")
  if ("HGNC_approved" %in% names(merged)) setnames(merged, "HGNC_approved", "HGNC_gene_name")
  # Drop helper columns used during matching
  drop_helpers <- intersect(c("ext_upper", "mapped_from_alias"), names(merged))
  if (length(drop_helpers) > 0) merged[, (drop_helpers) := NULL]
# Backup original
bak <- paste0(kin_file, ".bak.", format(now(), "%Y%m%dT%H%M%S"))
file.copy(kin_file, bak)
cat("Backup written to:", bak, "\n")
# Write updated CSV (preserve column order: original cols then new Manning cols if not already present)
# Ensure Manning columns exist in merged
for (c in c("Manning_Group","Manning_Family","Manning_Subfamily")) if (!c %in% names(merged)) merged[, (c) := NA]
fwrite(merged, kin_file)
# Report summary
n_total <- nrow(merged)
n_matched <- sum(!is.na(merged$Manning_Group) | !is.na(merged$Manning_Family) | !is.na(merged$Manning_Subfamily))
cat(sprintf("Merged Manning annotations: %d/%d matched (%.1f%%)\n", n_matched, n_total, 100* n_matched / n_total))
# List unmatched symbols (up to 20)
unmatched <- merged[is.na(Manning_Group) & is.na(Manning_Family) & is.na(Manning_Subfamily), external_gene_name]
if (length(unmatched)>0) {
  cat("Unmatched symbols (first 20):\n")
  cat(paste(head(unmatched,20), collapse=", "), "\n")
}
cat("Done.\n")
