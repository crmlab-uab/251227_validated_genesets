out_dt <- rbindlist(results, use.names = TRUE, fill = TRUE)
fwrite(out_dt, outfile)
cat(sprintf("Done. Output: %s (%d rows)\n", outfile, nrow(out_dt)), file=stderr())
headers <- add_headers(Accept = "application/json", `User-Agent` = "bRNA3F/1.0 (+https://github.com/bRNA3F)")

do_request <- function(url, tries = 3) {
	for (i in seq_len(tries)) {
		if (i > 1) Sys.sleep(1 * (i-1))
		res <- tryCatch(GET(url, headers, timeout(30)), error = function(e) e)
		if (inherits(res, "error")) next
		return(res)
	}
	return(res)
candidates <- list.files(input_dir, pattern = "human_kinome__.*\\.csv$", full.names = TRUE)
sym_col <- intersect(symbol_cols, names(kin))[1]
#!/usr/bin/env Rscript
# Robust HGNC metadata fetch for kinases baseline (cleaned and deduplicated)
outdir <- file.path("genesets","curated","kinases","val_sources")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
outfile <- file.path(outdir, "hgnc_kinase_groups.csv")
do_request <- function(url, tries = 3) {
	for (i in seq_len(tries)) {
		if (i > 1) Sys.sleep(1 * (i-1))

#!/usr/bin/env Rscript
# Robust HGNC metadata fetch for kinases baseline (CLEANED, DEDUPLICATED, SINGLE SCRIPT)
suppressWarnings(suppressMessages({
	library(httr)
	library(jsonlite)
	library(data.table)
}))

# Determine repo root and output locations
cmdArgs <- commandArgs(trailingOnly = FALSE)
fileArgIdx <- grep("--file=", cmdArgs)
fileArg <- if (length(fileArgIdx) > 0) sub("--file=", "", cmdArgs[fileArgIdx[1]]) else NA_character_
script_dir <- if (!is.na(fileArg) && nzchar(fileArg)) dirname(normalizePath(fileArg)) else normalizePath(".")
repo_root <- normalizePath(file.path(script_dir, "..", "..", ".."))

input_dir <- file.path(repo_root, "genesets","curated","kinases","inputs")
outdir <- file.path(repo_root, "genesets","curated","kinases","val_sources")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
outfile <- file.path(outdir, "hgnc_kinase_groups.csv")

# Robust input file selection: prefer kinases_human.csv, else latest kinases_human__*.csv or human_kinome__*.csv
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

# Detect symbol column
symbol_cols <- c("external_gene_name","hgnc_symbol","symbol","gene","Gene","gene_symbol")
possible <- setdiff(names(kin), c("ensembl_gene_id","hgnc_id","description","go_id","mgi_id"))
sym_col <- intersect(symbol_cols, names(kin))[1]
if (is.null(sym_col) || is.na(sym_col)) {
	possible <- setdiff(names(kin), c("ensembl_gene_id","hgnc_id","description","go_id","mgi_id"))
	if (length(possible) == 0) stop("Cannot detect gene symbol column in baseline kinome CSV")
	sym_col <- possible[1]
}

# Prepare HGNC REST API
headers <- add_headers(Accept = "application/json", `User-Agent` = "bRNA3F/1.0 (+https://github.com/bRNA3F)")
base_url <- "https://rest.genenames.org/fetch/symbol/"

do_request <- function(url, tries = 3) {
	for (i in seq_len(tries)) {
		if (i > 1) Sys.sleep(1 * (i-1))
		res <- tryCatch(GET(url, headers, timeout(30)), error = function(e) e)
		if (inherits(res, "error")) next
		if (inherits(res, "response") && status_code(res) == 200) return(res)
	}
	return(NULL)
}

# Fetch HGNC metadata for each symbol
symbols <- unique(na.omit(trimws(kin[[sym_col]])))
cat(sprintf("Fetching HGNC metadata for %d unique symbols...\n", length(symbols)), file=stderr())
results <- vector("list", length(symbols))
for (i in seq_along(symbols)) {
	sym <- symbols[i]
	url <- paste0(base_url, URLencode(sym, reserved=TRUE))
	res <- do_request(url)
	if (!is.null(res)) {
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

# Combine and write output
out_dt <- rbindlist(results, use.names=TRUE, fill=TRUE)
fwrite(out_dt, outfile)
cat(sprintf("Done. Output: %s (%d rows)\n", outfile, nrow(out_dt)), file=stderr())
		res <- tryCatch(GET(url, headers, timeout(30)), error = function(e) e)
		if (inherits(res, "error")) next
		return(res)
	}
	return(res)
cat("Using baseline kinome input:", latest, "\n", file=stderr())
kin <- tryCatch(fread(latest), error = function(e) stop("Failed to read baseline kinome CSV: ", e$message))

# detect symbol column
symbol_cols <- c("external_gene_name","hgnc_symbol","symbol","gene","Gene","gene_symbol")
sym_col <- intersect(symbol_cols, names(kin))[1]
if (is.null(sym_col) || is.na(sym_col)) {
	# fallback: pick first non-id like column
	possible <- setdiff(names(kin), c("ensembl_gene_id","hgnc_id","description","go_id","mgi_id"))
	if (length(possible) == 0) stop("Cannot detect gene symbol column in baseline kinome CSV")
	## Updated HGNC fetch: query HGNC per-symbol using baseline kinome input
	## This avoids relying on deprecated/unsupported 'group' query fields.

	library(httr)
	library(jsonlite)
	library(data.table)

	outdir <- file.path("genesets","curated","kinases","val_sources")
	dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
	outfile <- file.path(outdir, "hgnc_kinase_groups.csv")

	headers <- add_headers(Accept = "application/json", `User-Agent` = "bRNA3F/1.0 (+https://github.com/bRNA3F)")

	do_request <- function(url, tries = 3) {
		for (i in seq_len(tries)) {
			if (i > 1) Sys.sleep(1 * (i-1))
			res <- tryCatch(GET(url, headers, timeout(30)), error = function(e) e)
			if (inherits(res, "error")) next
			return(res)
		}
		return(res)
	}


	symbols <- unique(na.omit(kin[[sym_col]]))
	cat(sprintf("Querying HGNC for %d unique symbols (sample)...\n", length(symbols)), file=stderr())

	results <- list()
	base_url <- "https://rest.genenames.org/fetch/symbol/"
	count <- 0
	for (s in symbols) {
		count <- count + 1

	#!/usr/bin/env Rscript
	# Robust HGNC metadata fetch for kinases baseline (rewritten)
	suppressWarnings(suppressMessages({
		library(httr)
		library(jsonlite)
		library(data.table)
	}))

	# --- Setup paths ---
	cmdArgs <- commandArgs(trailingOnly = FALSE)
	fileArgIdx <- grep("--file=", cmdArgs)
	fileArg <- if (length(fileArgIdx) > 0) sub("--file=", "", cmdArgs[fileArgIdx[1]]) else NA_character_
	script_dir <- if (!is.na(fileArg) && nzchar(fileArg)) dirname(normalizePath(fileArg)) else normalizePath(".")
	repo_root <- normalizePath(file.path(script_dir, "..", "..", ".."))
	input_dir <- file.path(repo_root, "genesets","curated","kinases","inputs")
	outdir <- file.path(repo_root, "genesets","curated","kinases","val_sources")
	dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
	outfile <- file.path(outdir, "hgnc_kinase_groups.csv")

	# --- Robust input file selection ---
	preferred <- file.path(input_dir, "kinases_human.csv")
	if (file.exists(preferred)) {
		kin_file <- preferred
	} else {
		# Accept kinases_human__*.csv or human_kinome__*.csv, pick latest
		candidates <- list.files(input_dir, pattern = "^(kinases_human|human_kinome)__.*\\.csv$", full.names = TRUE)
		if (length(candidates) == 0) stop("No kinases_human.csv or kinases_human__*.csv or human_kinome__*.csv found in ", input_dir)
		candidate_info <- file.info(candidates)
		kin_file <- rownames(candidate_info)[which.max(candidate_info$mtime)]
	}
	cat("Using baseline kinome input:", kin_file, "\n", file=stderr())
	kin <- tryCatch(fread(kin_file), error = function(e) stop("Failed to read baseline kinome CSV: ", e$message))

	# --- Detect symbol column ---
	symbol_cols <- c("external_gene_name","hgnc_symbol","symbol","gene","Gene","gene_symbol")
	sym_col <- intersect(symbol_cols, names(kin))[1]
	if (is.null(sym_col) || is.na(sym_col)) {
		possible <- setdiff(names(kin), c("ensembl_gene_id","hgnc_id","description","go_id","mgi_id"))
		if (length(possible) == 0) stop("Cannot detect gene symbol column in baseline kinome CSV")
		sym_col <- possible[1]
	}
	cat(sprintf("Detected symbol column: %s (n=%d)\n", sym_col, nrow(kin)), file=stderr())
	symbols <- unique(na.omit(kin[[sym_col]]))
	cat(sprintf("Querying HGNC for %d unique symbols...\n", length(symbols)), file=stderr())

	# --- HGNC API fetch ---
	headers <- add_headers(Accept = "application/json", `User-Agent` = "bRNA3F/1.0 (+https://github.com/bRNA3F)")
	do_request <- function(url, tries = 3) {
		for (i in seq_len(tries)) {
			if (i > 1) Sys.sleep(1 * (i-1))
			res <- tryCatch(GET(url, headers, timeout(30)), error = function(e) e)
			if (inherits(res, "error")) next
			return(res)
		}
		return(res)
	}

	results <- vector("list", length(symbols))
	names(results) <- as.character(symbols)
	base_url <- "https://rest.genenames.org/fetch/symbol/"
	count <- 0
	for (s in symbols) {
		count <- count + 1
		enc <- URLencode(as.character(s), reserved = TRUE)
		url <- paste0(base_url, enc)
		res <- do_request(url)
		if (inherits(res, "error")) {
			results[[as.character(s)]] <- data.table(query_symbol = s, status = NA_character_)
			next
		}
		st <- status_code(res)
		if (st == 404) {
			results[[as.character(s)]] <- data.table(query_symbol = s, status = "not_found")
			next
		}
		if (st >= 400) {
			# fallback to search endpoint
			search_url <- paste0("https://rest.genenames.org/search/", enc)
			res2 <- tryCatch(GET(search_url, headers, timeout(30)), error = function(e) e)
			if (!inherits(res2, "error") && status_code(res2) >=200 && status_code(res2) < 300) res <- res2 else {
				results[[as.character(s)]] <- data.table(query_symbol = s, status = sprintf("http_%d", st))
				next
			}
		}
		txt <- tryCatch(content(res, as = "text", encoding = "UTF-8"), error = function(e) NA_character_)
		if (is.na(txt) || nchar(txt) == 0) { results[[as.character(s)]] <- data.table(query_symbol = s, status = "empty_response"); next }
		parsed <- tryCatch(fromJSON(txt, simplifyVector = FALSE), error = function(e) NULL)
		if (is.null(parsed) || is.null(parsed$response) || is.null(parsed$response$docs) || length(parsed$response$docs)<1) {
			results[[as.character(s)]] <- data.table(query_symbol = s, status = "no_docs"); next }
		doc <- parsed$response$docs[[1]]
		row <- data.table(
			query_symbol = s,
			symbol = if (!is.null(doc$symbol)) as.character(doc$symbol) else NA_character_,
			name = if (!is.null(doc$name)) as.character(doc$name) else NA_character_,
			hgnc_id = if (!is.null(doc$hgnc_id)) as.character(doc$hgnc_id) else NA_character_,
			status = if (!is.null(doc$status)) as.character(doc$status) else NA_character_,
			locus_group = if (!is.null(doc$locus_group)) as.character(doc$locus_group) else NA_character_,
			gene_group = if (!is.null(doc$gene_group)) as.character(doc$gene_group) else NA_character_,
			gene_group_id = if (!is.null(doc$gene_group_id)) as.character(doc$gene_group_id) else NA_character_,
			ensembl_gene_id = if (!is.null(doc$ensembl_gene_id)) as.character(doc$ensembl_gene_id) else NA_character_
		)
		results[[as.character(s)]] <- row
		if (count %% 10 == 0) Sys.sleep(0.5)
	}

	out_dt <- rbindlist(results, use.names = TRUE, fill = TRUE)
	fwrite(out_dt, outfile)
	cat(sprintf("Done. Output: %s (%d rows)\n", outfile, nrow(out_dt)), file=stderr())
