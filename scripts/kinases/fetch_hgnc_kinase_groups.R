# Author: C. Ryan Miller
# Created: 2025-12-28 02:17 CST
# Commit: 26ec675324c74c530b55664519f79329a3d427d8

# fetch_hgnc_kinase_groups.R
# Fetch HGNC kinase group annotation using the HGNC REST API and save as hgnc_kinase_groups.csv

library(httr)
library(jsonlite)
library(data.table)

# Unified HGNC kinase groups fetcher
# Uses the HGNC REST API. Two query styles are supported:
#  - mode = "fetch" (default): uses /fetch/group/kinase
#  - mode = "search": uses /search/status:Approved+AND+group:kinase
# This script writes `hgnc_kinase_groups.csv` in the working directory.

mode <- Sys.getenv("HGNC_FETCH_MODE", "fetch")
if (!(mode %in% c("fetch", "search"))) mode <- "fetch"

if (mode == "fetch") {
	url <- "https://rest.genenames.org/fetch/group/kinase"
} else {
	url <- "https://rest.genenames.org/search/status:Approved+AND+group:kinase"
}

headers <- add_headers(Accept = "application/json")

cat(sprintf("Querying HGNC REST API (mode=%s) for kinase group annotation...\n", mode), file=stderr())
res <- GET(url, headers)
stop_for_status(res)
json <- content(res, as = "text", encoding = "UTF-8")
data <- fromJSON(json)

kinases <- data$response$docs
kinases <- as.data.table(kinases)

cols <- intersect(c("symbol", "name", "group", "subgroup", "family", "subfamily", "hgnc_id", "ensembl_gene_id"), names(kinases))
kinases <- kinases[, ..cols]

fwrite(kinases, "hgnc_kinase_groups.csv")
cat(sprintf("Done. Output: %s (%d kinases)\n", "hgnc_kinase_groups.csv", nrow(kinases)), file=stderr())
