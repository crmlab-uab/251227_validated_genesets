# fetch_hgnc_kinase_groups.R
# Fetch HGNC kinase group annotation using the HGNC REST API and save as hgnc_kinase_groups.csv

library(httr)
library(jsonlite)
library(data.table)

# HGNC REST API endpoint for kinases
url <- "https://rest.genenames.org/fetch/group/kinase"

# Set headers for HGNC API
headers <- add_headers(Accept = "application/json")

# Fetch data
cat("Querying HGNC REST API for kinase group annotation...\n", file=stderr())
res <- GET(url, headers)
stop_for_status(res)
json <- content(res, as = "text", encoding = "UTF-8")
data <- fromJSON(json)

# Parse results
kinases <- data$response$docs
kinases <- as.data.table(kinases)

# Select relevant columns (HGNC symbol, name, family, subfamily, group, etc.)
cols <- intersect(c("symbol", "name", "group", "subgroup", "family", "subfamily", "hgnc_id", "ensembl_gene_id"), names(kinases))
kinases <- kinases[, ..cols]

# Save to CSV
fwrite(kinases, "hgnc_kinase_groups.csv")
cat(sprintf("Done. Output: %s (%d kinases)\n", "hgnc_kinase_groups.csv", nrow(kinases)), file=stderr())
