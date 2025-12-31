# fetch_kinhub_html_and_parse.R
# Author: bRNA3F AI Agent
# Created: 2025-12-27 20:00 CST
# Purpose: Download and parse KinHub HTML table to TSV for validation source pipeline

library(httr)
library(xml2)
library(rvest)
library(data.table)

kinhub_url <- "http://www.kinhub.org/kinases.html"
outfile <- commandArgs(trailingOnly=TRUE)[1]
if (is.na(outfile) || outfile == "") outfile <- "kinhub_mapping_raw.tsv"

cat(sprintf("Fetching KinHub HTML from %s...\n", kinhub_url), file=stderr())
page <- httr::GET(kinhub_url)
if (httr::status_code(page) != 200) stop("Failed to fetch KinHub HTML")
html <- httr::content(page, as = "text")
doc <- read_html(html)
tbl <- html_element(doc, "table")
tbl_df <- html_table(tbl, fill=TRUE)
# Remove empty columns and rows
tbl_df <- tbl_df[rowSums(is.na(tbl_df) | tbl_df == "") < ncol(tbl_df), ]
tbl_df <- tbl_df[,colSums(is.na(tbl_df) | tbl_df == "") < nrow(tbl_df)]
# Remove rows with no UniProt ID (should be last column)
tbl_df <- tbl_df[grepl("P[0-9]|Q[0-9]", tbl_df[[ncol(tbl_df)]]), ]
# Standardize column names
colnames(tbl_df)[1:4] <- c("HGNC","Group","Family","SubFamily")
# Write as TSV (HGNC,Group,Family,SubFamily)
fwrite(tbl_df[,1:4], file=outfile, sep="\t", col.names=FALSE)
cat(sprintf("Wrote KinHub mapping TSV: %s\n", outfile), file=stderr())
