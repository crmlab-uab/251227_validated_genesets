# fetch_mouse_kinome_from_kinhub.R
# Author: bRNA3F AI Agent
# Created: 2025-12-27 19:07 CST  # All timestamps must use Central Time (America/Chicago)
# Purpose: Fetch the definitive mouse kinome by mapping the KinHub human kinase list to mouse orthologs, then annotate with all required IDs (Symbol, Ensembl, Entrez, UniProt, MGI, Description).
# NOTE: Per project standards, all timestamps must use Central Time (America/Chicago), not UTC.

# --- SETUP ---
library(data.table)
library(biomaRt)
library(httr)

# --- 1. Download KinHub human kinase list (parse HTML table) ---
library(xml2)
library(rvest)
kinhub_url <- "http://www.kinhub.org/kinases.html"
cat("Fetching KinHub human kinase table...\n", file=stderr())
page <- httr::GET(kinhub_url)
html <- httr::content(page, as = "text")
doc <- read_html(html)
# Find the first table (the kinase table)
tbl <- html_element(doc, "table")
tbl_df <- html_table(tbl, fill=TRUE)
# Remove empty columns and rows
tbl_df <- tbl_df[rowSums(is.na(tbl_df) | tbl_df == "") < ncol(tbl_df), ]
tbl_df <- tbl_df[,colSums(is.na(tbl_df) | tbl_df == "") < nrow(tbl_df)]
# Remove rows with no UniProt ID (should be last column)
tbl_df <- tbl_df[grepl("P[0-9]|Q[0-9]", tbl_df[[ncol(tbl_df)]]), ]
# Standardize column names
colnames(tbl_df) <- c("Alias1","Alias2","HGNC_Symbol","Description","Group","Family","Subfamily","UniProt")
human_kinases <- as.data.table(tbl_df)
cat(sprintf("Parsed %d human kinases from KinHub\n", nrow(human_kinases)), file=stderr())
## Use HGNC_Symbol and UniProt for mapping

# --- 2. Map human kinases to mouse orthologs using Ensembl BioMart ---
cat("Mapping to mouse orthologs via BioMart...\n", file=stderr())
human = useMart("ensembl", dataset="hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset="mmusculus_gene_ensembl")
# Get orthologs for all human kinases, with retry on HTTP 500
get_orthos_safe <- function() {
    tryCatch({
        getLDS(attributes = c("hgnc_symbol", "uniprotswissprot", "ensembl_gene_id"),
                     filters = "hgnc_symbol",
                     values = human_kinases$HGNC_Symbol,
                     mart = human,
                     attributesL = c("mgi_symbol", "ensembl_gene_id", "uniprotswissprot", "mgi_id"),
                     martL = mouse)
    }, error = function(e) {
        cat("BioMart error: ", conditionMessage(e), "\nRetrying in 30 seconds...\n", file=stderr())
        Sys.sleep(30)
        tryCatch({
            getLDS(attributes = c("hgnc_symbol", "uniprotswissprot", "ensembl_gene_id"),
                         filters = "hgnc_symbol",
                         values = human_kinases$HGNC_Symbol,
                         mart = human,
                         attributesL = c("mgi_symbol", "ensembl_gene_id", "uniprotswissprot", "mgi_id"),
                         martL = mouse)
        }, error = function(e2) {
            cat("BioMart failed after retry: ", conditionMessage(e2), "\n", file=stderr())
            return(data.table())
        })
    })
}
orthos = get_orthos_safe()
if (nrow(orthos) == 0) {
    cat("BioMart mapping failed. Attempting to use local ortholog file...\n", file=stderr())
    local_orth_file <- "human_mouse_orthologs.csv"
    if (file.exists(local_orth_file)) {
        orthos <- fread(local_orth_file)
        cat(sprintf("Loaded %d orthologs from %s\n", nrow(orthos), local_orth_file), file=stderr())
    } else {
        stop("BioMart mapping failed and no local ortholog file found.")
    }
}
setDT(orthos)
if (!all(c("HGNC_Symbol","Human_UniProt","Human_Ensembl","Mouse_Symbol","Mouse_Ensembl","Mouse_UniProt","MGI_ID") %in% names(orthos))) {
    stop("Local ortholog file does not have required columns.")
}

# --- 3. Clean and deduplicate ---
cat("Cleaning and deduplicating...\n", file=stderr())
kinome_mouse <- unique(orthos[, .(Mouse_Symbol, Mouse_Ensembl, Mouse_Entrez, Mouse_UniProt, MGI_ID, Mouse_Desc, HGNC_Symbol, Human_UniProt)])
kinome_mouse <- kinome_mouse[!is.na(Mouse_Symbol) & Mouse_Symbol != ""]

# --- 4. Save output ---
outfile <- "mouse_kinome_from_kinhub.csv"
fwrite(kinome_mouse, outfile)
cat(sprintf("Done. Output: %s (%d kinases)\n", outfile, nrow(kinome_mouse)), file=stderr())
