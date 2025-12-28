## Read repo-level genesets config if present (allows modular runs)
if (requireNamespace("yaml", quietly=TRUE)) {
  cfg_candidates <- c("../genesets_config.yaml", "../config/genesets_config.yaml", "genesets_config.yaml")
  cfgf <- cfg_candidates[file.exists(cfg_candidates)][1]
  if (!is.na(cfgf) && nzchar(cfgf)) {
    cfg <- yaml::read_yaml(cfgf)
  } else {
    cfg <- list()
  }
} else {
  cfg <- list()
}

# helpers to get values with defaults
cfg_get <- function(path, default) {
  parts <- strsplit(path, "\\.")[[1]]
  cur <- cfg
  for (p in parts) {
    if (is.null(cur[[p]])) return(default)
    cur <- cur[[p]]
  }
  return(cur)
}

# default inputs/outputs
man_f <- cfg_get("input_files.manning", "manning_2002_TableS1.csv")
base_gene_f <- cfg_get("input_files.base_gene_list", "kinases_human.csv")
out_human <- cfg_get("outputs.human", "kinases_human.csv")
out_mouse <- cfg_get("outputs.mouse", "kinases_mouse.csv")

# Step guard: skip build if disabled
if (!cfg_get("steps.build", TRUE)) {
  cat("Skipping build step per genesets_config.yaml\n")
  quit(status=0)
}

# Parse species argument (default: human)
if (requireNamespace("optparse", quietly=TRUE)) {
  library(optparse)
  option_list <- list(
    make_option(c("-s","--species"), type="character", default="human", help="Species: mouse or human", metavar="character")
  )
  opt <- parse_args(OptionParser(option_list=option_list))
  species <- tolower(opt$species)
} else {
  species <- "human"
}

# Determine dataset, outfile, orgdb, symbol_col, kegg code
if (species == "mouse") {
  dataset <- "mmusculus_gene_ensembl"
  outfile <- out_mouse
  orgdb <- NULL
  kegg_code <- "mmu"
} else {
  dataset <- "hsapiens_gene_ensembl"
  outfile <- out_human
  orgdb <- NULL
  kegg_code <- "hsa"
}

# Load required libraries
suppressWarnings(suppressMessages({
  library(data.table)
  library(biomaRt)
  library(KEGGREST)
  library(AnnotationDbi)
  if (species == "human") library(org.Hs.eg.db)
  if (species == "mouse") library(org.Mm.eg.db)
}))

# 1. Fetch all protein-coding genes with kinase activity (GO:0004672, GO:0004674, GO:0016301, GO:0016773)
mart <- useMart("ensembl", dataset=dataset)
kinase_go_terms <- c("GO:0004672", "GO:0004674", "GO:0016301", "GO:0016773")
cat("Querying BioMart for kinase genes...\n", file=stderr())
all_kinases <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description", "go_id"),
                     filters = "go",
                     values = kinase_go_terms,
                     mart = mart)
setDT(all_kinases)
all_kinases <- unique(all_kinases[, .(ensembl_gene_id, external_gene_name, description)])

# 2. Annotate kinase group/family from HGNC (human only)
if (species == "human") {
  cat("Fetching HGNC kinase group annotation...\n", file=stderr())
  url <- "https://rest.genenames.org/search/status:Approved+AND+group:kinase"
  headers <- add_headers(Accept = "application/json")
  res <- GET(url, headers)
  stop_for_status(res)
  # Handle raw response for older httr
  # Download and parse HGNC gene group DB table
  hgnc_url <- "https://www.genenames.org/cgi-bin/genegroup/download-all"
  cat("\nDownloading HGNC gene group DB table from: ", hgnc_url, "\n", file=stderr())
  tmpfile <- tempfile(fileext = ".txt")
  download.file(hgnc_url, tmpfile, quiet = TRUE)
  hgnc_groups <- fread(tmpfile, sep = "\t", header = TRUE, quote = '"')
  # Filter for any group name containing 'kinase' (case-insensitive)
  kinases_hgnc <- hgnc_groups[grepl("kinase", `Group name`, ignore.case=TRUE)]
  cat("\nHGNC kinase-related groups: ", nrow(kinases_hgnc), " genes\n", file=stderr())
  # Defensive: check for expected columns
  if (!"Approved symbol" %in% names(kinases_hgnc)) stop("HGNC gene group table missing 'Approved symbol' column")
  setnames(kinases_hgnc, "Approved symbol", "external_gene_name")
  # Only keep relevant columns if present
  keep_cols <- intersect(c("external_gene_name", "Group name", "Family name", "Subfamily name", "HGNC ID"), names(kinases_hgnc))
  kinases_hgnc <- kinases_hgnc[, ..keep_cols]
  # Defensive: ensure HGNC ID is character and propagate NA if missing
  if ("HGNC ID" %in% names(kinases_hgnc)) kinases_hgnc$`HGNC ID` <- as.character(kinases_hgnc$`HGNC ID`)
  all_kinases <- merge(all_kinases, kinases_hgnc, by="external_gene_name", all.x=TRUE)
  if (!"HGNC ID" %in% names(all_kinases)) all_kinases$`HGNC ID` <- NA_character_
  # Annotate 'Protein' column: Y if gene has a non-NA Group name (i.e., present in HGNC kinase group)
  all_kinases$Protein <- ifelse(!is.na(all_kinases$`Group name`) & all_kinases$`Group name` != "", "Y", "N")
}

# 3. Annotate metabolic and lipid kinases using KEGG
# --- Restore [ xyz ] tag parsing and ID extraction ---
if (species == "mouse") {
  # Extract MGI ID (e.g., MGI:123456) from description if present
  mgi_pat <- "MGI:[0-9]+"
  mgi_match <- regexpr(mgi_pat, all_kinases$description)
  all_kinases$MGI_ID <- ifelse(mgi_match > 0, regmatches(all_kinases$description, mgi_match), NA_character_)
  # Remove bracketed [ xyz ] text from description
  all_kinases$description <- gsub("\\s*\\[[^]]*]", "", all_kinases$description)
}
if (species == "human") {
  # Remove bracketed [ xyz ] text from description
  all_kinases$description <- gsub("\\s*\\[[^]]*]", "", all_kinases$description)
}
cat("Annotating metabolic and lipid kinases using KEGG...\n", file=stderr())
# Metabolic
metabolic_pathways <- keggList("pathway", kegg_code)
metabolic_ids <- names(metabolic_pathways)[grepl("Metabolic pathways|metabolism", metabolic_pathways, ignore.case = TRUE)]
metabolic_genes <- unique(unlist(lapply(metabolic_ids, function(pid) {
  kegg_genes <- keggGet(pid)[[1]]$GENE
  if (is.null(kegg_genes)) return(NULL)
  gsub(" .*", "", kegg_genes)
})))
# Lipid
lipid_pathways <- keggList("pathway", kegg_code)
lipid_ids <- names(lipid_pathways)[grepl("lipid", lipid_pathways, ignore.case = TRUE)]
lipid_genes <- unique(unlist(lapply(lipid_ids, function(pid) {
  kegg_genes <- keggGet(pid)[[1]]$GENE
  if (is.null(kegg_genes)) return(NULL)
  gsub(" .*", "", kegg_genes)
})))
# Map KEGG Entrez IDs to gene symbols
entrez_to_symbol <- AnnotationDbi::select(orgdb, keys=unique(c(metabolic_genes, lipid_genes)), keytype="ENTREZID", columns=c("SYMBOL"))
metabolic_syms <- unique(na.omit(entrez_to_symbol$SYMBOL[entrez_to_symbol$ENTREZID %in% metabolic_genes]))
lipid_syms <- unique(na.omit(entrez_to_symbol$SYMBOL[entrez_to_symbol$ENTREZID %in% lipid_genes]))
# Add columns
all_kinases$Metabolic <- ifelse(all_kinases$external_gene_name %in% metabolic_syms, "Y", "N")
all_kinases$Lipid <- ifelse(all_kinases$external_gene_name %in% lipid_syms, "Y", "N")

# 4. Save final annotated kinome table
if (species == "mouse") {
  # Filter out Gm* and *Rik genes
  before_n <- nrow(all_kinases)
  all_kinases <- all_kinases[!grepl("^Gm[0-9]+$", external_gene_name, ignore.case=TRUE) & !grepl("Rik$", external_gene_name, ignore.case=TRUE)]
  after_n <- nrow(all_kinases)
  cat(sprintf("Filtered out Gm/Rik genes: %d removed\n", before_n - after_n), file=stderr())
}
if (species == "human") {
  # Ensure Group name column is present for human
  if (!"Group name" %in% names(all_kinases)) all_kinases$`Group name` <- NA_character_
}
# Exclude rows with missing or blank gene symbol
before_n <- nrow(all_kinases)
all_kinases <- all_kinases[!is.na(external_gene_name) & external_gene_name != ""]
after_n <- nrow(all_kinases)
cat(sprintf("Filtered out missing/blank gene symbols: %d removed\n", before_n - after_n), file=stderr())
fwrite(all_kinases, outfile)
cat(sprintf("Done. Output: %s (%d kinases)\n", outfile, nrow(all_kinases)), file=stderr())
