#!/usr/bin/env Rscript
# Author: C. Ryan Miller
# Created: 2025-12-28 02:17 CST
# Commit: copied from build_kinome_annotation.R

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
man_f <- cfg_get("input_files.manning", "data/manning_2002_TableS1.csv")
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

# Provide AnnotationDb object for downstream select() calls
if (species == "human") {
  orgdb <- org.Hs.eg.db
} else {
  orgdb <- org.Mm.eg.db
}

# 1. Fetch all protein-coding genes with kinase activity (GO:0004672, GO:0004674, GO:0016301, GO:0016773)
mart <- useMart("ensembl", dataset=dataset)
kinase_go_terms <- c("GO:0004672", "GO:0004674", "GO:0016301", "GO:0016773")
cat("Querying BioMart for source genes...\n", file=stderr())
all_src <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description", "go_id"),
                 filters = "go",
                 values = kinase_go_terms,
                 mart = mart)
setDT(all_src)
all_src <- unique(all_src[, .(ensembl_gene_id, external_gene_name, description)])

# 2. Annotate group/family from HGNC (human only)
if (species == "human") {
  cat("Fetching HGNC group annotation...\n", file=stderr())
  hgnc_url <- "https://www.genenames.org/cgi-bin/genegroup/download-all"
  tmpfile <- tempfile(fileext = ".txt")
  download.file(hgnc_url, tmpfile, quiet = TRUE)
  hgnc_groups <- fread(tmpfile, sep = "\t", header = TRUE, quote = '"')
  kinases_hgnc <- hgnc_groups[grepl("kinase", `Group name`, ignore.case=TRUE)]
  if (!"Approved symbol" %in% names(kinases_hgnc)) stop("HGNC gene group table missing 'Approved symbol' column")
  setnames(kinases_hgnc, "Approved symbol", "external_gene_name")
  keep_cols <- intersect(c("external_gene_name", "Group name", "Family name", "Subfamily name", "HGNC ID"), names(kinases_hgnc))
  kinases_hgnc <- kinases_hgnc[, ..keep_cols]
  if ("HGNC ID" %in% names(kinases_hgnc)) kinases_hgnc$`HGNC ID` <- as.character(kinases_hgnc$`HGNC ID`)
  all_src <- merge(all_src, kinases_hgnc, by="external_gene_name", all.x=TRUE)
  if (!"HGNC ID" %in% names(all_src)) all_src$`HGNC ID` <- NA_character_
  all_src$Protein <- ifelse(!is.na(all_src$`Group name`) & all_src$`Group name` != "", "Y", "N")
}

# 3. Annotate metabolic and lipid pathways using KEGG
cat("Annotating metabolic and lipid genes using KEGG...\n", file=stderr())
metabolic_pathways <- keggList("pathway", kegg_code)
metabolic_ids <- names(metabolic_pathways)[grepl("Metabolic pathways|metabolism", metabolic_pathways, ignore.case = TRUE)]
metabolic_genes <- unique(unlist(lapply(metabolic_ids, function(pid) {
  kegg_genes <- keggGet(pid)[[1]]$GENE
  if (is.null(kegg_genes)) return(NULL)
  gsub(" .*", "", kegg_genes)
})))
lipid_ids <- names(metabolic_pathways)[grepl("lipid", metabolic_pathways, ignore.case = TRUE)]
lipid_genes <- unique(unlist(lapply(lipid_ids, function(pid) {
  kegg_genes <- keggGet(pid)[[1]]$GENE
  if (is.null(kegg_genes)) return(NULL)
  gsub(" .*", "", kegg_genes)
})))
entrez_to_symbol <- AnnotationDbi::select(orgdb, keys=unique(c(metabolic_genes, lipid_genes)), keytype="ENTREZID", columns=c("SYMBOL"))
metabolic_syms <- unique(na.omit(entrez_to_symbol$SYMBOL[entrez_to_symbol$ENTREZID %in% metabolic_genes]))
lipid_syms <- unique(na.omit(entrez_to_symbol$SYMBOL[entrez_to_symbol$ENTREZID %in% lipid_genes]))
all_src$Metabolic <- ifelse(all_src$external_gene_name %in% metabolic_syms, "Y", "N")
all_src$Lipid <- ifelse(all_src$external_gene_name %in% lipid_syms, "Y", "N")

# 4. Save final annotated source table
if (species == "mouse") {
  before_n <- nrow(all_src)
  all_src <- all_src[!grepl("^Gm[0-9]+$", external_gene_name, ignore.case=TRUE) & !grepl("Rik$", external_gene_name, ignore.case=TRUE)]
  after_n <- nrow(all_src)
  cat(sprintf("Filtered out Gm/Rik genes: %d removed\n", before_n - after_n), file=stderr())
}
before_n <- nrow(all_src)
all_src <- all_src[!is.na(external_gene_name) & external_gene_name != ""]
after_n <- nrow(all_src)
cat(sprintf("Filtered out missing/blank gene symbols: %d removed\n", before_n - after_n), file=stderr())
fwrite(all_src, outfile)
cat(sprintf("Done. Output: %s (%d genes)\n", outfile, nrow(all_src)), file=stderr())
