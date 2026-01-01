#!/usr/bin/env Rscript
# 03_build_annotations.R
# Author: C. Ryan Miller
# Created: 2025-12-28 02:17 CST
# Updated: 2025-12-31 (refactored to use config_loader.R, added column provenance)
# Purpose: Build comprehensive kinase annotations from BioMart, HGNC, KEGG
# Usage: Rscript 03_build_annotations.R [--species=human|mouse]
#
# Column Provenance Tracking:
#   - Columns are suffixed with source identifier (e.g., _biomart, _hgnc, _kegg)
#   - A column_provenance.yaml file is generated alongside output
#   - Core identifier columns (ensembl_gene_id, external_gene_name) have no suffix

# Load centralized config (provides repo_root, paths, cfg, helper functions)
.script_dir <- (function() {
 cmd_args <- commandArgs(trailingOnly = FALSE)
 file_arg <- grep("--file=", cmd_args, value = TRUE)
 if (length(file_arg) > 0) {
   return(dirname(normalizePath(sub("--file=", "", file_arg[1]))))
 }
 return(normalizePath("."))
})()
source(file.path(.script_dir, "lib", "config_loader.R"))

# Step guard: skip build if disabled in config
if (!cfg_get("steps.build", TRUE)) {
  cat("Skipping build step per genesets_config.yaml\n")
  quit(status = 0)
}

# Parse species argument
if (requireNamespace("optparse", quietly = TRUE)) {
  library(optparse)
  option_list <- list(
    make_option(c("-s", "--species"), type = "character", default = NULL,
                help = "Species: mouse or human", metavar = "character")
  )
  opt <- parse_args(OptionParser(option_list = option_list))
  species <- if (!is.null(opt$species)) tolower(opt$species) else cfg_get("species", "human")
} else {
  species <- cfg_get("species", "human")
}

# Determine dataset and file paths based on species
if (species == "mouse") {
  dataset <- "mmusculus_gene_ensembl"
  base_gene_file <- input_path("kinases_mouse_biomart.csv")
  kegg_code <- "mmu"
} else {
  dataset <- "hsapiens_gene_ensembl"
  base_gene_file <- input_path("kinases_human_biomart.csv")
  kegg_code <- "hsa"
}

# Manning table path from config
manning_file <- file.path(repo_root, cfg_get("input_files.manning", "curated/kinases/inputs/manning_TableS1__251229.csv"))

# Check required input files
check_file_nonempty(manning_file, paste0("Manning supplement missing or empty: ", manning_file))
check_file_nonempty(base_gene_file, paste0("Base gene list missing or empty: ", base_gene_file))

# Load required libraries
suppressWarnings(suppressMessages({
  library(data.table)
  library(biomaRt)
  library(KEGGREST)
  library(AnnotationDbi)
  if (species == "human") library(org.Hs.eg.db)
  if (species == "mouse") library(org.Mm.eg.db)
  if (!requireNamespace("httr", quietly = TRUE)) {
    stop("Package 'httr' is required for HGNC REST calls. Install with: install.packages('httr')")
  }
}))

# Set OrgDb object
orgdb <- if (species == "human") org.Hs.eg.db else if (species == "mouse") org.Mm.eg.db else NULL

# Initialize column provenance tracker
column_provenance <- list()
add_provenance <- function(columns, source, description = "") {
  for (col in columns) {
    column_provenance[[col]] <<- list(
      source = source,
      description = description,
      added = format(Sys.time(), "%Y-%m-%d")
    )
  }
}

# 1. Fetch all protein-coding genes with kinase activity from BioMart
mart <- useMart("ensembl", dataset = dataset)
kinase_go_terms <- c("GO:0004672", "GO:0004674", "GO:0016301", "GO:0016773")
cat("Querying BioMart for kinase genes...\n", file = stderr())

all_kinases <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "description", "go_id"),
  filters = "go",
  values = kinase_go_terms,
  mart = mart
)
setDT(all_kinases)
all_kinases <- unique(all_kinases[, .(ensembl_gene_id, external_gene_name, description)])

# Rename BioMart columns with suffix (except identifiers)
setnames(all_kinases, "description", "description_biomart")
add_provenance("ensembl_gene_id", "biomart", "Ensembl gene identifier from BioMart")
add_provenance("external_gene_name", "biomart", "Gene symbol from BioMart")
add_provenance("description_biomart", "biomart", "Gene description from BioMart")

# 2. Annotate kinase group/family from HGNC (human only)
if (species == "human") {
  cat("Fetching HGNC kinase group annotation...\n", file = stderr())
  url <- "https://rest.genenames.org/search/status:Approved+AND+group:kinase"
  headers <- httr::add_headers(Accept = "application/json")
  res <- httr::GET(url, headers)
  httr::stop_for_status(res)

  # Download and parse HGNC gene group DB table
  hgnc_url <- "https://www.genenames.org/cgi-bin/genegroup/download-all"
  cat("\nDownloading HGNC gene group DB table from:", hgnc_url, "\n", file = stderr())
  tmpfile <- tempfile(fileext = ".txt")
  download.file(hgnc_url, tmpfile, quiet = TRUE)
  hgnc_groups <- fread(tmpfile, sep = "\t", header = TRUE, quote = '"')

  # Filter for kinase-related groups
  kinases_hgnc <- hgnc_groups[grepl("kinase", `Group name`, ignore.case = TRUE)]
  cat("\nHGNC kinase-related groups:", nrow(kinases_hgnc), "genes\n", file = stderr())

  if (!"Approved symbol" %in% names(kinases_hgnc)) {
    stop("HGNC gene group table missing 'Approved symbol' column")
  }
  setnames(kinases_hgnc, "Approved symbol", "external_gene_name")

  keep_cols <- intersect(c("external_gene_name", "Group name", "Family name", "Subfamily name", "HGNC ID"), names(kinases_hgnc))
  kinases_hgnc <- kinases_hgnc[, ..keep_cols]

  if ("HGNC ID" %in% names(kinases_hgnc)) {
    kinases_hgnc$`HGNC ID` <- as.character(kinases_hgnc$`HGNC ID`)
  }

  # Rename HGNC columns with suffix before merge
  hgnc_rename <- c(
    "Group name" = "Group_hgnc",
    "Family name" = "Family_hgnc",
    "Subfamily name" = "SubFamily_hgnc",
    "HGNC ID" = "HGNC_ID_hgnc"
  )
  for (old_name in names(hgnc_rename)) {
    if (old_name %in% names(kinases_hgnc)) {
      setnames(kinases_hgnc, old_name, hgnc_rename[[old_name]])
    }
  }

  all_kinases <- merge(all_kinases, kinases_hgnc, by = "external_gene_name", all.x = TRUE)

  # Track HGNC provenance
  add_provenance("Group_hgnc", "hgnc", "Kinase group from HGNC gene groups database")
  add_provenance("Family_hgnc", "hgnc", "Kinase family from HGNC gene groups database")
  add_provenance("SubFamily_hgnc", "hgnc", "Kinase subfamily from HGNC gene groups database")
  add_provenance("HGNC_ID_hgnc", "hgnc", "HGNC identifier from HGNC gene groups database")

  # Annotate 'Protein' column (derived from HGNC group membership)
  all_kinases$Protein_hgnc <- ifelse(!is.na(all_kinases$Group_hgnc) & all_kinases$Group_hgnc != "", "Y", "N")
  add_provenance("Protein_hgnc", "hgnc", "Protein kinase flag (Y if has HGNC group annotation)")
}

# 3. Annotate metabolic and lipid kinases using KEGG
if (species == "mouse") {
  # Extract MGI ID from description if present
  mgi_pat <- "MGI:[0-9]+"
  mgi_match <- regexpr(mgi_pat, all_kinases$description_biomart)
  all_kinases$MGI_ID_biomart <- ifelse(mgi_match > 0, regmatches(all_kinases$description_biomart, mgi_match), NA_character_)
  all_kinases$description_biomart <- gsub("\\s*\\[[^]]*]", "", all_kinases$description_biomart)
  add_provenance("MGI_ID_biomart", "biomart", "MGI identifier extracted from BioMart description")
}
if (species == "human") {
  all_kinases$description_biomart <- gsub("\\s*\\[[^]]*]", "", all_kinases$description_biomart)
}

cat("Annotating metabolic and lipid kinases using KEGG...\n", file = stderr())

# Metabolic pathways
metabolic_pathways <- keggList("pathway", kegg_code)
metabolic_ids <- names(metabolic_pathways)[grepl("Metabolic pathways|metabolism", metabolic_pathways, ignore.case = TRUE)]
metabolic_genes <- unique(unlist(lapply(metabolic_ids, function(pid) {
  kegg_genes <- keggGet(pid)[[1]]$GENE
  if (is.null(kegg_genes)) return(NULL)
  gsub(" .*", "", kegg_genes)
})))

# Lipid pathways
lipid_pathways <- keggList("pathway", kegg_code)
lipid_ids <- names(lipid_pathways)[grepl("lipid", lipid_pathways, ignore.case = TRUE)]
lipid_genes <- unique(unlist(lapply(lipid_ids, function(pid) {
  kegg_genes <- keggGet(pid)[[1]]$GENE
  if (is.null(kegg_genes)) return(NULL)
  gsub(" .*", "", kegg_genes)
})))

# Map KEGG Entrez IDs to gene symbols
entrez_to_symbol <- AnnotationDbi::select(orgdb, keys = unique(c(metabolic_genes, lipid_genes)),
                                          keytype = "ENTREZID", columns = c("SYMBOL"))
metabolic_syms <- unique(na.omit(entrez_to_symbol$SYMBOL[entrez_to_symbol$ENTREZID %in% metabolic_genes]))
lipid_syms <- unique(na.omit(entrez_to_symbol$SYMBOL[entrez_to_symbol$ENTREZID %in% lipid_genes]))

# Integrate Mammalian Metabolic Enzyme Database (if available)
metabolic_csv <- input_path("Mammalian_Metabolic_Final.csv")
metabolic_verified_syms <- character(0)
if (file.exists(metabolic_csv)) {
  metabolic_dt <- tryCatch(fread(metabolic_csv), error = function(e) NULL)
  if (!is.null(metabolic_dt) && "Gene Symbol" %in% names(metabolic_dt)) {
    metabolic_verified_syms <- toupper(trimws(metabolic_dt[["Gene Symbol"]]))
    cat(sprintf("Loaded %d metabolic enzyme gene symbols from Mammalian Metabolic Enzyme Database.\n",
                length(metabolic_verified_syms)), file = stderr())
  }
}

# Add Metabolic and Lipid columns with KEGG suffix
all_kinases$Metabolic_kegg <- ifelse(
  toupper(all_kinases$external_gene_name) %in% metabolic_verified_syms,
  "Verified",
  ifelse(all_kinases$external_gene_name %in% metabolic_syms, "Y", "N")
)
all_kinases$Lipid_kegg <- ifelse(all_kinases$external_gene_name %in% lipid_syms, "Y", "N")
add_provenance("Metabolic_kegg", "kegg", "Metabolic pathway membership from KEGG (Verified = Mammalian Metabolic DB)")
add_provenance("Lipid_kegg", "kegg", "Lipid pathway membership from KEGG")

# 4. Merge Manning and KinHub Group/Family for cross-source validation
cat("Merging Manning and KinHub classifications for validation...\n", file = stderr())

# Load Manning table
manning_dt <- fread(manning_file)
if ("Group" %in% names(manning_dt)) {
  manning_cols <- manning_dt[, .(Name, Group, Family, Subfamily)]
  setnames(manning_cols, c("Name", "Group", "Family", "Subfamily"),
           c("manning_name", "Group_manning", "Family_manning", "SubFamily_manning"))
  # Match by gene symbol (Manning uses human symbols)
  all_kinases <- merge(all_kinases, manning_cols,
                       by.x = "external_gene_name", by.y = "manning_name",
                       all.x = TRUE)
  add_provenance("Group_manning", "manning", "Kinase group from Manning et al. 2002 Table S1")
  add_provenance("Family_manning", "manning", "Kinase family from Manning et al. 2002 Table S1")
  add_provenance("SubFamily_manning", "manning", "Kinase subfamily from Manning et al. 2002 Table S1")
  cat(sprintf("  Manning matched: %d kinases\n", sum(!is.na(all_kinases$Group_manning))), file = stderr())
}

# Load KinHub reconciliation if available
kinhub_file <- input_path("manning_kinhub_reconciliation.csv")
if (file.exists(kinhub_file)) {
  kinhub_dt <- fread(kinhub_file)
  if ("KinHub_Group" %in% names(kinhub_dt)) {
    kinhub_cols <- kinhub_dt[, .(KinHub_HGNC, KinHub_Group, KinHub_Family)]
    setnames(kinhub_cols, c("KinHub_HGNC", "KinHub_Group", "KinHub_Family"),
             c("hgnc_symbol", "Group_kinhub", "Family_kinhub"))
    all_kinases <- merge(all_kinases, kinhub_cols,
                         by.x = "external_gene_name", by.y = "hgnc_symbol",
                         all.x = TRUE)
    add_provenance("Group_kinhub", "kinhub", "Kinase group from KinHub database")
    add_provenance("Family_kinhub", "kinhub", "Kinase family from KinHub database")
    cat(sprintf("  KinHub matched: %d kinases\n", sum(!is.na(all_kinases$Group_kinhub))), file = stderr())
  }
}

# 4b. Cross-source validation of Group/Family
cat("\n=== Cross-Source Group/Family Validation ===\n", file = stderr())

# Compare HGNC vs Manning
if ("Group_hgnc" %in% names(all_kinases) && "Group_manning" %in% names(all_kinases)) {
  both_have_group <- !is.na(all_kinases$Group_hgnc) & !is.na(all_kinases$Group_manning)
  if (sum(both_have_group) > 0) {
    group_match <- toupper(all_kinases$Group_hgnc) == toupper(all_kinases$Group_manning)
    group_match[!both_have_group] <- NA
    n_match <- sum(group_match, na.rm = TRUE)
    n_mismatch <- sum(!group_match, na.rm = TRUE)
    cat(sprintf("  HGNC vs Manning Group: %d match, %d mismatch\n", n_match, n_mismatch), file = stderr())

    # Flag mismatches
    all_kinases$Group_hgnc_manning_match <- ifelse(both_have_group, group_match, NA)
    add_provenance("Group_hgnc_manning_match", "validation", "TRUE if HGNC and Manning group agree")
  }
}

# Compare HGNC vs KinHub
if ("Group_hgnc" %in% names(all_kinases) && "Group_kinhub" %in% names(all_kinases)) {
  both_have_group <- !is.na(all_kinases$Group_hgnc) & !is.na(all_kinases$Group_kinhub)
  if (sum(both_have_group) > 0) {
    group_match <- toupper(all_kinases$Group_hgnc) == toupper(all_kinases$Group_kinhub)
    group_match[!both_have_group] <- NA
    n_match <- sum(group_match, na.rm = TRUE)
    n_mismatch <- sum(!group_match, na.rm = TRUE)
    cat(sprintf("  HGNC vs KinHub Group: %d match, %d mismatch\n", n_match, n_mismatch), file = stderr())

    all_kinases$Group_hgnc_kinhub_match <- ifelse(both_have_group, group_match, NA)
    add_provenance("Group_hgnc_kinhub_match", "validation", "TRUE if HGNC and KinHub group agree")
  }
}

# Compare Manning vs KinHub
if ("Group_manning" %in% names(all_kinases) && "Group_kinhub" %in% names(all_kinases)) {
  both_have_group <- !is.na(all_kinases$Group_manning) & !is.na(all_kinases$Group_kinhub)
  if (sum(both_have_group) > 0) {
    group_match <- toupper(all_kinases$Group_manning) == toupper(all_kinases$Group_kinhub)
    group_match[!both_have_group] <- NA
    n_match <- sum(group_match, na.rm = TRUE)
    n_mismatch <- sum(!group_match, na.rm = TRUE)
    cat(sprintf("  Manning vs KinHub Group: %d match, %d mismatch\n", n_match, n_mismatch), file = stderr())

    all_kinases$Group_manning_kinhub_match <- ifelse(both_have_group, group_match, NA)
    add_provenance("Group_manning_kinhub_match", "validation", "TRUE if Manning and KinHub group agree")
  }
}

# Create consensus Group column (prefer HGNC > Manning > KinHub)
all_kinases$Group_consensus <- ifelse(
  !is.na(all_kinases$Group_hgnc) & all_kinases$Group_hgnc != "",
  all_kinases$Group_hgnc,
  ifelse(!is.na(all_kinases$Group_manning) & all_kinases$Group_manning != "",
         all_kinases$Group_manning,
         all_kinases$Group_kinhub)
)
add_provenance("Group_consensus", "derived", "Consensus group: HGNC > Manning > KinHub priority")

cat("\n", file = stderr())

# 4c. Clean up and filter
if (species == "mouse") {
  before_n <- nrow(all_kinases)
  all_kinases <- all_kinases[!grepl("^Gm[0-9]+$", external_gene_name, ignore.case = TRUE) &
                             !grepl("Rik$", external_gene_name, ignore.case = TRUE)]
  after_n <- nrow(all_kinases)
  cat(sprintf("Filtered out Gm/Rik genes: %d removed\n", before_n - after_n), file = stderr())
}

if (species == "human") {
  if (!"Group_hgnc" %in% names(all_kinases)) all_kinases$Group_hgnc <- NA_character_
}

# Exclude rows with missing or blank gene symbol
before_n <- nrow(all_kinases)
all_kinases <- all_kinases[!is.na(external_gene_name) & external_gene_name != ""]
after_n <- nrow(all_kinases)
cat(sprintf("Filtered out missing/blank gene symbols: %d removed\n", before_n - after_n), file = stderr())

# 5. Write output to canonical outputs directory
date_tag <- format(Sys.time(), "%y%m%d")
ensure_dir(paths$output_dir)

outfile <- output_path(paste0("kinases_", species, "_annotated.csv"))
fwrite(all_kinases, outfile)
check_file_nonempty(outfile, paste0("Output file missing or empty: ", outfile))

# Write MD5 checksum
if (requireNamespace("tools", quietly = TRUE)) {
  md5 <- tools::md5sum(outfile)
  md5file <- paste0(outfile, ".md5")
  cat(sprintf("%s  %s\n", unname(md5), basename(outfile)), file = md5file)
}

# Write column provenance YAML
provenance_file <- output_path(paste0("kinases_", species, "_annotated.provenance.yaml"))
if (requireNamespace("yaml", quietly = TRUE)) {
  # Build provenance structure
  provenance_yaml <- list(
    file = basename(outfile),
    generated = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    species = species,
    sources = list(
      biomart = list(
        name = "Ensembl BioMart",
        dataset = dataset,
        url = "https://www.ensembl.org/biomart"
      ),
      hgnc = list(
        name = "HGNC Gene Groups Database",
        url = "https://www.genenames.org/cgi-bin/genegroup/download-all"
      ),
      kegg = list(
        name = "KEGG Pathway Database",
        organism = kegg_code,
        url = "https://www.kegg.jp/"
      )
    ),
    columns = column_provenance
  )
  yaml::write_yaml(provenance_yaml, provenance_file)
  cat(sprintf("✓ Column provenance: %s\n", provenance_file), file = stderr())
} else {
  cat("Note: Install 'yaml' package for provenance YAML output\n", file = stderr())
}

cat(sprintf("Done. Output: %s (%d kinases)\n", outfile, nrow(all_kinases)), file = stderr())
