#!/usr/bin/env Rscript
# build_kinase_list.R
# Author: C. Ryan Miller, MD, PhD (rmiller@uab.edu)
# Created: 2025-12-31
# Purpose: Generate curated kinase list from HGNC kinase gene group (primary source)
#          Validates against Manning Table S1 and KinHub
#          All columns track provenance via suffix convention
#
# Usage: Rscript build_kinase_list.R [--skip-mouse]
#
# Outputs:
#   - kinases_human_curated.csv: Main kinase list with all annotations
#   - kinases_human_curated.provenance.yaml: Column provenance metadata

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
  library(httr)
  library(jsonlite)
})

# Load centralized config
.script_dir <- (function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("--file=", "", file_arg[1]))))
  }
  return(normalizePath("."))
})()
source(file.path(.script_dir, "lib", "config_loader.R"))

# Parse command-line arguments
option_list <- list(
  make_option(c("--skip-mouse"), action = "store_true", default = FALSE,
              help = "Skip mouse ortholog mapping (faster, no BioMart queries)")
)
opt <- parse_args(OptionParser(option_list = option_list))

cat("=== Kinase List Generator ===\n")
cat("Primary source: HGNC Kinase Gene Group\n")
cat("Validation sources: Manning Table S1, KinHub\n\n")

# Initialize column provenance tracker
column_provenance <- list()
add_provenance <- function(col, source, description) {
  column_provenance[[col]] <<- list(
    source = source,
    description = description,
    added = format(Sys.time(), "%Y-%m-%d")
  )
}

# ============================================================================
# STEP 1: Load HGNC Gene Groups and Filter for Kinases (Primary Source)
# ============================================================================
cat("Step 1: Loading HGNC gene groups file...\n")

# Check for local HGNC gene groups file, download if missing
hgnc_groups_file <- input_path("hgnc_gene_groups.csv")
if (!file.exists(hgnc_groups_file)) {
  cat("  Downloading HGNC gene groups from genenames.org...\n")
  download.file(
    "https://www.genenames.org/cgi-bin/genegroup/download-all/csv",
    hgnc_groups_file,
    quiet = TRUE
  )
}

check_file_nonempty(hgnc_groups_file, "HGNC gene groups file not found")

# Load HGNC gene groups
hgnc_all <- fread(hgnc_groups_file, header = TRUE, sep = "\t")
cat("  Loaded", nrow(hgnc_all), "total gene-group entries\n")

# Filter for groups containing "kinase" (case-insensitive)
# Also include "EPH receptors" which are receptor tyrosine kinases but named without "kinase"
kinase_groups <- hgnc_all[
  grepl("kinase", `Group name`, ignore.case = TRUE) |
  `Group name` == "EPH receptors"
]
cat("  Found", nrow(kinase_groups), "entries in kinase-related groups\n")

# Show group distribution
n_groups <- length(unique(kinase_groups$`Group name`))
cat("  Spanning", n_groups, "HGNC kinase groups\n")

# Select and rename columns
kinases <- kinase_groups[, .(
  HGNC_symbol = `Approved symbol`,
  HGNC_ID = `HGNC ID`,
  Name_hgnc = `Approved name`,
  Group_hgnc = `Group name`,
  Group_ID_hgnc = `Group ID`,
  Status_hgnc = Status,
  Locus_type_hgnc = `Locus type`
)]

# Remove duplicates (genes can be in multiple groups - keep first occurrence)
kinases <- kinases[!duplicated(HGNC_symbol)]
cat("  Unique kinases:", nrow(kinases), "\n")

# Add kinase_id primary key (K0001 format, ordered by HGNC_symbol)
kinases <- kinases[order(HGNC_symbol)]
kinases[, kinase_id := sprintf("K%04d", seq_len(.N))]
cat("  Assigned kinase_id: K0001 to K", sprintf("%04d", nrow(kinases)), "\n", sep = "")

# Derive substrate classification from HGNC group (not mutually exclusive)
# A kinase can phosphorylate multiple substrate types
# Also flag non-kinase entries that appear in kinase groups

# First, identify non-kinase entries (phosphatases, regulators, anchoring proteins)
# Note: PFKFB (6-phosphofructo-2-kinase/fructose-2,6-biphosphatase) are bifunctional
# enzymes with BOTH kinase and phosphatase activity - they are TRUE kinases
# MAP kinase phosphatases (MKPs/DUSPs) are phosphatases despite having "kinase" in name
kinases[, `:=`(
  Is_phosphatase = ifelse(
    (grepl("phosphatase", Group_hgnc, ignore.case = TRUE) &
     !grepl("kinase", Group_hgnc, ignore.case = TRUE)) |
    grepl("MAP kinase phosphatase", Group_hgnc, ignore.case = TRUE),  # MKPs are phosphatases
    "Y", "N"
  ),
  Is_kinase_regulator = ifelse(
    grepl("MOB kinase activator|kinase activator|kinase inhibitor", Group_hgnc, ignore.case = TRUE),
    "Y", "N"
  ),
  Is_anchoring_protein = ifelse(
    grepl("anchoring", Group_hgnc, ignore.case = TRUE),
    "Y", "N"
  ),
  Is_pseudokinase = ifelse(
    grepl("pseudokinase|PEAK family", Group_hgnc, ignore.case = TRUE),
    "Y", "N"
  )
)]

# Initial substrate classification based on HGNC group names only
# (Manning-based refinement happens after Manning data is loaded in Step 3)
kinases[, `:=`(
  # Protein substrates: protein kinases, receptor kinases, most "kinase" groups
  Substrate_protein = ifelse(
    Is_phosphatase == "N" & Is_kinase_regulator == "N" & Is_anchoring_protein == "N" & (
      grepl("protein kinase|tyrosine kinase|serine.threonine|receptor.kinase|CDK|MAPK|kinome|Src|JAK|NEK|aurora|polo|casein|phosphorylase kinase|EPH receptor|NIMA|histidine kinase|WNK|cyclin.dependent|PKC|PKA",
            Group_hgnc, ignore.case = TRUE) |
      # Generic protein kinase detection: has "kinase" but not substrate-specific
      (grepl("kinase", Group_hgnc, ignore.case = TRUE) &
       !grepl("lipid|diacylglycerol|sphingosine|phosphatidylinositol|PI3K|PI4K|adenylate|guanylate|nucleoside|thymidine|uridine|hexokinase|pyruvate|phosphofructo|creatine|choline|ethanolamine|glycerol|ribose|ribulose|pantothenate|mevalonate|phosphatase|activator|anchoring|membrane associated guanylate",
              Group_hgnc, ignore.case = TRUE))
    ),
    "Y", "N"
  ),

  # Lipid substrates: lipid kinases, PI kinases, diacylglycerol kinases, etc.
  Substrate_lipid = ifelse(
    Is_phosphatase == "N" & Is_kinase_regulator == "N" & Is_anchoring_protein == "N" &
    grepl("lipid|diacylglycerol|sphingosine|phosphatidylinositol|PI3K|PI4K|PI5K|choline|ethanolamine|ceramide",
          Group_hgnc, ignore.case = TRUE),
    "Y", "N"
  ),

  # Nucleotide substrates: adenylate, guanylate, nucleoside kinases
  Substrate_nucleotide = ifelse(
    Is_phosphatase == "N" & Is_kinase_regulator == "N" & Is_anchoring_protein == "N" &
    grepl("adenylate|guanylate|nucleoside|thymidine|uridine|cytidine|deoxycytidine",
          Group_hgnc, ignore.case = TRUE) &
    !grepl("membrane associated guanylate", Group_hgnc, ignore.case = TRUE),  # MAGUK are pseudokinases
    "Y", "N"
  ),

  # Carbohydrate substrates: hexokinases, pyruvate kinases, phosphofructokinases, glycerol kinases, ribose kinases
  # Note: Glycerol kinases phosphorylate glycerol (a sugar alcohol) - classified as carbohydrate
  # Diacylglycerol kinases phosphorylate diacylglycerol (a lipid) - classified as lipid (not here)
  Substrate_carbohydrate = ifelse(
    Is_phosphatase == "N" & Is_kinase_regulator == "N" & Is_anchoring_protein == "N" &
    grepl("hexokinase|pyruvate|phosphofructo|creatine|glycerol kinase|ribose|ribulose|fructokinase|galactokinase|glucokinase|phosphoglycerate|ketohexokinase",
          Group_hgnc, ignore.case = TRUE) &
    !grepl("diacylglycerol", Group_hgnc, ignore.case = TRUE),  # DGK are lipid kinases
    "Y", "N"
  ),

  # Other substrates: pantothenate (vitamin B5), mevalonate, etc.
  Substrate_other = ifelse(
    Is_phosphatase == "N" & Is_kinase_regulator == "N" & Is_anchoring_protein == "N" &
    grepl("pantothenate|mevalonate|shikimate|homoserine|nicotinamide|riboflavin",
          Group_hgnc, ignore.case = TRUE),
    "Y", "N"
  )
)]

# MAGUK proteins: membrane associated guanylate kinases have inactive GK domain
kinases[grepl("membrane associated guanylate kinase|MAGUK", Group_hgnc, ignore.case = TRUE),
        Is_pseudokinase := "Y"]

# Report initial substrate distribution (before Manning refinement)
n_protein <- sum(kinases$Substrate_protein == "Y")
n_lipid <- sum(kinases$Substrate_lipid == "Y")
n_nucleotide <- sum(kinases$Substrate_nucleotide == "Y")
n_carbohydrate <- sum(kinases$Substrate_carbohydrate == "Y")
n_other <- sum(kinases$Substrate_other == "Y")
n_phosphatase <- sum(kinases$Is_phosphatase == "Y")
n_regulator <- sum(kinases$Is_kinase_regulator == "Y")
n_anchoring <- sum(kinases$Is_anchoring_protein == "Y")
n_pseudokinase <- sum(kinases$Is_pseudokinase == "Y")
cat("  Substrate classification:\n")
cat("    Protein:", n_protein, "\n")
cat("    Lipid:", n_lipid, "\n")
cat("    Nucleotide:", n_nucleotide, "\n")
cat("    Carbohydrate:", n_carbohydrate, "\n")
cat("    Other (vitamin/small molecule):", n_other, "\n")
cat("  Non-kinase entries in kinase groups:\n")
cat("    Phosphatases:", n_phosphatase, "\n")
cat("    Kinase regulators:", n_regulator, "\n")
cat("    Anchoring proteins:", n_anchoring, "\n")
cat("    Pseudokinases:", n_pseudokinase, "\n")

add_provenance("kinase_id", "derived", "Primary key in K#### format (ordered by HGNC_symbol)")
add_provenance("Substrate_protein", "derived", "Y if phosphorylates protein substrates")
add_provenance("Substrate_lipid", "derived", "Y if phosphorylates lipid substrates")
add_provenance("Substrate_nucleotide", "derived", "Y if phosphorylates nucleotide substrates")
add_provenance("Substrate_carbohydrate", "derived", "Y if phosphorylates carbohydrate substrates")
add_provenance("Substrate_other", "derived", "Y if phosphorylates other substrates (vitamins, small molecules)")
add_provenance("Is_phosphatase", "derived", "Y if entry is a phosphatase (not a kinase)")
add_provenance("Is_kinase_regulator", "derived", "Y if entry is a kinase regulator (not a kinase itself)")
add_provenance("Is_anchoring_protein", "derived", "Y if A-kinase anchoring protein (not a kinase itself)")
add_provenance("Is_pseudokinase", "derived", "Y if catalytically inactive kinase domain")
add_provenance("HGNC_symbol", "hgnc", "Authoritative human gene symbol from HGNC gene groups")
add_provenance("HGNC_ID", "hgnc", "HGNC database identifier")
add_provenance("Name_hgnc", "hgnc", "Full gene name from HGNC")
add_provenance("Group_hgnc", "hgnc", "HGNC gene group name (contains 'kinase')")
add_provenance("Group_ID_hgnc", "hgnc", "HGNC gene group numeric ID")
add_provenance("Status_hgnc", "hgnc", "HGNC gene approval status")
add_provenance("Locus_type_hgnc", "hgnc", "HGNC locus type classification")

# ============================================================================
# STEP 2: Load KinHub for Manning → HGNC mapping validation
# ============================================================================
cat("\nStep 2: Loading KinHub for Manning → HGNC mapping...\n")

kinhub_file <- input_path("251231_KinHub_Human_Kinases.csv")
kinhub_mapping <- NULL

if (file.exists(kinhub_file)) {
  kinhub_raw <- fread(kinhub_file)
  cat("  Loaded KinHub with", nrow(kinhub_raw), "entries\n")

  # Clean column names (remove BOM, non-breaking spaces)
  old_names <- names(kinhub_raw)
  new_names <- gsub("\u00A0", " ", old_names)
  new_names <- gsub("^\uFEFF", "", new_names)
  new_names <- trimws(new_names)
  setnames(kinhub_raw, old_names, new_names)

  # Extract Manning → HGNC mapping from KinHub
  hgnc_col <- names(kinhub_raw)[1]  # First column is HGNC Name
  manning_col <- grep("Manning", names(kinhub_raw), value = TRUE)[1]

  if (!is.null(manning_col)) {
    kinhub_mapping <- kinhub_raw[, .(
      HGNC_kinhub = get(hgnc_col),
      Manning_kinhub = get(manning_col),
      Group_kinhub = Group,
      Family_kinhub = Family,
      SubFamily_kinhub = SubFamily,
      UniProt_kinhub = UniprotID
    )]
    kinhub_mapping <- kinhub_mapping[!is.na(Manning_kinhub) & Manning_kinhub != ""]
    kinhub_mapping <- kinhub_mapping[!duplicated(HGNC_kinhub)]
    cat("  KinHub Manning → HGNC mappings:", nrow(kinhub_mapping), "\n")
  }
} else {
  cat("  Warning: KinHub file not found\n")
}

# ============================================================================
# STEP 3: Add Manning Table S1 Validation Columns (using KinHub mapping)
# ============================================================================
cat("\nStep 3: Adding Manning Table S1 validation columns...\n")

manning_file <- file.path(repo_root, cfg_get("input_files.manning",
  "curated/kinases/inputs/Manning_tableS1__251228.csv"))

if (file.exists(manning_file)) {
  manning <- fread(manning_file, header = TRUE)
  cat("  Loaded", nrow(manning), "entries from Manning Table S1\n")

  # Select columns with _manning suffix
  manning_cols <- manning[, .(
    manning_name = Name,
    kinase_id_manning = Accession,
    Group_manning = Group,
    Family_manning = Family,
    SubFamily_manning = Subfamily,
    Pseudogene_manning = `Pseudogene?`,
    Novelty_manning = Novelty
  )]

  # Strategy 1: Direct match (HGNC symbol == Manning name, case-insensitive)
  manning_cols[, manning_name_upper := toupper(manning_name)]
  kinases[, HGNC_symbol_upper := toupper(HGNC_symbol)]

  kinases <- merge(kinases, manning_cols,
                   by.x = "HGNC_symbol_upper", by.y = "manning_name_upper",
                   all.x = TRUE)
  kinases[, HGNC_symbol_upper := NULL]

  n_direct <- sum(!is.na(kinases$Group_manning))
  cat("  Direct HGNC=Manning match:", n_direct, "\n")

  # Strategy 2: Use KinHub Manning → HGNC mapping for unmatched
  if (!is.null(kinhub_mapping) && sum(is.na(kinases$Group_manning)) > 0) {
    # For kinases without Manning match, try KinHub mapping
    unmatched <- kinases[is.na(Group_manning), .(HGNC_symbol)]

    # Find Manning entries that map to these HGNC symbols via KinHub
    kinhub_matched <- merge(unmatched, kinhub_mapping,
                            by.x = "HGNC_symbol", by.y = "HGNC_kinhub")

    if (nrow(kinhub_matched) > 0) {
      # Get Manning data for these via KinHub's Manning name
      manning_via_kinhub <- merge(kinhub_matched[, .(HGNC_symbol, Manning_kinhub)],
                                  manning_cols,
                                  by.x = "Manning_kinhub", by.y = "manning_name")

      if (nrow(manning_via_kinhub) > 0) {
        manning_via_kinhub[, Manning_kinhub := NULL]
        manning_via_kinhub[, manning_name_upper := NULL]

        # Update kinases with KinHub-mapped Manning data
        for (col in setdiff(names(manning_via_kinhub), "HGNC_symbol")) {
          kinases[manning_via_kinhub, (col) := get(paste0("i.", col)), on = "HGNC_symbol"]
        }

        n_kinhub_mapped <- sum(!is.na(manning_via_kinhub$Group_manning))
        cat("  KinHub-mapped Manning match:", n_kinhub_mapped, "\n")
      }
    }
  }

  # Add InManning flag
  kinases[, InManning := ifelse(!is.na(Group_manning), "Y", "N")]

  n_in_manning <- sum(kinases$InManning == "Y")
  cat("  Total Manning matched:", n_in_manning, "/", nrow(kinases), "\n")

  # Manning-based protein kinase refinement: if a gene is in Manning, it's a protein kinase
  # (Manning only contains protein kinases, not lipid/carbohydrate/nucleotide kinases)
  kinases[InManning == "Y" & Substrate_protein == "N", Substrate_protein := "Y"]
  n_manning_refined <- sum(kinases$InManning == "Y" & kinases$Substrate_protein == "Y")
  cat("  Protein kinases after Manning refinement:", n_manning_refined, "\n")

  add_provenance("kinase_id_manning", "manning", "Manning accession (SK###)")
  add_provenance("manning_name", "manning", "Kinase name from Manning et al. 2002")
  add_provenance("Group_manning", "manning", "Top-level kinase classification (10 groups)")
  add_provenance("Family_manning", "manning", "Mid-level kinase classification")
  add_provenance("SubFamily_manning", "manning", "Optional lower-level classification")
  add_provenance("Pseudogene_manning", "manning", "Pseudogene status: Y=yes, N=no, R=related")
  add_provenance("Novelty_manning", "manning", "Annotation status: Known vs Novel")
  add_provenance("InManning", "validation", "Present in Manning Table S1 (Y/N)")
} else {
  cat("  Warning: Manning file not found, skipping validation columns\n")
  kinases[, `:=`(kinase_id_manning = NA_character_, manning_name = NA_character_,
                 Group_manning = NA_character_, Family_manning = NA_character_,
                 SubFamily_manning = NA_character_, Pseudogene_manning = NA_character_,
                 Novelty_manning = NA_character_, InManning = "N")]
}

# ============================================================================
# STEP 4: Add KinHub Validation Columns (reuse already-loaded data)
# ============================================================================
cat("\nStep 4: Adding KinHub validation columns...\n")

if (!is.null(kinhub_mapping)) {
  # Merge KinHub columns with kinases
  kinhub_cols <- kinhub_mapping[, .(
    HGNC_kinhub,
    Group_kinhub,
    Family_kinhub,
    SubFamily_kinhub,
    UniProt_kinhub
  )]

  kinases <- merge(kinases, kinhub_cols,
                   by.x = "HGNC_symbol", by.y = "HGNC_kinhub",
                   all.x = TRUE)

  # Add InKinHub flag
  kinases[, InKinHub := ifelse(!is.na(Group_kinhub), "Y", "N")]

  n_in_kinhub <- sum(kinases$InKinHub == "Y")
  cat("  Matched in KinHub:", n_in_kinhub, "/", nrow(kinases), "\n")
  cat("  With UniProt ID:", sum(!is.na(kinases$UniProt_kinhub) & kinases$UniProt_kinhub != ""), "\n")

  add_provenance("Group_kinhub", "kinhub", "Kinase group from KinHub database")
  add_provenance("Family_kinhub", "kinhub", "Kinase family from KinHub database")
  add_provenance("SubFamily_kinhub", "kinhub", "Kinase subfamily from KinHub database")
  add_provenance("UniProt_kinhub", "kinhub", "UniProt accession from KinHub")
  add_provenance("InKinHub", "validation", "Present in KinHub database (Y/N)")
} else {
  cat("  Warning: KinHub file not found, skipping validation columns\n")
  kinases[, `:=`(Group_kinhub = NA_character_, Family_kinhub = NA_character_,
                 SubFamily_kinhub = NA_character_, UniProt_kinhub = NA_character_,
                 InKinHub = "N")]
}

# ============================================================================
# STEP 5: Cross-Source Validation
# ============================================================================
cat("\nStep 5: Cross-source validation...\n")

# Compare Manning vs KinHub Group classifications
kinases[, Group_match_manning_kinhub := ifelse(
  !is.na(Group_manning) & !is.na(Group_kinhub),
  toupper(Group_manning) == toupper(Group_kinhub),
  NA
)]

# Compare Manning vs KinHub Family classifications
kinases[, Family_match_manning_kinhub := ifelse(
  !is.na(Family_manning) & !is.na(Family_kinhub),
  toupper(Family_manning) == toupper(Family_kinhub),
  NA
)]

# Report validation stats
n_group_match <- sum(kinases$Group_match_manning_kinhub == TRUE, na.rm = TRUE)
n_group_mismatch <- sum(kinases$Group_match_manning_kinhub == FALSE, na.rm = TRUE)
n_family_match <- sum(kinases$Family_match_manning_kinhub == TRUE, na.rm = TRUE)
n_family_mismatch <- sum(kinases$Family_match_manning_kinhub == FALSE, na.rm = TRUE)

cat("  Group validation (Manning vs KinHub):", n_group_match, "match,", n_group_mismatch, "mismatch\n")
cat("  Family validation (Manning vs KinHub):", n_family_match, "match,", n_family_mismatch, "mismatch\n")

# Create consensus columns (prefer Manning > KinHub > HGNC group)
# For non-protein kinases, use HGNC group as the classification
kinases[, Group_consensus := ifelse(
  !is.na(Group_manning), Group_manning,
  ifelse(!is.na(Group_kinhub), Group_kinhub, Group_hgnc)
)]
kinases[, Family_consensus := ifelse(
  !is.na(Family_manning), Family_manning,
  ifelse(!is.na(Family_kinhub), Family_kinhub, Group_hgnc)  # Use HGNC group as family fallback
)]

add_provenance("Group_match_manning_kinhub", "validation", "TRUE if Manning and KinHub groups agree")
add_provenance("Family_match_manning_kinhub", "validation", "TRUE if Manning and KinHub families agree")
add_provenance("Group_consensus", "derived", "Consensus group (Manning > KinHub > HGNC priority)")
add_provenance("Family_consensus", "derived", "Consensus family (Manning > KinHub > HGNC priority)")

# ============================================================================
# STEP 5b: Add KEGG Metabolic/Lipid Pathway Annotations
# ============================================================================
cat("\nStep 5b: Adding KEGG metabolic/lipid pathway annotations...\n")

if (requireNamespace("KEGGREST", quietly = TRUE)) {
  library(KEGGREST)

  tryCatch({
    # Get metabolic pathway genes
    metabolic_pathways <- keggList("pathway", "hsa")
    metabolic_ids <- names(metabolic_pathways)[grepl("Metabolic pathways|metabolism",
                                                      metabolic_pathways, ignore.case = TRUE)]

    # KEGG returns alternating: EntrezID, "SYMBOL; description", EntrezID, "SYMBOL; description"...
    # Extract even-indexed elements (the gene symbol+description), then parse symbol before ";"
    metabolic_genes <- unique(unlist(lapply(metabolic_ids, function(pid) {
      genes <- tryCatch(keggGet(pid)[[1]]$GENE, error = function(e) NULL)
      if (!is.null(genes) && length(genes) >= 2) {
        genes[seq(2, length(genes), by = 2)]  # Even indices have symbols
      }
    })))
    metabolic_syms <- gsub(";.*", "", metabolic_genes)
    cat("  KEGG metabolic pathway genes:", length(metabolic_syms), "\n")

    # Get lipid pathway genes
    lipid_ids <- names(metabolic_pathways)[grepl("lipid|fatty acid|sphingolipid|glycerolipid|phospholipid",
                                                  metabolic_pathways, ignore.case = TRUE)]
    lipid_genes <- unique(unlist(lapply(lipid_ids, function(pid) {
      genes <- tryCatch(keggGet(pid)[[1]]$GENE, error = function(e) NULL)
      if (!is.null(genes) && length(genes) >= 2) {
        genes[seq(2, length(genes), by = 2)]  # Even indices have symbols
      }
    })))
    lipid_syms <- gsub(";.*", "", lipid_genes)
    cat("  KEGG lipid pathway genes:", length(lipid_syms), "\n")

    # Check for Mammalian Metabolic Enzyme Database file for "Verified" status
    metabolic_verified_syms <- character(0)
    metabolic_csv <- input_path("mammalian_metabolic_enzymes.csv")
    if (file.exists(metabolic_csv)) {
      metabolic_dt <- tryCatch(fread(metabolic_csv), error = function(e) NULL)
      if (!is.null(metabolic_dt) && "Gene Symbol" %in% names(metabolic_dt)) {
        metabolic_verified_syms <- toupper(trimws(metabolic_dt[["Gene Symbol"]]))
        cat("  Mammalian Metabolic Enzyme DB genes:", length(metabolic_verified_syms), "\n")
      }
    }

    # Add Metabolic_kegg column
    kinases[, Metabolic_kegg := ifelse(
      toupper(HGNC_symbol) %in% metabolic_verified_syms,
      "Verified",
      ifelse(HGNC_symbol %in% metabolic_syms, "Y", "N")
    )]

    # Add Lipid_kegg column
    kinases[, Lipid_kegg := ifelse(HGNC_symbol %in% lipid_syms, "Y", "N")]

    n_metabolic <- sum(kinases$Metabolic_kegg != "N", na.rm = TRUE)
    n_lipid <- sum(kinases$Lipid_kegg == "Y", na.rm = TRUE)
    cat("  Kinases in metabolic pathways:", n_metabolic, "\n")
    cat("  Kinases in lipid pathways:", n_lipid, "\n")

    add_provenance("Metabolic_kegg", "kegg", "Metabolic pathway membership from KEGG (Verified = Mammalian Metabolic DB)")
    add_provenance("Lipid_kegg", "kegg", "Lipid pathway membership from KEGG")

  }, error = function(e) {
    cat("  Warning: KEGG query failed:", e$message, "\n")
    kinases[, `:=`(Metabolic_kegg = NA_character_, Lipid_kegg = NA_character_)]
  })
} else {
  cat("  Warning: KEGGREST package not available, skipping KEGG annotations\n")
  kinases[, `:=`(Metabolic_kegg = NA_character_, Lipid_kegg = NA_character_)]
}

# ============================================================================
# STEP 6: Add Mouse Ortholog Mapping
# ============================================================================
if (!opt$`skip-mouse`) {
  cat("\nStep 6: Mapping mouse orthologs via BioMart...\n")

  if (requireNamespace("biomaRt", quietly = TRUE)) {
    library(biomaRt)

    tryCatch({
      human_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
      hgnc_symbols <- unique(kinases$HGNC_symbol)

      orthologs <- getBM(
        attributes = c("external_gene_name", "mmusculus_homolog_associated_gene_name",
                       "mmusculus_homolog_ensembl_gene"),
        filters = "external_gene_name",
        values = hgnc_symbols,
        mart = human_mart
      )
      setDT(orthologs)

      orthologs <- orthologs[mmusculus_homolog_associated_gene_name != ""]
      orthologs <- orthologs[!duplicated(external_gene_name)]

      setnames(orthologs, c("external_gene_name", "mmusculus_homolog_associated_gene_name",
                            "mmusculus_homolog_ensembl_gene"),
               c("HGNC_symbol", "Mouse_Symbol", "Ensembl_Mouse"))

      kinases <- merge(kinases, orthologs, by = "HGNC_symbol", all.x = TRUE)

      n_mouse <- sum(!is.na(kinases$Mouse_Symbol))
      cat("  Mouse orthologs found:", n_mouse, "/", nrow(kinases), "\n")

      add_provenance("Mouse_Symbol", "biomart", "Mouse ortholog gene symbol from Ensembl BioMart")
      add_provenance("Ensembl_Mouse", "biomart", "Mouse Ensembl gene ID from BioMart")

    }, error = function(e) {
      cat("  Warning: BioMart query failed:", e$message, "\n")
      kinases[, `:=`(Mouse_Symbol = NA_character_, Ensembl_Mouse = NA_character_)]
    })
  } else {
    cat("  Warning: biomaRt package not available, skipping mouse mapping\n")
    kinases[, `:=`(Mouse_Symbol = NA_character_, Ensembl_Mouse = NA_character_)]
  }
} else {
  cat("\nStep 6: Skipped mouse ortholog mapping (--skip-mouse)\n")
  kinases[, `:=`(Mouse_Symbol = NA_character_, Ensembl_Mouse = NA_character_)]
}

# ============================================================================
# STEP 7: Reorder columns and write output
# ============================================================================
cat("\nStep 7: Writing output...\n")

# Define column order
col_order <- c(
  # Primary key and core identifiers
  "kinase_id", "HGNC_symbol", "HGNC_ID", "Name_hgnc",
  # Entry type flags (non-kinase entries in kinase groups)
  "Is_phosphatase", "Is_kinase_regulator", "Is_anchoring_protein", "Is_pseudokinase",
  # Substrate classification (not mutually exclusive, only for true kinases)
  "Substrate_protein", "Substrate_lipid", "Substrate_nucleotide",
  "Substrate_carbohydrate", "Substrate_other",
  # HGNC classification
  "Group_hgnc", "Group_ID_hgnc", "Status_hgnc", "Locus_type_hgnc",
  # Manning validation
  "kinase_id_manning", "manning_name", "Group_manning", "Family_manning", "SubFamily_manning",
  "Pseudogene_manning", "Novelty_manning", "InManning",
  # KinHub validation
  "Group_kinhub", "Family_kinhub", "SubFamily_kinhub", "UniProt_kinhub", "InKinHub",
  # Cross-source validation
  "Group_match_manning_kinhub", "Family_match_manning_kinhub",
  # Consensus
  "Group_consensus", "Family_consensus",
  # KEGG pathway annotations
  "Metabolic_kegg", "Lipid_kegg",
  # Mouse
  "Mouse_Symbol", "Ensembl_Mouse"
)

# Select columns that exist
existing_cols <- intersect(col_order, names(kinases))
kinases <- kinases[, ..existing_cols]

# Create output directory
ensure_dir(paths$output_dir)

# Write main output
output_file <- output_path("kinases_human_curated.csv")
fwrite(kinases, output_file)
cat("  ✓ Kinase list:", output_file, "\n")
cat("    (", nrow(kinases), "kinases,", ncol(kinases), "columns)\n")

# Write provenance YAML
if (requireNamespace("yaml", quietly = TRUE)) {
  provenance_file <- output_path("kinases_human_curated.provenance.yaml")

  provenance_yaml <- list(
    file = basename(output_file),
    generated = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    primary_source = "HGNC Kinase Gene Group",
    validation_sources = c("Manning Table S1 (2002)", "KinHub"),
    total_kinases = nrow(kinases),
    sources = list(
      hgnc = list(
        name = "HGNC Gene Names",
        url = "https://www.genenames.org/",
        role = "Primary source"
      ),
      manning = list(
        name = "Manning et al. 2002",
        file = basename(manning_file),
        url = "https://doi.org/10.1126/science.1075762",
        role = "Validation source"
      ),
      kinhub = list(
        name = "KinHub Database",
        file = "251231_KinHub_Human_Kinases.csv",
        url = "http://www.kinhub.org/",
        role = "Validation source"
      ),
      kegg = list(
        name = "KEGG Pathway Database",
        url = "https://www.kegg.jp/",
        role = "Metabolic and lipid pathway annotations"
      ),
      biomart = list(
        name = "Ensembl BioMart",
        url = "https://www.ensembl.org/biomart",
        role = "Mouse ortholog mapping"
      )
    ),
    columns = column_provenance
  )

  yaml::write_yaml(provenance_yaml, provenance_file)
  cat("  ✓ Provenance:", provenance_file, "\n")
}

# Summary statistics
cat("\n=== Summary ===\n")
cat("Total entries (HGNC kinase groups):", nrow(kinases), "\n")

# Count true active kinases (exclude non-kinase entries and pseudokinases)
n_true_kinases <- nrow(kinases[Is_phosphatase == "N" & Is_kinase_regulator == "N" &
                                Is_anchoring_protein == "N" & Is_pseudokinase == "N"])
n_including_pseudo <- nrow(kinases[Is_phosphatase == "N" & Is_kinase_regulator == "N" & Is_anchoring_protein == "N"])
cat("True active kinases:", n_true_kinases, "\n")
cat("Including pseudokinases:", n_including_pseudo, "\n")

cat("\nNon-kinase entries in kinase-related groups:\n")
cat("  Phosphatases:", sum(kinases$Is_phosphatase == "Y", na.rm = TRUE), "\n")
cat("  Kinase regulators:", sum(kinases$Is_kinase_regulator == "Y", na.rm = TRUE), "\n")
cat("  Anchoring proteins:", sum(kinases$Is_anchoring_protein == "Y", na.rm = TRUE), "\n")
cat("  Pseudokinases:", sum(kinases$Is_pseudokinase == "Y", na.rm = TRUE), "\n")

cat("\nSubstrate classification (not mutually exclusive):\n")
cat("  Protein substrate:", sum(kinases$Substrate_protein == "Y", na.rm = TRUE), "\n")
cat("  Lipid substrate:", sum(kinases$Substrate_lipid == "Y", na.rm = TRUE), "\n")
cat("  Nucleotide substrate:", sum(kinases$Substrate_nucleotide == "Y", na.rm = TRUE), "\n")
cat("  Carbohydrate substrate:", sum(kinases$Substrate_carbohydrate == "Y", na.rm = TRUE), "\n")
cat("  Other substrate:", sum(kinases$Substrate_other == "Y", na.rm = TRUE), "\n")

cat("\nValidation sources:\n")
cat("  In Manning Table S1:", sum(kinases$InManning == "Y", na.rm = TRUE), "\n")
cat("  In KinHub:", sum(kinases$InKinHub == "Y", na.rm = TRUE), "\n")
cat("  In metabolic pathways:", sum(kinases$Metabolic_kegg != "N", na.rm = TRUE), "\n")
cat("  In lipid pathways:", sum(kinases$Lipid_kegg == "Y", na.rm = TRUE), "\n")
cat("  With mouse ortholog:", sum(!is.na(kinases$Mouse_Symbol)), "\n")

# Group distribution (consensus)
cat("\nGroup distribution (consensus):\n")
group_dist <- kinases[!is.na(Group_consensus), .N, by = Group_consensus][order(-N)]
print(group_dist)

n_no_group <- sum(is.na(kinases$Group_consensus))
if (n_no_group > 0) {
  cat("(", n_no_group, "kinases without group classification)\n")
}

cat("\n✓ Kinase list generation complete.\n")
