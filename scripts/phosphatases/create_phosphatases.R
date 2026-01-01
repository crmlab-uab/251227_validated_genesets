#!/usr/bin/env Rscript
# Create curated human phosphatase geneset from HGNC gene groups
# Author: Claude Code
# Created: 2026-01-01
#
# Usage: Rscript create_phosphatases.R [--skip-mouse]
#
# Retrieves ALL phosphatases from HGNC, then annotates by substrate type

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(optparse)
})

# Parse command-line arguments
option_list <- list(
  make_option(c("--skip-mouse"), action = "store_true", default = FALSE,
              help = "Skip mouse ortholog mapping (faster, no BioMart queries)")
)
opt <- parse_args(OptionParser(option_list = option_list))

# Paths
input_file <- "curated/phosphatases/inputs/hgnc_gene_groups.csv"
output_file <- "curated/phosphatases/outputs/phosphatases_human_curated.csv"

# Read HGNC gene groups
hgnc <- read_tsv(input_file, show_col_types = FALSE)

# Define ALL phosphatase groups from HGNC
# Approach: include any group with "phosphatase" in name
phosphatase_groups <- c(

  # ===== PROTEIN PHOSPHATASES (catalytic) =====
  # Ser/Thr phosphatases
  "Protein phosphatase catalytic subunits",
  "Protein phosphatases",
  "Protein phosphatases, Mg2+/Mn2+ dependent",
  "Serine/threonine phosphatases",
  "Calcineurin subunits",


  # Tyrosine phosphatases
  "Protein tyrosine phosphatases receptor type",
  "Protein tyrosine phosphatases non-receptor type",
  "Protein tyrosine phosphatase 4A family",
  "LAR protein receptor tyrosine phosphatase family",

  # Dual specificity phosphatases
  "Atypical dual specificity phosphatases",
  "MAP kinase phosphatases",
  "CDC14 phosphatases",
  "Slingshot protein phosphatases",
  "CTD family phosphatases",
  "PTEN protein phosphatases",
  "Class II Cys-based phosphatases",
  "Class III Cys-based CDC25 phosphatases",

  # HAD phosphatases
  "HAD Asp-based protein phosphatases",
  "HAD Asp-based non-protein phosphatases",

  # Special protein phosphatases
  "EYA transcriptional coactivator and phosphatases",

  # ===== PROTEIN PHOSPHATASE REGULATORY SUBUNITS =====
  "Protein phosphatase 1 regulatory subunits",
  "Protein phosphatase 2 regulatory subunits",
  "Protein phosphatase 2 modulatory subunits",
  "Protein phosphatase 2 scaffold subunits",
  "Protein phosphatase 3 regulatory subunits",
  "Protein phosphatase 4 regulatory subunits",
  "Protein phosphatase 6 regulatory subunits",
  "Myosin phosphatase targeting family",
  "Phosphatase and actin regulators",

  # ===== LIPID PHOSPHATASES =====
  "Lipid phosphatases",
  "Phosphoinositide phosphatases",
  "Phospholipid phosphatases",
  "Phospholipid phosphatase related",
  "Sphingosine-1-phosphate phosphatases",

  # ===== CARBOHYDRATE/METABOLIC PHOSPHATASES =====
  "Acid phosphatases",
  "Alkaline phosphatases",
  "Bisphosphoglycerate phosphatases",
  "Fructose-1,6-bisphosphatases",
  "Glucose 6-phosphatases, catalytic",
  "Sugar phosphatases",
  "Phosphoglycerate mutases",  # have phosphatase activity

  # ===== NUCLEOTIDE PHOSPHATASES =====
  "Ectonucleotide pyrophosphatase/phosphodiesterase family",

  # ===== BIFUNCTIONAL (kinase + phosphatase) =====
  "6-phosphofructo-2-kinase/fructose-2,6-biphosphatase family"
)

cat("=== Phosphatase List Generator ===\n")
cat("Primary source: HGNC Gene Groups\n\n")

# Filter to phosphatase groups
phos <- hgnc %>%
  filter(`Group name` %in% phosphatase_groups) %>%
  filter(Status == "Approved") %>%
  filter(`Locus type` == "gene with protein product")

cat("Found", nrow(phos), "gene-group entries in", length(phosphatase_groups), "phosphatase groups\n")

# Get unique genes (some appear in multiple groups)
phos_unique <- phos %>%
  group_by(`Approved symbol`) %>%
  summarise(
    HGNC_ID = first(`HGNC ID`),
    Name_hgnc = first(`Approved name`),
    Status_hgnc = first(Status),
    Locus_type_hgnc = first(`Locus type`),
    Chromosome = first(Chromosome),
    NCBI_Gene_ID = first(`NCBI Gene ID`),
    Ensembl_ID = first(`Ensembl gene ID`),
    Group_hgnc = paste(unique(`Group name`), collapse = "; "),
    Group_ID_hgnc = paste(unique(`Group ID`), collapse = "; "),
    .groups = "drop"
  ) %>%
  rename(HGNC_symbol = `Approved symbol`)

cat("Unique phosphatases:", nrow(phos_unique), "\n")

# Classify phosphatases by substrate and type
phos_final <- phos_unique %>%
  mutate(
    # ===== SUBSTRATE CLASSIFICATION =====
    # Protein substrates
    Substrate_protein = case_when(
      grepl("Protein phosphatase|tyrosine phosphatase|Serine/threonine|Calcineurin|dual specificity|CDC14|CDC25|Slingshot|CTD|PTEN|EYA|MAP kinase|HAD.*protein|Class II Cys|Class III Cys|LAR|4A family",
            Group_hgnc, ignore.case = TRUE) ~ "Y",
      TRUE ~ "N"
    ),

    # Lipid substrates
    Substrate_lipid = case_when(
      grepl("Lipid phosphatase|Phosphoinositide|Phospholipid|Sphingosine|PTEN",
            Group_hgnc, ignore.case = TRUE) ~ "Y",
      TRUE ~ "N"
    ),

    # Nucleotide substrates
    Substrate_nucleotide = case_when(
      grepl("nucleotide|pyrophosphatase", Group_hgnc, ignore.case = TRUE) ~ "Y",
      TRUE ~ "N"
    ),

    # Carbohydrate/metabolic substrates
    Substrate_carbohydrate = case_when(
      grepl("Acid phosphatase|Alkaline phosphatase|Bisphosphoglycerate|Fructose|Glucose|Sugar|Phosphoglycerate|phosphofructo",
            Group_hgnc, ignore.case = TRUE) ~ "Y",
      TRUE ~ "N"
    ),

    # Other substrates (HAD non-protein, etc.)
    Substrate_other = case_when(
      grepl("HAD.*non-protein", Group_hgnc, ignore.case = TRUE) ~ "Y",
      TRUE ~ "N"
    ),

    # ===== CATALYTIC VS REGULATORY =====
    Is_catalytic = case_when(
      grepl("regulatory|modulatory|scaffold|targeting|Phosphatase and actin", Group_hgnc) ~ "N",
      TRUE ~ "Y"
    ),

    Is_regulatory = case_when(
      grepl("regulatory|modulatory|scaffold|targeting|Phosphatase and actin", Group_hgnc) ~ "Y",
      TRUE ~ "N"
    ),

    # ===== RECEPTOR VS NON-RECEPTOR (for PTPs) =====
    Is_receptor_type = case_when(
      grepl("receptor type|LAR", Group_hgnc) ~ "Y",
      TRUE ~ "N"
    ),

    # ===== PRIMARY CLASSIFICATION =====
    Class_primary = case_when(
      # Protein phosphatases - specific classes first
      grepl("tyrosine phosphatase.*non-receptor|non-receptor.*tyrosine", Group_hgnc) ~ "Non-receptor PTP",
      grepl("tyrosine phosphatase.*receptor type|LAR", Group_hgnc) ~ "Receptor PTP",
      grepl("4A family", Group_hgnc) ~ "PTP4A family",
      grepl("MAP kinase", Group_hgnc) ~ "MAP kinase phosphatase",
      grepl("CDC14|Slingshot", Group_hgnc) ~ "Dual specificity phosphatase",
      grepl("Atypical dual specificity", Group_hgnc) ~ "Atypical dual specificity",
      grepl("CDC25|Class III Cys", Group_hgnc) ~ "CDC25 phosphatase",
      grepl("CTD", Group_hgnc) ~ "CTD phosphatase",
      grepl("EYA", Group_hgnc) ~ "EYA phosphatase",
      grepl("Calcineurin", Group_hgnc) ~ "Calcineurin",
      grepl("PTEN", Group_hgnc) ~ "PTEN phosphatase",
      grepl("Protein phosphatase catalytic|Mg2\\+/Mn2\\+|Serine/threonine", Group_hgnc) ~ "Ser/Thr phosphatase",
      grepl("HAD.*protein", Group_hgnc) ~ "HAD protein phosphatase",
      grepl("HAD.*non-protein", Group_hgnc) ~ "HAD non-protein phosphatase",
      grepl("Class II Cys", Group_hgnc) ~ "Class II Cys-based",

      # Lipid phosphatases
      grepl("Phosphoinositide", Group_hgnc) ~ "Phosphoinositide phosphatase",
      grepl("Phospholipid phosphatase[^s]|Phospholipid phosphatases$", Group_hgnc) ~ "Phospholipid phosphatase",
      grepl("Sphingosine", Group_hgnc) ~ "Sphingosine phosphatase",
      grepl("Lipid phosphatase", Group_hgnc) ~ "Lipid phosphatase",

      # Metabolic phosphatases
      grepl("Acid phosphatase", Group_hgnc) ~ "Acid phosphatase",
      grepl("Alkaline phosphatase", Group_hgnc) ~ "Alkaline phosphatase",
      grepl("Glucose 6-phosphatase", Group_hgnc) ~ "Glucose-6-phosphatase",
      grepl("Fructose-1,6-bisphosphatase", Group_hgnc) ~ "Fructose-1,6-bisphosphatase",
      grepl("Bisphosphoglycerate", Group_hgnc) ~ "Bisphosphoglycerate phosphatase",
      grepl("Sugar phosphatase", Group_hgnc) ~ "Sugar phosphatase",
      grepl("Phosphoglycerate mutase", Group_hgnc) ~ "Phosphoglycerate mutase",

      # Nucleotide phosphatases
      grepl("Ectonucleotide", Group_hgnc) ~ "Ectonucleotide phosphatase",

      # Bifunctional
      grepl("phosphofructo.*biphosphatase", Group_hgnc) ~ "Bifunctional kinase/phosphatase",

      # Regulatory subunits
      grepl("PP1.*regulatory|Protein phosphatase 1 regulatory", Group_hgnc) ~ "PP1 regulatory",
      grepl("PP2.*regulatory|PP2.*modulatory|PP2.*scaffold|Protein phosphatase 2", Group_hgnc) ~ "PP2 regulatory",
      grepl("PP3.*regulatory|Protein phosphatase 3 regulatory", Group_hgnc) ~ "PP3 regulatory",
      grepl("PP4.*regulatory|Protein phosphatase 4 regulatory", Group_hgnc) ~ "PP4 regulatory",
      grepl("PP6.*regulatory|Protein phosphatase 6 regulatory", Group_hgnc) ~ "PP6 regulatory",
      grepl("Myosin|Phosphatase and actin", Group_hgnc) ~ "Cytoskeletal regulatory",

      TRUE ~ "Other"
    )
  ) %>%
  # Add phosphatase ID
  arrange(HGNC_symbol) %>%
  mutate(phosphatase_id = sprintf("P%04d", row_number()))

# Add mouse orthologs via BioMart
if (!opt$`skip-mouse`) {
  cat("\nMapping mouse orthologs via BioMart...\n")

  if (requireNamespace("biomaRt", quietly = TRUE)) {
    library(biomaRt)

    tryCatch({
      ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

      orthologs <- getBM(
        attributes = c("external_gene_name", "mmusculus_homolog_associated_gene_name",
                       "mmusculus_homolog_ensembl_gene"),
        filters = "external_gene_name",
        values = phos_final$HGNC_symbol,
        mart = ensembl
      )

      # Remove empty mappings and duplicates
      orthologs <- orthologs[orthologs$mmusculus_homolog_associated_gene_name != "", ]
      orthologs <- orthologs[!duplicated(orthologs$external_gene_name), ]

      colnames(orthologs) <- c("HGNC_symbol", "Mouse_Symbol", "Ensembl_Mouse")

      phos_final <- phos_final %>%
        left_join(orthologs, by = "HGNC_symbol")

      n_mouse <- sum(!is.na(phos_final$Mouse_Symbol))
      cat("  Mouse orthologs found:", n_mouse, "/", nrow(phos_final), "\n")

    }, error = function(e) {
      cat("  Warning: BioMart query failed:", e$message, "\n")
      phos_final <<- phos_final %>%
        mutate(Mouse_Symbol = NA_character_, Ensembl_Mouse = NA_character_)
    })
  } else {
    cat("  Warning: biomaRt package not available, skipping mouse mapping\n")
    phos_final <- phos_final %>%
      mutate(Mouse_Symbol = NA_character_, Ensembl_Mouse = NA_character_)
  }
} else {
  cat("\nSkipped mouse ortholog mapping (--skip-mouse)\n")
  phos_final <- phos_final %>%
    mutate(Mouse_Symbol = NA_character_, Ensembl_Mouse = NA_character_)
}

# Final column selection
phos_final <- phos_final %>%
  dplyr::select(
    phosphatase_id, HGNC_symbol, HGNC_ID, Name_hgnc,
    # Substrate classification
    Substrate_protein, Substrate_lipid, Substrate_nucleotide,
    Substrate_carbohydrate, Substrate_other,
    # Type classification
    Is_catalytic, Is_regulatory, Is_receptor_type,
    Class_primary,
    # HGNC info
    Group_hgnc, Group_ID_hgnc,
    Status_hgnc, Locus_type_hgnc,
    Chromosome, NCBI_Gene_ID, Ensembl_ID,
    # Mouse orthologs
    Mouse_Symbol, Ensembl_Mouse
  )

# Summary stats
cat("\n=== Phosphatase Summary ===\n")
cat("Total unique phosphatases:", nrow(phos_final), "\n")

cat("\nSubstrate classification (not mutually exclusive):\n")
cat("  Protein substrate:", sum(phos_final$Substrate_protein == "Y"), "\n")
cat("  Lipid substrate:", sum(phos_final$Substrate_lipid == "Y"), "\n")
cat("  Nucleotide substrate:", sum(phos_final$Substrate_nucleotide == "Y"), "\n")
cat("  Carbohydrate substrate:", sum(phos_final$Substrate_carbohydrate == "Y"), "\n")
cat("  Other substrate:", sum(phos_final$Substrate_other == "Y"), "\n")

cat("\nCatalytic vs Regulatory:\n")
cat("  Catalytic:", sum(phos_final$Is_catalytic == "Y"), "\n")
cat("  Regulatory:", sum(phos_final$Is_regulatory == "Y"), "\n")

cat("\nBy primary class:\n")
print(sort(table(phos_final$Class_primary), decreasing = TRUE))

cat("\nWith mouse ortholog:", sum(!is.na(phos_final$Mouse_Symbol)), "\n")

# Write output
write_csv(phos_final, output_file)
cat("\nWrote:", nrow(phos_final), "phosphatases to", output_file, "\n")
