## Index script: grouped validation and export helpers for kinases
# Run from repository root: Rscript scripts/kinases/validations_index.R

message("Sourcing kinase validation scripts (relative to repo root)")

source("scripts/kinases/comprehensive_kinase_validation.R")
source("scripts/kinases/kinase_validation.R")
source("scripts/kinases/validate_against_genome.R")
source("scripts/kinases/scan_invalid_entrez.R")
source("scripts/kinases/clean_comprehensive_mouse_kinome.R")
source("scripts/kinases/export_kinase_gmt.R")
source("scripts/kinases/build_kinome_annotation.R")

message("Validation scripts sourced. Use the functions defined in those files.")
