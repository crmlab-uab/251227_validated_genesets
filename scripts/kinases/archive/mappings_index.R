## Index script: grouped mapping and merging helpers for kinases
# Run from repository root: Rscript scripts/kinases/mappings_index.R

message("Sourcing kinase mapping scripts (relative to repo root)")

source("scripts/kinases/map_human_to_mouse_uniprot.R")
source("scripts/kinases/map_mouse_uniprot_biomart.R")
source("scripts/kinases/merge_kinase_uniprot_validation.R")
source("scripts/kinases/augment_matching_with_aliases.R")

message("Mapping scripts sourced. Use the functions defined in those files.")
