## Index script: grouped mapping and merging helpers for kinases (moved to lib)
message("Sourcing kinase mapping scripts (relative to repo root)")

source("scripts/kinases/lib/map_human_to_mouse_uniprot.R")
source("scripts/kinases/lib/map_mouse_uniprot_biomart.R")
source("scripts/kinases/lib/merge_kinase_uniprot_validation.R")
source("scripts/kinases/lib/augment_matching_with_aliases.R")

message("Mapping scripts sourced. Use the functions defined in those files.")
