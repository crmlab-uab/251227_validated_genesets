## Index script: grouped annotation helpers for kinases
# Run from repository root: Rscript scripts/kinases/annotations_index.R

message("Sourcing kinase annotation scripts (relative to repo root)")

source("scripts/kinases/annotate_lipid_kinases_from_kegg.R")
source("scripts/kinases/annotate_metabolic_kinases_from_kegg.R")
source("scripts/kinases/add_manning_annotation.R")
source("scripts/kinases/summarize_kinase_groups.R")

message("Annotation scripts sourced. Use the functions defined in those files.")
