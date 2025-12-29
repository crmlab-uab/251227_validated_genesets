## Index script: grouped fetcher entrypoints for kinases
# Run from repository root: Rscript scripts/kinases/fetchers_index.R

message("Sourcing kinase fetcher scripts (relative to repo root)")

source("scripts/kinases/fetch_kinome.R")
source("scripts/kinases/fetch_kinhub_and_merge.R")
source("scripts/kinases/fetch_mouse_kinome_from_kinhub.R")
source("scripts/kinases/fetch_comprehensive_mouse_kinome_biomart.R")
source("scripts/kinases/fetch_definitive_mouse_kinome.R")
source("scripts/kinases/fetch_uniprot_kinase_mapping_api.R")
source("scripts/kinases/fetch_uniprot_mouse_mapping_api.R")
source("scripts/kinases/fetch_val_sources_and_merge.R")
source("scripts/kinases/fetch_hgnc_kinase_groups.R")
## consolidated fetch_hgnc_kinase_groups.R replaces the v2 variant

message("Fetcher scripts sourced. Use the functions defined in those files.")
