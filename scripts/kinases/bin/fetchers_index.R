## Index script: grouped fetcher scripts (relative to repo root) - moved to bin
message("Sourcing kinase fetcher scripts (relative to repo root)")

source("scripts/kinases/bin/fetch_kinome.R")
source("scripts/kinases/bin/fetch_kinhub_and_merge.R")
source("scripts/kinases/bin/fetch_mouse_kinome_from_kinhub.R")
source("scripts/kinases/bin/fetch_comprehensive_mouse_kinome_biomart.R")
source("scripts/kinases/bin/fetch_definitive_mouse_kinome.R")
source("scripts/kinases/bin/fetch_uniprot_kinase_mapping_api.R")
source("scripts/kinases/bin/fetch_uniprot_mouse_mapping_api.R")
source("scripts/kinases/bin/fetch_val_sources_and_merge.R")
source("scripts/kinases/bin/fetch_hgnc_kinase_groups.R")

message("Fetcher scripts sourced. Use the functions defined in those files.")
