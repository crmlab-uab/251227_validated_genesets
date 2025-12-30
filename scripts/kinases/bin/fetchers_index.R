## Index script: grouped fetcher scripts (relative to repo root) - moved to bin
message("Sourcing kinase fetcher scripts (relative to repo root)")

# Primary BioMart fetcher
if (file.exists('scripts/kinases/bin/01_fetch_geneset_BioMart.R')) source('scripts/kinases/bin/01_fetch_geneset_BioMart.R') else if (file.exists('scripts/kinases/bin/fetch_kinome.R')) source('scripts/kinases/bin/fetch_kinome.R')

# Validation sources
if (file.exists('scripts/kinases/bin/02_fetch_validation_sources.R')) source('scripts/kinases/bin/02_fetch_validation_sources.R')
if (file.exists('scripts/kinases/bin/02_fetch_validation_sources_hgnc.R')) source('scripts/kinases/bin/02_fetch_validation_sources_hgnc.R')

if (file.exists('scripts/kinases/bin/fetch_mouse_kinome_from_kinhub.R')) source('scripts/kinases/bin/fetch_mouse_kinome_from_kinhub.R')

# Optional comprehensive mouse biomart fetcher (legacy)
if (file.exists('scripts/kinases/bin/fetch_comprehensive_mouse_kinome_biomart.R')) source('scripts/kinases/bin/fetch_comprehensive_mouse_kinome_biomart.R')

# Mapping / Uniprot fetchers
if (file.exists('scripts/kinases/bin/06_map_human_to_mouse.R')) source('scripts/kinases/bin/06_map_human_to_mouse.R') else if (file.exists('scripts/kinases/bin/fetch_definitive_mouse_kinome.R')) source('scripts/kinases/bin/fetch_definitive_mouse_kinome.R')
if (file.exists('scripts/kinases/bin/fetch_uniprot_kinase_mapping_api.R')) source('scripts/kinases/bin/fetch_uniprot_kinase_mapping_api.R')
if (file.exists('scripts/kinases/bin/fetch_uniprot_mouse_mapping_api.R')) source('scripts/kinases/bin/fetch_uniprot_mouse_mapping_api.R')

# Fallback merge script
if (file.exists('scripts/kinases/bin/fetch_val_sources_and_merge.R')) source('scripts/kinases/bin/fetch_val_sources_and_merge.R')

message("Fetcher scripts sourced. Use the functions defined in those files.")
