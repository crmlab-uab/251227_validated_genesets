## Index script: grouped validation entrypoints for kinases (moved to bin)
message("Sourcing kinase validation scripts (canonical lib/ helpers preferred)")

# Source canonical helpers from `lib/` where available
source("scripts/kinases/lib/comprehensive_kinase_validation.R")
source("scripts/kinases/bin/kinase_validation.R")
source("scripts/kinases/lib/scan_invalid_entrez.R")

message("Validation scripts sourced. Use the functions defined in those files.")
