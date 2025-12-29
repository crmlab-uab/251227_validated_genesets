# Author: C. Ryan Miller
# Created: 2025-12-28 02:25 CST
# Commit: 85c49284e5027d5940219124b28324bb8eba3205

#!/usr/bin/env Rscript
## Orchestrator: run functions in order (alias matching before external annotations)
library(data.table)
cfg <- if (file.exists("../genesets_config.yaml")) yaml::read_yaml("../genesets_config.yaml") else list()
cfg_get <- function(path, default) {
  parts <- strsplit(path, "\\.")[[1]]
  cur <- cfg
  for (p in parts) { if (is.null(cur[[p]])) return(default); cur <- cur[[p]] }
  cur
}

source("scripts/utilities/scripts/utilities/functions_genesets.R")

species <- ifelse(!is.null(cfg$species), cfg$species, "human")
cache_file <- cfg_get("hgnc.cache", "kinases/hgnc_lookup_cache.rds")

cat("Step: fetch_genes\n")
genes <- fetch_genes(species=species)
fwrite(genes, cfg_get("outputs.raw_fetch", "kinases/kinases_raw_fetch.csv"))

cat("Step: map_hgnc\n")
genes <- map_hgnc(genes, cache_file=cache_file)
fwrite(genes, cfg_get("outputs.mapped_hgnc", "kinases/kinases_mapped_hgnc.csv"))

cat("Step: augment_aliases (before Manning)\n")
# Manning CSV now located in kinases/data/
manning_file <- cfg_get("input_files.manning", "kinases/data/manning_2002_TableS1.csv")
if (file.exists(manning_file)) {
  man <- fread(manning_file)
} else man <- data.table()
genes <- augment_aliases(genes, man, cache_file=cache_file)
fwrite(genes, cfg_get("outputs.augmented_aliases", "kinases/kinases_augmented_aliases.csv"))

cat("Step: apply_manning\n")
genes <- apply_manning(genes, man)
fwrite(genes, cfg_get("outputs.with_manning", "kinases/kinases_with_manning.csv"))

cat("Step: dedupe\n")
genes <- dedupe_genes(genes)
fwrite(genes, cfg_get("outputs.final", "kinases/kinases_human.csv"))

if (isTRUE(cfg$steps$merge_val_sources)) {
  cat("Step: merge_val_sources\n")
  # default val_sources dir
  val_dir <- cfg_get("resources.val_sources_dir", "kinases/val_sources")
  if (dir.exists(val_dir)) {
    # delegate to fetch_val_sources_and_merge.R which handles CSV/GMT/HTML inputs
    source(file.path("kinases", "fetch_val_sources_and_merge.R"))
  } else {
    cat("Validation sources directory missing; skipping merge_val_sources\n")
  }
}

cat("Exporting GMT\n")
export_gmt(genes, cfg_get("outputs.gmt", "kinases/kinases_human.gmt"))

cat("Pipeline finished. Final output:", cfg_get("outputs.final", "kinases/kinases_human.csv"), "\n")
