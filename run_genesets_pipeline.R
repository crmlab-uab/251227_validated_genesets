# Author: C. Ryan Miller
# Created: 2025-12-28 02:25 CST
# Commit: 85c49284e5027d5940219124b28324bb8eba3205

#!/usr/bin/env Rscript
# Runner for genesets pipeline: reads genesets_config.yaml and runs steps in order
library(yaml)
library(glue)

cfgf <- if (file.exists("genesets_config.yaml")) "genesets_config.yaml" else if (file.exists("config/genesets_config.yaml")) "config/genesets_config.yaml" else NULL
if (is.null(cfgf)) stop("genesets_config.yaml not found; create one or use defaults")
cfg <- yaml::read_yaml(cfgf)

run_step <- function(script, args = character()) {
  cmd <- c(script, args)
  cat(glue("Running: Rscript {script} {paste(args, collapse=' ')}\n"))
  res <- system2("Rscript", cmd)
  if (res != 0) stop(glue("Step failed: {script} (exit {res})"))
}

# default species
species <- ifelse(!is.null(cfg$species), cfg$species, "human")

if (isTRUE(cfg$steps$build) || isTRUE(cfg$steps$annotate_manning) || isTRUE(cfg$steps$augment_aliases) || isTRUE(cfg$steps$merge_kinhub)) {
  # prefer consolidated function-driven orchestrator that implements desired ordering
  run_step("kinases/run_pipeline_functions.R")
} else {
  cat("No steps enabled in genesets_config.yaml; nothing to run.\n")
}

cat("Pipeline completed.\n")
