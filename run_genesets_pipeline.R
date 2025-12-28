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

if (isTRUE(cfg$steps$build)) {
  run_step("kinases/build_kinome_annotation.R", c("--species", species))
}
if (isTRUE(cfg$steps$annotate_manning)) {
  run_step("kinases/add_manning_annotation.R")
}
if (isTRUE(cfg$steps$augment_aliases)) {
  run_step("kinases/augment_matching_with_aliases.R")
}
if (isTRUE(cfg$steps$merge_kinhub)) {
  run_step("kinases/fetch_kinhub_and_merge.R")
}

cat("Pipeline completed.\n")
