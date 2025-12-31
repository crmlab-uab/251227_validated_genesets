#!/usr/bin/env Rscript

# Wrapper to export GMT using lib/export_kinase_gmt.R
# Add file existence check for the main input file used in export_kinase_gmt.R
input_file <- 'mouse_kinome_definitive.csv'
output_dir <- '../../../curated/kinases/outputs'
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
if (!file.exists(input_file) || file.info(input_file)$size == 0) {
  stop(paste0('Required input file missing or empty: ', input_file, '\nCheck that previous pipeline steps completed successfully.'))
}
source('../../../scripts/kinases/lib/export_kinase_gmt.R')
