# Author: C. Ryan Miller
# Created: 2025-12-28 02:17 CST
# Commit: 26ec675324c74c530b55664519f79329a3d427d8

# Summarize kinases by group
library(data.table)
library(yaml)
# Load config and set input/output dirs from YAML if available
config_file <- Sys.getenv('KINASES_CONFIG', unset = 'genesets_config.yaml')
if (file.exists(config_file)) {
	cfg <- yaml::read_yaml(config_file)
	input_dir <- if (!is.null(cfg$input_dir)) cfg$input_dir else 'curated/kinases/inputs'
	output_dir <- if (!is.null(cfg$output_dir)) cfg$output_dir else 'curated/kinases/outputs'
} else {
	input_dir <- 'curated/kinases/inputs'
	output_dir <- 'curated/kinases/outputs'
}
input_file <- file.path(input_dir, '251227_curated_kinases_uniprot_validated.csv')
output_file <- file.path(output_dir, 'kinase_group_summary.csv')
kinases <- fread(input_file)
group_summary <- kinases[, .N, by=Group][order(-N)]
print(group_summary)
write.csv(group_summary, output_file, row.names=FALSE)
cat('✓ Group summary written to', output_file, '\n')
