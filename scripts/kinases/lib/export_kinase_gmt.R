# Export mouse kinases as GMT for GSEA, with clarified human symbol columns
# Author: C. Ryan Miller, MD, PhD (rmiller@uab.edu)

library(data.table)

# Column mapping based on curated file:
# Human_Symbol: Manning (canonical)
# Human_Symbol2: KinHub
# Human_Symbol3: Coral
# Human_Symbol4: (spare/legacy, often matches canonical)



# Load config and set output_dir from YAML if available
config_file <- Sys.getenv('KINASES_CONFIG', unset = 'genesets_config.yaml')
if (file.exists(config_file)) {
    cfg <- yaml::read_yaml(config_file)
    output_dir <- if (!is.null(cfg$output_dir)) cfg$output_dir else 'curated/kinases/outputs'
} else {
    output_dir <- 'curated/kinases/outputs'
}
gmt_mouse_file <- file.path(output_dir, 'mouse_kinome_all_sources.gmt')

# Load data
kinases <- fread(input_file)

# Only export mouse kinome GMT set (all mouse symbols)
cat('Exporting mouse kinome GMT...\n')
cat('MOUSE_KINOME\tAll curated mouse kinases\t',
    paste(kinases$Mouse_symbol, collapse='\t'),
    '\n',
    file=gmt_mouse_file)

cat('✓ GMT export complete.\n')
