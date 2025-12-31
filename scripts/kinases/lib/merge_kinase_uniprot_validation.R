# Merge UniProt mapping with kinase list and biomart mapping
#
# NOTE: This is the canonical script for merging UniProt and BioMart mappings in the kinases pipeline.
# All merging logic should be consolidated here. Do not duplicate merging logic elsewhere.
# Author: C. Ryan Miller, MD, PhD (rmiller@uab.edu)

library(data.table)


# Load config and set input/output dirs from YAML if available
library(yaml)
config_file <- Sys.getenv('KINASES_CONFIG', unset = 'genesets_config.yaml')
if (file.exists(config_file)) {
	cfg <- yaml::read_yaml(config_file)
	input_dir <- if (!is.null(cfg$input_dir)) cfg$input_dir else 'curated/kinases/inputs'
	output_dir <- if (!is.null(cfg$output_dir)) cfg$output_dir else 'curated/kinases/outputs'
} else {
	input_dir <- 'curated/kinases/inputs'
	output_dir <- 'curated/kinases/outputs'
}
candidates <- list.files(input_dir, pattern='201006_composite_kinases_curated.*\\.csv$', full.names=TRUE, ignore.case=TRUE)
if (length(candidates) == 0) stop('Missing input snapshot: please place 201006_composite_kinases_curated__YYMMDD.csv in ', input_dir)
kinase_file <- sort(candidates, decreasing=TRUE)[1]
biomart_file <- file.path(output_dir, '251227_curated_kinases_biomart.csv')
uniprot_file <- file.path(output_dir, 'uniprot_mouse_kinase_idmapping.tab')
output_file <- file.path(output_dir, '251227_curated_kinases_uniprot_validated.csv')

# Load data
kinases <- fread(kinase_file, header=TRUE)
biomart <- fread(biomart_file, header=TRUE)
uniprot <- fread(uniprot_file, header=TRUE)

# Collapse UniProt entries per gene symbol
uniprot_collapsed <- uniprot[, .(UniProt_Mouse_UniProtAPI = paste(unique(Entry), collapse=';')), by=Mouse_Symbol]

# Collapse biomart mapping per gene symbol
biomart_collapsed <- biomart[, .(UniProt_Mouse_biomart = paste(na.omit(unique(UniProt_Mouse_biomart)), collapse=';')), by=Mouse_Symbol]

# Merge all sources to kinases
merged <- merge(kinases, biomart_collapsed, by='Mouse_Symbol', all.x=TRUE)
merged <- merge(merged, uniprot_collapsed, by='Mouse_Symbol', all.x=TRUE)

# Annotate source columns
merged[, biomart_source := ifelse(UniProt_Mouse_biomart != '' & !is.na(UniProt_Mouse_biomart), 'biomart', NA)]
merged[, uniprotapi_source := ifelse(UniProt_Mouse_UniProtAPI != '' & !is.na(UniProt_Mouse_UniProtAPI), 'uniprot_api', NA)]

# Collapse to one row per Mouse_Symbol, keeping all annotation columns

# Rename human symbol columns for kinome tree clarity
collapsed <- merged[, .(
	Human_Symbol_Manning = first(Human_Symbol),
	Human_Symbol_KinHub = first(Human_Symbol2),
	Human_Symbol_Coral = first(Human_Symbol3),
	Human_Symbol_Legacy = first(Human_Symbol4),
	Description = first(Description),
	Group = first(Group),
	Family = first(Family),
	Subfamily = first(Subfamily),
	UniProt_Human = first(UniProt_Human),
	MGI_Mouse = first(MGI_Mouse),
	Entrez_Mouse = first(Entrez_Mouse),
	Source = first(Source),
	Date = first(Date),
	UniProt_Mouse_biomart = first(UniProt_Mouse_biomart),
	UniProt_Mouse_UniProtAPI = first(UniProt_Mouse_UniProtAPI),
	biomart_source = first(biomart_source),
	uniprotapi_source = first(uniprotapi_source)
), by=Mouse_Symbol]

fwrite(collapsed, output_file)
cat('✓ Collapsed kinase UniProt validation table saved to', output_file, '\n')
