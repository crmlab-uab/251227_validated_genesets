# Merge UniProt mapping with kinase list and biomart mapping
# Author: C. Ryan Miller, MD, PhD (rmiller@uab.edu)

library(data.table)

# Input files
kinase_file <- 'genesets/curated/kinases/201006_composite_kinases_curated.csv'
biomart_file <- '251227_curated_kinases_biomart.csv'
uniprot_file <- 'uniprot_mouse_kinase_idmapping.tab'
output_file <- '251227_curated_kinases_uniprot_validated.csv'

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
