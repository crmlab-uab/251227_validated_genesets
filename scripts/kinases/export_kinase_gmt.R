# Export mouse kinases as GMT for GSEA, with clarified human symbol columns
# Author: C. Ryan Miller, MD, PhD (rmiller@uab.edu)

library(data.table)

# Column mapping based on curated file:
# Human_Symbol: Manning (canonical)
# Human_Symbol2: KinHub
# Human_Symbol3: Coral
# Human_Symbol4: (spare/legacy, often matches canonical)

input_file <- '251227_curated_kinases_uniprot_validated.csv'
gmt_mouse_file <- 'mouse_kinome_all_sources.gmt'
gmt_human_manning_file <- 'human_kinome_manning.gmt'
gmt_human_kinhub_file <- 'human_kinome_kinhub.gmt'
gmt_human_coral_file <- 'human_kinome_coral.gmt'

# Load data
kinases <- fread(input_file)

# Clarify column names for output (not changing file, just for reference)
setnames(kinases, c('Human_Symbol','Human_Symbol2','Human_Symbol3','Human_Symbol4'),
         c('Human_Symbol_Manning','Human_Symbol_KinHub','Human_Symbol_Coral','Human_Symbol_Legacy'))

# 1. Export mouse kinome as a single GMT set (all mouse symbols)
cat('Exporting mouse kinome GMT...\n')
cat('MOUSE_KINOME\tAll curated mouse kinases\t',
    paste(kinases$Mouse_Symbol, collapse='\t'),
    '\n',
    file=gmt_mouse_file)

# 2. Export human kinome sets for each source
cat('Exporting human kinome GMTs by source...\n')
# Manning
cat('HUMAN_KINOME_MANNING\tHuman kinome (Manning)\t',
    paste(na.omit(unique(kinases$Human_Symbol_Manning)), collapse='\t'),
    '\n',
    file=gmt_human_manning_file)
# KinHub
cat('HUMAN_KINOME_KINHUB\tHuman kinome (KinHub)\t',
    paste(na.omit(unique(kinases$Human_Symbol_KinHub)), collapse='\t'),
    '\n',
    file=gmt_human_kinhub_file)
# Coral
cat('HUMAN_KINOME_CORAL\tHuman kinome (Coral)\t',
    paste(na.omit(unique(kinases$Human_Symbol_Coral)), collapse='\t'),
    '\n',
    file=gmt_human_coral_file)

cat('✓ GMT export complete.\n')
