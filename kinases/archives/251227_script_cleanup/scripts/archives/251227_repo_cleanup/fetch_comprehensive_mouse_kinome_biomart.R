# Fetch a comprehensive list of mouse kinases from Ensembl BioMart
# Includes all genes annotated with GO:0004672 (protein kinase activity) and related terms
# Author: C. Ryan Miller, MD, PhD (rmiller@uab.edu)

library(biomaRt)
library(data.table)

# Connect to Ensembl BioMart
ensembl <- useMart('ensembl', dataset='mmusculus_gene_ensembl')

# Fetch all genes with protein kinase activity or related GO terms
kinase_go_terms <- c('GO:0004672', # protein kinase activity
                     'GO:0004674', # protein serine/threonine kinase activity
                     'GO:0004713', # protein tyrosine kinase activity
                     'GO:0008289', # lipid kinase activity
                     'GO:0016301', # kinase activity
                     'GO:0004679', # AMP-activated protein kinase activity
                     'GO:0004697', # cyclin-dependent protein serine/threonine kinase activity
                     'GO:0004702', # receptor signaling protein serine/threonine kinase activity
                     'GO:0004714', # transmembrane receptor protein tyrosine kinase activity
                     'GO:0004715', # non-membrane spanning protein tyrosine kinase activity
                     'GO:0001727', # sphingosine kinase activity
                     'GO:0008440', # ceramide kinase activity
                     'GO:0004673'  # protein histidine kinase activity
)

cat('Querying BioMart for mouse kinases...\n')

kinases <- getBM(
  attributes = c('mgi_symbol', 'ensembl_gene_id', 'entrezgene_id', 'uniprotswissprot', 'description', 'go_id', 'name_1006'),
  filters = 'go',
  values = kinase_go_terms,
  mart = ensembl
)

kinases <- as.data.table(kinases)
setnames(kinases, c('mgi_symbol','ensembl_gene_id','entrezgene_id','uniprotswissprot','description','go_id','name_1006'),
         c('Mouse_Symbol','Ensembl_Gene_ID','Entrez_ID','UniProt_ID','Description','GO_ID','GO_Term'))

fwrite(kinases, 'comprehensive_mouse_kinome_biomart.csv')
cat('✓ Comprehensive mouse kinome table saved to comprehensive_mouse_kinome_biomart.csv\n')
