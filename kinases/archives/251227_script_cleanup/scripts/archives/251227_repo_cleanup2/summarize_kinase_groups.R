# Summarize kinases by group
library(data.table)
kinases <- fread('251227_curated_kinases_uniprot_validated.csv')
group_summary <- kinases[, .N, by=Group][order(-N)]
print(group_summary)
write.csv(group_summary, 'kinase_group_summary.csv', row.names=FALSE)
cat('✓ Group summary written to kinase_group_summary.csv\n')
