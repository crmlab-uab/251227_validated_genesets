# Author: C. Ryan Miller
# Created: 2025-12-28 02:17 CST
# Commit: 26ec675324c74c530b55664519f79329a3d427d8

# Summarize kinases by group
library(data.table)
kinases <- fread('251227_curated_kinases_uniprot_validated.csv')
group_summary <- kinases[, .N, by=Group][order(-N)]
print(group_summary)
write.csv(group_summary, 'kinase_group_summary.csv', row.names=FALSE)
cat('✓ Group summary written to kinase_group_summary.csv\n')
