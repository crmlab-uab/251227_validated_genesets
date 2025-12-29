# Author: C. Ryan Miller
# Created: 2025-12-28 02:17 CST
# Commit: 26ec675324c74c530b55664519f79329a3d427d8

library(data.table)
library(org.Mm.eg.db)
inputs_dir <- 'genesets/curated/kinases/inputs'
candidates <- list.files(inputs_dir, pattern='201006_composite_kinases_curated.*\\.csv$', full.names=TRUE, ignore.case=TRUE)
if (length(candidates) == 0) stop('Missing input snapshot: please place 201006_composite_kinases_curated__YYMMDD.csv in ', inputs_dir)
kinases_file <- sort(candidates, decreasing=TRUE)[1]
kinases <- fread(kinases_file, header=TRUE)
all_entrez <- keys(org.Mm.eg.db, keytype='ENTREZID')
bad <- kinases[!is.na(Entrez_Mouse) & Entrez_Mouse != '' & !(as.character(Entrez_Mouse) %in% all_entrez)]
if (nrow(bad) > 0) {
  print(bad[, .(Mouse_Symbol, Entrez_Mouse)])
} else {
  cat('All Entrez_Mouse IDs are valid.\n')
}
