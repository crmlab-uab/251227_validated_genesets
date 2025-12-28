library(data.table)
library(org.Mm.eg.db)
kinases <- fread('201006_composite_kinases_curated.csv', header=TRUE)
all_entrez <- keys(org.Mm.eg.db, keytype='ENTREZID')
bad <- kinases[!is.na(Entrez_Mouse) & Entrez_Mouse != '' & !(as.character(Entrez_Mouse) %in% all_entrez)]
if (nrow(bad) > 0) {
  print(bad[, .(Mouse_Symbol, Entrez_Mouse)])
} else {
  cat('All Entrez_Mouse IDs are valid.\n')
}
