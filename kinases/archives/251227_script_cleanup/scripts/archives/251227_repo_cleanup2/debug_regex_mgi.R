# Debug regex for extracting MGI ID from Description
library(data.table)
kinases <- fread('comprehensive_mouse_kinome_biomart.csv')
cat('First 5 Description values:\n')
print(head(kinases$Description, 5))
cat('First 5 regex attempts:\n')
print(sub('.*Acc:([^]]+)]', '\\1', head(kinases$Description, 5)))
print(sub('.*Acc:([^\]]+)\]', '\\1', head(kinases$Description, 5), perl=TRUE))
print(sub('.*Acc:([^\]]+)\]', '\1', head(kinases$Description, 5), perl=TRUE))
