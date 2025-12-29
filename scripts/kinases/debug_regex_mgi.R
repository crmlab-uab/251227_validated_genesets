# Author: C. Ryan Miller
# Created: 2025-12-28 02:17 CST
# Commit: 26ec675324c74c530b55664519f79329a3d427d8

# Debug regex for extracting MGI ID from Description
library(data.table)
kinases <- fread('comprehensive_mouse_kinome_biomart.csv')
cat('First 5 Description values:\n')
print(head(kinases$Description, 5))
cat('First 5 regex attempts:\n')
print(sub('.*Acc:([^]]+)]', '\\1', head(kinases$Description, 5)))
print(sub('.*Acc:([^\]]+)\]', '\\1', head(kinases$Description, 5), perl=TRUE))
print(sub('.*Acc:([^\]]+)\]', '\1', head(kinases$Description, 5), perl=TRUE))
