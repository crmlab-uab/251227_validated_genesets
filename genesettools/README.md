# genesettools

Minimal utilities to download and load gene set GMT files (MSigDB and custom GMTs).

Usage:

```r
# install via remotes::install_local("/path/to/genesettools")
library(genesettools)
sets <- load_msigdb(c("H","C2"), species = "Hs")
length(sets)
```
