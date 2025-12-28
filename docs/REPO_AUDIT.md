 # Repository Genesets Audit
 Date: 2025-12-28

 ## Inventory
 See `docs/file_inventory.txt` for the list of gene set files discovered under `data/`.

 ### Files summarized

 - `data/genesets/sources/kinases/201006_composite_kinases_curated.csv` — curated kinases composite (539 genes)
 - `data/genesets/sources/kinases/manning_2002_TableS1.csv` — Manning Table S1 (raw)

 ## Findings

 - Repository contains both canonical source artifacts under `data/genesets/sources/kinases` and a minimal `data/genesets/curated/kinases` symlink for compatibility.
 - No other gene set files detected under `data/` in this repo.

 ## Recommendations

 1. Keep the canonical copies under `data/genesets/sources/` and move any other curated or project-specific sets into `data/genesets/custom/`.
 2. Maintain per-file metadata (`*.meta.yaml`) and checksums (already added for the curated kinases file).
 3. Add manifest entries for all canonical files (manifest exists and includes the kinases file).
 4. Add CI checks to validate manifest and checksums (deferred per request).

 ## Actions taken

 - Consolidated curated kinases files into `data/genesets/sources/kinases/`.
 - Added `data/genesets/manifest.yaml` with the kinases entry.
 - Added `data/genesets/sources/kinases/201006_composite_kinases_curated.meta.yaml` and checksum file.
 - Created a symlink at `data/genesets/curated/kinases/201006_composite_kinases_curated.csv` pointing to the canonical file for backward compatibility.
