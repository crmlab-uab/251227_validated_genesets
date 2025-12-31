# Kinases Pipeline (ordered entrypoints)

This file documents the canonical, ordered pipeline entrypoints under `scripts/kinases/bin/`.



Order (run sequentially):

- `01_fetch_geneset_BioMart.R`: fetch baseline gene lists from BioMart (species). Produces `kinases_human.csv` / `kinases_mouse.csv` in working dir.
- `02_fetch_validation_sources.R`: fetch external validation sources (KinHub, HGNC, etc) and produce raw validation files under `val_sources/`.
- `03_build_annotations.R`: assemble annotations (Manning, HGNC merges, alias augmentation, metabolic/lipid/provenance) into an annotated kinases table. This is the canonical annotation step.
- `04_map_human_to_mouse.R`: map canonical human kinases to mouse via BioMart/UniProt (produces `mouse_kinome_definitive.csv`).
- `05_export_gmt.R`: export final GMT files from validated canonical outputs into `genesets/curated/kinases/outputs/`.

**Notes:**
- All mapping and validation logic is now consolidated in canonical scripts under `scripts/kinases/lib/` (see `map_human_to_mouse_uniprot.R`, `merge_kinase_uniprot_validation.R`, `comprehensive_kinase_validation.R`).
- Deprecated scripts (`05_generate_from_biomart.R`, `05_merge_and_flag_validation_sources.R`, `06_merge_and_validate.R`) have been removed. Use only the canonical entrypoints above.

Notes:

- `lib/` contains reusable helper functions (mapping, merging, validators) and intentionally is not numbered — helpers are modular and do not imply execution order.
- `archive/` preserves original scripts and is not part of the runnable ordered pipeline.
- Checksums: canonical outputs should have colocated `.md5` files next to their CSVs.
- Example run (from repo root):




```
Rscript scripts/kinases/bin/01_fetch_geneset_BioMart.R --species human
Rscript scripts/kinases/bin/02_fetch_validation_sources.R --source=kinhub
Rscript scripts/kinases/bin/02_fetch_validation_sources.R --source=hgnc
Rscript scripts/kinases/bin/03_build_annotations.R
Rscript scripts/kinases/bin/04_map_human_to_mouse.R
Rscript scripts/kinases/bin/05_export_gmt.R
```

Maintainers: follow README.md and INVENTORY.md for more details.
