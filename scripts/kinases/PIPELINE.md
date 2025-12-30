# Kinases Pipeline (ordered entrypoints)

This file documents the canonical, ordered pipeline entrypoints under `scripts/kinases/bin/`.

Order (run sequentially):

- `01_fetch_geneset_BioMart.R`: fetch baseline gene lists from BioMart (species). Produces `kinases_human.csv` / `kinases_mouse.csv` in working dir.
- `02_fetch_validation_sources*.R`: fetch external validation sources (KinHub, HGNC groups) and produce raw validation files under `val_sources/` or `kinases/`.
- `03_build_annotations.R`: assemble annotations (Manning, HGNC merges, alias augmentation) into an annotated kinases table.
- `04_generate_from_biomart.R`: run the generic BioMart generator (InterPro + GO) that writes intermediates to `output/temp/` and union to `genesets/curated/kinases/outputs/kinases_{species}_union.csv`.
- `05_merge_and_validate.R`: merge validation sources into the union and run validation checks (produces merged validation CSVs and reports).
- `06_map_human_to_mouse.R`: map canonical human kinases to mouse via BioMart/UniProt (produces `mouse_kinome_definitive.csv`).
- `07_export_gmt.R`: export final GMT files from validated canonical outputs into `genesets/curated/kinases/outputs/`.

Notes:

- `lib/` contains reusable helper functions (mapping, merging, validators) and intentionally is not numbered — helpers are modular and do not imply execution order.
- `archive/` preserves original scripts and is not part of the runnable ordered pipeline.
- Checksums: canonical outputs should have colocated `.md5` files next to their CSVs.
- Example run (from repo root):

```
Rscript scripts/kinases/bin/01_fetch_geneset_BioMart.R --species human
Rscript scripts/kinases/bin/04_generate_from_biomart.R --species human
Rscript scripts/kinases/bin/05_merge_and_validate.R
Rscript scripts/kinases/bin/07_export_gmt.R
```

Maintainers: follow README.md and INVENTORY.md for more details.
