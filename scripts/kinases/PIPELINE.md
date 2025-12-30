# Kinases Pipeline (ordered entrypoints)

This file documents the canonical, ordered pipeline entrypoints under `scripts/kinases/bin/`.

Order (run sequentially):

- `01_fetch_geneset_BioMart.R`: fetch baseline gene lists from BioMart (species). Produces `kinases_human.csv` / `kinases_mouse.csv` in working dir.
- `02_fetch_validation_sources.R`: fetch external validation sources (KinHub, HGNC, etc) and produce raw validation files under `val_sources/`.
- `04_build_annotations.R`: assemble annotations (Manning, HGNC merges, alias augmentation) into an annotated kinases table.
- `05_generate_from_biomart.R`: run the generic BioMart generator (InterPro + GO) that writes intermediates to `output/temp/` and union to `genesets/curated/kinases/outputs/kinases_{species}_union.csv`.
- `05_merge_and_flag_validation_sources.R`: merge validation sources into the union and run validation checks (produces merged validation CSVs and reports).
- `06_merge_and_validate.R`: merge validation sources and run validations for the kinases pipeline.
- `07_map_human_to_mouse.R`: map canonical human kinases to mouse via BioMart/UniProt (produces `mouse_kinome_definitive.csv`).
- `08_export_gmt.R`: export final GMT files from validated canonical outputs into `genesets/curated/kinases/outputs/`.

Notes:

- `lib/` contains reusable helper functions (mapping, merging, validators) and intentionally is not numbered — helpers are modular and do not imply execution order.
- `archive/` preserves original scripts and is not part of the runnable ordered pipeline.
- Checksums: canonical outputs should have colocated `.md5` files next to their CSVs.
- Example run (from repo root):


```
Rscript scripts/kinases/bin/01_fetch_geneset_BioMart.R --species human
Rscript scripts/kinases/bin/02_fetch_validation_sources.R --source=kinhub
Rscript scripts/kinases/bin/02_fetch_validation_sources.R --source=hgnc
Rscript scripts/kinases/bin/04_build_annotations.R
Rscript scripts/kinases/bin/05_generate_from_biomart.R --species human
Rscript scripts/kinases/bin/05_merge_and_flag_validation_sources.R
Rscript scripts/kinases/bin/06_merge_and_validate.R
Rscript scripts/kinases/bin/07_map_human_to_mouse.R
Rscript scripts/kinases/bin/08_export_gmt.R
```

Maintainers: follow README.md and INVENTORY.md for more details.
