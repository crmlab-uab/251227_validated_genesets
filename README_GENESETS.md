# Genesets pipeline (config-driven)

This README describes the small, config-driven genesets pipeline used to build and annotate kinases gene sets.

Files added:
- `genesets_config.yaml` — repository-level config (defaults and step flags).
- `run_genesets_pipeline.R` — runner that executes steps in order using `Rscript`.

Primary scripts (in `kinases/`):
- `build_kinome_annotation.R` — query BioMart/KEGG/HGNC to build base kinome table. Use `--species` to set `human` or `mouse`.
- `add_manning_annotation.R` — merge Manning Table S1 into the kinome table; uses HGNC REST for canonical symbols.
	- Manning supplement CSV is located at `kinases/data/manning_2002_TableS1.csv` (default); you can override via `genesets_config.yaml`.
- `augment_matching_with_aliases.R` — expand matching using HGNC aliases and Ensembl crosswalks (cached lookups).
- `fetch_kinhub_and_merge.R` — optional KinHub scrape and merge.

Usage:
1. Edit `genesets_config.yaml` to set inputs/outputs and enable/disable steps.
2. Run the pipeline:
```
Rscript run_genesets_pipeline.R
```

Notes:
- Each script still has sensible defaults so they can be run individually without the config file.
- HGNC REST calls are cached in `kinases/hgnc_lookup_cache.rds` when created.

Validation sources:
- Place validation files (CSV, GMT, or HTML) in `genesets/curated/` (preferred).
- The pipeline will also accept `val_sources/` or `kinases/val_sources/` if present.
- CSVs will be merged automatically by Ensembl ID or gene symbol when possible.
- GMTs will be parsed and merged by gene symbol (adds `val_sources` column).
- HTML pages with KinHub content will be delegated to the KinHub parser if detected (filename contains `kinhub`), otherwise the first HTML table is attempted.
