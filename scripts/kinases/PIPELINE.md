# Kinases Pipeline

This file documents the canonical, ordered pipeline entrypoints under `scripts/kinases/`.

## Configuration

All scripts use centralized configuration from `lib/config_loader.R`, which:
- Automatically finds repo root by searching for `genesets_config.yaml`
- Loads paths from config file
- Provides helper functions: `input_path()`, `output_path()`, `val_sources_path()`
- Ensures all outputs go to canonical locations (no hardcoded paths)

## Pipeline Steps (run sequentially)

| Step | Script | Purpose | Output |
|------|--------|---------|--------|
| 01 | `01_fetch_geneset_BioMart.R` | Fetch baseline kinase genes from BioMart | `inputs/kinases_human_biomart.csv` or `inputs/kinases_mouse_biomart.csv` |
| 02 | `02_fetch_validation_sources.R` | Fetch external validation sources (KinHub, HGNC) | `outputs/kinases_human_kinhub.csv`, `outputs/kinases_human_hgnc.csv` |
| 03 | `03_build_annotations.R` | Build comprehensive annotations (HGNC, KEGG metabolic/lipid) | `outputs/kinases_human_annotated.csv` or `outputs/kinases_mouse_annotated.csv` |
| 04 | `04_map_human_to_mouse.R` | Map human kinases to mouse orthologs via BioMart | `outputs/kinases_mouse_orthologs.csv` |
| 05 | `05_export_gmt.R` | Export GMT files for GSEA | `outputs/kinases_human_allsources.gmt`, `outputs/kinases_mouse_allsources.gmt` |

## Example Run (from repo root)

```bash
cd /data/251227_validated_genesets

# Human kinases with mouse orthologs
Rscript scripts/kinases/01_fetch_geneset_BioMart.R --species=human
Rscript scripts/kinases/02_fetch_validation_sources.R --source=kinhub
Rscript scripts/kinases/03_build_annotations.R --species=human
Rscript scripts/kinases/04_map_human_to_mouse.R
Rscript scripts/kinases/05_export_gmt.R

# Verify outputs
ls -la curated/kinases/outputs/
```

## Output Locations

All outputs use config-driven canonical paths:
- **Inputs**: `curated/kinases/inputs/`
- **Outputs**: `curated/kinases/outputs/`

## Notes

- `lib/` contains reusable helper functions and `config_loader.R`
- `archives/` preserves original scripts (not part of runnable pipeline)
- All canonical outputs have colocated `.md5` checksum files
- Species can be set via `--species=human|mouse` argument or in `genesets_config.yaml`

## Maintainers

See README.md and INVENTORY.md for more details.
