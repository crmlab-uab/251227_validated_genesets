# Genesets Pipeline (config-driven)

This README describes the config-driven genesets pipeline used to build and annotate kinases gene sets.

## Configuration

- **`genesets_config.yaml`** - Repository-level config (paths, species, step flags)
- **`scripts/kinases/lib/config_loader.R`** - Centralized config loader used by all scripts

## Pipeline Scripts

Located in `scripts/kinases/`:

| Script | Purpose |
|--------|---------|
| `01_fetch_geneset_BioMart.R` | Query BioMart for kinase genes by GO terms |
| `02_fetch_validation_sources.R` | Fetch external validation (KinHub, HGNC) |
| `03_build_annotations.R` | Build comprehensive annotations (HGNC groups, KEGG metabolic/lipid) |
| `04_map_human_to_mouse.R` | Map human kinases to mouse orthologs |
| `05_export_gmt.R` | Export GMT files for GSEA |

## Usage

### Option 1: Run individual scripts

```bash
cd /data/251227_validated_genesets

# Human kinases with mouse orthologs
Rscript scripts/kinases/01_fetch_geneset_BioMart.R --species=human
Rscript scripts/kinases/03_build_annotations.R --species=human
Rscript scripts/kinases/04_map_human_to_mouse.R
Rscript scripts/kinases/05_export_gmt.R
```

### Option 2: Configure via YAML

Edit `genesets_config.yaml`:

```yaml
species: "human"  # or "mouse"
input_dir: "curated/kinases/inputs"
output_dir: "curated/kinases/outputs"

steps:
  build: true
  annotate_manning: true
  augment_aliases: true
  merge_val_sources: true
```

## Output Locations

All outputs go to canonical config-driven paths:
- **Inputs**: `curated/kinases/inputs/`
  - `kinases_human_biomart.csv` / `kinases_mouse_biomart.csv`
- **Outputs**: `curated/kinases/outputs/`
  - `kinases_human_annotated.csv` / `kinases_mouse_annotated.csv`
  - `kinases_mouse_orthologs.csv`
  - `kinases_human_allsources.gmt` / `kinases_mouse_allsources.gmt`
  - `kinases_human_kinhub.csv`, `kinases_human_hgnc.csv` (validation sources)

## Notes

- All scripts use `lib/config_loader.R` for consistent path resolution
- No hardcoded relative paths - all paths derived from config
- Each script can be run individually or as part of the full pipeline
- MD5 checksums generated for all output files

## Validation Sources

Place validation files in `curated/kinases/val_sources/`:
- CSVs merged by Ensembl ID or gene symbol
- GMTs parsed and merged by gene symbol
- HTML pages with KinHub content parsed automatically

See `scripts/kinases/PIPELINE.md` for detailed pipeline documentation.
