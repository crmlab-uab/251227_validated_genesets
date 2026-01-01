# Session Notes: Phosphatases & GMT Export
**Date:** 2026-01-01
**Branch:** master

## Summary

Added comprehensive phosphatase gene set with substrate classification (matching kinases approach). Created GMT export utility for all curated gene sets. Updated TF dataset with mouse orthologs.

## Changes Made

### New Files

1. **curated/phosphatases/** - Phosphatase curation pipeline
   - `outputs/phosphatases_human_curated.csv` - 437 phosphatases with annotations
   - `outputs/phosphatases_human_curated.provenance.yaml` - Data provenance

2. **scripts/phosphatases/create_phosphatases.R** - Phosphatase extraction script
   - Extracts all phosphatase groups from HGNC (43 groups)
   - Classifies by substrate: protein, lipid, nucleotide, carbohydrate, other
   - Annotates catalytic vs regulatory subunits
   - Adds mouse orthologs via BioMart

3. **scripts/utilities/export_curated_gmt.R** - GMT export utility
   - Generates GMT files for all curated sets
   - Creates human and mouse versions
   - Outputs combined GMT files

4. **curated/gmt/** - GMT output files
   - kinases_human.gmt, kinases_mouse.gmt
   - phosphatases_human.gmt, phosphatases_mouse.gmt
   - tf_human.gmt, tf_mouse.gmt
   - curated_genesets_human.gmt, curated_genesets_mouse.gmt

### Modified Files

1. **curated/tf/outputs/tf_human_curated.csv** - Added mouse orthologs (1,429/1,852)

## Phosphatase Statistics

| Category | Count |
|----------|-------|
| Total phosphatases | 437 |
| Catalytic | 222 |
| Regulatory | 215 |
| Protein substrate | 357 |
| Lipid substrate | 44 |
| With mouse orthologs | 427 |

## GMT Gene Set Contents

### Human
- KINASES_HUMAN_ALL: 546 genes
- KINASES_HUMAN_PROTEIN: 388 genes
- PHOSPHATASES_HUMAN_ALL: 437 genes
- PHOSPHATASES_HUMAN_PROTEIN: 357 genes
- PHOSPHATASES_HUMAN_CATALYTIC: 222 genes
- TF_HUMAN_ALL: 1,852 genes

### Mouse
- KINASES_MOUSE_ALL: 530 genes
- PHOSPHATASES_MOUSE_ALL: 427 genes
- TF_MOUSE_ALL: 1,429 genes

## Data Sources

- HGNC Gene Groups (primary)
- Ensembl BioMart (mouse orthologs)

## Next Steps

- Add validation against PhosphoSitePlus
- Consider adding phosphatase complexes
- Sync GMT files to bRNA3F template repo (done this session)
