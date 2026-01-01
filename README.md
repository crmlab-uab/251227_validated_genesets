# Validated Gene Sets for RNA-seq Analysis

**Version:** 1.1.0
**Date:** December 31, 2025
**Author:** C. Ryan Miller, MD, PhD (rmiller@uab.edu)

## Purpose

This repository contains **independently validated gene sets** for use in RNA-seq differential expression and pathway analysis. All gene sets have been curated from published literature and validated against authoritative databases (UniProt, Gene Ontology, Bioconductor annotations).

## Repository Organization

```
251227_validated_genesets/
├── curated/
│   └── kinases/
│       ├── inputs/           # Source gene lists and validation inputs
│       ├── outputs/          # Final validated gene sets (CSV, GMT)
│       └── val_sources/      # Intermediate validation files
├── scripts/
│   └── kinases/
│       ├── 01_fetch_geneset_BioMart.R
│       ├── 02_fetch_validation_sources.R
│       ├── 03_build_annotations.R
│       ├── 04_map_human_to_mouse.R
│       ├── 05_export_gmt.R
│       ├── lib/              # Reusable helpers & config_loader.R
│       └── PIPELINE.md       # Canonical pipeline documentation
├── pathways/                 # Curated pathway gene sets
│   └── VERHAAK_GBM_Subtypes.gmt
├── docs/                     # Documentation
│   ├── VALIDATION_METHODS.md
│   └── CHANGELOG.md
└── genesets_config.yaml      # Pipeline configuration
```

## Quick Start

```bash
cd /data/251227_validated_genesets

# Run full human kinase pipeline with mouse orthologs
Rscript scripts/kinases/01_fetch_geneset_BioMart.R --species=human
Rscript scripts/kinases/03_build_annotations.R --species=human
Rscript scripts/kinases/04_map_human_to_mouse.R
Rscript scripts/kinases/05_export_gmt.R

# Check outputs
ls -la curated/kinases/outputs/
```

## Pipeline Scripts

All scripts use centralized configuration from `lib/config_loader.R`:
- Automatic repo root detection
- Config-driven paths (no hardcoded relative paths)
- Consistent output locations

| Script | Purpose |
|--------|---------|
| `01_fetch_geneset_BioMart.R` | Fetch kinase genes from Ensembl BioMart |
| `02_fetch_validation_sources.R` | Fetch KinHub/HGNC validation data |
| `03_build_annotations.R` | Build comprehensive annotations (HGNC, KEGG) |
| `04_map_human_to_mouse.R` | Map human kinases to mouse orthologs |
| `05_export_gmt.R` | Export GMT files for GSEA |

See `scripts/kinases/PIPELINE.md` for detailed documentation.

## Validated Gene Sets

### Human Kinome (v1.1)
- **Source:** BioMart (GO kinase activity terms), HGNC groups, KEGG pathways
- **Total genes:** 810 annotated kinases
- **Annotations:** HGNC kinase groups, metabolic kinases, lipid kinases
- **Date validated:** December 31, 2025
- **File:** `curated/kinases/outputs/03_kinases_human__YYMMDD.csv`

### Mouse Kinome (v1.1)
- **Source:** Human-to-mouse ortholog mapping via BioMart
- **Total genes:** 538 kinases
- **Validation:** Cross-referenced with curated kinase lists
- **Date validated:** December 31, 2025
- **Files:**
  - `curated/kinases/outputs/04_mouse_kinome_definitive__YYMMDD.csv`
  - `curated/kinases/outputs/05_mouse_kinome_all_sources__YYMMDD.gmt`

### GBM Subtype Signatures (Verhaak)
- **Source:** Verhaak et al. Cancer Cell 2010
- **Subtypes:** Classical, Mesenchymal, Proneural
- **Species:** Human to Mouse orthologs
- **File:** `pathways/VERHAAK_GBM_Subtypes.gmt`

## Usage in Analysis Projects

### 1. Copy gene sets to your analysis repository

```bash
# From your analysis repo (e.g., 251207_mNSC)
cp /data/251227_validated_genesets/curated/kinases/outputs/*.csv genesets/
cp /data/251227_validated_genesets/curated/kinases/outputs/*.gmt genesets/custom/
cp /data/251227_validated_genesets/pathways/*.gmt genesets/custom/
```

### 2. Reference in your config

```yaml
# config.yaml
custom_gmt_files:
  - "genesets/custom/VERHAAK_GBM_Subtypes.gmt"
  - "genesets/custom/05_mouse_kinome_all_sources__251231.gmt"

kinase_list: "genesets/03_kinases_human__251231.csv"
```

### 3. Verify checksums before analysis

```bash
# All outputs have .md5 checksum files
md5sum -c curated/kinases/outputs/*.md5
```

## Validation Methodology

All gene sets undergo multi-database validation:

1. **Ensembl BioMart** - GO kinase activity terms (GO:0004672, GO:0004674, GO:0016301, GO:0016773)
2. **HGNC** - Kinase gene groups and families
3. **KEGG** - Metabolic and lipid pathway annotations
4. **UniProt** - Protein function validation
5. **Cross-species mapping** - BioMart ortholog queries

See `docs/VALIDATION_METHODS.md` for detailed protocols.

## Version History

- **v1.1.0** (2025-12-31): Config-driven IO refactoring, updated kinase annotations
- **v1.0.0** (2025-12-27): Initial release with mouse kinome (540 genes, 528 validated)

## Quality Standards

All gene sets must meet:
- **>=95% validation rate** against authoritative databases
- **Traceable provenance** to published literature
- **Species annotation** clearly documented
- **Validation script** included for reproducibility
- **MD5 checksums** for all output files

## Contact

Questions? Open an issue or contact: rmiller@uab.edu

## License

MIT License - See LICENSE file for details
