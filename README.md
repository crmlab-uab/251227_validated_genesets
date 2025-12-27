# Validated Gene Sets for RNA-seq Analysis

**Version:** 1.0.0  
**Date:** December 27, 2025, 02:40 CST  
**Author:** C. Ryan Miller, MD, PhD (rmiller@uab.edu)

## Purpose

This repository contains **independently validated gene sets** for use in RNA-seq differential expression and pathway analysis. All gene sets have been curated from published literature and validated against authoritative databases (UniProt, Gene Ontology, Bioconductor annotations).

## Repository Organization

```
251227_validated_genesets/
├── kinases/                    # Mouse kinome (540 genes, 97.8% validated)
│   ├── 201006_composite_kinases_curated.csv
│   ├── mouse_kinome_comprehensive_validation.csv
│   ├── VALIDATION_REPORT.md
│   └── comprehensive_kinase_validation.R
├── pathways/                   # Curated pathway gene sets
│   └── VERHAAK_GBM_Subtypes.gmt
├── docs/                       # Documentation
│   ├── VALIDATION_METHODS.md
│   └── CHANGELOG.md
└── scripts/                    # Sync and verification tools
    └── sync_to_analysis_repo.sh
```

## Validated Gene Sets

### Mouse Kinome (v1.0)
- **Source:** Composite from KinHub, Coral, UniProt databases
- **Total genes:** 540 kinases
- **Validation:** 528/540 (97.8%) confirmed via GO kinase activity terms
- **Date curated:** October 6, 2020
- **Date validated:** December 27, 2025
- **Species:** Mouse (Mus musculus)
- **File:** `kinases/201006_composite_kinases_curated.csv`

### GBM Subtype Signatures (Verhaak)
- **Source:** Verhaak et al. Cancer Cell 2010
- **Subtypes:** Classical, Mesenchymal, Proneural
- **Species:** Human → Mouse orthologs
- **File:** `pathways/VERHAAK_GBM_Subtypes.gmt`

## Usage in Analysis Projects

### 1. Copy gene sets to your analysis repository

```bash
# From your analysis repo (e.g., 251207_mNSC)
cp /data/251227_validated_genesets/kinases/*.csv data/genesets/
cp /data/251227_validated_genesets/pathways/*.gmt data/genesets/custom/

# Generate checksums for verification
cd data/genesets
md5sum *.csv custom/*.gmt > ../docs/geneset_checksums.txt
```

### 2. Reference in your config

```yaml
# config.yaml
custom_gmt_files:
  - "data/genesets/custom/VERHAAK_GBM_Subtypes.gmt"

kinase_list: "data/genesets/201006_composite_kinases_curated.csv"
```

### 3. Verify checksums before analysis

```bash
# Compare checksums to ensure files are unmodified
md5sum -c docs/geneset_checksums.txt
```

## Validation Methodology

All gene sets undergo multi-database validation:

1. **Gene Ontology (GO)** - Functional annotation from Bioconductor
2. **UniProt** - Protein function and keywords
3. **Bioconductor AnnotationDbi** - Cross-species ortholog mapping
4. **Literature curation** - Manual verification against source publications

See `docs/VALIDATION_METHODS.md` for detailed protocols.

## Citation Policy

If you use these gene sets in publications, please cite:

### Kinome Validation
- **This repository:** `251227_validated_genesets v1.0` (DOI: pending)
- **Original sources:** KinHub, Coral, Manning kinase classification

### Pathway Gene Sets
- **Verhaak GBM subtypes:** Verhaak et al. Cancer Cell 2010 (PMID: 20129251)

## Version History

- **v1.0.0** (2025-12-27): Initial release with mouse kinome (540 genes, 528 validated)

## Contributing

To add new validated gene sets:
1. Follow validation protocol in `docs/VALIDATION_METHODS.md`
2. Generate comprehensive validation report
3. Update this README with new gene set metadata
4. Update `CHANGELOG.md`

## Quality Standards

All gene sets must meet:
- ✅ **≥95% validation rate** against authoritative databases
- ✅ **Traceable provenance** to published literature
- ✅ **Species annotation** clearly documented
- ✅ **Validation script** included for reproducibility

## Contact

Questions? Open an issue or contact: rmiller@uab.edu

## License

MIT License - See LICENSE file for details
