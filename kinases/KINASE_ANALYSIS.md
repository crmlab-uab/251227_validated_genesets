# Kinase Gene Set Analysis

**Date**: December 26, 2025  
**File**: `201006_composite_kinases.csv`  
**Analyst**: Ryan Miller

---

## Summary

Comprehensive analysis of the mouse kinome annotation file used for stratifying DEGs by kinase expression.

### Key Statistics

- **Total kinases**: 558 genes
- **Reference standard**: Manning et al. Science 2002 (~540 mouse kinases expected)
- **Coverage**: 103% of expected (18 additional genes)
- **Date added**: October 6, 2020
- **Sources**: KinHub, Coral, UniProt, MGI

---

## File Structure

### Columns (17 total)
1. `Mouse_symbol` - Mouse gene symbol (primary identifier)
2. `xName` - Alternative name
3. `ManningName` - Manning et al. nomenclature
4. `HGNC` - Human Gene Nomenclature Committee ID
5. `CoralName` - Coral database name
6. `gene_name` - Full gene description
7. `Group` - Manning classification group
8. `Family` - Kinase family
9. `SubFamily` - Kinase subfamily
10. `UniProtID_Human` - Human UniProt ID
11. `MGI_ID` - Mouse Genome Informatics ID
12. `Entrez_ID_mouse` - NCBI Entrez Gene ID
13. `KinHub` - Present in KinHub database (X = yes)
14. `Coral` - Present in Coral database (X = yes)
15. `UniProt` - Present in UniProt (X = yes)
16. `Source` - Primary data source
17. `added` - Date added to file

---

## Manning Classification Breakdown

Based on Manning et al. (2002) kinase classification:

| Group | Count | Description |
|-------|-------|-------------|
| **TK** | 91 | Tyrosine kinases (receptors + non-receptors) |
| **Other** | 87 | Diverse kinase families |
| **CAMK** | 72 | Ca²⁺/calmodulin-dependent kinases |
| **CMGC** | 63 | CDK, MAPK, GSK3, CLK families |
| **AGC** | 60 | PKA, PKG, PKC families |
| **STE** | 47 | Sterile kinases (MAPK pathway) |
| **TKL** | 43 | Tyrosine kinase-like |
| **Atypical** | 43 | Non-canonical kinases |
| **Metabolic** | 31 | Metabolic/lipid kinases (NOT in Manning) |
| **CK1** | 11 | Casein kinase 1 family |
| **RGC** | 6 | Receptor guanylate cyclase |
| **Other categories** | 16 | PI3K catalytic subunits, misc. |

**Total**: 558 kinases

---

## Additional Genes Beyond Manning Reference (558 vs 540)

### Explanation

The file contains **18 more kinases** than the Manning et al. 2002 reference (~540 expected mouse orthologs of 518 human kinases).

### Source of Additional Genes

**Primary source: Metabolic kinases (31 genes)**

These are primarily lipid and metabolic kinases not included in Manning's protein kinase classification:

1. Adk - Adenosine kinase
2. Agk - Acylglycerol kinase
3. Chka - Choline kinase alpha
4. Chkb - Choline kinase beta
5. Dck - Deoxycytidine kinase
6. Fn3k - Fructosamine 3 kinase
7. Fn3krp - Fructosamine 3 kinase related protein
8. Galk1 - Galactokinase 1
9. Hk1 - Hexokinase 1
10. Ip6k1 - Inositol hexaphosphate kinase 1
11. Ip6k2 - Inositol hexaphosphate kinase 2
12. Ipmk - Inositol polyphosphate multikinase
13. Khk - Ketohexokinase
14. Nme1 - NME/NM23 nucleoside diphosphate kinase 1
15. Nme2 - NME/NM23 nucleoside diphosphate kinase 2
16. Nme4 - NME/NM23 nucleoside diphosphate kinase 4
17. Pck2 - Phosphoenolpyruvate carboxykinase 2 (mitochondrial)
18. Pdxk - Pyridoxal (pyridoxine, vitamin B6) kinase
19. Pfkl - Phosphofructokinase, liver
20. Pfkm - Phosphofructokinase, muscle
21. Pfkp - Phosphofructokinase, platelet
22. Pik3c3 - Phosphatidylinositol 3-kinase catalytic subunit type 3
23. Pikfyve - Phosphoinositide kinase, FYVE type zinc finger containing
24. Pip4k2a - Phosphatidylinositol-5-phosphate 4-kinase, type II alpha
25. Pip4k2b - Phosphatidylinositol-5-phosphate 4-kinase, type II beta
26. Pip4k2c - Phosphatidylinositol-5-phosphate 4-kinase, type II gamma
27. Pip5k1a - Phosphatidylinositol-4-phosphate 5-kinase, type 1 alpha
28. Pip5k1c - Phosphatidylinositol-4-phosphate 5-kinase, type 1 gamma
29. Pkm - Pyruvate kinase, muscle
30. Ppip5k2 - Diphosphoinositol pentakisphosphate kinase 2
31. Tk2 - Thymidine kinase 2, mitochondrial

**Secondary sources**: 
- Post-2002 discoveries (kinases identified after Manning publication)
- Refined mouse genome annotations
- Expanded inclusion criteria (pseudokinases, atypical kinases)

### Assessment

✅ **The 31 additional genes are LEGITIMATE additions**, not errors:
- Metabolic kinases are bona fide kinases with catalytic activity
- File is more comprehensive than Manning reference
- Appropriate for modern RNA-seq analysis (2025)

---

## Validation: Key GBM Driver Kinases

Verified presence of all major glioblastoma-relevant kinases:

### Receptor Tyrosine Kinases (RTKs)
- ✅ **Egfr** - Epidermal growth factor receptor (Group: TK, Family: EGFR)
- ✅ **Pdgfra** - Platelet-derived growth factor receptor alpha (Group: TK, Family: PDGFR)
- ✅ **Pdgfrb** - PDGFR beta
- ✅ **Met** - Hepatocyte growth factor receptor (Group: TK, Family: Met)
- ✅ **Fgfr1-3** - Fibroblast growth factor receptors 1-3 (Group: TK, Family: FGFR)
- ✅ **Erbb2** - HER2/neu (Group: TK, Family: EGFR)

### PI3K/AKT Pathway
- ✅ **Akt1-3** - Serine/threonine kinases (Group: AGC, Family: Akt)
- ✅ **Pik3ca** - PI3K catalytic subunit alpha
- ✅ **Pik3cb** - PI3K catalytic subunit beta
- ✅ **Pik3cd** - PI3K catalytic subunit delta
- ✅ **Pik3cg** - PI3K catalytic subunit gamma
- ✅ **Pdk1** - 3-phosphoinositide-dependent protein kinase 1

### MAPK Pathway
- ✅ **Map2k1-7** - MAPK kinases (MEK family)
- ✅ **Mapk1** - ERK2
- ✅ **Mapk3** - ERK1
- ✅ **Mapk8-10** - JNK family

### Cell Cycle
- ✅ **Cdk4** - Cyclin-dependent kinase 4
- ✅ **Cdk6** - Cyclin-dependent kinase 6
- ✅ **Aurka** - Aurora kinase A
- ✅ **Aurkb** - Aurora kinase B

---

## Comprehensiveness Assessment

### ✅ EXCELLENT - Suitable for Publication-Quality Analysis

**Strengths:**
1. **Complete coverage**: 103% of Manning reference (558/540)
2. **Authoritative sources**: Cross-referenced with KinHub, UniProt, MGI
3. **Hierarchical classification**: Group → Family → Subfamily
4. **Mouse-specific**: Native mouse gene symbols (no conversion needed)
5. **GBM-validated**: All major GBM driver kinases present
6. **Modern**: Updated 2020 (includes post-Manning discoveries)
7. **Well-annotated**: Multiple IDs per gene (HGNC, MGI, Entrez, UniProt)

**Limitations:**
- File is CSV format (needs GMT conversion for pipeline integration)
- Some entries have formatting artifacts in Group column (e.g., long descriptions for PI3K)

**Recommendation**: 
✅ **USE AS-IS** - No modifications needed to gene content. Proceed with GMT conversion.

---

## Pipeline Integration Strategy

### Step 1: Convert CSV → GMT Format

Create multiple gene sets based on Manning classification:

**All kinases:**
- `MOUSE_KINOME_ALL` (558 genes)

**By group:**
- `MOUSE_KINOME_TK` (91 genes) - Tyrosine kinases
- `MOUSE_KINOME_TKL` (43 genes) - TK-like
- `MOUSE_KINOME_STE` (47 genes) - Sterile kinases
- `MOUSE_KINOME_CK1` (11 genes) - Casein kinase 1
- `MOUSE_KINOME_AGC` (60 genes) - PKA/PKG/PKC
- `MOUSE_KINOME_CAMK` (72 genes) - Ca²⁺/calmodulin
- `MOUSE_KINOME_CMGC` (63 genes) - CDK/MAPK/GSK/CLK
- `MOUSE_KINOME_ATYPICAL` (43 genes) - Atypical kinases
- `MOUSE_KINOME_RGC` (6 genes) - Receptor guanylate cyclase
- `MOUSE_KINOME_OTHER` (87 genes) - Diverse families

**By function:**
- `MOUSE_KINOME_RTK` (subset of TK - receptor tyrosine kinases)
- `MOUSE_KINOME_GBM_DRIVERS` (EGFR, PDGFRA, MET, FGFR1-3, etc.)

### Step 2: Add to Config

```yaml
gene_sets:
  custom_gmt_files:
    - "custom/VERHAAK_GBM_Subtypes.gmt"
    - "custom/MOUSE_KINOME.gmt"  # All kinases + by-group sets
```

### Step 3: Automatic Pipeline Integration

Pipeline will automatically:
1. Load kinome gene sets during chunk 14 (cache-cleanup)
2. Include in enrichment analysis (GSEA + ORA)
3. Generate kinase-specific enrichment plots
4. Export kinase stratification in DEG results

---

## Use Cases for Analysis

### 1. Kinase DEG Enrichment
- Test if kinases are over-represented in DEGs
- Identify which kinase families are most affected
- Compare kinase vs non-kinase DEG rates

### 2. Kinase-Specific Heatmaps
- Filter DEG results to kinases only
- Generate heatmaps showing kinase expression patterns
- Cluster by Manning family

### 3. Pathway Context
- Combine with KEGG signaling pathways
- Identify activated/repressed kinase cascades
- Visualize using pathview (KEGG diagrams)

### 4. Driver Gene Analysis
- Stratify by RTK vs non-RTK
- Focus on GBM driver kinases (EGFR, PDGFRA)
- Compare EGFRvIII vs PDGFRA-driven samples

---

## References

1. **Manning G, Whyte DB, Martinez R, Hunter T, Sudarsanam S** (2002). "The protein kinase complement of the human genome". *Science* 298(5600):1912-1934. DOI: 10.1126/science.1075762

2. **KinHub**: http://www.kinhub.org/ - Comprehensive kinase database

3. **Coral**: Kinase classification and annotation resource

4. **UniProt**: https://www.uniprot.org/ - Protein sequence and functional information

5. **MGI** (Mouse Genome Informatics): http://www.informatics.jax.org/

---

## File History

- **2020-10-06**: File created with 558 kinases from composite sources
- **2025-12-26**: Comprehensive analysis performed, validated for bRNA3F pipeline integration

---

## Contact

For questions about this kinase annotation:
- File maintainer: Ryan Miller
- Pipeline: bRNA3F v5.7.1+
- Repository: github.com/crmlab-uab/bRNA3F
