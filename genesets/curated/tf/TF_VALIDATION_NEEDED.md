# Transcription Factors - Validation Status

**File**: 221003_TranscriptionFactors_v_1.01.csv  
**Status**: ⚠️ NOT YET VALIDATED  
**Species**: Human (ENSG IDs, HGNC symbols)  
**Count**: 1,640 transcription factors

## Data Structure
- **Ensembl_ID**: ENSG identifiers (human)
- **HGNC_symbol**: Human gene symbol
- **DBD**: DNA-binding domain classification
- **TF**: Transcription factor name
- **TF_assessment**: Quality/confidence score

## Validation Needed

### 1. Ortholog Conversion (Human→Mouse)
- Use biomaRt or cached ortholog mapping
- Many TFs may lack 1:1 orthologs
- Expected ~70-80% conversion rate

### 2. GO Term Validation
GO terms for transcription factor activity:
- GO:0003700 (DNA-binding transcription factor activity)
- GO:0000981 (DNA-binding transcription factor activity, RNA polymerase II-specific)
- GO:0001228 (DNA-binding transcription factor activity, RNA polymerase II-specific activating)
- GO:0001227 (DNA-binding transcription factor activity, RNA polymerase II-specific repressing)

### 3. UniProt Cross-Reference
- Verify functional annotations
- Check DNA-binding domain classifications

## Recommended Approach

Create `comprehensive_tf_validation.R` similar to kinase validation:
1. Load human TF list
2. Convert to mouse orthologs
3. Query GO.db for TF-specific terms
4. Query UniProt.ws for DNA-binding domains
5. Generate validation report with success rate

## Priority
**Medium** - TFs are valuable for pathway analysis but less critical than kinases for immediate RNA-seq analysis.

---
**Created**: 2025-12-27 02:55 CST  
**Author**: C. Ryan Miller, MD, PhD (rmiller@uab.edu)
