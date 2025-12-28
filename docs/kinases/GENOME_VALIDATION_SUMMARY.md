# Genome Validation Summary

**Created**: 2025-12-27 03:20 CST  
**Author**: C. Ryan Miller, MD, PhD (rmiller@uab.edu)  
**Genome**: GENCODE vM37 (Mouse primary assembly)

## Results

**Total sources (kinases)**: 539  
**Found in genome**: 538 (99.8%)  
**Protein-coding**: 538 (100% of found)  
**Not found**: 1

## Fixes Applied

### Cilk1 (ICK) - MGI and Entrez IDs Added
- **Before**: Missing MGI_Mouse and Entrez_Mouse  
- **After**: MGI:1934157, Entrez:269548  
- **Validation**: ✓ ENSMUSG00000009828, protein_coding, chr9

### Cdk11a - Not in Mouse Genome
- **Status**: Human gene CDC2L1/CDK11A has no mouse ortholog  
- **Note**: Mouse has Cdk11b (MGI:88353) but NOT Cdk11a  
- **Action**: Should be REMOVED from mouse source list

## Unvalidated Sources (GO Terms)

From previous validation, 11 sources (kinases) lack GO kinase evidence:

1. **Abr** - In genome ✓, needs literature validation
2. **Cdk11a** - NOT in genome ✗, remove  
3. **Cilk1** - In genome ✓, IDs now fixed, re-run GO validation
4. **Mlkl** - In genome ✓, pseudokinase (necroptosis)
5. **Plk5** - In genome ✓, pseudokinase (tumor suppressor)
6. **Ptk7** - In genome ✓, pseudokinase (Wnt signaling)
7. **Stk19** - In genome ✓, needs validation
8. **Tbck** - In genome ✓, questionable kinase activity
9. **Trib1** - In genome ✓, pseudokinase (adaptor)
10. **Trib2** - In genome ✓, pseudokinase (adaptor)
11. **Trib3** - In genome ✓, pseudokinase (adaptor)
12. **Trrap** - In genome ✓, pseudokinase (PIKK family, chromatin remodeling)

## Recommendations

1. **Remove Cdk11a** (not in mouse genome)
2. **Re-run GO validation** for Cilk1 with correct Entrez ID
3. **Flag 6 confirmed pseudokinases** (Mlkl, Plk5, Ptk7, Trib1/2/3, Trrap)
4. **Manually review 3 questionable** (Abr, Stk19, Tbck)

