# Pseudokinase and Unvalidated Kinase Annotations

**Created**: 2025-12-27 03:05 CST  
**Author**: C. Ryan Miller, MD, PhD (rmiller@uab.edu)

## Known Pseudokinases (Removed from Active Kinome)

### hEGFR
- **Status**: REMOVED from mouse kinome list
- **Reason**: Human gene, not mouse ortholog
- **Action**: Deleted from 201006_composite_kinases_curated.csv

### Mlkl (MLKL)
- **Name**: Mixed lineage kinase domain-like protein
- **Classification**: TKL family pseudokinase
- **Evidence**: Lacks catalytic activity; functions as necroptosis executioner via oligomerization
- **GO Validation**: None (expected - not a kinase)
- **Disposition**: **PSEUDOKINASE** - retain for necroptosis research

### Plk5 (PLK5)
- **Name**: Polo-like kinase 5
- **Classification**: PLK family pseudokinase
- **Evidence**: Lacks kinase activity; functions as tumor suppressor through kinase-independent mechanisms
- **GO Validation**: None (expected - not a kinase)
- **Disposition**: **PSEUDOKINASE** - retain for cell cycle research

### Ptk7 (CCK4)
- **Name**: Inactive tyrosine-protein kinase 7
- **Classification**: TK family pseudokinase (explicitly labeled "Inactive")
- **Evidence**: Lacks catalytic activity; functions as co-receptor in Wnt signaling
- **GO Validation**: None (expected - not a kinase)
- **Disposition**: **PSEUDOKINASE** - retain for Wnt pathway research

### Trib1, Trib2, Trib3 (Tribbles homologs)
- **Name**: Tribbles pseudokinase family
- **Classification**: CAMK family pseudokinases
- **Evidence**: Lack ATP-binding capacity; function as adaptor proteins for E3 ubiquitin ligases
- **GO Validation**: None (expected - not kinases)
- **Disposition**: **PSEUDOKINASES** - retain for signaling scaffold research

## Kinases Requiring Manual Review

### Abr (ABR)
- **Name**: Active breakpoint cluster region-related protein
- **Classification**: Atypical/BCR family
- **GO Validation**: FAILED (0 kinase GO terms)
- **UniProt**: Q12979
- **Notes**: Related to BCR (breakpoint cluster region protein); may have RhoGEF activity but unclear kinase function
- **Action Needed**: Literature review for kinase activity evidence

### Cilk1 (ICK)
- **Name**: Serine/threonine-protein kinase ICK (Intestinal cell kinase)
- **Classification**: CMGC/RCK family
- **GO Validation**: FAILED (0 kinase GO terms)
- **Missing Data**: No Entrez ID in dataset
- **Notes**: Should have kinase activity (CMGC family member); validation failure may be data issue
- **Action Needed**: Re-query with correct IDs

### Stk19 (G11)
- **Name**: Serine/threonine-protein kinase 19
- **Classification**: Atypical/G11 family
- **GO Validation**: FAILED (0 kinase GO terms)
- **UniProt**: P49842
- **Notes**: Nuclear kinase; literature suggests functional kinase activity
- **Action Needed**: Manual GO term check or literature validation

### Tbck (TBCK)
- **Name**: TBC domain-containing protein kinase-like protein
- **Classification**: Other/TBCK family
- **GO Validation**: FAILED (0 kinase GO terms)
- **UniProt**: Q8TEA7
- **Notes**: Name suggests "kinase-like" - may be pseudokinase or misannotated
- **Action Needed**: Literature review for catalytic activity

### Trrap (TRRAP)
- **Name**: Transformation/transcription domain-associated protein
- **Classification**: Atypical/PIKK family
- **GO Validation**: FAILED (0 kinase GO terms)
- **UniProt**: Q9Y4A5
- **Notes**: PIKK family member; known pseudokinase - lacks kinase activity, functions in chromatin remodeling
- **Action Needed**: **Likely PSEUDOKINASE** - confirm and flag

## Summary Statistics

**Original list**: 540 kinases  
**hEGFR removed**: 539 remaining  
**Confirmed pseudokinases**: 6 (Mlkl, Plk5, Ptk7, Trib1, Trib2, Trib3)  
**Under review**: 5 (Abr, Cilk1, Stk19, Tbck, Trrap)  
**GO validated active kinases**: 528/539 (98.0%)  

## Recommendations

1. **Add Pseudokinase column** to CSV with values: "ACTIVE", "PSEUDO", "REVIEW"
2. **Flag 6 confirmed pseudokinases** for separate analysis tracks
3. **Manually validate 5 under-review** genes via literature search
4. **Re-run validation** after fixing Cilk1 Entrez ID issue
