# RNA-seq Analysis Verification Report

**Generated:**  2025-11-25 03:46:28 

---

## 1. Calcium Pathway Enrichment Verification

**Total calcium-related pathway hits:** 1134

**Significant (FDR < 0.05):** 92

### Summary by Contrast

| Contrast | Database | N Pathways | Avg NES | Best FDR |
|----------|----------|------------|---------|----------|
| G32A_vs_Ctrl_D35 | cgp | 5 | -1.40 | 0.0000 |
| G32A_vs_Ctrl_D35 | reactome | 2 | 1.74 | 0.0030 |
| G32A_vs_Ctrl_D35 | gobp | 1 | 1.91 | 0.0191 |
| Maturation_G32A_specific | cgp | 3 | 0.90 | 0.0000 |
| Maturation_G32A_specific | reactome | 1 | -1.52 | 0.0354 |
| Maturation_R403C_specific | cgp | 2 | 2.32 | 0.0000 |
| Maturation_R403C_specific | reactome | 1 | -1.76 | 0.0026 |
| Maturation_R403C_specific | gomf | 1 | 1.76 | 0.0358 |
| Maturation_R403C_specific | wiki | 1 | 1.54 | 0.0404 |
| Maturation_R403C_specific | gobp | 1 | -1.79 | 0.0436 |
| R403C_vs_Ctrl_D35 | cgp | 4 | -1.02 | 0.0000 |
| R403C_vs_Ctrl_D35 | reactome | 3 | 0.76 | 0.0000 |
| R403C_vs_Ctrl_D35 | gobp | 4 | 0.17 | 0.0002 |
| R403C_vs_Ctrl_D35 | wiki | 1 | -1.72 | 0.0057 |
| R403C_vs_Ctrl_D65 | cgp | 2 | 1.54 | 0.0038 |
| Time_Ctrl | cgp | 4 | -2.56 | 0.0000 |
| Time_Ctrl | gobp | 8 | 1.63 | 0.0017 |
| Time_Ctrl | gomf | 2 | 1.65 | 0.0282 |
| Time_G32A | cgp | 3 | -2.43 | 0.0000 |
| Time_G32A | reactome | 2 | 1.99 | 0.0007 |

### Interpretation

✅ **Calcium pathway enrichment is supported by GSEA**

Calcium-related pathways show significant enrichment in the following contexts:

- **Upregulated** (67 pathways): Maturation_R403C_specific, Maturation_G32A_specific, R403C_vs_Ctrl_D35, Time_R403C, Time_G32A, Time_Ctrl, G32A_vs_Ctrl_D35, R403C_vs_Ctrl_D65 
- **Downregulated** (25 pathways): Time_Ctrl, G32A_vs_Ctrl_D35, R403C_vs_Ctrl_D35, Time_G32A, Time_R403C, Maturation_R403C_specific, Maturation_G32A_specific 

---

## 2. Ribosome Pathway Statistics

**Total ribosome pathway enrichments:** 18

**Significant (FDR < 0.05):** 14

### KEY FINDING: Ribosome Pathway Enrichment

| Contrast | Pathway | NES | FDR | Genes |
|----------|---------|-----|-----|-------|
| G32A_vs_Ctrl_D35 | postsynaptic ribosome | 2.46 | 0.0000 | 65 |
| G32A_vs_Ctrl_D35 | presynaptic ribosome | 2.33 | 0.0000 | 51 |
| R403C_vs_Ctrl_D35 | postsynaptic ribosome | 2.29 | 0.0000 | 65 |
| R403C_vs_Ctrl_D35 | presynaptic ribosome | 2.23 | 0.0000 | 51 |
| G32A_vs_Ctrl_D65 | presynaptic ribosome | -2.49 | 0.0000 | 51 |
| G32A_vs_Ctrl_D65 | postsynaptic ribosome | -2.49 | 0.0000 | 65 |
| R403C_vs_Ctrl_D65 | postsynaptic ribosome | -2.58 | 0.0000 | 65 |
| R403C_vs_Ctrl_D65 | presynaptic ribosome | -2.48 | 0.0000 | 51 |
| Time_Ctrl | postsynaptic ribosome | 2.41 | 0.0000 | 65 |
| Time_Ctrl | presynaptic ribosome | 2.42 | 0.0000 | 51 |
| Maturation_G32A_specific | postsynaptic ribosome | -3.02 | 0.0000 | 65 |
| Maturation_G32A_specific | presynaptic ribosome | -2.90 | 0.0000 | 51 |
| Maturation_R403C_specific | postsynaptic ribosome | -2.89 | 0.0000 | 65 |
| Maturation_R403C_specific | presynaptic ribosome | -2.71 | 0.0000 | 51 |

---

## 3. Calcium Gene Differential Expression

**Calcium genes analyzed:** 15
Genes: NNAT, CACNG3, CACNA1C, CACNA1S, ATP2A1, RYR1, MYLK3, CASR, VDR, STIM1, STIM2, ORAI1, CALB1, CALR, PNPO

**Total DE events:** 6

### Most Frequently DE Calcium Genes

- **NNAT**: 4 contrasts
- **PNPO**: 2 contrasts


