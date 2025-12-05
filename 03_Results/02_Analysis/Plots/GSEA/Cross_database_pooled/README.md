# Cross-Database Pooled GSEA Dotplots

## Overview

This folder contains **cross-database pooled dotplots** that visualize GSEA (Gene Set Enrichment Analysis) results by combining significant pathways across multiple pathway databases into a single unified view. This approach provides a comprehensive snapshot of enriched biological processes for each contrast.

**Generating Script:** `02_Analysis/3.9.viz_pooled_dotplots.R`

## What These Plots Show

Each dotplot displays:

| Visual Element | Meaning |
|----------------|---------|
| **X-axis** | Gene Ratio = (Leading Edge genes) / (Set Size) |
| **Y-axis** | Pathway descriptions, sorted by Gene Ratio |
| **Dot Color** | NES (Normalized Enrichment Score): Blue (negative/downregulated) to Orange (positive/upregulated) |
| **Dot Size** | -log10(FDR): Larger dots = more significant |
| **Black Outline** | Pathways with FDR < 0.05 (highly significant) |

### Color Scheme (Colorblind-Safe)

The NES color gradient uses the project's global colorblind-safe palette from `01_Scripts/R_scripts/color_config.R`:

- **Blue (#2166AC)**: Negative NES (downregulated pathways)
- **White (#F7F7F7)**: NES near zero (neutral)
- **Orange (#B35806)**: Positive NES (upregulated pathways)

## Contrasts Included

All 9 contrasts from the differential expression analysis are included:

### D35 Mutation Effects (Early)
- `G32A_vs_Ctrl_D35` - G32A mutation vs Control at Day 35
- `R403C_vs_Ctrl_D35` - R403C mutation vs Control at Day 35

### D65 Mutation Effects (Late)
- `G32A_vs_Ctrl_D65` - G32A mutation vs Control at Day 65
- `R403C_vs_Ctrl_D65` - R403C mutation vs Control at Day 65

### Time/Maturation Effects
- `Time_Ctrl` - Maturation effect in Control (D65 vs D35)
- `Time_G32A` - Maturation effect in G32A mutants
- `Time_R403C` - Maturation effect in R403C mutants

### Interaction Contrasts (Mutation-Specific Maturation)
- `Maturation_G32A_specific` - G32A-specific maturation changes (interaction)
- `Maturation_R403C_specific` - R403C-specific maturation changes (interaction)

## Pathway Pooling Strategy

### How It Works

1. **Extract significant pathways** from GSEA results (FDR < 0.05)
2. **Pool across databases** - combine results from multiple pathway databases
3. **Apply neuronal filtering** - remove non-neuronal tissue-specific pathways
4. **Select top pathways per database** - take top N per database by FDR
5. **Generate unified dotplot** - show all selected pathways in one figure

### Thresholds & Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| FDR cutoff | < 0.05 | Only significant pathways included |
| Significance outline | p.adjust < 0.05 | Black border for highly significant |
| Top N (focused) | 10 per database | For manuscript-ready figures |
| Top N (comprehensive) | 5 per database | For overview with more databases |
| NES limit | Auto-calculated (max 4) | Caps color scale at reasonable bounds |

## Output Structure

### Subdirectories

```
Cross_database_pooled/
├── README.md                    # This file
├── pooled_dotplots_summary.txt  # Legacy summary (deprecated)
├── focused/                     # Manuscript-ready plots
│   └── {contrast}_pooled_focused.pdf
├── comprehensive/               # Full database coverage plots
│   └── {contrast}_pooled_comprehensive.pdf
└── csv_data/                    # Intermediate data files
    ├── {contrast}_raw_comprehensive.csv
    ├── {contrast}_raw_focused.csv
    ├── {contrast}_filtered_comprehensive.csv
    └── {contrast}_filtered_focused.csv
```

### Database Collections

**Focused Version** (4 databases - for publication):
- `kegg` - KEGG pathways (curated metabolic/signaling pathways)
- `reactome` - Reactome pathways (curated biological processes)
- `syngo` - SynGO synaptic ontology (100% neuronal-relevant)
- `mitocarta` - MitoCarta mitochondrial pathways (critical for DRP1 study)

*Rationale: Focused on curated pathway databases plus domain-specific databases (synaptic, mitochondrial) most relevant to iPSC-derived cortical neuron DRP1 mutations.*

**Comprehensive Version** (11 databases - for exploration):
- All focused databases (kegg, reactome, syngo, mitocarta)
- `hallmark` - MSigDB Hallmark gene sets (broad biological processes)
- `gobp` - GO Biological Process
- `gocc` - GO Cellular Component
- `gomf` - GO Molecular Function
- `wiki` - WikiPathways
- `canon` - MSigDB Canonical Pathways
- `tf` - Transcription Factor targets

## Neuronal Pathway Filtering

Since this study uses iPSC-derived cortical neurons, tissue-specific pathways from other organs are filtered out:

### Excluded (Non-Neuronal)
- Pancreatic (PANCREA, BETA_CELL, ISLET)
- Cardiac (CARDIAC, HEART, MYOCARDI)
- Kidney (KIDNEY, RENAL, NEPHRON)
- Intestinal (INTESTIN, COLON, GUT)
- Lung (LUNG, ALVEOL, BRONCH)
- Liver (LIVER, HEPAT, BILE_ACID)
- Breast/Prostate (BREAST, MAMMARY, PROSTAT)
- Other (ADIPOGEN, XENOBIOTIC, COAGULATION)

### Always Included (Neuronal-Relevant)
- Neural (NEURO, NEURAL, SYNAP, AXON, DENDRIT)
- Brain regions (BRAIN, CORTEX, HIPPOCAM, GLIAL)
- Core biology (RIBOSOM, MITOCHOND, APOPTOSIS)
- Cell cycle (E2F_TARGET, G2M_CHECKPOINT, MYC_TARGET)
- All SynGO pathways (automatically included)

## How to Read the Plots

### Interpreting Gene Ratio (X-axis)
- **Higher Gene Ratio** = Larger fraction of pathway genes in the leading edge
- Leading edge genes are those contributing most to the enrichment signal
- Gene Ratio = count(core_enrichment genes) / setSize

### Interpreting NES (Color)
- **Orange/Positive NES**: Pathway genes tend to be upregulated in the condition
- **Blue/Negative NES**: Pathway genes tend to be downregulated
- **Magnitude**: Stronger colors indicate larger effect sizes

### Interpreting Significance (Size & Outline)
- **Larger dots**: More statistically significant (lower FDR)
- **Black outline**: Highly significant (FDR < 0.05)
- **No outline**: Significant but FDR between 0.05-0.10

### Example Interpretation
> "In the G32A_vs_Ctrl_D35 plot, Oxidative Phosphorylation (orange, large dot, black outline) has a positive NES with high gene ratio and significance - indicating robust upregulation of mitochondrial energy production genes in G32A mutants at Day 35."

## Regenerating Plots

```bash
Rscript 02_Analysis/3.9.viz_pooled_dotplots.R
```

### Dependencies
- R packages: `here`, `ggplot2`, `dplyr`, `stringr`
- Required checkpoints:
  - `03_Results/02_Analysis/checkpoints/all_gsea_results.rds` (MSigDB databases)
  - `03_Results/02_Analysis/checkpoints/syngo_gsea_results.rds` (SynGO)
  - `03_Results/02_Analysis/checkpoints/mitocarta_gsea_results.rds` (MitoCarta)
- Helper functions: `01_Scripts/R_scripts/gsea_dotplot_helpers.R`
- Color configuration: `01_Scripts/R_scripts/color_config.R`

## Related Files

- **Color Config**: `01_Scripts/R_scripts/color_config.R` - Global color palette
- **Dotplot Helpers**: `01_Scripts/R_scripts/gsea_dotplot_helpers.R` - Shared functions
- **Main GSEA Pipeline**: `02_Analysis/1.1.main_pipeline.R` - Generates GSEA checkpoints
- **Master GSEA Table**: `03_Results/02_Analysis/master_gsea_table.csv` - All GSEA results

---
*Last updated: 2025-12-05*
*Generated by: `02_Analysis/3.9.viz_pooled_dotplots.R`*
