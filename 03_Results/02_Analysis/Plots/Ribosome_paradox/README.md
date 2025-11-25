# Ribosome Paradox Plots

## Purpose

This folder contains the core finding of the DRP1 mutation analysis: a **translation paradox** where cytoplasmic ribosome biogenesis increases during maturation while synaptic ribosome programs sharply decline. This paradox suggests an ATP/translation bottleneck at synapses despite increased ribosome production, revealing a fundamental energy-translation coupling failure in DRP1 mutant neurons.

## Generating Scripts

- **R Script**: `02_Analysis/viz_ribosome_paradox.R`
- **Python Script**: `02_Analysis/7.ribosome_upset_plot.py` (generates UpSet plot, if present)
- **Runtime**: ~3-5 minutes
- **Dependencies**:
  - `03_Results/02_Analysis/checkpoints/all_gsea_results.rds` - Complete GSEA results
  - `03_Results/02_Analysis/checkpoints/syngo_gsea_results.rds` - SynGO-specific enrichment
  - GSEA results from main pipeline (`1a.Main_pipeline.R`)

## The Ribosome Paradox

### Core Finding

During neuronal maturation (D35 → D65), DRP1 mutations trigger a compensatory response where:

1. **Cytoplasmic ribosome biogenesis** ↑ (NES ≈ +2.25, FDR < 1e-11)
2. **Presynaptic ribosomes** ↓ (NES ≈ -2.9, FDR < 1e-12)
3. **Postsynaptic ribosomes** ↓ (NES ≈ -3.0, FDR < 1e-12)

### Interpretation

The neuron responds to DRP1-mediated mitochondrial dysfunction by **increasing ribosome production**, but these ribosomes **fail to support synaptic translation** due to:
- Local ATP depletion at synapses (mitochondrial positioning failure)
- Inability to transport/anchor ribosomes to synaptic compartments
- Energy bottleneck preventing translation initiation

This represents a **failed compensation**: the cell attempts to solve a translation problem by making more ribosomes, but the root cause (ATP delivery) remains unresolved.

## Three Ribosomal Pools

The paradox involves three distinct ribosome pools with opposite trajectories:

| Pool | Location | Early (D35) | Maturation (TrajDev) | Late (D65) | Pattern |
|------|----------|-------------|----------------------|------------|---------|
| **1. Cytoplasmic Biogenesis** | Nucleolus/cytoplasm | DOWN (NES: -2.6) | **UP** (NES: +2.25) | Normal | **Compensation** |
| **2. Synaptic Ribosomes (pre/post)** | Synaptic compartments | UP (NES: +2.3-2.4) | **DOWN** (NES: -2.9 to -3.0) | DOWN (NES: -2.5) | **Collapse** |
| **3. Mitochondrial Ribosomes** | Mitochondrial matrix | DOWN (NES: -2.2) | UP (NES: +1.9) | Normal/UP | **Compensation** |

**Key insight**: The cell compensates for ribosome deficits by increasing production (#1 and #3) but fails to maintain functional synaptic translation (#2).

## Plots Generated

### Ribosome_Paradox_Three_Pools.pdf

**Purpose**: Shows the paradoxical trajectories of three ribosomal pools across developmental time

**Plot Structure**:
- **X-axis**: Developmental stage
  - `Early (D35)` = Mutation vs Control at Day 35
  - `TrajDev (Maturation)` = Mutation-specific maturation changes
  - `Late (D65)` = Mutation vs Control at Day 65
- **Y-axis**: Normalized Enrichment Score (NES)
  - Positive = pathway upregulated
  - Negative = pathway downregulated
- **Color/Shape**: Three ribosome pools
  - **Pool 1** (red): Cytoplasmic ribosome biogenesis (GO:BP)
  - **Pool 2** (blue): Synaptic ribosomes - pre/postsynaptic (SynGO)
  - **Pool 3** (green): Mitochondrial ribosomes (MitoCarta)
- **Lines**: Connect same pathway across time points
- **Facets**: Separate panels for G32A and R403C mutations

**Statistical annotations**:
- Points with FDR < 0.05 are solid
- Non-significant points (FDR > 0.05) are hollow
- Error bars represent 95% confidence intervals (if shown)

**Key observations**:
1. **Divergent trajectories**: Pool 1 and 3 go UP during TrajDev, Pool 2 goes DOWN
2. **Strong effect sizes**: |NES| > 2.0 for most pathways indicates robust effects
3. **Mutation consistency**: Both G32A and R403C show the same paradox pattern
4. **Compensation at TrajDev**: The maturation phase (TrajDev) shows the strongest compensatory upregulation

### Ribosome_Temporal_Trajectory.pdf

**Purpose**: Focuses specifically on the temporal trajectory of ribosome biogenesis compensation

**Plot Structure**:
- **X-axis**: Developmental trajectory (Early → TrajDev → Late)
- **Y-axis**: NES for ribosome biogenesis pathway
- **Lines**: Smooth trajectory connecting timepoints
- **Shaded region**: 95% confidence interval or FDR significance threshold
- **Color**: Mutation type (G32A vs R403C)

**Key features**:
- Shows the "V-shape" or "U-shape" pattern: DOWN at Early, UP at TrajDev, Normal at Late
- Highlights the transient nature of the compensatory response
- Demonstrates that compensation is strongest during the critical period (D35-D65)

**Interpretation**:
- Early disruption triggers compensatory response
- Compensation peaks during maturation
- By Late stage, biogenesis normalizes but synaptic translation remains impaired

## Data Files

### Ribosome_Paradox_Data.csv (12 KB)

Complete pathway-level data for all three ribosome pools across all contrasts.

**Key columns**:
- `Pool` - Ribosome pool category (1, 2, or 3)
- `Pathway_Short` - Short pathway name for plotting
- `Description` - Full pathway identifier (e.g., `GOBP_RIBOSOME_BIOGENESIS`)
- `Contrast` - Experimental contrast (G32A_vs_Ctrl_D35, Maturation_G32A_specific, etc.)
- `Mutation` - G32A or R403C
- `Timepoint` - Early (D35), TrajDev (Maturation), or Late (D65)
- `NES` - Normalized Enrichment Score (effect size)
- `p.adjust` - FDR-corrected p-value
- `pvalue` - Raw p-value
- `setSize` - Number of genes in pathway

**Pool definitions** (from CSV data):
- **Pool 1**: Cytoplasmic ribosome biogenesis
  - `GOBP_RIBOSOME_BIOGENESIS` (287 genes)
  - `GOBP_RIBOSOMAL_SMALL_SUBUNIT_BIOGENESIS` (98 genes)
  - `GOBP_RIBOSOMAL_LARGE_SUBUNIT_BIOGENESIS` (62 genes)

- **Pool 2**: Synaptic ribosomes (SynGO)
  - `presynaptic ribosome` (51 genes)
  - `postsynaptic ribosome` (65 genes)

- **Pool 3**: Mitochondrial ribosomes (MitoCarta)
  - `Mitochondrial_central_dogma` (216 genes)
  - `Translation` (145 genes)
  - `Mitochondrial_ribosome` (77 genes)
  - `Mitochondrial_ribosome_assembly` (24 genes)
  - `Translation_factors` (15 genes, less consistent effects)

## Methods

**GSEA Parameters**:
- **Software**: fgsea (R package) via clusterProfiler wrapper
- **Gene ranking**: t-statistic from limma-voom differential expression
- **Databases**:
  - Gene Ontology Biological Process (GO:BP) - cytoplasmic ribosome pathways
  - SynGO - synaptic ribosome compartments
  - MitoCarta 3.0 - mitochondrial ribosome/translation pathways
- **Thresholds**: FDR < 0.05 for significance
- **Permutations**: 10,000
- **Min/max gene set size**: 15-500 genes

**Contrast Framework**:
- **Early (D35)**: Mutation vs Control at Day 35 (baseline mutation effect)
- **TrajDev (Maturation)**: Interaction term capturing mutation-specific developmental changes (D35 → D65)
- **Late (D65)**: Mutation vs Control at Day 65 (mature neuron effect)

**Statistical Model**:
- limma-voom with TMM normalization
- Design matrix includes genotype, timepoint, and genotype×timepoint interaction
- Contrasts extract mutation-specific maturation effects

## Key Statistics

### G32A Mutation

| Pool | Pathway | Early NES | TrajDev NES | Late NES | TrajDev FDR |
|------|---------|-----------|-------------|----------|-------------|
| 1. Cytoplasmic | Ribosome Biogenesis | -2.60 | **+2.25** | -1.22 | **5.5e-12** |
| 2. Synaptic (Pre) | Presynaptic Ribosome | +2.33 | **-2.90** | -2.49 | **3.2e-12** |
| 2. Synaptic (Post) | Postsynaptic Ribosome | +2.46 | **-3.02** | -2.49 | **1.9e-15** |
| 3. Mitochondrial | Mito Ribosome | -2.18 | **+1.89** | -1.43 | **7.9e-4** |

### R403C Mutation

| Pool | Pathway | Early NES | TrajDev NES | Late NES | TrajDev FDR |
|------|---------|-----------|-------------|----------|-------------|
| 1. Cytoplasmic | Ribosome Biogenesis | -1.72 | **+1.59** | +0.61 | **3.7e-3** |
| 2. Synaptic (Pre) | Presynaptic Ribosome | +2.23 | **-2.71** | -2.48 | **3.6e-10** |
| 2. Synaptic (Post) | Postsynaptic Ribosome | +2.29 | **-2.89** | -2.58 | **2.5e-12** |
| 3. Mitochondrial | Mito Ribosome | -1.63 | **+1.75** | +1.05 | **1.7e-2** |

**Common pattern**: Both mutations show TrajDev compensation in Pools 1 & 3, but collapse in Pool 2

## Biological Interpretation

### The Energy-Translation Bottleneck Model

```
DRP1 Mutation
    ↓
Mitochondrial Positioning Failure
    ↓
Synaptic ATP Depletion
    ↓
┌─────────────────────────┬────────────────────────┐
│ COMPENSATION            │ FAILURE                │
│ (Pools 1 & 3)          │ (Pool 2)               │
├─────────────────────────┼────────────────────────┤
│ ↑ Ribosome biogenesis  │ ↓ Synaptic translation │
│ ↑ Mito ribosomes       │ ↓ Pre/postsynaptic     │
│                         │   ribosome programs    │
└─────────────────────────┴────────────────────────┘
    ↓
Translation Crisis at Synapses
    ↓
Synaptic/Network Dysfunction
```

### Why Compensation Fails

1. **Wrong target**: Increasing ribosome production doesn't solve ATP depletion
2. **Compartmentalization**: Newly synthesized ribosomes can't reach or function at ATP-depleted synapses
3. **Energy priority**: Limited ATP diverted to somatic ribosome biogenesis, further starving synapses
4. **Critical period vulnerability**: D35-D65 maturation window requires high synaptic translation for synaptogenesis

### Clinical Relevance

This paradox explains the **developmental epilepsy phenotype**:
- Synapse formation requires local protein synthesis
- ATP depletion prevents synaptic translation
- Compensation fails to restore synaptic function
- Result: Abnormal circuit development → network hyperexcitability → seizures

## Interpretation Guide

### Reading the Three Pools Plot

1. **Look for divergence at TrajDev**: This is where the paradox manifests—some pools go up, others go down
2. **Compare Early vs Late**: Shows whether effects are transient (resolved by Late) or persistent
3. **Focus on Pool 2 (synaptic)**: This is the critical failure point
4. **Assess significance**: Only solid points (FDR < 0.05) represent reliable effects

### Understanding NES Values

- **|NES| < 1.5**: Weak effect, may not be biologically meaningful
- **|NES| = 1.5-2.0**: Moderate effect, likely biologically relevant
- **|NES| > 2.0**: Strong effect, highly biologically significant
- **|NES| > 2.5**: Very strong effect, central to phenotype

### Compensation vs Collapse

- **Compensation pattern**: Early DOWN → TrajDev UP → Late Normal
  - System attempts to correct deficits
  - Seen in Pools 1 & 3

- **Collapse pattern**: Early UP → TrajDev DOWN → Late DOWN
  - Initial stress response followed by failure
  - Seen in Pool 2 (synaptic ribosomes)

## Related Analyses

- **Mito_translation_cascade/**: Mechanistic details of the energy → translation → synapse cascade
- **Synaptic_ribosomes/**: Detailed analysis of pre vs postsynaptic ribosome gene sets
- **Cross_database_validation/**: Validates compensation pattern across 10 independent databases
- **Critical_period_trajectories/**: GSVA-based temporal modeling of pathway dynamics
- **Publication_Figures/**: Integrated multi-panel figures for manuscript

## Notes

- The paradox is strongest in the **TrajDev** (maturation) contrast, indicating this is a developmental phenomenon
- Both mutations (G32A and R403C) show the same paradox pattern, suggesting a common downstream mechanism
- Mitochondrial ribosome compensation (Pool 3) is distinct from cytoplasmic (Pool 1), indicating organelle-specific regulation
- Synaptic translation failure (Pool 2) is spatially and mechanistically different from somatic ribosome biogenesis
- The critical period (D35-D65) corresponds to peak synaptogenesis in cortical neurons, making this window especially vulnerable

---

**Last Updated**: 2025-11-25
**Analysis Version**: Post-cleanup, validated across databases
**Key Finding**: Translation paradox—ribosome biogenesis ↑, synaptic translation ↓
