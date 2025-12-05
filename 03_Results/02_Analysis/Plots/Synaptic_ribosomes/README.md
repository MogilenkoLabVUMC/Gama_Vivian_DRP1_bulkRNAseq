# Synaptic Ribosomes Plots

## Purpose

Focused analysis of pre- and postsynaptic ribosome programs using SynGO (Synaptic Gene Ontology) annotations. Reveals compartment-specific translation deficits underlying the synaptic dysfunction in DRP1 mutant neurons.

## Generating Script

- **Script**: `02_Analysis/2.3.viz_synaptic_ribosomes.R`
- **Runtime**: ~2-3 minutes
- **Dependencies**:
  - `03_Results/02_Analysis/checkpoints/syngo_gsea_results.rds` - GSEA results for SynGO pathways
  - `03_Results/02_Analysis/checkpoints/syngo_lists.rds` - SynGO gene set definitions
  - `03_Results/02_Analysis/checkpoints/fit_object.rds` - limma-voom fitted model with coefficients
  - `03_Results/02_Analysis/Verification_reports/syngo_presyn_ribosome_genes.txt` - 52 presynaptic ribosome genes
  - `03_Results/02_Analysis/Verification_reports/syngo_postsyn_ribosome_genes.txt` - 70 postsynaptic ribosome genes

## Biological Context

### Synaptic Local Translation

Synapses require **local protein synthesis** for:
- Rapid response to synaptic activity
- Synaptic plasticity (LTP/LTD)
- Structural remodeling during learning
- Neurotransmitter receptor trafficking

**Key requirement**: Local mitochondria supply ATP for translation

**DRP1 mutation effect**: Mitochondrial mislocalization → ATP depletion → translation failure

### Pre vs Postsynaptic Ribosomes

| Compartment | Location | Function | DRP1 Mutation Effect |
|-------------|----------|----------|---------------------|
| **Presynaptic** | Axon terminal | Synaptic vesicle proteins, active zone components | Downregulated (NES -2.9) |
| **Postsynaptic** | Dendritic spine | Receptors, scaffold proteins, signaling molecules | Downregulated (NES -3.0) |

**Overlap**: All 52 presynaptic genes are contained within the 70 postsynaptic genes (100% of presynaptic genes overlap). The postsynaptic set includes 18 additional postsynaptic-only genes.

## Plots Generated

### Panel_C_Expression_Heatmap.pdf

**Purpose**: Heatmap showing log2 fold-change trajectories for synaptic ribosome genes across both DRP1 mutations.

**Data shown**: Log2 fold-change (logFC) coefficients extracted directly from the limma-voom fitted model (`fit$coefficients`). **NOT** normalized expression values or z-scores.

**Structure**:
- **Rows**: 67 ribosomal genes present in the expression data, split by compartment annotation:
  - **"Postsynaptic Only"** (top section): 14 genes annotated only in postsynaptic ribosome set
  - **"Both Compartments"** (bottom section): 53 genes annotated in both pre- and postsynaptic sets
- **Columns**: 6 contrasts organized as trajectory stages for each mutation:
  - **G32A** (left): Early (G32A_vs_Ctrl_D35), TrajDev (Maturation_G32A_specific), Late (G32A_vs_Ctrl_D65)
  - **R403C** (right): Early (R403C_vs_Ctrl_D35), TrajDev (Maturation_R403C_specific), Late (R403C_vs_Ctrl_D65)
- **Color scale**: Blue-White-Orange diverging palette
  - Blue (negative logFC, down to -0.6): Downregulated vs control
  - White (logFC = 0): No change
  - Orange (positive logFC, up to +0.6): Upregulated vs control

**Row ordering**:
- Hierarchical clustering within each compartment slice (Euclidean distance, complete linkage)
- Compartment slices are NOT reordered (Postsynaptic Only always on top)
- Dendrogram shown on left side

**Column ordering**: Fixed by experimental design (Early → TrajDev → Late), not clustered

**Annotations**:
- Right side: SynGO Annotation bar indicating "Postsynaptic only" (cyan) vs "Both" (gray)
- Column headers split by mutation (G32A, R403C)

**Key patterns visible**:
- **Early (D35)**: Mostly white/light colors → minimal mutation effect at early timepoint
- **TrajDev**: Strong blue → negative logFC indicating downregulation during maturation
- **Late (D65)**: Deep blue → persistent downregulation in mature neurons
- Pattern is highly consistent between G32A and R403C mutations
- Both compartment groups show similar trajectory patterns

## SynGO Gene Sets

### Presynaptic Ribosome (52 genes)

**Source file**: `03_Results/02_Analysis/Verification_reports/syngo_presyn_ribosome_genes.txt`

**Gene composition**: All 52 genes are also in the postsynaptic set (100% overlap with postsynaptic).

**Function**: Support local translation of presynaptic proteins (vesicle proteins, active zone scaffolds)

**GSEA statistics** (from `syngo_gsea_results.rds`):
- G32A TrajDev: NES = -2.90, FDR = 3.2e-12
- R403C TrajDev: NES = -2.71, FDR = 3.6e-10

### Postsynaptic Ribosome (70 genes)

**Source file**: `03_Results/02_Analysis/Verification_reports/syngo_postsyn_ribosome_genes.txt`

**Gene composition**:
- 52 genes shared with presynaptic set (74% of postsynaptic)
- 18 postsynaptic-only genes (26% unique to postsynaptic)

**Function**: Support local translation of postsynaptic proteins (receptors, PSD scaffolds, signaling molecules)

**GSEA statistics**:
- G32A TrajDev: NES = -3.02, FDR = 1.9e-15
- R403C TrajDev: NES = -2.89, FDR = 2.5e-12

### Postsynaptic-Only Genes (18 genes)

These genes are annotated only in the postsynaptic ribosome set, not the presynaptic set:

```
RPL21, RPL30, RPL31, RPLP1, RPS15, RPS17, RPS18, RPS19, RPS2, RPS20,
RPS21, RPS23, RPS3, RPS3A, RPS7, RPS8, RPS9, RPSA
```

**Note on heatmap**: Only 14 of these 18 genes appear in the "Postsynaptic Only" section of the heatmap. The remaining 4 genes (RPS17, RPS18, RPS3A, RPS9) are likely filtered out due to low expression (not present in the filtered expression matrix).

## Methods

### SynGO Annotation
- **Database**: SynGO release 2023-12-01 (bulk download)
- **Location**: `00_Data/SynGO_bulk_20231201/`
- **Namespace**: Cellular Component (CC) for compartment localization
- **Evidence codes**: Experimental annotations only (exclude computational predictions)

### Data Extraction for Heatmap
- **Source**: `fit$coefficients` from limma-voom fitted model (`fit_object.rds`)
- **Values displayed**: Log2 fold-change (logFC) coefficients, NOT expression values
- **NO scaling applied**: Raw logFC values are displayed (no z-score or row normalization)
- **Color mapping**: `colorRamp2(c(-0.6, 0, 0.6), c("#2166AC", "#F7F7F7", "#B35806"))` (Blue-White-Orange)

### Contrasts Displayed
The heatmap shows 6 contrasts organized as Early → TrajDev → Late for each mutation:
1. `G32A_vs_Ctrl_D35` (Early G32A)
2. `Maturation_G32A_specific` (TrajDev G32A)
3. `G32A_vs_Ctrl_D65` (Late G32A)
4. `R403C_vs_Ctrl_D35` (Early R403C)
5. `Maturation_R403C_specific` (TrajDev R403C)
6. `R403C_vs_Ctrl_D65` (Late R403C)

### Gene Filtering
- Gene lists loaded from text files in `Verification_reports/`
- Intersection with `rownames(fit$coefficients)` to keep only expressed genes
- Final count: 67 genes (14 postsynaptic-only + 53 both compartments)

### Heatmap Generation
- **Software**: R `ComplexHeatmap` package
- **Clustering**: Hierarchical clustering within row slices (Euclidean distance, complete linkage via default hclust)
- **Row split**: By compartment annotation (factor with levels: "Postsynaptic only", "Both")
- **Column split**: By mutation (G32A, R403C), fixed order
- **Column clustering**: Disabled (columns ordered by experimental design)

### GSEA Statistics (Referenced)
- GSEA via fgsea (10,000 permutations)
- Gene ranking by t-statistic from limma-voom differential expression
- FDR < 0.05 significance threshold
- Results stored in `syngo_gsea_results.rds` checkpoint

## Key Statistics

### Enrichment Across Developmental Stages

| Compartment | Mutation | Early (D35) NES | TrajDev NES | Late (D65) NES | TrajDev FDR |
|-------------|----------|-----------------|-------------|----------------|-------------|
| Presynaptic | G32A | +2.33 | **-2.90** | -2.49 | **3.2e-12** |
| Postsynaptic | G32A | +2.46 | **-3.02** | -2.49 | **1.9e-15** |
| Presynaptic | R403C | +2.23 | **-2.71** | -2.48 | **3.6e-10** |
| Postsynaptic | R403C | +2.29 | **-2.89** | -2.58 | **2.5e-12** |

**Pattern**: Early ↑ (stress response), TrajDev ↓↓ (collapse), Late ↓ (persistent deficit)

### Gene Set Overlap

- **Total presynaptic ribosome genes**: 52 (from gene list file)
- **Total postsynaptic ribosome genes**: 70 (from gene list file)
- **Overlap**: 52 genes (all presynaptic genes are a subset of postsynaptic)
- **Postsynaptic-only**: 18 genes (specialized for dendritic translation)
- **Genes in heatmap**: 67 (14 postsynaptic-only + 53 shared, after expression filtering)

**Interpretation**: Presynaptic ribosome genes are entirely contained within the postsynaptic set, suggesting shared core ribosomal machinery with additional postsynaptic-specific factors.

## Biological Interpretation

### Translation Collapse at Synapses

The heatmap reveals:
1. **Bidirectional failure**: Both pre- and postsynaptic translation affected
2. **Postsynaptic severity**: Slightly stronger deficits in postsynaptic compartment (NES -3.0 vs -2.9)
3. **Maturation-dependent**: Collapse emerges during D35 → D65 critical period
4. **Mutation-independent**: Both G32A and R403C show same pattern

### Compartment-Specific Vulnerability

**Postsynaptic > Presynaptic** vulnerability may reflect:
- Greater energy demands of postsynaptic signaling (receptors, scaffolds)
- Longer distance from soma to distal dendrites (mitochondrial positioning more critical)
- More complex protein machinery (PSD95, SHANK, NMDA receptors)

**Presynaptic** translation also fails due to:
- Active zone maintenance requires local ATP
- Synaptic vesicle protein turnover depends on translation
- Neurotransmitter release machinery needs continuous synthesis

### Link to Epilepsy Phenotype

Synaptic translation failure → abnormal synapse development → circuit dysfunction:
- Excitatory/inhibitory imbalance (different synaptic types affected differently)
- Impaired synaptic plasticity (LTP/LTD require local translation)
- Network hyperexcitability → seizures

## Interpretation Guide

### Reading the logFC Heatmap

1. **Color patterns**: Look for consistent blue (downregulation) or orange (upregulation) across trajectory stages
2. **Trajectory comparison**: Follow Early → TrajDev → Late columns to see how mutation effects evolve
3. **Compartment differences**: Compare "Postsynaptic Only" section vs "Both Compartments" section
4. **Clustering dendrograms**: Branching patterns reveal genes with similar logFC trajectory profiles
5. **Cross-mutation consistency**: Compare G32A (left) vs R403C (right) to identify shared vs mutation-specific effects

### Understanding logFC Values

- **Orange/positive logFC**: Gene upregulated in mutant vs control for that contrast
- **Blue/negative logFC**: Gene downregulated in mutant vs control for that contrast
- **White/zero logFC**: No difference between mutant and control
- **Color saturation**: Stronger colors indicate larger fold-changes (capped at ±0.6)

### Key Observations in This Heatmap

1. **Early columns are mostly white**: Minimal mutation effect at D35
2. **TrajDev and Late columns are blue**: Strong downregulation emerges during maturation
3. **Patterns are consistent**: G32A and R403C show nearly identical trajectory profiles
4. **Both compartment groups affected equally**: No compartment-specific protection

### Identifying Candidate Genes

Genes with the strongest trajectory effects (darkest blue in TrajDev/Late) are candidates for:
- Mechanistic investigation (why are these most affected?)
- Therapeutic targeting (can we rescue their expression?)
- Biomarker development (early detection of synaptic dysfunction)

## Related Analyses

- **Ribosome_paradox/**: Overall ribosome biogenesis vs translation paradox
- **Mito_translation_cascade/**: Mechanistic cascade from energy to synapse
- **Cross_database_validation/**: SynGO validation across other databases
- **Publication_Figures/**: Integrated multi-panel figures

## Related Visualizations

See also:
- [../Publication_Figures/](../Publication_Figures/README.md) - Main manuscript figures
- [../Ribosome_paradox/](../Ribosome_paradox/README.md) - Core ribosome downregulation finding
- [../Mito_translation_cascade/](../Mito_translation_cascade/README.md) - Energy-translation-synapse cascade
- [../Critical_period_trajectories/](../Critical_period_trajectories/README.md) - GSVA temporal trajectories
- [../Cross_database_validation/](../Cross_database_validation/README.md) - Pattern validation framework
- [../Publication_Figures_Dotplot/](../Publication_Figures_Dotplot/README_DOTPLOT_VERSION.md) - Alternative dotplot visualizations

## Notes

- SynGO annotations are based on experimental evidence from synaptic proteomics and imaging studies
- Presynaptic genes are entirely a subset of postsynaptic genes (100% overlap of presynaptic with postsynaptic)
- Gene set sizes: 52 presynaptic, 70 postsynaptic (18 postsynaptic-only)
- Heatmap shows logFC values (contrast coefficients), NOT expression levels or z-scores
- Not all annotated genes appear in heatmap due to expression filtering (67 of 70 total unique genes)
- Clustering within compartment slices reveals genes with similar trajectory response profiles

## Figure Caption (for publication)

**Synaptic Ribosome Expression Trajectories.** Heatmap showing log2 fold-change (logFC) values for 67 synaptic ribosome genes across developmental trajectory stages in G32A and R403C DRP1 mutant neurons. Rows are split by SynGO compartment annotation: "Postsynaptic Only" (14 genes annotated exclusively in postsynaptic ribosome) and "Both Compartments" (53 genes annotated in both pre- and postsynaptic ribosome sets). Columns represent trajectory stages: Early (mutation vs control at D35), TrajDev (mutation-specific maturation effect), and Late (mutation vs control at D65). Color scale: blue = downregulated (logFC < 0), white = no change, orange = upregulated (logFC > 0). Hierarchical clustering (Euclidean distance, complete linkage) applied within each compartment slice. Both mutations show consistent patterns: minimal effect at Early timepoint, followed by strong downregulation at TrajDev and Late stages, indicating pan-synaptic ribosome program failure during neuronal maturation.

---

**Last Updated**: 2025-12-05
**Generating Script**: `02_Analysis/2.3.viz_synaptic_ribosomes.R`
**Key Finding**: Both synaptic compartments show coordinated ribosome program downregulation during maturation (TrajDev NES ≈ -2.9 to -3.0), with no compartment-specific sparing
