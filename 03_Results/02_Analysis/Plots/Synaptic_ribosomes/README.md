# Synaptic Ribosomes Plots

## Purpose

Focused analysis of pre- and postsynaptic ribosome programs using SynGO (Synaptic Gene Ontology) annotations. Reveals compartment-specific translation deficits underlying the synaptic dysfunction in DRP1 mutant neurons.

## Generating Script

- **Script**: `02_Analysis/viz_synaptic_ribosomes.R`
- **Runtime**: ~2-3 minutes
- **Dependencies**:
  - `03_Results/02_Analysis/checkpoints/syngo_gsea_results.rds`
  - `03_Results/02_Analysis/checkpoints/fit_object.rds`
  - `03_Results/02_Analysis/checkpoints/gene_intersections.rds` (if available)

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

**Overlap**: 52 presynaptic genes are subset of 70 postsynaptic genes (74% overlap)

## Plots Generated

### Panel_C_Expression_Heatmap.pdf

**Purpose**: Expression heatmap showing pre- and postsynaptic ribosome gene expression across samples

**Structure**:
- **Rows**: Ribosomal genes (grouped by compartment)
  - Presynaptic-specific ribosome genes (52 genes)
  - Postsynaptic-only ribosome genes (18 genes)
  - Shared genes may be marked separately
- **Columns**: Individual samples (grouped by genotype and timepoint)
  - Ctrl_D35, G32A_D35, R403C_D35 (early samples)
  - Ctrl_D65, G32A_D65, R403C_D65 (late samples)
- **Color**: Expression level (log2 normalized counts)
  - Blue/Purple: Low expression
  - Red/Yellow: High expression
  - Color scale typically centered at median or zero

**Annotations**:
- Hierarchical clustering of rows (genes) to reveal co-expression modules
- Columns may be ordered by experimental design or clustered by similarity
- Compartment labels (pre vs post) on left side
- Sample group labels on top

**Key features**:
- Shows gene-level detail underlying pathway enrichment
- Reveals whether pre vs post genes have different expression patterns
- Identifies leading-edge genes (most contributing to enrichment)

**Key findings**:
- Postsynaptic ribosome genes show stronger downregulation than presynaptic
- Expression patterns consistent across both mutations
- Maturation exacerbates the expression deficits

## SynGO Gene Sets

### Presynaptic Ribosome (51-52 genes)

**Annotated genes** (examples):
- RPL* (large ribosomal subunit proteins, e.g., RPL5, RPL7, RPL10)
- RPS* (small ribosomal subunit proteins, e.g., RPS3, RPS6, RPS27)

**Function**: Support local translation of presynaptic proteins (vesicle proteins, active zone scaffolds)

**GSEA statistics**:
- G32A TrajDev: NES = -2.90, FDR = 3.2e-12
- R403C TrajDev: NES = -2.71, FDR = 3.6e-10

### Postsynaptic Ribosome (65-70 genes)

**Annotated genes** (examples):
- All presynaptic ribosome genes (52 genes, 74% overlap)
- Additional postsynaptic-enriched ribosome genes (18 genes, 26% unique)

**Function**: Support local translation of postsynaptic proteins (receptors, PSD scaffolds, signaling molecules)

**GSEA statistics**:
- G32A TrajDev: NES = -3.02, FDR = 1.9e-15
- R403C TrajDev: NES = -2.89, FDR = 2.5e-12

**Postsynaptic-only genes** (18 genes): More specialized for dendritic spine translation

## Methods

**SynGO Annotation**:
- **Database**: SynGO release 2023-12-01 (bulk download)
- **Namespace**: Cellular Component (CC) for compartment localization
- **Evidence codes**: Experimental annotations only (exclude computational predictions)

**Expression Data**:
- **Normalization**: TMM-normalized log2 CPM
- **Source**: limma-voom processed data from main pipeline
- **Filtering**: Genes with mean expression > 1 CPM retained

**Heatmap Generation**:
- **Software**: R with ComplexHeatmap or pheatmap
- **Clustering**: Hierarchical clustering using Euclidean distance and complete linkage
- **Scaling**: Row-wise z-score scaling (mean = 0, SD = 1) for visualization

**Statistical Testing**:
- GSEA via fgsea (10,000 permutations)
- Gene ranking by t-statistic from limma-voom differential expression
- FDR < 0.05 significance threshold

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

- **Total presynaptic ribosome genes**: 51-52
- **Total postsynaptic ribosome genes**: 65-70
- **Overlap**: 52 genes (all presynaptic genes included in postsynaptic set)
- **Postsynaptic-only**: 18 genes (specialized for dendritic translation)

**Interpretation**: Presynaptic ribosome genes are largely a subset of postsynaptic genes, suggesting shared core ribosomal machinery with additional postsynaptic-specific factors.

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

### Reading the Expression Heatmap

1. **Color patterns**: Look for blocks of similar color indicating co-expressed gene modules
2. **Sample groups**: Compare Ctrl vs Mutants at each timepoint to identify when deficits emerge
3. **Compartment differences**: Compare presynaptic vs postsynaptic gene sections
4. **Clustering dendrograms**: Branching patterns reveal which genes/samples are most similar

### Gene Expression Levels

- **High expression (red/yellow)**: Active translation in these samples
- **Low expression (blue/purple)**: Reduced translation capacity
- **Intermediate (white)**: Baseline or moderate expression

### Identifying Leading-Edge Genes

Leading-edge genes (most contributing to GSEA enrichment) often show:
- Consistent direction of change across all mutation samples
- Large fold-changes
- High statistical significance

These are candidate therapeutic targets or biomarkers.

## Related Analyses

- **Ribosome_paradox/**: Overall ribosome biogenesis vs translation paradox
- **Mito_translation_cascade/**: Mechanistic cascade from energy to synapse
- **Cross_database_validation/**: SynGO validation across other databases
- **Publication_Figures/**: Integrated multi-panel figures

## Notes

- SynGO annotations are based on experimental evidence from synaptic proteomics and imaging studies
- Presynaptic and postsynaptic sets have substantial overlap (not mutually exclusive)
- Gene set sizes may vary slightly depending on SynGO release version and filtering
- Expression heatmap shows all genes in the sets, not just differentially expressed genes
- Clustering reveals co-expression patterns that may suggest functional submodules

---

**Last Updated**: 2025-11-25
**Analysis Version**: Post-Session 5 (streamlined to key panels, postsynaptic-specific view added)
**Key Finding**: Both synaptic compartments fail during maturation (NES ≈ -3.0), postsynaptic slightly worse
