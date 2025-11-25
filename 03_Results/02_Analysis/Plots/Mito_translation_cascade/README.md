# Mitochondrial Translation Cascade Plots

## Purpose

Mechanistic visualization of the energy → translation → synapse cascade that underlies the DRP1 mutation phenotype. Shows how mitochondrial dysfunction propagates through translation systems to disrupt synaptic function.

## Generating Script

- **Script**: `02_Analysis/viz_mito_translation_cascade.R`
- **Runtime**: ~2-3 minutes
- **Dependencies**:
  - `03_Results/02_Analysis/checkpoints/all_gsea_results.rds`
  - `03_Results/02_Analysis/checkpoints/fit_object.rds`
  - GSEA results from main pipeline

## Cascade Model

```
Mitochondrial Dysfunction (DRP1 mutation)
    ↓
Energy Crisis (ATP depletion, OXPHOS disruption)
    ↓
Translation Failure (Ribosome dysfunction, protein synthesis collapse)
    ↓
Synaptic Dysfunction (Synaptogenesis failure, neurotransmission defects)
    ↓
Network Hyperexcitability (Epilepsy phenotype)
```

## Plots Generated

### Mechanistic_Cascade_Heatmap.pdf

**Purpose**: Detailed heatmap showing pathway enrichment across the mechanistic cascade

**Structure**:
- **Rows**: Pathways grouped into functional modules:
  - **Energy Crisis**: OXPHOS, ATP synthesis, mitochondrial membrane potential
  - **Mitochondrial Positioning**: Fission/fusion, transport, anchoring
  - **Translation Systems**: Ribosome biogenesis, translation initiation/elongation, tRNA
  - **Synaptic Function**: Synaptogenesis, neurotransmitter release, postsynaptic signaling
- **Columns**: Developmental contrasts (Early, TrajDev, Late) × Mutations (G32A, R403C)
- **Color**: Normalized Enrichment Score (NES)
  - Blue: Downregulated (pathway suppressed)
  - Red: Upregulated (pathway enhanced)
  - White: Not significant (FDR > 0.05)

**Annotations**:
- Hierarchical clustering of rows (pathways) to reveal co-regulated modules
- Dendrogram shows pathway relationships
- Significant cells (FDR < 0.05) may have asterisks or bold borders

**Key features**:
- **Mechanistic arrows**: Shows causal relationships between modules
- **Row clustering enabled**: Reveals coordinated regulation within modules
- **Large fonts**: Enhanced readability for presentations

**Key findings**:
- Energy crisis pathways downregulated early, compensate during maturation
- Translation systems show paradoxical pattern (biogenesis up, function down)
- Synaptic pathways persistently disrupted across all stages

### Module_Summary_Heatmap.pdf

**Purpose**: Simplified summary heatmap showing module-level (averaged) enrichment

**Structure**:
- **Rows**: Four major modules (Energy, Translation, Synapse, Signaling)
- **Columns**: Same as detailed heatmap (contrasts × mutations)
- **Values**: Average NES across all pathways within each module

**Advantages**:
- Cleaner visualization for high-level overview
- Easier to identify module-level patterns
- Better for presentations where detail isn't needed

**Interpretation**:
- **Energy module**: Shows compensation pattern (Early down, TrajDev up)
- **Translation module**: Mixed signals (biogenesis vs function paradox)
- **Synapse module**: Persistent dysfunction (no compensation)

## Module Definitions

### 1. Energy Crisis Module

**Pathways included**:
- OXPHOS (Complexes I-V)
- ATP synthesis (ATP5F1* genes)
- Mitochondrial membrane potential
- Electron transport chain
- TCA cycle

**Biological interpretation**:
- Measures mitochondrial energetic capacity
- Downregulation indicates ATP depletion
- Compensation indicates adaptive upregulation of energy production

### 2. Mitochondrial Positioning Module (if shown)

**Pathways included**:
- Mitochondrial fission/fusion
- Mitochondrial transport along axons
- Mitochondrial anchoring at synapses

**Biological interpretation**:
- Measures mitochondrial dynamics and localization
- DRP1 is a key fission factor
- Dysfunction causes mitochondrial mislocalization

### 3. Translation Systems Module

**Pathways included**:
- Ribosome biogenesis (cytoplasmic, mitochondrial)
- Translation initiation (eIF factors)
- Translation elongation (eEF factors)
- tRNA charging and processing

**Biological interpretation**:
- Measures protein synthesis capacity
- Paradox: biogenesis up, functional translation down
- Reflects ATP bottleneck preventing translation despite ribosome availability

### 4. Synaptic Function Module

**Pathways included**:
- Synaptogenesis and synaptic assembly
- Presynaptic neurotransmitter release
- Postsynaptic receptor trafficking
- Synaptic vesicle cycle
- Synaptic plasticity (LTP/LTD)

**Biological interpretation**:
- Measures synaptic integrity and function
- Persistent downregulation indicates failure to compensate
- Links molecular deficits to circuit-level phenotype (epilepsy)

## Methods

**Pathway Selection**:
- Curated from GO:BP, KEGG, Reactome, MitoCarta
- Filtered to include only pathways relevant to cascade model
- Manually assigned to modules based on biological function

**Clustering**:
- Hierarchical clustering using Euclidean distance and complete linkage
- Applied to rows (pathways) to reveal co-regulated gene sets
- Dendrogram shows pathway relationships

**Normalization**:
- NES values used directly (already normalized by GSEA)
- Color scale centered at 0 (neutral) for symmetric interpretation

**Visualization**:
- Software: R with ComplexHeatmap package
- Color palette: RdBu (red-white-blue) for diverging scale
- Annotations: Module boundaries marked with gaps or lines

## Key Findings

### Cascade Propagation

The heatmap reveals how dysfunction propagates through the cascade:

1. **Stage 1 (Early)**: Energy crisis initiated
   - OXPHOS and ATP synthesis downregulated
   - Triggers compensatory response

2. **Stage 2 (TrajDev)**: Compensation attempt
   - Energy pathways upregulate (partial rescue)
   - Translation biogenesis increases (failed compensation)
   - Synaptic pathways begin to fail

3. **Stage 3 (Late)**: New equilibrium
   - Energy normalizes (compensation successful)
   - Translation remains paradoxical
   - Synaptic dysfunction persists (compensation failure)

### Module Coupling

**Strong coupling** (co-regulated):
- Energy ↔ Mitochondrial positioning (expected, both depend on DRP1)
- Translation biogenesis ↔ Mitochondrial ribosomes (coordinate upregulation)

**Weak coupling** (decoupled):
- Energy ↔ Synaptic function (energy recovers, synapse doesn't)
- Translation biogenesis ↔ Synaptic translation (paradox)

**Interpretation**: The cascade is not strictly linear; compensation succeeds upstream but fails downstream.

## Interpretation Guide

### Reading the Heatmap

1. **Identify modules**: Look for groups of pathways with similar patterns
2. **Follow the cascade**: Read from top (energy) to bottom (synapse)
3. **Compare mutations**: G32A vs R403C columns reveal mutation-specific effects
4. **Focus on TrajDev**: This column shows the compensatory response

### Color Intensity

- **Dark blue (NES < -2)**: Strong downregulation, major deficit
- **Light blue (NES -1 to -2)**: Moderate downregulation
- **White (NES ≈ 0)**: No change or not significant
- **Light red (NES +1 to +2)**: Moderate upregulation, compensation
- **Dark red (NES > +2)**: Strong upregulation, robust compensation

### Pattern Recognition

- **Horizontal stripes**: Pathways with consistent effects across all conditions (robust)
- **Vertical stripes**: Conditions with similar effects across pathways (mutation-specific)
- **Checkerboard**: Interaction effects (pathway response depends on context)

## Related Analyses

- **Ribosome_paradox/**: Focus on translation module paradox
- **Cross_database_validation/**: Validates modules across databases
- **Synaptic_ribosomes/**: Detailed synaptic translation analysis
- **Publication_Figures/**: Integrated figures for manuscript

## Notes

- Module assignments are based on expert curation, not clustering
- Pathways may belong to multiple modules (e.g., mitochondrial ribosome in both Energy and Translation)
- Row clustering within modules reveals functional sub-modules
- Missing values (white cells) indicate pathways not significant in that contrast (FDR > 0.05)

---

**Last Updated**: 2025-11-25
**Analysis Version**: Post-Session 5 improvements (color conflicts resolved, clustering enabled)
**Key Insight**: Cascade shows sequential failure—energy compensates, translation paradoxical, synapse fails
