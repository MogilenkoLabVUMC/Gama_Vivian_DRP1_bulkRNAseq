# Publication Figures

## Purpose

Publication-ready integrated figures combining multiple analyses into cohesive multi-panel layouts for manuscript submission. These figures synthesize the key findings from cross-database validation, ribosome paradox, and trajectory pattern analysis.

## Generating Script

- **Script**: `02_Analysis/6.publication_figures.py`
- **Runtime**: ~5-10 minutes
- **Dependencies**:
  - `03_Results/02_Analysis/Python_exports/gsea_trajectory.csv`
  - `03_Results/02_Analysis/checkpoints/syngo_gsea_results.rds`
  - Pattern classification data from earlier analyses

## Figures Generated

### Fig1_Ribosome_Paradox.pdf/.png (Main Figure 1)

**Purpose**: Core ribosome paradox showing opposing trajectories of ribosome biogenesis vs synaptic translation

**Panels**:
- Three-pool trajectory plot (cytoplasmic, synaptic, mitochondrial ribosomes)
- Separate facets for G32A and R403C mutations
- Shows Early (D35) → TrajDev (Maturation) → Late (D65) progression

**Key finding**: Ribosome biogenesis ↑ (NES +2.25) while synaptic ribosomes ↓ (NES -3.0) during maturation

**Formats**: PDF (vector, for publication) and PNG (raster, for presentations)

### Fig1b_Ribosome_Detailed_Breakdown.pdf/.png (Supplementary)

**Purpose**: Detailed breakdown of ribosome pathway components

**Shows**:
- Individual ribosome subpathways (small subunit, large subunit, assembly)
- Separate trajectories for each component
- Statistical significance markers

**Use case**: Supplementary figure providing granular detail on ribosome biogenesis compensation

### Fig1b_Ribosome_Gene_Overlap_UpSet.pdf/.png (Supplementary)

**Purpose**: UpSet plot showing gene set overlaps between ribosome pools

**Plot type**: UpSet plot (alternative to Venn diagram for multi-way intersections)

**Shows**:
- Overlap between cytoplasmic, synaptic (pre/post), and mitochondrial ribosome genes
- Set sizes for each compartment
- Intersection sizes for all combinations

**Key insight**: Identifies compartment-specific vs shared ribosomal genes

### Fig2_MitoCarta_Trajectory_Patterns.pdf/.png (Main Figure 2)

**Purpose**: Temporal trajectory patterns for mitochondrial pathways (MitoCarta database)

**Visualization**:
- Heatmap showing pathway enrichment across developmental stages
- Rows: MitoCarta pathways (grouped by pattern: Compensation, Progressive, etc.)
- Columns: Early, TrajDev, Late for both G32A and R403C
- Color: NES (blue = down, red = up)

**Key finding**: 45 G32A and 29 R403C pathways show compensation pattern in MitoCarta

### Fig3_Pattern_Classification_Summary.pdf/.png (Main Figure 3)

**Purpose**: Summary bar plot of temporal pattern distributions across databases

**Axes**:
- X-axis: Pattern type (Compensation, Progressive, Natural_worsening, etc.)
- Y-axis: Number of pathways
- Colors: G32A (blue) vs R403C (red)
- Facets: One panel per database (MitoCarta, SynGO, GO:BP, Hallmark, etc.)

**Key finding**: Compensation is the dominant pattern across all databases; G32A shows stronger compensation than R403C

### Fig3b_SynGO_Trajectory_Patterns.pdf/.png (Main Figure 3, panel b)

**Purpose**: Temporal trajectory patterns for synaptic pathways (SynGO database)

**Similar to Fig2 but focused on**:
- Synaptic-specific pathways
- Pre vs postsynaptic compartment differences
- Synaptic function categories (transmission, plasticity, organization)

**Key finding**: Synaptic pathways show mixed patterns with substantial compensation alongside collapse

### Fig4_Semantic_Pathway_Overview.pdf/.png (Main Figure 4)

**Purpose**: High-level semantic overview of pathway categories and their enrichment patterns

**Visualization**:
- Large-scale heatmap or network diagram
- Groups pathways by biological function (energy, translation, synapse, signaling, etc.)
- Shows cross-database consistency

**Key finding**: Energy-translation-synapse axis shows coordinated disruption and compensation

### Supplementary_Color_Palette_Reference.pdf/.png

**Purpose**: Reference guide for color schemes used across all figures

**Contents**:
- Color mappings for mutations (G32A, R403C)
- Color scales for NES values
- Pattern category colors
- Database-specific colors
- Accessibility notes (colorblind-safe palettes)

**Use case**: Methods supplement; ensures reproducibility and consistency

## Methods

**Visualization Framework**:
- **Software**: Python 3.x with matplotlib, seaborn, UpSetPlot
- **Color palettes**: Colorblind-safe (viridis, RdBu_r)
- **Figure dimensions**: Optimized for journal submission (typically 7-10 inches wide)
- **Resolution**: PDF (vector) for publication, PNG at 300 dpi for presentations

**Data Integration**:
- Combines GSEA results from 10 databases
- Pattern classifications from trajectory analysis
- Statistical thresholds: FDR < 0.05 for significance

**Export Settings**:
- Dual format: PDF (editable vector graphics) + PNG (high-resolution raster)
- Consistent fonts across all panels (typically Arial or Helvetica)
- Colorblind-accessible palettes where possible

## Figure Organization by Manuscript Section

| Figure | Section | Key Message |
|--------|---------|-------------|
| Fig1 | Results 1 | Ribosome paradox—core finding |
| Fig1b (Supp) | Results 1 | Detailed ribosome component analysis |
| Fig2 | Results 2 | Mitochondrial pathway compensation |
| Fig3 | Results 3 | Cross-database pattern validation |
| Fig3b | Results 3 | Synaptic-specific patterns |
| Fig4 | Discussion | Integrated biological model |

## Interpretation Guide

### Multi-Panel Integration

These figures are designed to build a coherent narrative:
1. **Fig1**: Establishes the paradox (what)
2. **Fig2**: Shows mitochondrial compensation mechanism (how)
3. **Fig3**: Validates findings across databases (robustness)
4. **Fig4**: Integrates into biological model (why)

### Color Coding Consistency

Across all figures:
- **Blue**: Downregulated pathways (negative NES)
- **Red**: Upregulated pathways (positive NES)
- **G32A**: Typically shown in left panels or blue/green in comparative plots
- **R403C**: Typically shown in right panels or orange/red in comparative plots

### Statistical Rigor

All figures include:
- FDR-corrected p-values (threshold: 0.05)
- Effect sizes (NES) for biological magnitude
- Sample sizes in figure legends
- Confidence intervals or error bars where appropriate

## Related Analyses

- **Ribosome_paradox/**: Source data and detailed analysis
- **Cross_database_validation/**: Pattern classification methodology
- **Mito_translation_cascade/**: Mechanistic cascade details
- **Critical_period_trajectories/**: Temporal trajectory modeling

## File Naming Convention

Format: `[Main/Supplementary]_Fig[Number][Panel]_[Description].[pdf/png]`

Examples:
- `Fig1_Ribosome_Paradox.pdf` - Main figure 1
- `Fig1b_Ribosome_Detailed_Breakdown.pdf` - Supplementary to Fig1
- `Supplementary_*` - Supplementary materials

---

**Last Updated**: 2025-11-25
**Status**: Publication-ready figures
**Formats**: PDF (vector) + PNG (300 dpi raster)
