# Publication Figures - Dotplot Version

## Overview

This directory contains the **dotplot versions** of the publication figures. The original heatmap versions are preserved in the parent `Publication_Figures/` directory.

## What Changed

### Visualization Method
- **Original**: Heatmaps with colored cells showing NES values
- **Dotplot**: Scatter plots where each cell is represented by a dot

### Dotplot Encoding
Each dot encodes three pieces of information:

1. **Color** = NES (Normalized Enrichment Score)
   - Uses the same Blue-White-Orange diverging colormap
   - Blue = downregulated pathways
   - Orange = upregulated pathways
   - White = no enrichment

2. **Size** = Statistical significance (-log10(padj))
   - Larger dots = more significant (lower padj)
   - Smaller dots = less significant (higher padj)
   - Size range: 20-300 points
   - **Visual size legend included** showing example dots at different significance levels

3. **Edge** = Significance threshold indicator
   - **Black outline** = padj < 0.05 (significant)
   - **Gray outline** = padj ≥ 0.05 (not significant)

### Size Legend
All figures include a visual size legend at the bottom showing:
- Example dots at padj = 0.001, 0.01, 0.05, and 0.1
- Labels explaining significance levels (highly sig, very sig, threshold, not sig)
- Clear indication of which dots have black edges (significant) vs gray edges (not significant)

### Advantages of Dotplot Version
- **More intuitive**: Size immediately shows significance level
- **Better visual clarity**: Dots don't bleed together like heatmap cells
- **Easier to spot patterns**: Large black-outlined dots stand out clearly
- **Publication-ready**: Commonly used in pathway enrichment visualizations

## Generated Files

All files have `_dotplot` suffix to distinguish from original heatmaps:

### Main Figures
- `Fig1_Ribosome_Paradox_dotplot.{pdf,png}` - 3-panel ribosome trajectory comparison
- `Fig2_MitoCarta_Trajectory_Patterns_dotplot.{pdf,png}` - All mitochondrial pathways
- `Fig3b_SynGO_Trajectory_Patterns_dotplot.{pdf,png}` - All synaptic pathways
- `Fig3_Pattern_Classification_Summary.{pdf,png}` - Bar chart (unchanged from original)
- `Fig4_Semantic_Pathway_Overview_dotplot.{pdf,png}` - Comprehensive semantic categories

### Pattern Classification
Figure 3 is identical to the original because it's a bar chart, not a heatmap/dotplot.

## Technical Implementation

### New Components Created
1. **DotplotRenderer** class (`01_Scripts/RNAseq-toolkit/scripts/GSEA/GSEA_plotting_python/dotplot_renderer.py`)
   - Handles dot size calculation based on padj values
   - Manages edge color assignment based on significance
   - Provides consistent rendering across all figures

2. **Modified Script** (`02_Analysis/3.2.publication_figures_dotplot.py`)
   - All plotting functions updated to use scatter plots
   - Data filtering logic unchanged
   - Output directory separate from original

### Data Filtering (Unchanged)
- CGP database excluded (cancer-focused)
- Pathways with ≥3 trajectory data points
- At least one significant result required
- Same semantic categories and exclusion criteria

## How to Regenerate

Run the dotplot script:
```bash
python3 02_Analysis/3.2.publication_figures_dotplot.py
```

This will regenerate all dotplot figures without affecting the original heatmaps.

## Files Created

```
Publication_Figures_Dotplot/
├── Fig1_Ribosome_Paradox_dotplot.pdf (85K)
├── Fig1_Ribosome_Paradox_dotplot.png (460K)
├── Fig2_MitoCarta_Trajectory_Patterns_dotplot.pdf (73K)
├── Fig2_MitoCarta_Trajectory_Patterns_dotplot.png (280K)
├── Fig3b_SynGO_Trajectory_Patterns_dotplot.pdf (67K)
├── Fig3b_SynGO_Trajectory_Patterns_dotplot.png (295K)
├── Fig3_Pattern_Classification_Summary.pdf (49K)
├── Fig3_Pattern_Classification_Summary.png (182K)
├── Fig4_Semantic_Pathway_Overview_dotplot.pdf (152K)
├── Fig4_Semantic_Pathway_Overview_dotplot.png (1.2M)
└── README_DOTPLOT_VERSION.md
```

**Note**: All figures include visual size legends explaining the padj-to-dot-size relationship.

## Original Files (Preserved)

The original heatmap versions remain unchanged in:
```
Publication_Figures/
├── Fig1_Ribosome_Paradox.pdf
├── Fig1_Ribosome_Paradox.png
├── Fig2_MitoCarta_Trajectory_Patterns.pdf
├── Fig2_MitoCarta_Trajectory_Patterns.png
├── Fig3b_SynGO_Trajectory_Patterns.pdf
├── Fig3b_SynGO_Trajectory_Patterns.png
├── Fig3_Pattern_Classification_Summary.pdf
├── Fig3_Pattern_Classification_Summary.png
├── Fig4_Semantic_Pathway_Overview.pdf
└── Fig4_Semantic_Pathway_Overview.png
```

## Comparison Guide

### When to Use Heatmaps (Original)
- When you want to emphasize the gradient of enrichment scores
- For pattern recognition across many pathways
- When color continuity is important

### When to Use Dotplots (New Version)
- When statistical significance is a key message
- For presentations where size is easier to interpret than color intensity
- When you want to highlight specific significant pathways
- Standard in pathway enrichment publications

## Related Visualizations

See also:
- [../Publication_Figures/](../Publication_Figures/README.md) - Original heatmap versions (main manuscript)
- [../Ribosome_paradox/](../Ribosome_paradox/README.md) - Core ribosome downregulation finding
- [../Mito_translation_cascade/](../Mito_translation_cascade/README.md) - Mechanistic cascade visualization
- [../Synaptic_ribosomes/](../Synaptic_ribosomes/README.md) - Synaptic translation deep-dive
- [../Critical_period_trajectories/](../Critical_period_trajectories/README.md) - GSVA temporal trajectories
- [../Cross_database_validation/](../Cross_database_validation/README.md) - Pattern validation framework

## Notes

- Both versions use the same colorblind-safe Blue-White-Orange palette
- Filtering criteria are identical between versions
- Pattern classification (Fig3) is the same in both directories
- All figures maintain the trajectory framework: Early → TrajDev → Late

Generated: 2025-11-25
Last updated: 2025-11-25 (added visual size legends to all figures)
