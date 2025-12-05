# Trajectory Flow Visualizations

## Overview

This folder contains visualizations showing **pathway trajectory dynamics** from early (D35)
through maturation to late (D65) timepoints for DRP1 mutations. Two main visualization types:

1. **Alluvial Diagrams** - Flow/Sankey plots showing how early defects branch into outcomes
2. **Bump Charts** - Trajectory plots showing individual pathway NES changes over time

## Folder Structure

```
Trajectory_Flow/
├── README.md                               # This file
├── interactive_bump_dashboard.html         # KEY: Interactive explorer (HTML)
├── bump_focused_curved_nes.pdf/png         # KEY: Paper figure - focused curved
├── bump_focused_FINAL_paper_combined.pdf/png  # KEY: Paper figure - final combined
├── bump_curved_nes_significant.pdf/png     # KEY: Paper figure - all significant curved
├── alluvial/                               # All alluvial/Sankey diagrams
│   ├── README.md
│   ├── alluvial_binary_*.html/png          # Interactive binary flow
│   ├── alluvial_graded_*.html/png          # Interactive graded severity flow
│   ├── alluvial_ggalluvial_*.pdf/png       # Classical ggalluvial plots
│   └── alluvial_combined.pdf/png           # Side-by-side comparison
└── bump/                                   # Non-key bump chart variants
    ├── bump_*_uniform_*.pdf/png            # Uniform line width variants
    ├── bump_*_weighted_*.pdf/png           # Weighted line width variants
    ├── bump_*_highlight_*.pdf/png          # With pathway labels
    └── interactive_bump_*.html             # Interactive variant explorers
```

## Key Paper Figures (Root Folder)

### interactive_bump_dashboard.html
**Interactive pathway trajectory explorer** - Comprehensive dashboard for exploring all
pathway trajectories with filtering, toggling, and tooltips.

**Features:**
- Filter by pattern, database, significance
- Toggle between NES and Rank views
- Enable/disable trajectory curves (TrajDev visualization)
- Color by pattern or NES value
- Search and highlight specific pathways
- Hover for detailed stats (NES, p-values, pattern interpretation)

### bump_focused_curved_nes.pdf/png
**Focused modules (MitoCarta + SynGO) with trajectory curves** - Shows the curvature
representing TrajDev (trajectory deviation) magnitude and direction.

**How to read:**
- X-axis: Early (D35) to Late (D65)
- Y-axis: NES (Normalized Enrichment Score)
- Line curvature: TrajDev magnitude (bulge = strong trajectory effect)
- Curve direction: TrajDev direction (upward bulge = upregulated during maturation)

### bump_focused_FINAL_paper_combined.pdf/png
**Final paper figure** combining weighted lines, trajectory curves, and key pathway labels.
Best representation of the trajectory analysis for publication.

### bump_curved_nes_significant.pdf/png
**All significant pathways** with trajectory curves - comprehensive view of all pathways
showing meaningful patterns (excludes Complex).

## Bump Chart Module

The bump chart visualization system is implemented in `01_Scripts/Python/viz_bump_charts.py`.

### Key Concepts

#### Scopes
- **focused**: MitoCarta + SynGO databases (~100-200 pathways)
- **significant**: All meaningful patterns except Complex (~2-4k pathways)
- **all**: All pathways (~12k) - rarely used due to density

#### Y-Axis Types
- **NES**: Normalized Enrichment Score - shows effect strength and direction
- **Rank**: Relative ordering by NES - shows position changes

#### Visual Modes
- **uniform**: All lines same width (baseline)
- **weighted**: Line width inversely proportional to pattern frequency
  - Dominant (>30%): Thin, low opacity (background)
  - Common (>10%): Medium width
  - Uncommon (>1%): Thicker
  - Rare (<1%): Thickest, full opacity (foreground)

### Layering Algorithm

Lines are rendered in specific order to ensure rare patterns are visible:

1. **Sort patterns** by weight category (dominant → rare)
2. **Within category**, sort by frequency (more frequent first)
3. **Render in order**: Frequent patterns draw first (bottom layer)
4. **Result**: Rare patterns appear on top, visible against the mass

### Filtering Logic

Pathways pass through multiple filters:
1. **Scope filter**: Database membership (focused) or pattern significance (significant)
2. **Pattern filter**: Must have meaningful pattern (not Complex, not Insufficient_data)
3. **NES completeness**: Must have valid NES values for both Early and Late

### Curving Algorithm (TrajDev Visualization)

When `show_curves=True`, lines are rendered as quadratic Bezier curves:

```
P0 = (0, Early_NES)           # Start point
P2 = (1, Late_NES)            # End point
P1 = (0.5, Y_mid + offset)    # Control point (determines curve)

Where:
  Y_mid = (Early_NES + Late_NES) / 2
  offset = TrajDev_NES × curve_strength  (for NES mode)
  offset = -TrajDev_NES × scale × curve_strength  (for Rank mode)
```

**Curve interpretation:**
- Upward bulge: Positive TrajDev (pathway upregulated during maturation)
- Downward bulge: Negative TrajDev (pathway downregulated during maturation)
- Straight line: No significant TrajDev or TrajDev ≈ 0

**Note:** Curves only appear for pathways with significant TrajDev (p.adjust < 0.05).

### Highlight Selection

When `show_highlights=True`, pathways are labeled using priority ranking:

1. **Filter**: Remove irrelevant pathways (viral, bacterial, non-neuronal)
2. **Prioritize**: Semantic categories (Mitochondrial > Synaptic > Ribosomal > Other)
3. **Sort**: By priority, then by |NES change| magnitude
4. **Select**: Top N per pattern (default 5)
5. **Position**: Use adjustText library to prevent label overlap

## Scripts

| Script | Output | Description |
|--------|--------|-------------|
| `3.5.viz_trajectory_flow.py` | `alluvial/*.html/png` | Interactive Plotly alluvial diagrams |
| `3.6.viz_alluvial_ggalluvial.R` | `alluvial/*.pdf/png` | Classical ggalluvial diagrams |
| `3.7.viz_bump_chart.py` | Root + `bump/*` | Static bump charts (all variants) |
| `3.8.viz_interactive_bump.py` | `bump/interactive_*.html` | Interactive bump chart variants |
| `3.8.viz_interactive_bump_dashboard.py` | `interactive_bump_dashboard.html` | Main interactive dashboard |

## Running the Scripts

```bash
# Generate alluvial diagrams (both Python and R)
python3 02_Analysis/3.5.viz_trajectory_flow.py
Rscript 02_Analysis/3.6.viz_alluvial_ggalluvial.R

# Generate all bump chart variants
python3 02_Analysis/3.7.viz_bump_chart.py

# Generate interactive bump charts
python3 02_Analysis/3.8.viz_interactive_bump.py
python3 02_Analysis/3.8.viz_interactive_bump_dashboard.py
```

## Biological Interpretation

### Key Finding: 99.9% of Early Defects Improve

During neuronal maturation (D35 → D65):
- **~60%** improve through **Active Compensation** (significant TrajDev opposing the defect)
- **~40%** improve through **Passive Buffering** (normal developmental processes)
- **~0.1%** worsen (Progressive pattern)
- **46-61 pathways** show **Late_onset** (new defects emerging without early problems)

### Pattern Trajectories on Bump Charts

| Pattern | Early NES | Curve Direction | Late NES |
|---------|-----------|-----------------|----------|
| Compensation | Strong defect | Opposing (bulge toward 0) | Reduced |
| Progressive | Strong defect | Amplifying (bulge away from 0) | Worsened |
| Natural_improvement | Strong defect | Minimal curve | Reduced |
| Natural_worsening | Strong defect | Minimal curve | Worsened |
| Late_onset | Near zero | Variable | New defect |
| Transient | Strong defect | Variable | Resolved |

## Related Documentation

- `docs/PATTERN_CLASSIFICATION.md` - Full pattern classification framework
- `01_Scripts/Python/pattern_definitions.py` - Canonical pattern classification code
- `01_Scripts/Python/viz_bump_charts.py` - Bump chart visualization module
- `03_Results/02_Analysis/master_gsea_table.csv` - Source data

---
Generated: 2024-12-04
