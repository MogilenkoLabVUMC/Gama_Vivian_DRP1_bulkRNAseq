# Bump Chart Variants

## Purpose

This folder contains **non-key bump chart variants** - alternative visualizations used
during analysis but not selected as primary paper figures. The key figures remain in
the parent folder.

## File Naming Convention

```
bump_{scope}_{style}_{metric}.{ext}
```

- **scope**: `focused` (MitoCarta+SynGO) or `significant` (all meaningful patterns)
- **style**: `uniform`, `weighted`, `highlight`, `curved`, `FINAL_paper_combined`
- **metric**: `nes` (effect size) or `rank` (relative position)
- **ext**: `pdf`, `png`, or `html`

## Variant Types

### Uniform (`bump_*_uniform_*`)
All lines same width and opacity - baseline visualization without weighting.

### Weighted (`bump_*_weighted_*`)
Line width inversely proportional to pattern frequency:
- Dominant patterns (>30%): Thin, low opacity (background)
- Rare patterns (<1%): Thick, high opacity (foreground)

### Highlight (`bump_*_highlight_*`)
Weighted lines with key pathway labels added.

### Curved (`bump_*_curved_*`)
Lines rendered as Bezier curves showing TrajDev magnitude/direction.
**Note:** Key curved figures are in parent folder.

### FINAL Paper Combined
Weighted + Curved + Highlights combined.
**Note:** Focused version is in parent folder.

## Interactive Variants

- `interactive_bump_*.html` - Plotly-based interactive versions
- Toggle patterns, filter by significance, hover for details

## Key Files (in Parent Folder)

The following key files remain in the parent `Trajectory_Flow/` folder:
- `interactive_bump_dashboard.html` - Main interactive explorer
- `bump_focused_curved_nes.pdf/png` - Focused modules with curves
- `bump_focused_FINAL_paper_combined.pdf/png` - Final paper figure
- `bump_curved_nes_significant.pdf/png` - All significant with curves

## Scripts

| Script | Description |
|--------|-------------|
| `02_Analysis/3.7.viz_bump_chart.py` | Static bump chart generator |
| `02_Analysis/3.8.viz_interactive_bump.py` | Interactive variant generator |

## Running

```bash
python3 02_Analysis/3.7.viz_bump_chart.py
python3 02_Analysis/3.8.viz_interactive_bump.py
```

---
Generated: 2024-12-04
