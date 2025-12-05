# Alluvial Diagrams (Sankey Flows)

## Purpose

Alluvial/Sankey diagrams showing the **FLOW** of pathway dynamics from early defects
through maturation mechanisms to late outcomes.

## Files

### Interactive Alluvial (Plotly - Python)
- `alluvial_binary_G32A.html` - Binary view: Has Early Defect vs No Defect (Interactive)
- `alluvial_binary_R403C.html` - Same for R403C mutation
- `alluvial_graded_G32A.html` - Graded view: Strong/Moderate/No defect
- `alluvial_graded_R403C.html` - Same for R403C mutation
- Corresponding `.png` files for static viewing

### Classical Alluvial (ggalluvial - R)
- `alluvial_ggalluvial_G32A.pdf/png` - Classical alluvial for G32A
- `alluvial_ggalluvial_R403C.pdf/png` - Classical alluvial for R403C
- `alluvial_combined.pdf/png` - Side-by-side comparison

## Structure

```
         Early Status          Mechanism           Late Outcome
         ────────────         ─────────           ────────────
        ┌─────────────┐      ┌─────────┐         ┌──────────┐
        │Early Defect │ ───► │ Active  │ ───────►│ Improved │
        └─────────────┘      │(TrajDev)│         └──────────┘
              │              └─────────┘              ▲
              │                   │                   │
              │              ┌─────────┐              │
              └─────────────►│ Passive │─────────────┘
                             │(Buffer) │
                             └─────────┘

        ┌─────────────┐      ┌──────────┐        ┌───────────┐
        │No Early Def │ ───► │Late-onset│ ──────►│ New Defect│
        └─────────────┘      └──────────┘        └───────────┘
```

- **LEFT**: Early defect status/severity
- **MIDDLE**: Mechanism (Active TrajDev vs Passive buffering)
- **RIGHT**: Late outcome (Improved, Resolved, Worsened)
- **Late_onset**: Separate stream with no left-side input

## Scripts

| Script | Output | Description |
|--------|--------|-------------|
| `02_Analysis/3.5.viz_trajectory_flow.py` | `alluvial_binary_*.html/png`, `alluvial_graded_*.html/png` | Interactive Plotly alluvial |
| `02_Analysis/3.6.viz_alluvial_ggalluvial.R` | `alluvial_ggalluvial_*.pdf/png`, `alluvial_combined.*` | Classical R ggalluvial |

## Running

```bash
# Python (interactive Plotly)
python3 02_Analysis/3.5.viz_trajectory_flow.py

# R (classical ggalluvial)
Rscript 02_Analysis/3.6.viz_alluvial_ggalluvial.R
```

## See Also

- Parent folder README for bump charts and overall documentation
- `docs/PATTERN_CLASSIFICATION.md` for pattern definitions

---
Generated: 2024-12-04
