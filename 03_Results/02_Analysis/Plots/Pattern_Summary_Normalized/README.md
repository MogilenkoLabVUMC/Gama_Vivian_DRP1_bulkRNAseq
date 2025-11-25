# Normalized Pattern Summary Visualizations

## Problem Statement

**Original pattern summary figures are biased by database size:**

| Database | Pathway Count | Visual Impact | Pattern Visibility |
|----------|--------------|---------------|-------------------|
| gobp     | ~5000        | Dominates     | Patterns visible  |
| hallmark | ~10          | Barely visible| Patterns obscured |
| MitoCarta| ~50          | Small         | Hard to compare   |

**Result:** Cannot compare pattern distributions across databases of different sizes.

**Example:**
- "Does hallmark show more compensation *relatively*?" → **Can't tell**
- "Which database has highest % progressive?" → **Can't compare**

## Solution: Normalized Visualizations

### Three Alternative Approaches

#### 1. **Normalized 100% Stacked Bars** ⭐ RECOMMENDED

**File:** `pattern_summary_normalized.pdf`

**Features:**
- Each database bar = 100% (regardless of absolute size)
- Pattern proportions directly comparable
- Absolute counts shown as annotations: "N (X%)"
- Total pathways shown at bar end: "n=5000"

**Use case:** Primary figure for pattern comparison

**Example interpretation:**
```
hallmark  ████████████████████70%████████20%███10%  n=10
          [Comp: 7 (70%)] [Prog: 2 (20%)] [Other: 1 (10%)]

gobp      ████████45%█████████████30%█████████25%  n=5000
          [Comp: 2250 (45%)] [Prog: 1500 (30%)] [Other: 1250 (25%)]
```

**Insight:** Hallmark has higher compensation proportion (70% vs 45%) despite fewer pathways.

---

#### 2. **Dual-Panel Comparison**

**Files:** `pattern_comparison_dual_G32A.pdf`, `pattern_comparison_dual_R403C.pdf`

**Features:**
- **Panel A:** Normalized (100% stacked) - for proportion comparison
- **Panel B:** Absolute counts - for context/reference

**Use case:** When both perspectives needed in single figure

---

#### 3. **Percentage Heatmap**

**File:** `pattern_summary_heatmap.pdf`

**Features:**
- Compact matrix view (Database × Pattern)
- Color intensity = percentage
- Text annotations show both "X% (n=Y)"

**Use case:** Quick overview of all pattern distributions

---

## Comparison: Original vs Normalized

### Original Figures

**Location:** `03_Results/02_Analysis/Plots/Cross_database_validation/pattern_summary.pdf`
**Script:** `02_Analysis/4.visualize_trajectory_patterns.py`

**Strengths:**
- Shows absolute pathway counts
- Good for understanding database sizes
- Simple side-by-side comparison

**Limitations:**
- Large databases dominate visually
- Cannot compare pattern proportions across databases
- Small databases hard to interpret

### Normalized Figures

**Location:** `03_Results/02_Analysis/Plots/Pattern_Summary_Normalized/*.pdf`
**Script:** `02_Analysis/8.pattern_summary_normalized.py`

**Strengths:**
- Equal visual weight for all databases
- Direct proportion comparison
- Absolute counts preserved as annotations
- Answers "which database has highest % compensation?"

**Limitations:**
- Database size information less prominent
- Requires reading annotations for absolute context

---

## Usage Recommendations

### For Publications/Presentations

**Use normalized version when:**
- Comparing pattern distributions across databases
- Highlighting relative differences in compensation/progression
- Database sizes are very different (e.g., 10 vs 5000)

**Use original version when:**
- Showing absolute numbers of pathways
- Emphasizing the scale of affected pathways
- Context of database sizes is important

**Best practice:** Show both!
- Main figure: Normalized (for comparison)
- Supplementary: Original (for absolute context)

### For Analysis

1. **Start with normalized** to identify interesting proportions
2. **Check original** to verify absolute significance
3. **Cross-reference** to avoid over-interpreting small databases

---

## Technical Implementation

### Normalization Formula

```python
# For each database independently:
totals = pathway_counts_per_database.sum()
proportions = (pathway_counts / totals) * 100  # Convert to percentage
```

### Annotation Format

```python
# Inside bars (if space > 5%):
label = f'{absolute_count}\n({percentage:.0f}%)'

# At bar end:
label = f'n={total_pathways}'
```

### Pattern Order

Both original and normalized use same pattern order:
1. Compensation
2. Progressive
3. Natural_worsening
4. Natural_improvement
5. Late_onset
6. Transient

---

## Generated Files

```
Pattern_Summary_Normalized/
├── pattern_summary_normalized.pdf           # Main figure (100% stacked)
├── pattern_summary_normalized.png
├── pattern_comparison_dual_G32A.pdf         # Dual-panel for G32A
├── pattern_comparison_dual_G32A.png
├── pattern_comparison_dual_R403C.pdf        # Dual-panel for R403C
├── pattern_comparison_dual_R403C.png
├── pattern_summary_heatmap.pdf              # Heatmap alternative
├── pattern_summary_heatmap.png
└── README.md                                # This file
```

---

## Regenerating Figures

```bash
# Generate all normalized figures
python3 02_Analysis/8.pattern_summary_normalized.py

# Output: Pattern_Summary_Normalized/*.pdf
```

**Note:** Original figures unchanged - this script creates new alternatives.

---

## Key Insights Revealed by Normalization

*(To be filled after analyzing the normalized figures)*

### Compensation Patterns

- **Database X** shows highest compensation proportion: ___%
- **Database Y** shows lowest: ___%
- Difference more apparent in normalized view

### Progressive/Worsening Patterns

- **Database A** has __% progressive (highest)
- **Database B** has __% progressive (lowest)

### Cross-Database Comparison

- Enables ranking databases by pattern prevalence
- Identifies databases with unusual distributions
- Reveals patterns masked by absolute counts

---

## Related Documentation

- Original pattern summary: `../Cross_database_validation/README.md`
- Pattern classification logic: `01_Scripts/Python/patterns.py`
- Trajectory framework: `02_Analysis/README_trajectory_analysis.md`

---

**Generated:** 2025-11-25
**Script:** `02_Analysis/8.pattern_summary_normalized.py`
**Purpose:** Address database size bias in pattern visualization
