# Pattern Classification Framework for DRP1 Trajectory Analysis

**Canonical Reference:** `01_Scripts/Python/pattern_definitions.py`

This document describes the pattern classification system used to characterize pathway enrichment trajectories across neuronal maturation in DRP1 mutants. The framework distinguishes between **active** adaptive responses (requiring significant trajectory deviation) and **passive** changes (following normal developmental buffering).

## Biological Framework

### Trajectory Components

The analysis uses three trajectory stages to characterize pathway dynamics:

| Stage | GSEA Contrast | Biological Meaning |
|-------|---------------|-------------------|
| **Early** | `{Mutation}_vs_Ctrl_D35` | Initial mutation effect at Day 35 (immature neurons) |
| **TrajDev** | `Maturation_{Mutation}_specific` | Mutation-specific deviation from normal maturation |
| **Late** | `{Mutation}_vs_Ctrl_D65` | Mature state mutation effect at Day 65 |

### Key Insight: TrajDev as Active Plasticity

TrajDev (Trajectory Deviation) represents how mutant neurons deviate from normal developmental trajectories:

```
TrajDev = (D65_mut - D35_mut) - (D65_ctrl - D35_ctrl)
```

A **significant TrajDev** indicates that the mutant's maturation actively differs from control maturation—this represents transcriptional plasticity, not just passive developmental changes. This distinction is critical for identifying:

- **Active Compensation:** The system actively adjusts to counteract defects
- **Active Progression:** Defects compound through active maladaptive responses
- **Passive Changes:** Defects follow normal developmental buffering without active response

## Super-Category System (Simplified)

For main text interpretation, patterns are grouped into 5 super-categories:

| Super-Category | Includes | Use Case |
|----------------|----------|----------|
| **Active_Compensation** | Compensation | Main narrative - active adaptive responses |
| **Active_Progression** | Progressive | Main narrative - active worsening (rare in data) |
| **Passive** | Natural_improvement, Natural_worsening | Main narrative - developmental buffering |
| **Late_onset** | Late_onset | Maturation-dependent effects |
| **Other** | Transient, Complex | Requires individual inspection |

**Usage:**
- Main figures and text: Use super-categories for clearer narrative
- Detailed analysis: Use full 7-pattern taxonomy
- Methods/Supplements: Document both levels

**Implementation:** See `pattern_definitions.py` for `SUPER_CATEGORY_MAP` and `add_super_category_columns()`.

---

## Pattern Definitions

### The 7-Pattern Classification System

| Pattern | Active? | Criteria | Biological Interpretation |
|---------|---------|----------|---------------------------|
| **Compensation** | Active | Early defect + TrajDev opposes + Late improved | Adaptive plasticity; system actively compensates for mutation effects |
| **Progressive** | Active | Early defect + TrajDev amplifies + Late worsened | Cumulative damage; active maladaptive response |
| **Natural_worsening** | Passive | Early defect + TrajDev NS + Late worsened | Passive deterioration; no adaptive capacity |
| **Natural_improvement** | Passive | Early defect + TrajDev NS + Late improved | Passive recovery; developmental buffering |
| **Late_onset** | - | No Early defect + Late defect emerges | Maturation-dependent dysfunction |
| **Transient** | - | Strong Early defect + Late resolved | Developmental delay that recovers |
| **Complex** | - | Inconsistent or multiphasic | Requires individual inspection |

### Classification Criteria Details

**Early Defect Requirements:**
- p.adjust < 0.05 AND |NES| > 0.5 (High confidence)
- p.adjust < 0.10 AND |NES| > 0.5 (Medium confidence)

**TrajDev Significance (Active vs Passive):**
- Significant: p.adjust < 0.05 AND |NES| > 0.5
- Direction must be biologically meaningful (opposes or amplifies Early)

**Late Outcome Assessment:**
- Improved: |Late|/|Early| < 0.7 (≥30% reduction) OR |Late| < 0.5
- Worsened: |Late|/|Early| > 1.3 (≥30% increase)
- Resolved: |Late| < 0.5

**Special Patterns:**
- Late_onset: Early p.adjust ≥ 0.10 OR |Early| ≤ 0.5; Late p.adjust < 0.05 AND |Late| > 1.0
- Transient: Early p.adjust < 0.05 AND |Early| > 1.0; |Late| < 0.5

## Thresholds

### GSEA (NES-based)

| Threshold | Value | Description |
|-----------|-------|-------------|
| `PADJ_SIGNIFICANT` | 0.05 | Primary significance cutoff |
| `PADJ_TRENDING` | 0.10 | Medium confidence cutoff |
| `NES_EFFECT` | 0.5 | Minimum effect size |
| `NES_STRONG` | 1.0 | Strong effect (for Transient, Late_onset) |
| `IMPROVEMENT_RATIO` | 0.7 | |Late|/|Early| threshold for improvement |
| `WORSENING_RATIO` | 1.3 | |Late|/|Early| threshold for worsening |

### GSVA (Enrichment Score-based)

GSVA scores typically range from -1 to 1. Thresholds are scaled accordingly:

| Threshold | Value | Equivalent to |
|-----------|-------|---------------|
| `GSVA_EFFECT` | 0.15 | NES 0.5 |
| `GSVA_STRONG` | 0.30 | NES 1.0 |

## Confidence Levels

Classifications include confidence levels to distinguish robust from borderline findings:

| Confidence | Early Defect | Interpretation |
|------------|-------------|----------------|
| **High** | p.adjust < 0.05 AND |NES| > 0.5 | Robust classification |
| **Medium** | p.adjust < 0.10 AND |NES| > 0.5 | Potential, requires validation |

**Reporting Recommendations:**
- **Primary results:** Report High confidence counts
- **Sensitivity analysis:** Report High + Medium (potential) counts
- **Individual pathways:** Note confidence level when discussing specific pathways

## Algorithm Pseudocode

```python
def classify_pattern(early_nes, early_padj, trajdev_nes, trajdev_padj, late_nes, late_padj):
    # Step 1: Assess Early defect
    early_sig_defect = (early_padj < 0.05) and (|early_nes| > 0.5)
    early_trending = (early_padj < 0.10) and (|early_nes| > 0.5)
    early_strong = (early_padj < 0.05) and (|early_nes| > 1.0)
    early_no_defect = (early_padj >= 0.10) or (|early_nes| <= 0.5)

    # Step 2: Assess Late outcome
    late_sig_defect = (late_padj < 0.05) and (|late_nes| > 1.0)
    late_resolved = |late_nes| < 0.5
    improved = (|late|/|early| < 0.7) or late_resolved
    worsened = |late|/|early| > 1.3

    # Step 3: Assess TrajDev (Active vs Passive)
    trajdev_sig = (trajdev_padj < 0.05) and (|trajdev_nes| > 0.5)
    trajdev_opposes = sign(trajdev_nes) != sign(early_nes)
    trajdev_amplifies = sign(trajdev_nes) == sign(early_nes)

    # Step 4: Classify (order matters!)
    # Late_onset: no early defect, significant late defect
    if early_no_defect and late_sig_defect:
        return "Late_onset", "High"

    if early_sig_defect or early_trending:
        confidence = "High" if early_sig_defect else "Medium"

        # Active patterns FIRST - significant TrajDev indicates active compensation,
        # not passive transient recovery
        if trajdev_sig and trajdev_opposes and improved:
            return "Compensation", confidence

        if trajdev_sig and trajdev_amplifies and worsened:
            return "Progressive", confidence

        # Transient AFTER active patterns - only when TrajDev is NOT significant opposing
        if early_strong and late_resolved:
            return "Transient", "High"

        # Passive patterns (TrajDev not significant)
        if not trajdev_sig:
            if improved:
                return "Natural_improvement", confidence
            if worsened:
                return "Natural_worsening", confidence

    return "Complex", None
```

**Important:** Active patterns (Compensation, Progressive) are checked BEFORE Transient because significant opposing TrajDev indicates active compensatory mechanisms, not passive developmental recovery.

## Implementation Files

| File | Purpose |
|------|---------|
| `01_Scripts/Python/pattern_definitions.py` | **Canonical source** - all thresholds and classify_pattern() |
| `01_Scripts/Python/patterns.py` | User-facing wrapper with backward compatibility |
| `01_Scripts/Python/semantic_categories.py` | PATTERN_COLORS for visualization |
| `02_Analysis/9b.create_master_gsva_table.R` | R implementation for GSVA analysis |

## Usage Examples

### Python (GSEA Analysis)

```python
from Python.pattern_definitions import classify_pattern, add_pattern_classification

# Single pathway classification
pattern, confidence = classify_pattern(
    early_nes=-1.5, early_padj=0.01,
    trajdev_nes=1.2, trajdev_padj=0.02,
    late_nes=-0.3, late_padj=0.5
)
# Returns: ('Compensation', 'High')

# Batch classification
df = add_pattern_classification(gsea_wide_df, mutations=['G32A', 'R403C'])
# Adds Pattern_G32A, Confidence_G32A, Pattern_R403C, Confidence_R403C columns
```

### R (GSVA Analysis)

```r
# Run the master GSVA table generator
source("02_Analysis/9b.create_master_gsva_table.R")

# Outputs include Pattern_G32A, Pattern_R403C, Confidence_* columns
```

## Validation Checklist

Before finalizing pattern classifications:

- [ ] Confirm Early, TrajDev, Late contrasts are correctly mapped
- [ ] Verify NES and p.adjust columns are named consistently
- [ ] Check that TrajDev direction is biologically meaningful
- [ ] Compare High vs High+Medium counts for sensitivity
- [ ] Cross-validate key findings across multiple databases

## Important Note: Descriptive Classification

This pattern classification system is **descriptive**, not inferential. Patterns summarize trajectory dynamics using pre-defined criteria—they do not represent formal statistical hypothesis tests comparing pattern frequencies or testing whether a pathway "belongs to" one pattern versus another.

Claims about pathway trajectories will need validation, and would need to show the underlying trajectory data to allow visual verification.

---

## Methods Section Template

> Pathway enrichment patterns were classified using a significance-based trajectory framework that distinguishes active adaptive responses from passive developmental changes. Classification required both statistical significance (p.adjust < 0.05 for High confidence, p.adjust < 0.10 for Medium confidence) and biological effect size (|NES| > 0.5). Seven mutually exclusive patterns were defined based on three trajectory stages: Early (mutation effect at Day 35), TrajDev (mutation-specific deviation from control maturation trajectory), and Late (mutation effect at Day 65).
>
> **Compensation** was assigned to pathways with an Early defect (p.adjust < 0.05, |NES| > 0.5) that showed significant trajectory deviation opposing the defect direction (TrajDev p.adjust < 0.05, |NES| > 0.5, opposite sign to Early) resulting in Late improvement (≥30% reduction in |NES| or |Late NES| < 0.5). **Progressive** patterns required the same Early criteria but with trajectory deviation amplifying the defect direction (same sign as Early) and Late worsening (≥30% increase in |NES|).
>
> Passive patterns (**Natural_improvement**, **Natural_worsening**) were assigned when Early defects improved or worsened without significant trajectory deviation (TrajDev p.adjust ≥ 0.05 or |NES| ≤ 0.5), indicating reliance on normal developmental buffering rather than active transcriptional adaptation.
>
> **Late_onset** patterns captured pathways with no significant Early defect (p.adjust ≥ 0.10 or |NES| ≤ 0.5) but emergence of strong Late dysfunction (p.adjust < 0.05, |NES| > 1.0). **Transient** patterns identified strong Early defects (p.adjust < 0.05, |NES| > 1.0) that fully resolved by Late stage (|NES| < 0.5). Pathways not fitting these criteria were classified as **Complex**.

## Change Log

| Date | Version | Changes |
|------|---------|---------|
| 2025-11-26 | 1.1 | Added super-category system; added descriptive classification note |
| 2025-11-26 | 1.0 | Initial canonical definitions with significance requirements |

---

**Canonical Source:** `01_Scripts/Python/pattern_definitions.py`
**Last Updated:** 2025-11-26
