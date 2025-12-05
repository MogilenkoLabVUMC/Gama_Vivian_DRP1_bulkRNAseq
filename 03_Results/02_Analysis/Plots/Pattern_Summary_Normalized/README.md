# Pattern_Summary_Normalized

## Overview

Normalized pattern distribution visualizations showing trajectory pattern classifications for **ALL 12,221 pathways** tested in GSEA analysis.

**CRITICAL CLARIFICATION:** These figures show pattern proportions for **ALL tested pathways**, not just those with p.adjust < 0.05. Every pathway tested in GSEA receives a pattern classification based on its trajectory dynamics across Early, TrajDev, and Late stages.

## Pathway Inclusion Criteria

### What Is Included

- **ALL 12,221 unique pathways** from the master GSEA table
- Pathways are included regardless of statistical significance
- Pattern classification is applied to every pathway based on trajectory metrics (NES values and p-values)

### Database Breakdown

| Database | Pathways | Description |
|----------|----------|-------------|
| gobp | 4,785 | Gene Ontology Biological Process |
| cgp | 2,985 | Chemical and Genetic Perturbations |
| reactome | 1,295 | Reactome pathways |
| gomf | 1,025 | Gene Ontology Molecular Function |
| wiki | 706 | WikiPathways |
| gocc | 627 | Gene Ontology Cellular Component |
| kegg | 355 | KEGG pathways |
| tf | 272 | Transcription Factor targets |
| MitoCarta | 72 | Mitochondrial pathways |
| hallmark | 50 | MSigDB Hallmark collection |
| SynGO | 32 | Synaptic Gene Ontology |
| canon | 17 | Canonical pathways |

### Pattern Distribution Reality

Most pathways (~78%) are classified as "Complex" because they do not meet the criteria for defined trajectory patterns:

| Pattern | G32A | R403C | Description |
|---------|------|-------|-------------|
| Complex | 9,566 (78%) | 8,851 (72%) | Does not fit defined trajectory criteria |
| Compensation | 1,462 (12%) | 1,612 (13%) | Active adaptation |
| Natural_improvement | 851 (7%) | 1,283 (10%) | Passive recovery |
| Sign_reversal | 244 (2%) | 386 (3%) | Trajectory reversal |
| Late_onset | 94 (<1%) | 73 (<1%) | Maturation-dependent |
| Transient | 2 (<1%) | 12 (<1%) | Temporary defect |
| Natural_worsening | 2 (<1%) | 4 (<1%) | Passive deterioration |

## What the Figures Show

### Visualization Choice: "Complex" Excluded

The figures **exclude the "Complex" category** from visualization to focus on pathways with interpretable trajectory patterns. This means:

- **Bars represent only ~22% of pathways for G32A** (2,655 of 12,221)
- **Bars represent only ~28% of pathways for R403C** (3,370 of 12,221)
- The "n=X" annotation at right shows total pathways per database **after excluding Complex**

### How to Interpret

1. **Proportions are within "meaningful patterns" only** - not the full database
2. **Each bar normalized to 100%** - enables comparison across databases of different sizes
3. **Absolute counts shown as "N (%)"** inside bars for reference
4. **"n=X" totals** at right = pathways in that database with non-Complex patterns

### Example Interpretation

If hallmark shows 50% Compensation with n=25:
- 25 hallmark pathways have non-Complex patterns
- 12-13 of those 25 show Compensation
- The remaining 25 hallmark pathways (50 - 25 = 25) are classified as Complex

## Generating Script

**`02_Analysis/3.4.pattern_summary_normalized.py`**

## Contents

| File | Description |
|------|-------------|
| `pattern_summary_normalized.pdf/png` | Main figure: 100% stacked bars for both mutations (G32A and R403C) |
| `pattern_comparison_dual_G32A.pdf/png` | Dual-panel: normalized (left) vs absolute counts (right) for G32A |
| `pattern_comparison_dual_R403C.pdf/png` | Dual-panel: normalized (left) vs absolute counts (right) for R403C |

## Data Source

- **Input:** `03_Results/02_Analysis/master_gsea_table.csv`
- **Classification:** Pattern assignments from `01_Scripts/Python/pattern_definitions.py`
- **Filtering:** Complex pattern excluded; all other patterns (Compensation, Sign_reversal, Progressive, Natural_worsening, Natural_improvement, Late_onset, Transient) shown

## Pattern Classification Logic

Patterns are classified based on:
1. **Early effect:** NES and p.adjust for D35 mutation vs control
2. **TrajDev effect:** NES and p.adjust for mutation-specific maturation deviation
3. **Late effect:** NES and p.adjust for D65 mutation vs control

Key thresholds (from `pattern_definitions.py`):
- p.adjust < 0.05 = significant
- |NES| > 0.5 = meaningful effect
- |Late|/|Early| < 0.7 = improvement
- |Late|/|Early| > 1.3 = worsening

See `docs/PATTERN_CLASSIFICATION.md` for full specification.

## Regeneration

```bash
python3 02_Analysis/3.4.pattern_summary_normalized.py
```

**Runtime:** < 1 minute

## Usage Guidance

**Appropriate uses:**
- Comparing pattern proportions across databases
- Showing relative prevalence of trajectory response types
- Supporting claims about mutation-specific response profiles

**Inappropriate uses:**
- Claiming "X% of pathways show compensation" (must specify "of non-Complex pathways")
- Direct comparison of bar heights between databases (use dual-panel figures instead)

**Recommended caption template:**
> "Pattern distribution for pathways with classifiable trajectory dynamics (excluding Complex patterns; see Methods). Proportions normalized per database. N=pathways per database with defined patterns."
