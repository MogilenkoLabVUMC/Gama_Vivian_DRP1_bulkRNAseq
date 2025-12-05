# Cross-Database Validation Plots

## Purpose

This folder contains cross-database GSEA validation analysis that classifies pathway enrichment patterns across the developmental trajectory (D35 → D65) for both DRP1 mutations (G32A and R403C). The analysis validates findings across 10 independent pathway databases to identify robust, replicated biological signals versus database-specific artifacts.

## Generating Script

- **Script**: `02_Analysis/3.4.pattern_summary_normalized.py`
- **Runtime**: ~2-5 minutes
- **Dependencies**:
  - `03_Results/02_Analysis/master_gsea_table.csv` (authoritative source with significance-based patterns)
  - GSEA results from main pipeline checkpoints

**Note**: As of 2025-11-27, this folder uses the authoritative `master_gsea_table.csv` with significance-based pattern classifications. Old files (`pathways_classified.csv`, `pattern_summary.csv`) using magnitude-only classifications have been deprecated and moved to `.deprecated/` folder.

## Pattern Classification Framework

**Canonical Reference:** `01_Scripts/Python/pattern_definitions.py`
**Full Documentation:** `docs/PATTERN_CLASSIFICATION.md`

Pathways are classified into temporal patterns based on their enrichment across three developmental stages (Early, TrajDev, Late). The key distinction is between **active** patterns (requiring significant trajectory deviation) and **passive** patterns (following normal developmental buffering):

| Pattern | Active? | Criteria | Biological Interpretation |
|---------|---------|----------|---------------------------|
| **Compensation** | Active | Early defect + TrajDev significant opposing + Late improved | Adaptive plasticity; system actively compensates |
| **Progressive** | Active | Early defect + TrajDev significant amplifying + Late worsened | Active maladaptive transcriptional response |
| **Natural_worsening** | Passive | Early defect + TrajDev NS + Late worsened | Passive deterioration; lacks adaptive capacity |
| **Natural_improvement** | Passive | Early defect + TrajDev NS + Late improved | Passive recovery; normal developmental buffering |
| **Late_onset** | - | No Early defect + Late defect emerges | Maturation-dependent vulnerability |
| **Transient** | - | Strong Early defect + Late fully resolved | Developmental delay that recovers |
| **Complex** | - | Does not fit other patterns | Non-linear or multiphasic dynamics |

**Key thresholds:**
- Significance: p.adjust < 0.05 (High confidence), < 0.10 (Medium confidence)
- Effect size: |NES| > 0.5 (minimum), |NES| > 1.0 (strong)
- Improvement: |Late|/|Early| < 0.7; Worsening: |Late|/|Early| > 1.3

## Plots Generated

### trajectory_comparative_[database].pdf (10 files)

**Databases analyzed**:
- `gobp` - Gene Ontology Biological Process
- `gocc` - Gene Ontology Cellular Component
- `gomf` - Gene Ontology Molecular Function
- `hallmark` - MSigDB Hallmark gene sets
- `kegg` - KEGG pathways
- `reactome` - Reactome pathways
- `tf` - Transcription factor targets
- `wiki` - WikiPathways
- `MitoCarta` - Mitochondrial pathways (custom)
- `SynGO` - Synaptic gene ontology (custom)

**Plot Structure**:
- **X-axis**: Developmental contrast
  - `Early` = Mutation vs Control at D35 (G32A_vs_Ctrl_D35, R403C_vs_Ctrl_D35)
  - `TrajDev` = Mutation-specific maturation (Maturation_G32A_specific, Maturation_R403C_specific)
  - `Late` = Mutation vs Control at D65 (G32A_vs_Ctrl_D65, R403C_vs_Ctrl_D65)
- **Y-axis**: Pathway name (top pathways by significance)
- **Color**: Normalized Enrichment Score (NES)
  - Red: Upregulated (positive NES)
  - Blue: Downregulated (negative NES)
  - Gray: Not significant (FDR > 0.05)
- **Facets**: Separate panels for G32A and R403C mutations
- **Pattern annotation**: Each pathway labeled with its temporal pattern classification

**Key Features**:
- Only shows pathways significant in at least one contrast (FDR < 0.05)
- NES values normalized for color scale consistency
- Pattern labels enable quick identification of trajectory types
- Ordered by hierarchical clustering of NES profiles

### pattern_summary.pdf

**Purpose**: Summarizes the distribution of temporal patterns across all databases and both mutations

**Plot Structure**:
- **X-axis**: Temporal pattern category
- **Y-axis**: Number of pathways
- **Colors**: Separate bars for G32A (blue) and R403C (red)
- **Facets**: One panel per database

**Interpretation**:
- **Compensation** is the most common pattern across databases, especially in MitoCarta and GO categories
- G32A shows more **Compensation** pathways (45 MitoCarta) vs R403C (29 MitoCarta)
- R403C shows more **Natural_worsening** in MitoCarta (26 vs 14 for G32A)
- Pattern distributions validate major findings across independent databases

## Data Files

### Authoritative Data Source

**Primary**: `03_Results/02_Analysis/master_gsea_table.csv` (109,990 rows)
- Complete pathway-level data with significance-based pattern classifications
- Long format: one row per pathway per contrast
- Pattern columns: `Pattern_G32A`, `Pattern_R403C`, `Confidence_G32A`, `Confidence_R403C`
- Generated by: `4.1.create_master_pathway_table.py`

### Deprecated Files (Moved to `.deprecated/` folder)

⚠️ **DO NOT USE**: `pathways_classified.csv` and `pattern_summary.csv` have been deprecated (2025-11-27)
- **Reason**: Used old magnitude-only pattern classifier
- **Issue**: R403C showed 980 Progressive pathways vs 0 in authoritative table
- **Replacement**: Use `master_gsea_table.csv` or `load_classified_pathways()` from `data_loader.py`
- See `.deprecated/DEPRECATION_NOTICE.md` for migration guide

## Methods

**GSEA Parameters**:
- **Software**: fgsea (R package)
- **Gene ranking**: t-statistic from limma-voom differential expression
- **Thresholds**: FDR < 0.05 for significance
- **Permutations**: 10,000
- **Min/max gene set size**: 15-500 genes

**Pattern Classification Algorithm** (see `docs/PATTERN_CLASSIFICATION.md` for full details):

The algorithm distinguishes **active** patterns (significant TrajDev) from **passive** patterns:

1. Extract NES and p.adjust for each pathway across Early, TrajDev, and Late contrasts
2. Assess defect significance (requires BOTH p.adjust < 0.05 AND |NES| > 0.5)
3. Assess outcome change (improvement: |Late|/|Early| < 0.7; worsening: > 1.3)
4. Assess TrajDev activity (significant: p.adjust < 0.05 AND |NES| > 0.5)
5. Apply classification rules:
   - Early defect + TrajDev sig opposing + Late improved → **Compensation** (active)
   - Early defect + TrajDev sig amplifying + Late worsened → **Progressive** (active)
   - Early defect + TrajDev NS + Late improved → **Natural_improvement** (passive)
   - Early defect + TrajDev NS + Late worsened → **Natural_worsening** (passive)
   - No Early defect + Late strong defect → **Late_onset**
   - Strong Early defect + Late resolved → **Transient**
   - Inconsistent patterns → **Complex**

Classifications include confidence levels: High (p < 0.05) and Medium (0.05 ≤ p < 0.10).

**Statistical Testing**:
- Individual pathway significance: FDR-corrected p-values from GSEA
- Sample size: N=3 biological replicates per group
- Differential expression: limma-voom with TMM normalization

## Key Findings

### Cross-Database Validation

1. **Compensation is widespread**: All databases show substantial compensation pathways
   - G32A: 45 MitoCarta, 36 Hallmark, 824 Reactome
   - R403C: 29 MitoCarta, 42 Hallmark, 893 Reactome

2. **Mitochondrial pathways robustly compensate**: MitoCarta shows strongest compensation signal
   - 45/72 significant pathways (62%) show compensation in G32A
   - Includes: OXPHOS, ribosome assembly, mtDNA maintenance

3. **Natural worsening more common in R403C**:
   - MitoCarta: 26 pathways (R403C) vs 14 (G32A)
   - Suggests R403C has broader, more persistent mitochondrial dysfunction

4. **Synaptic pathways show mixed patterns**: SynGO analysis reveals
   - Both mutations: ~20 compensation pathways (pre/postsynaptic compartments)
   - Progressive pathways minimal (2 each)
   - Transient effects rare

5. **GO:BP captures most pathways**: 3111 (G32A) and 3267 (R403C) compensation pathways
   - Reflects broader biological process coverage
   - Validates findings from more specific databases

### Biological Interpretation

The cross-database consistency supports three major conclusions:

1. **Adaptive compensation is the primary response** to DRP1 mutations during neuronal maturation
2. **Mitochondrial pathways are central** to both mutation and compensation mechanisms
3. **G32A shows stronger compensation** while **R403C shows more persistent dysfunction**, suggesting domain-specific effects (GTPase vs stalk domain)

## Interpretation Guide

### Reading Trajectory Comparative Plots

1. **Look for consistent patterns**: Pathways with similar color patterns across Early → TrajDev → Late indicate robust temporal trends
2. **Compare mutations**: Pathways with different patterns between G32A and R403C panels indicate mutation-specific effects
3. **Focus on pattern labels**: Quickly identify whether pathways compensate, worsen, or show complex dynamics
4. **Color intensity**: Darker colors (red/blue) indicate stronger enrichment effects

### Pattern Summary Bar Plot

- **High compensation bars**: Indicate successful adaptive responses
- **High progressive/natural_worsening bars**: Indicate failure to compensate
- **Database consistency**: Similar patterns across databases validate findings
- **Mutation differences**: Different bar heights reveal mutation-specific trajectory profiles

## Related Analyses

- **Ribosome_paradox/**: Deep-dive into ribosome biogenesis compensation vs synaptic translation failure
- **Mito_translation_cascade/**: Mechanistic analysis of mitochondria → translation → synapse cascade
- **Critical_period_trajectories/**: GSVA-based temporal trajectory modeling
- **Publication_Figures/**: Integrated multi-panel figures for manuscript

## Related Visualizations

See also:
- [../Publication_Figures/](../Publication_Figures/README.md) - Main manuscript figures
- [../Ribosome_paradox/](../Ribosome_paradox/README.md) - Core translation paradox finding
- [../Mito_translation_cascade/](../Mito_translation_cascade/README.md) - Mechanistic cascade visualization
- [../Synaptic_ribosomes/](../Synaptic_ribosomes/README.md) - Synaptic translation deep-dive
- [../Critical_period_trajectories/](../Critical_period_trajectories/README.md) - GSVA temporal analysis
- [../Pattern_Summary_Normalized/](../Pattern_Summary_Normalized/README.md) - Normalized pattern visualizations

## Notes

- Pattern classifications are deterministic based on NES direction and FDR thresholds
- Some pathways may be unclassifiable if they show inconsistent or contradictory enrichment patterns
- Database-specific pathway definitions can lead to different classifications for overlapping gene sets
- MitoCarta and SynGO are custom databases specifically curated for this analysis domain

---

**Last Updated**: 2025-11-25
**Analysis Version**: Post-cleanup, validated patterns
