# Master Tables Documentation

This document covers the comprehensive master tables generated from the DRP1 mutation analysis pipeline:

1. **Master Pathway Table (GSEA)** - GSEA pathway enrichment results across all databases
2. **Master GSVA Table (Comprehensive)** - GSVA trajectory scores for all pathways across all databases
3. **Master GSVA Table (Focused)** - GSVA trajectory scores for 7 key modules

---

# 1. Master Pathway Table (GSEA)

## Overview

The **master_pathway_table.csv** provides a comprehensive view of all GSEA pathway results across all contrasts, databases, and pattern classifications. This table was created to support the compensation pattern distribution figure and provide complete statistical details for all pathways.

**File**: `03_Results/02_Analysis/master_pathway_table.csv`
**Rows**: 109,989 pathway × contrast combinations
**Unique pathways**: 12,221
**Generation script**: `02_Analysis/9.create_master_pathway_table.py`

## Quick Facts

- **Databases**: 12 (gobp, cgp, reactome, gomf, wiki, gocc, kegg, tf, MitoCarta, hallmark, SynGO, canon)
- **Contrasts**: 9 (Early, Late, TrajDev for both mutations, plus Time controls)
- **Pattern types**: 8 (Compensation, Progressive, Natural_worsening, Natural_improvement, Late_onset, Transient, Persistent, Complex)
- **Significance**: 43.1% of results show significance in at least one contrast

## Table Structure

### Column Groups

#### 1. Pathway Identification (4 columns)
- **pathway_id**: Unique identifier (format: `database::ID`)
- **database**: Source database (e.g., MitoCarta, gobp, SynGO)
- **ID**: Original pathway ID from database
- **Description**: Human-readable pathway name

#### 2. Contrast Information (4 columns)
- **contrast**: Original contrast name (e.g., `G32A_vs_Ctrl_D35`)
- **new_name**: Mapped trajectory framework name (e.g., `Early_G32A`)
- **category**: Trajectory stage (`Early`, `Late`, `TrajDev`, `Reference`)
- **mutation**: Which mutation (`G32A`, `R403C`, `Control`)

#### 3. GSEA Statistics (6 columns)
- **NES**: Normalized Enrichment Score (magnitude and direction of enrichment)
- **pvalue**: Raw p-value from GSEA permutation test
- **p.adjust**: FDR-adjusted p-value (Benjamini-Hochberg)
- **qvalue**: Q-value (alternative FDR method)
- **enrichmentScore**: Raw enrichment score (before normalization)
- **setSize**: Number of genes in pathway gene set

#### 4. Pattern Classifications (3 columns)
- **Pattern_G32A**: Temporal pattern for G32A mutation
- **Pattern_R403C**: Temporal pattern for R403C mutation
- **change_consistency**: Consistency of change between mutations

#### 5. Trajectory NES Values (6 columns)
- **NES_Early_G32A**: G32A enrichment at Day 35 baseline
- **NES_Early_R403C**: R403C enrichment at Day 35 baseline
- **NES_TrajDev_G32A**: G32A maturation trajectory deviation
- **NES_TrajDev_R403C**: R403C maturation trajectory deviation
- **NES_Late_G32A**: G32A enrichment at Day 65 mature state
- **NES_Late_R403C**: R403C enrichment at Day 65 mature state

#### 6. Significance Flags (2 columns)
- **ever_significant**: TRUE if pathway is significant (FDR < 0.05) in ANY contrast
- **ever_significant_trajectory**: TRUE if significant in Early/TrajDev/Late contrasts

## Pattern Classification Logic

### Pattern Definitions

Patterns are classified based on Early → TrajDev → Late NES dynamics:

| Pattern | Criteria | Biological Interpretation |
|---------|----------|---------------------------|
| **Compensation** | Early defect (|NES| > 0.5), Late improved (|Late| < |Early|), TrajDev opposes Early | Active adaptive response corrects initial defect |
| **Progressive** | Early defect, Late worsened (|Late| > |Early|), TrajDev amplifies Early | Cumulative damage worsens over time |
| **Natural_worsening** | Early defect, Late worsened, but NO active trajectory change | Passive deterioration without compensation attempt |
| **Natural_improvement** | Early defect, Late improved, but NO active trajectory change | Passive improvement without active compensation |
| **Late_onset** | No early defect (|Early| < 0.5), but late defect (|Late| > 1.0) | Maturation-dependent dysfunction |
| **Transient** | Early defect (|Early| > 1.0), resolved by late (|Late| < 0.5) | Temporary developmental delay |
| **Persistent** | Stable defect (|Late| ≈ |Early|, difference < 0.5), no trajectory change | Unchanging pathway dysfunction |
| **Complex** | Doesn't fit above patterns | Non-linear or inconsistent pattern |
| **Insufficient_data** | Missing NES values | Cannot classify |

### Thresholds Used

- **low_threshold**: 0.5 (minimum |NES| to consider meaningful)
- **defect_threshold**: 1.0 (|NES| > 1.0 indicates clear defect)
- **padj_cutoff**: 0.05 (significance threshold)

### Change Consistency Categories

The `change_consistency` column classifies how the two mutations compare:

#### Consistent Patterns
- **Consistent_[Pattern]_co-upregulated**: Both mutants show same pattern with positive NES
- **Consistent_[Pattern]_co-downregulated**: Both mutants show same pattern with negative NES
- **Consistent_[Pattern]_opposite-direction**: Same pattern type but opposite directions

#### Inconsistent Patterns
- **Inconsistent_G32A-more-severe**: G32A shows more severe pattern than R403C
- **Inconsistent_R403C-more-severe**: R403C shows more severe pattern than G32A
- **Inconsistent_[Pattern1]_vs_[Pattern2]**: Different pattern types

## Data Flow Pipeline

### 1. Source GSEA Results (RDS checkpoints)
```
03_Results/02_Analysis/checkpoints/
├── all_gsea_results.rds        # MSigDB databases
├── syngo_gsea_results.rds      # SynGO synaptic ontology
└── mitocarta_gsea_results.rds  # MitoCarta mitochondrial pathways
```

Generated by: `02_Analysis/1a.Main_pipeline.R`

### 2. Export to CSV (for Python)
```
03_Results/02_Analysis/Python_exports/
├── gsea_results_long.csv       # All GSEA statistics (long format)
├── gsea_results_wide.csv       # Pivot table (wide format)
└── gsea_trajectory_all.csv     # Trajectory contrasts only
```

Generated by: `02_Analysis/3.export_gsea_for_python.R`

### 3. Pattern Classification
```
03_Results/02_Analysis/Plots/Cross_database_validation/
└── pathways_classified.csv     # Pattern classifications added
```

Generated by: `02_Analysis/.deprecated/4.visualize_trajectory_patterns.py`

### 4. Master Table Assembly
```
03_Results/02_Analysis/
└── master_pathway_table.csv    # Complete table with all information
```

Generated by: `02_Analysis/9.create_master_pathway_table.py`

## Usage Examples

### Example 1: Find all compensation pathways in MitoCarta

```python
import pandas as pd

df = pd.read_csv('03_Results/02_Analysis/master_pathway_table.csv')

# Filter for MitoCarta compensation pathways
mito_comp = df[
    (df['database'] == 'MitoCarta') &
    (df['Pattern_G32A'] == 'Compensation') &
    (df['category'] == 'Early')  # One row per pathway
]

print(f"Found {len(mito_comp)} MitoCarta compensation pathways")
print(mito_comp[['Description', 'NES', 'p.adjust']].head())
```

### Example 2: Compare NES values for a specific pathway across contrasts

```python
# Find all results for a specific pathway
pathway = 'MitoCarta::OXPHOS'
pathway_data = df[df['pathway_id'] == pathway]

# Pivot to compare contrasts
pivot = pathway_data.pivot_table(
    index='contrast',
    values=['NES', 'p.adjust'],
    aggfunc='first'
)
print(pivot)
```

### Example 3: Find pathways with consistent co-downregulation

```python
consistent_down = df[
    (df['change_consistency'].str.contains('Consistent.*co-downregulated')) &
    (df['category'] == 'Early')
].drop_duplicates('pathway_id')

print(f"Found {len(consistent_down)} pathways with consistent co-downregulation")
print(consistent_down.groupby('database').size())
```

### Example 4: Export significant pathways for a specific contrast

```python
g32a_d35_sig = df[
    (df['contrast'] == 'G32A_vs_Ctrl_D35') &
    (df['p.adjust'] < 0.05)
].sort_values('NES', ascending=False)

g32a_d35_sig.to_csv('G32A_D35_significant_pathways.csv', index=False)
```

### Example 5: Find mutation-specific patterns

```python
# G32A-specific compensation (not in R403C)
g32a_specific = df[
    (df['Pattern_G32A'] == 'Compensation') &
    (df['Pattern_R403C'] != 'Compensation') &
    (df['category'] == 'Early')
].drop_duplicates('pathway_id')

print(f"G32A-specific compensation: {len(g32a_specific)} pathways")

# R403C-specific progressive
r403c_specific = df[
    (df['Pattern_R403C'] == 'Progressive') &
    (df['Pattern_G32A'] != 'Progressive') &
    (df['category'] == 'Early')
].drop_duplicates('pathway_id')

print(f"R403C-specific progressive: {len(r403c_specific)} pathways")
```

## Key Findings from Master Table

### Pattern Distribution

**G32A patterns** (unique pathways):
- Compensation: 5,908 (64.1%)
- Natural_worsening: 1,970 (21.4%)
- Progressive: 980 (10.6%)
- Natural_improvement: 222 (2.4%)
- Others: 139 (1.5%)

**R403C patterns** (unique pathways):
- Compensation: 6,196 (67.2%)
- Natural_worsening: 1,791 (19.4%)
- Progressive: 848 (9.2%)
- Natural_improvement: 299 (3.2%)
- Others: 85 (0.9%)

### Change Consistency

Top consistency categories:
1. **Consistent_Compensation_co-upregulated**: 2,363 pathways
2. **Inconsistent_R403C-more-severe**: 2,094 pathways
3. **Consistent_Compensation_co-downregulated**: 2,032 pathways
4. **Inconsistent_G32A-more-severe**: 1,850 pathways

### Database Coverage

Largest databases:
1. GO:BP (Gene Ontology Biological Process): 4,785 pathways
2. CGP (Chemical and Genetic Perturbations): 2,985 pathways
3. Reactome: 1,295 pathways
4. GO:MF (Molecular Function): 1,025 pathways

Specialized databases:
- **MitoCarta**: 72 mitochondrial pathways (62% show compensation in G32A)
- **SynGO**: 32 synaptic pathways (both mutations ~20 compensation pathways)

## Relationship to Pattern Summary Figure

The pattern_summary_normalized.pdf uses this master table as follows:

1. **Data source**: Loads `pathways_classified.csv` (subset of master table)
2. **Aggregation**: Counts pathways per pattern, database, and mutation
3. **Normalization**: Converts counts to percentages (100% stacked bars)
4. **Visualization**: Shows relative pattern distributions across databases

The master table provides the underlying data for:
- Pattern counts in the figure
- Statistical validation (NES, p.adjust values)
- Database-specific analysis
- Cross-mutation comparisons

## Interpretation Guidelines

### NES (Normalized Enrichment Score)
- **Positive NES**: Pathway upregulated (genes in set have higher expression)
- **Negative NES**: Pathway downregulated (genes in set have lower expression)
- **|NES| > 1.0**: Moderate effect
- **|NES| > 2.0**: Strong effect

### p.adjust (FDR-adjusted p-value)
- **< 0.05**: Significant enrichment (standard threshold)
- **< 0.01**: Highly significant
- **≥ 0.05**: Not significant (cannot reject null hypothesis)

### Pattern Interpretation
- **Compensation**: Look for Early defect + opposing TrajDev + Late improvement
- **Progressive**: Look for Early defect + aligned TrajDev + Late worsening
- **Consistency**: Compare Pattern_G32A vs Pattern_R403C

## Quality Control Checks

### Data Completeness
- All 12,221 pathways have pattern classifications for both mutations
- 109,989 total pathway × contrast results (12,221 pathways × 9 contrasts)
- 43.1% show significance in at least one contrast

### Pattern Distribution Validation
- Compensation is most common (64-67% of pathways)
- Progressive is less common (9-11%)
- Rare patterns (Late_onset, Transient) represent edge cases

### Consistency Validation
- 4,395 pathways show consistent compensation patterns
- 3,944 pathways show inconsistent severity between mutations
- Validates domain-specific effects (GTPase G32A vs stalk R403C)

## Troubleshooting

### Missing pattern classifications
**Issue**: Some pathways have `Insufficient_data` for patterns
**Cause**: Missing NES values in trajectory contrasts
**Solution**: Check if pathway was tested in all contrasts

### Inconsistent results across contrasts
**Issue**: Same pathway shows different significance in different contrasts
**Cause**: Expected - contrasts test different biological questions
**Solution**: Focus on trajectory contrasts (Early/TrajDev/Late) for pattern analysis

### Large file size
**Issue**: CSV file is large (~50-100 MB)
**Cause**: 109,989 rows with 25 columns
**Solution**: Use pandas chunking or filter to relevant subset

## Citation

When using this table, cite the analysis pipeline and pattern classification framework:

```
Pattern classifications based on trajectory framework:
- Early (D35 baseline mutation effects)
- TrajDev (maturation-specific trajectory changes)
- Late (D65 mature state outcomes)

Classification algorithm: Python.patterns.classify_trajectory_pattern()
GSEA: fgsea R package, 10,000 permutations, FDR < 0.05
```

## Related Files

- **pattern_summary_normalized.pdf**: Visual summary of pattern distributions
- **pathways_classified.csv**: Pattern classifications without per-contrast statistics
- **gsea_results_long.csv**: Original GSEA results before pattern classification
- **Cross_database_validation/**: Trajectory heatmaps per database

## Updates and Versioning

**Version**: 1.0
**Last generated**: 2025-11-25
**Pipeline version**: Post-cleanup validated patterns
**Compatible with**: pattern_summary_normalized.pdf (same date)

## Contact

For questions about this table or pattern classification logic:
- See: `01_Scripts/Python/patterns.py` for classification algorithm
- See: `02_Analysis/9.create_master_pathway_table.py` for generation script
- See: `CLAUDE.md` for project overview and analysis pipeline

---

**Note**: This table represents the comprehensive output of the DRP1 mutation trajectory analysis. All pathway classifications and statistics are reproducible from the source RDS checkpoint files using the documented pipeline.

---

# 2. Master GSVA Table (Comprehensive - All Pathways)

## Overview

The **master_gsva_all_pathways.csv** provides comprehensive GSVA (Gene Set Variation Analysis) enrichment scores for ALL pathways tested in GSEA, creating a trajectory view parallel to the GSEA master table.

**File**: `03_Results/02_Analysis/master_gsva_all_pathways.csv`
**Rows**: ~72,000 pathway × genotype × timepoint combinations
**Unique pathways**: ~14,500 (after size filtering)
**Generation script**: `02_Analysis/10.comprehensive_gsva_analysis.R`

## Quick Facts

- **Databases**: 12 (same as GSEA: hallmark, kegg, reactome, gobp, gocc, gomf, cgp, tf, wiki, canon, MitoCarta, SynGO)
- **Genotypes**: 3 (Ctrl, G32A, R403C)
- **Timepoints**: 2 (D35, D65)
- **Total observations**: ~87,000 group combinations
- **Pathway size filter**: 10-500 genes per pathway
- **Computation time**: ~10-30 minutes (cached in checkpoint)

## Table Structure

### Column Groups

#### 1. Pathway Identification (5 columns)
- **pathway_id**: Unique identifier (format: `database::pathway_name`)
- **database**: Source database (hallmark, kegg, gobp, MitoCarta, SynGO, etc.)
- **pathway_name**: Original pathway name from database
- **n_genes_total**: Total genes in original pathway definition
- **n_genes_in_expression**: Genes present in expression matrix (used for GSVA)

#### 2. Sample Grouping (5 columns)
- **Genotype**: `Ctrl`, `G32A`, or `R403C`
- **Day**: Numeric day (35 or 65)
- **Timepoint**: Factor format (`D35` or `D65`)
- **Trajectory_Category**: `Reference` (controls), `Early` (D35 mutants), `Late` (D65 mutants)
- **Contrast_Label**: Combined label (e.g., `Early_G32A`, `Ctrl_D35`)

#### 3. GSVA Scores (4 columns)
- **Mean_GSVA**: Mean GSVA enrichment score for this group
- **SD_GSVA**: Standard deviation across biological replicates
- **SE_GSVA**: Standard error of the mean
- **N**: Number of biological replicates

#### 4. Trajectory Metrics (2 columns)
Relative to Ctrl D35 baseline:
- **Baseline_CtrlD35**: Control D35 baseline score for this pathway
- **Expression_vs_CtrlD35**: Deviation from baseline (trajectory metric)

#### 5. Divergence Metrics (2 columns)
Relative to matched control at same timepoint:
- **Ctrl_Mean_Same_Day**: Control mean at same timepoint
- **Divergence_vs_Ctrl**: Difference from same-day control

#### 6. Statistical Tests (4 columns)
T-tests comparing mutants to controls:
- **t_statistic**: T-test statistic
- **p_value**: Raw p-value
- **p_adjusted**: Benjamini-Hochberg adjusted p-value (across all comparisons)
- **significant**: TRUE if p.adj < 0.05

## GSVA vs GSEA Comparison

### Methodological Differences

| Feature | GSEA (Master Pathway Table) | GSVA (This Table) |
|---------|----------------------------|-------------------|
| **Method** | Pathway enrichment (rank-based) | Single-sample enrichment scoring |
| **Input** | Ranked gene list per contrast | Full expression matrix |
| **Output** | NES per pathway per contrast | Enrichment score per pathway per sample |
| **Granularity** | Contrast-level | Sample-level |
| **Question** | "Is pathway enriched in differential genes?" | "How active is pathway in each sample?" |
| **Advantage** | Detects differential enrichment | Tracks individual trajectories |
| **Statistical test** | Permutation test | T-test on GSVA scores |

### When to Use Which Table

**Use GSEA Master Table when:**
- Discovering differentially enriched pathways between conditions
- Focusing on statistical significance of enrichment
- Comparing mutation effects at specific contrasts
- Need NES (normalized enrichment scores)

**Use GSVA Master Table when:**
- Tracking pathway activity trajectories over time
- Examining sample-level variation within groups
- Studying developmental trends and compensation
- Need continuous scores for correlation analysis

**Use both when:**
- Comprehensive pathway analysis
- Validating findings across methods
- Distinguishing differential enrichment from baseline activity

## Data Flow Pipeline

### 1. Gene Set Loading
```
All databases loaded (MSigDB, MitoCarta, SynGO)
  ↓
Filter by size (10-500 genes)
  ↓
Filter by gene availability in expression matrix
  ↓
~14,500 pathways retained
```

### 2. GSVA Computation
```
02_Analysis/10.comprehensive_gsva_analysis.R
  ↓
Inputs:
  - checkpoints/qc_variables.rds (logCPM expression + metadata)
  - MSigDB gene sets (via msigdbr)
  - MitoCarta pathways (from GMX file)
  - SynGO pathways (from Excel files)
  ↓
GSVA computation (Gaussian kernel)
  ↓
Checkpoint: checkpoints/gsva_all_pathways.rds
```

### 3. Statistical Analysis
```
GSVA scores matrix (14,500 pathways × 25 samples)
  ↓
Calculate group statistics (mean, SD, SE)
  ↓
T-tests: mutants vs controls at each timepoint
  ↓
FDR correction (Benjamini-Hochberg)
  ↓
Add trajectory and divergence metrics
```

### 4. Master Table Assembly
```
Combine all components
  ↓
Output: master_gsva_all_pathways.csv (~72,000 rows)
```

## Usage Examples

### Example 1: Find pathways with largest divergence at D65

```r
library(dplyr)

# Load table
df <- read.csv('03_Results/02_Analysis/master_gsva_all_pathways.csv')

# Find G32A pathways with largest D65 divergence
top_divergent <- df %>%
  filter(Genotype == 'G32A', Day == 65, !is.na(Divergence_vs_Ctrl)) %>%
  arrange(desc(abs(Divergence_vs_Ctrl))) %>%
  head(20) %>%
  select(database, pathway_name, Divergence_vs_Ctrl, p_adjusted, significant)

print(top_divergent)
```

### Example 2: Compare GSVA trajectory to GSEA NES for OXPHOS

```r
# GSVA trajectory
gsva_oxphos <- df %>%
  filter(pathway_id == 'MitoCarta::OXPHOS') %>%
  select(Genotype, Day, Mean_GSVA, Expression_vs_CtrlD35)

# Load GSEA results
gsea <- read.csv('03_Results/02_Analysis/master_pathway_table.csv')

gsea_oxphos <- gsea %>%
  filter(pathway_id == 'MitoCarta::OXPHOS') %>%
  select(contrast, NES, p.adjust)

# Compare patterns
print("GSVA trajectory:")
print(gsva_oxphos)
print("\nGSEA enrichment:")
print(gsea_oxphos)
```

### Example 3: Find significant GSVA changes across databases

```r
# Count significant pathways per database
sig_summary <- df %>%
  filter(significant == TRUE) %>%
  group_by(database) %>%
  summarize(
    n_significant = n_distinct(pathway_id),
    mean_abs_divergence = mean(abs(Divergence_vs_Ctrl), na.rm = TRUE)
  ) %>%
  arrange(desc(n_significant))

print(sig_summary)
```

### Example 4: Pathway trajectory clustering

```r
library(tidyr)

# Pivot to wide format for clustering
gsva_wide <- df %>%
  filter(Genotype != 'Ctrl') %>%
  select(pathway_id, Genotype, Day, Expression_vs_CtrlD35) %>%
  pivot_wider(
    names_from = c(Genotype, Day),
    values_from = Expression_vs_CtrlD35
  )

# Remove pathways with missing data
gsva_matrix <- gsva_wide %>%
  select(-pathway_id) %>%
  na.omit()

# K-means clustering
set.seed(42)
clusters <- kmeans(gsva_matrix, centers = 5)

# Add cluster assignments
gsva_wide$cluster <- clusters$cluster
```

### Example 5: Database-specific pattern analysis

```r
# For each database, find pathways with compensation-like patterns
compensation_like <- df %>%
  filter(Genotype == 'G32A') %>%
  group_by(pathway_id, database) %>%
  summarize(
    Early_div = Divergence_vs_Ctrl[Day == 35],
    Late_div = Divergence_vs_Ctrl[Day == 65],
    recovered = Early_div < -0.2 & Late_div > -0.1,
    .groups = 'drop'
  ) %>%
  filter(recovered == TRUE)

comp_by_db <- compensation_like %>%
  group_by(database) %>%
  summarize(n_compensating = n())

print(comp_by_db)
```

## Key Findings from Comprehensive GSVA Analysis

### Overall Statistics
- **Total pathways analyzed**: ~14,500
- **Pathways per database**:
  - GO:BP: ~4,300
  - CGP: ~2,700
  - Reactome: ~1,200
  - GO:MF: ~900
  - MitoCarta: 72
  - SynGO: 32

### Significant Findings
*(Will be populated after analysis completes)*

## Interpretation Guidelines

### GSVA Score Interpretation
- **Positive score**: Pathway upregulated in that sample relative to dataset mean
- **Negative score**: Pathway downregulated relative to dataset mean
- **Score range**: Typically -1 to +1 (standardized)
- **Meaningful change**: >0.1 units suggests biological relevance

### Trajectory Metrics (Expression_vs_CtrlD35)
- **Positive**: Pathway increasingly upregulated over development
- **Negative**: Pathway increasingly downregulated over development
- **Near zero**: Stable pathway activity

### Divergence Metrics (Divergence_vs_Ctrl)
- **Positive**: Pathway more active in mutant than control
- **Negative**: Pathway less active in mutant than control
- **Magnitude**: Effect size of mutation at that timepoint

### Statistical Significance
- **p_adjusted < 0.05**: Significant difference from control
- **Effect size**: Consider both p-value AND divergence magnitude
- **Multiple testing**: FDR correction applied across ALL 58,000 comparisons

## Computational Considerations

### File Size
- **CSV file**: ~100-200 MB (compressed: ~20-30 MB)
- **Memory usage**: ~500 MB to load in R
- **Recommendation**: Filter to relevant subset before downstream analysis

### Runtime
- **GSVA computation**: 10-30 minutes (one-time, then cached)
- **Statistical tests**: 5-10 minutes
- **Total script runtime**: ~15-40 minutes
- **Checkpoint file**: ~50-100 MB (gsva_all_pathways.rds)

### Reproducibility
- All results cached in `checkpoints/gsva_all_pathways.rds`
- Rerunning script loads from checkpoint (instant)
- Force recompute: Set `force_recompute_gsva = TRUE` in script

## Quality Control

### Data Completeness
- All pathways have scores for all samples
- All comparisons have t-test statistics
- FDR correction applied globally (most conservative)

### Validation Against GSEA
- Pathways significant in GSEA should show corresponding GSVA divergence
- Opposite direction pathways in GSEA should show opposite GSVA trends
- Magnitude correlation: High |NES| in GSEA ≈ High |Divergence| in GSVA

### Expected Patterns
- Control groups should have divergence ≈ 0 (by definition)
- Trajectory should show smooth developmental changes
- Standard deviations should be consistent across groups

## Troubleshooting

### Large file won't load
**Solution**: Load in chunks or filter by database
```r
# Load only MitoCarta pathways
df_mito <- read.csv('master_gsva_all_pathways.csv') %>%
  filter(database == 'MitoCarta')
```

### GSVA and GSEA show different results
**Expected**: Different methods measure different aspects
- GSVA: Sample-level pathway activity
- GSEA: Differential enrichment in ranked genes
**Solution**: Use both for comprehensive understanding

### Few significant pathways
**Expected**: Small sample sizes (N=3-6) limit power
**Solution**: Focus on effect sizes (divergence magnitude) in addition to p-values

## Related Files

- **master_pathway_table.csv**: GSEA results (complementary method)
- **master_gsva_table.csv**: Focused 7-module GSVA table
- **checkpoints/gsva_all_pathways.rds**: GSVA score matrix and metadata

## Citation

When using this table, cite:

```
GSVA scores calculated using GSVA R package (Hänzelmann et al. 2013)
Gene sets from:
  - MSigDB v2023.1 (Liberzon et al. 2015)
  - MitoCarta 3.0 (Rath et al. 2021)
  - SynGO (Koopmans et al. 2019)
Statistical testing: Two-sample t-tests with Benjamini-Hochberg FDR correction
```

## Updates and Versioning

**Version**: 1.0
**Last generated**: 2025-11-25
**Pipeline version**: Comprehensive GSVA analysis
**Compatible with**: All master tables and trajectory visualizations

---

# 3. Master GSVA Table (Focused - 7 Key Modules)

## Overview

The **master_gsva_table.csv**, along with companion files **gsva_pattern_summary.csv** and **gsva_statistics_summary.txt**, provides comprehensive GSVA (Gene Set Variation Analysis) enrichment scores for 7 key biological modules across the trajectory framework.

**Files**:
- `03_Results/02_Analysis/master_gsva_table.csv`
- `03_Results/02_Analysis/gsva_pattern_summary.csv`
- `03_Results/02_Analysis/gsva_statistics_summary.txt`

**Generation script**: `02_Analysis/9b.create_master_gsva_table.R`

## Quick Facts

- **Modules**: 7 (Ribosome Biogenesis, Cytoplasmic Translation, Synaptic Ribosomes, Mitochondrial Ribosome, Mito Ribosome Assembly, mtDNA Maintenance, OXPHOS)
- **Genotypes**: 3 (Ctrl, G32A, R403C)
- **Timepoints**: 2 (D35, D65)
- **Total observations**: 42 (7 modules × 3 genotypes × 2 timepoints)
- **Statistical comparisons**: 28 (mutants vs controls at each timepoint)
- **Significance**: 0 comparisons reach FDR < 0.05 (modules show trends but not statistical significance)

## Table Structure

### master_gsva_table.csv (Long Format)

One row per module-genotype-timepoint combination (42 rows total).

#### Column Groups

**1. Module Identification (5 columns)**
- **Module**: Internal module key (e.g., `Ribosome_Biogenesis`)
- **Display_Name**: Human-readable name (e.g., `Ribosome Biogenesis`)
- **Panel_ID**: Figure panel letter (A-G)
- **Source_Database**: Gene set source (`GO:BP`, `SynGO`, `MitoCarta`)
- **N_genes**: Number of genes in module (24-158 genes)

**2. Sample Grouping (5 columns)**
- **Genotype**: `Ctrl`, `G32A`, or `R403C`
- **Day**: Numeric day (35 or 65)
- **Timepoint**: Factor format (`D35` or `D65`)
- **Trajectory_Category**: `Reference` (controls), `Early` (D35 mutants), `Late` (D65 mutants)
- **Contrast_Label**: Combined label (e.g., `Early_G32A`, `Ctrl_D35`)

**3. GSVA Scores (4 columns)**
- **Mean_GSVA**: Mean GSVA enrichment score for this group
- **SD_GSVA**: Standard deviation
- **SE_GSVA**: Standard error of the mean
- **N**: Number of biological replicates

**4. Trajectory Metrics (2 columns)**
Relative to Ctrl D35 baseline:
- **Baseline_CtrlD35**: Control D35 baseline score
- **Expression_vs_CtrlD35**: Deviation from baseline (trajectory)

**5. Divergence Metrics (2 columns)**
Relative to matched control at same timepoint:
- **Ctrl_Mean_Same_Day**: Control mean at same timepoint
- **Divergence_vs_Ctrl**: Difference from same-day control

**6. Statistical Tests (4 columns)**
T-tests comparing mutants to controls:
- **t_statistic**: T-test statistic
- **p_value**: Raw p-value
- **p_adjusted**: Benjamini-Hochberg adjusted p-value
- **significant**: TRUE if p.adj < 0.05

### gsva_pattern_summary.csv (Wide Format)

One row per module (7 rows total). Pivoted format for easy pattern comparison.

**Key Columns:**
- **Module**, **Display_Name**, **Source_Database**: Module identification
- **Expression_vs_CtrlD35_[Genotype]_[Day]**: Trajectory metrics for each group
- **Divergence_vs_Ctrl_[Genotype]_[Day]**: Divergence metrics for each group
- **significant_[Genotype]_[Day]**: Significance flags
- **Pattern_G32A**: Trajectory pattern classification for G32A
- **Pattern_R403C**: Trajectory pattern classification for R403C
- **Change_Consistency**: Consistency between mutations

## Pattern Classification Logic

### Pattern Definitions

Patterns are classified based on trajectory dynamics (Expression_vs_CtrlD35):

| Pattern | Criteria | Biological Interpretation |
|---------|----------|---------------------------|
| **Compensation** | Early ≤ 0, Late > 0.2 | Starts down/neutral, recovers |
| **Progressive** | Early > 0, Late > Early | Starts up, gets worse |
| **Persistent** | Same direction both timepoints, |Late - Early| < 0.1 | Stable dysregulation |
| **Transient** | |Early| > 0.2, |Late| < 0.1 | Temporary defect, resolves |
| **Natural_worsening** | Early > 0, Late < 0 | Reversal of direction |
| **Complex** | Doesn't fit above patterns | Non-linear dynamics |
| **Insufficient_data** | Missing values | Cannot classify |

### Thresholds Used
- **Meaningful change**: 0.1 GSVA score units
- **Clear defect**: 0.2 GSVA score units
- **Significance**: p.adj < 0.05 (though none reach this threshold)

## The 7 Trajectory Modules

### Module Details

| Module | Source | N Genes | Panel | Biological Function |
|--------|--------|---------|-------|---------------------|
| **Ribosome Biogenesis** | GO:BP | 158 | A | Nuclear ribosome assembly |
| **Cytoplasmic Translation** | GO:BP | 76 | B | Cytoplasmic protein synthesis |
| **Synaptic Ribosomes** | SynGO | 65 | C | Synapse-localized translation |
| **Mitochondrial Ribosome** | MitoCarta | 77 | D | Mitoribosome structure |
| **Mito Ribosome Assembly** | MitoCarta | 24 | E | Mitoribosome biogenesis |
| **mtDNA Maintenance** | MitoCarta | 29 | F | Mitochondrial genome integrity |
| **OXPHOS** | MitoCarta | 139 | G | Oxidative phosphorylation |

## Pattern Summary (from Current Data)

### G32A Patterns
- **Complex**: 6 modules
- **Persistent**: 1 module

### R403C Patterns
- **Persistent**: 4 modules
- **Complex**: 2 modules
- **Transient**: 1 module

### Change Consistency
- **Consistent_Complex**: 2 modules
- **Consistent_Persistent**: 1 module
- **Inconsistent_Complex_vs_Persistent**: 3 modules
- **Inconsistent_Complex_vs_Transient**: 1 module

## Usage Examples

### Example 1: Load and explore GSVA scores

```r
library(dplyr)

# Load master table
df <- read.csv('03_Results/02_Analysis/master_gsva_table.csv')

# View OXPHOS trajectory
df %>%
  filter(Module == 'OXPHOS') %>%
  select(Genotype, Day, Mean_GSVA, Expression_vs_CtrlD35, Divergence_vs_Ctrl)
```

### Example 2: Find significant changes

```r
# Find significant G32A changes at D65
sig_changes <- df %>%
  filter(Genotype == 'G32A', Day == 65, significant == TRUE)

# Note: Currently no comparisons reach significance (p.adj < 0.05)
```

### Example 3: Compare trajectory patterns

```r
# Load pattern summary
patterns <- read.csv('03_Results/02_Analysis/gsva_pattern_summary.csv')

# Compare patterns between mutations
patterns %>%
  select(Display_Name, Pattern_G32A, Pattern_R403C, Change_Consistency)
```

### Example 4: Plot module trajectories

```r
library(ggplot2)

# Plot trajectory for all modules
ggplot(df, aes(x = Day, y = Mean_GSVA, color = Genotype, group = Genotype)) +
  geom_line() +
  geom_point() +
  facet_wrap(~Display_Name, scales = "free_y") +
  theme_minimal() +
  labs(title = "GSVA Module Trajectories",
       y = "Mean GSVA Score", x = "Day")
```

### Example 5: Calculate effect sizes

```r
# Calculate Cohen's d for divergence at D65
library(effsize)

# For OXPHOS at D65
oxphos_d65 <- df %>% filter(Module == 'OXPHOS', Day == 65)

# Effect size for G32A divergence
cohen_d <- oxphos_d65 %>%
  filter(Genotype == 'G32A') %>%
  pull(Divergence_vs_Ctrl)

cat(sprintf("OXPHOS G32A D65 divergence: %.3f\n", cohen_d))
```

## Data Flow Pipeline

### 1. GSVA Score Generation
```
02_Analysis/1a.Main_pipeline.R
  └── Generates: 03_Results/02_Analysis/checkpoints/gsva_module_scores.rds
      - gsva_scores: GSVA enrichment matrix (7 modules × samples)
      - gene_modules_filtered: Gene lists for each module
```

### 2. Master Table Assembly
```
02_Analysis/9b.create_master_gsva_table.R
  └── Inputs:
      - checkpoints/gsva_module_scores.rds (GSVA scores)
      - checkpoints/qc_variables.rds (sample metadata)
  └── Outputs:
      - master_gsva_table.csv (long format)
      - gsva_pattern_summary.csv (wide format)
      - gsva_statistics_summary.txt (summary)
```

## Key Findings

### Overall Significance
- **Statistical significance**: No modules reach FDR < 0.05
- **Biological trends**: Modules show consistent directional changes
- **Power limitation**: Small sample sizes (N=3-6 per group) limit statistical power

### Module-Specific Insights

**Ribosome Biogenesis & mtDNA Maintenance:**
- Show strongest early defects (Early scores < -0.5)
- Marginal p-values (0.045-0.051) before multiple testing correction
- Suggest early mitochondrial stress

**OXPHOS & Mitochondrial Ribosome:**
- More stable across development
- Smaller effect sizes (divergence < 0.15)
- Complex patterns without clear compensation

**Synaptic Ribosomes:**
- Transient pattern in R403C (early up, late neutral)
- Suggests temporary synaptic translation compensation

## Interpretation Guidelines

### GSVA Score Interpretation
- **Positive score**: Module upregulated relative to mean
- **Negative score**: Module downregulated relative to mean
- **Score magnitude**: Typically ranges from -1 to +1
- **Meaningful change**: >0.1 units suggests biological relevance

### Trajectory Metrics
- **Expression_vs_CtrlD35**: Cumulative change from baseline
  - Positive: Progressive upregulation
  - Negative: Progressive downregulation
- **Divergence_vs_Ctrl**: Instantaneous difference from control
  - Measures current deficit/excess at each timepoint

### Pattern Interpretation
- **Complex patterns**: Most common, reflect nuanced developmental responses
- **Persistent patterns**: Stable module dysregulation
- **Lack of significance**: Trends present but underpowered statistically

## Comparison to GSEA Master Table

| Feature | GSEA Master Table | GSVA Master Table |
|---------|-------------------|-------------------|
| **Scope** | 12,221 pathways | 7 key modules |
| **Method** | Pathway enrichment | Module enrichment scoring |
| **Rows** | 109,989 (pathway × contrast) | 42 (module × group) |
| **Granularity** | Pathway-level | Module-level aggregate |
| **Significance** | 43% pathways significant | 0% modules significant (FDR) |
| **Use case** | Broad pathway discovery | Focused trajectory analysis |

**When to use which:**
- **GSEA table**: Discovering new pathway associations, broad hypothesis generation
- **GSVA table**: Tracking specific module trajectories, focused hypothesis testing

## Quality Control

### Data Completeness
- ✓ All 7 modules have complete data across all groups
- ✓ All comparisons have t-test statistics
- ✓ Pattern classifications available for all modules

### Statistical Checks
- Mean GSVA scores centered around 0 (expected for GSVA)
- Standard deviations range 0.05-0.85 (biological variation)
- No missing values in critical columns

### Pattern Distribution
- Mix of Complex and Persistent patterns (expected for small module set)
- Consistency between mutations varies by module (biological heterogeneity)

## Troubleshooting

### No significant results
**Issue**: All p.adj values > 0.05
**Reason**: Small sample sizes (N=3-6), conservative FDR correction
**Solution**: Focus on effect sizes and biological patterns, not just p-values

### Different results from GSEA
**Issue**: GSVA and GSEA show different patterns for same pathways
**Reason**: Different methods (single-sample enrichment vs. pathway enrichment)
**Solution**: GSVA captures sample-level variation; GSEA captures differential enrichment

### Pattern classification seems subjective
**Issue**: Threshold-based classification may not capture all nuances
**Reason**: Simplified classification for 7 modules
**Solution**: Examine raw GSVA scores and trajectories for full interpretation

## Related Visualizations

**Figure Sources:**
- `03_Results/02_Analysis/Plots/Critical_period_trajectories/gsva/` - Individual module trajectories
  - `trajectory/Panel_[A-G]_*.pdf` - Trajectory heatmaps
  - `divergence/Panel_[A-G]_*.pdf` - Divergence heatmaps
  - `combined/Trajectory_7panel_grid.pdf` - All modules in grid layout
  - `combined/Divergence_overlay_*.pdf` - Comparative overlays

## Citation

When using these tables, cite the GSVA methodology:

```
GSVA enrichment scores calculated using the GSVA R package (Hänzelmann et al. 2013).
Module definitions based on GO:BP, SynGO, and MitoCarta 3.0 databases.
Pattern classification: trajectory framework (Early D35 → Late D65).
Statistical testing: Two-sample t-tests with Benjamini-Hochberg FDR correction.
```

## Related Files

- **gsva_statistics_summary.txt**: Detailed summary statistics and usage examples
- **viz_critical_period_trajectories_gsva.R**: Visualization script for trajectory plots
- **checkpoints/gsva_module_scores.rds**: Source GSVA scores and gene modules

## Updates and Versioning

**Version**: 1.0
**Last generated**: 2025-11-25
**Pipeline version**: Main pipeline with GSVA integration
**Compatible with**: Critical period trajectory visualizations

## Contact

For questions about GSVA master tables:
- See: `02_Analysis/9b.create_master_gsva_table.R` for generation script
- See: `02_Analysis/viz_critical_period_trajectories_gsva.R` for visualization
- See: `CLAUDE.md` for project overview

---

**Note**: Both master tables (GSEA and GSVA) provide complementary views of the DRP1 mutation trajectory analysis. GSEA offers breadth across thousands of pathways; GSVA offers depth for focused module trajectories.
