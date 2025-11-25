# Master Pathway Table Documentation

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
