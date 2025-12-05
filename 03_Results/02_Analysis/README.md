# Differential Expression Analysis Results

**Project:** DRP1 Mutation Bulk RNA-seq Analysis
**Last Updated:** 2025-12-04
**Pipeline Version:** Post-migration validated (v1.0+)
**Pattern System Version:** 2.0 (8-pattern taxonomy with Sign_reversal)

---

## Overview

This directory contains the complete differential expression and pathway enrichment analysis of iPSC-derived cortical neurons carrying DRP1 mutations (G32A and R403C) compared to isogenic controls at two developmental timepoints (Day 35 and Day 65).

### Experimental Design

| Factor | Levels | Description |
|--------|--------|-------------|
| **Genotype** | Control, G32A, R403C | DRP1 mutation status |
| **Timepoint** | D35, D65 | Days of neuronal maturation |
| **Replicates** | n=3 per condition | Biological replicates |
| **Total samples** | 18 | 3 genotypes × 2 timepoints × 3 replicates |

### Mutation Biology

- **G32A**: GTPase domain mutation affecting DRP1 enzymatic activity
- **R403C**: Stalk domain mutation affecting oligomerization and membrane binding

---

## Directory Structure

```
03_Results/02_Analysis/
├── README.md                        # This file
│
├── DE_results/                      # Differential expression statistics (9 contrasts)
├── checkpoints/                     # Cached R objects for reproducibility
├── Calcium_genes/                   # Calcium signaling gene analysis
├── Python_exports/                  # Data exports for Python visualization
├── Verification_reports/            # QC reports and gene lists
├── Summary/                         # DE gene count summaries
│
├── Plots/                           # All visualizations
│   ├── README.md                   # Comprehensive figure overview
│   ├── General/                    # QC plots (MDS, correlation, UpSet)
│   ├── Volcano/                    # Volcano plots (8 variants)
│   ├── GSEA/                       # Pathway enrichment (12 databases × 9 contrasts)
│   ├── Publication_Figures/        # Final manuscript figures (Python)
│   ├── Ribosome_paradox/           # Core finding visualization
│   ├── Mito_translation_cascade/
│   ├── Synaptic_ribosomes/
│   ├── Critical_period_trajectories/
│   ├── Cross_database_validation/
│   └── Pattern_Summary_Normalized/
│
├── master_gsea_table.csv            # Comprehensive GSEA results (110K rows)
├── master_gsva_focused_table.csv    # Focused GSVA (7 key modules, 42 rows)
├── master_gsva_all_table.csv        # Comprehensive GSVA (all pathways, 87K rows)
├── gsva_pattern_summary.csv         # GSVA pattern classifications (7 modules)
└── gsva_statistics_summary.txt      # GSVA summary statistics
```

---

## Analysis Pipeline

### Statistical Methods

#### 1. Differential Expression Analysis

**Software:** limma-voom (v3.56+) with edgeR normalization

**Steps:**
1. **Normalization:** TMM (Trimmed Mean of M-values) via edgeR
2. **Transformation:** voom log2-CPM with precision weights
3. **Linear modeling:** `voomLmFit()` with factorial design matrix
4. **Empirical Bayes:** `eBayes(robust = TRUE)` for moderated statistics
5. **Multiple testing:** Benjamini-Hochberg FDR correction

**Key Parameters:**

| Parameter | Value | Usage |
|-----------|-------|-------|
| fdr_cutoff | 0.05 | DEG classification via limma::decideTests (BH-adjusted) |
| p_cutoff | 0.05 | Volcano plots (p mode) raw p-value threshold |
| fc_cutoff | 2 (|log2FC| >= 2) | Volcano plot 4-fold change visualization threshold |

**Note:** DEG counts (bar charts, UpSet plots) use FDR < 0.05 only (no FC cutoff). Volcano plots use FDR ≤ 0.1 (fdr mode) for coloring decisions.

#### 2. Contrast Definitions

Nine contrasts organized into three biological questions:

**Mutation Effects (Baseline):**
- `G32A_vs_Ctrl_D35`: G32A mutation effect at early timepoint
- `R403C_vs_Ctrl_D35`: R403C mutation effect at early timepoint
- `G32A_vs_Ctrl_D65`: G32A mutation effect at late timepoint
- `R403C_vs_Ctrl_D65`: R403C mutation effect at late timepoint

**Maturation Effects (Time):**
- `Time_Ctrl`: Control maturation (D65 - D35)
- `Time_G32A`: G32A maturation trajectory
- `Time_R403C`: R403C maturation trajectory

**Interaction Effects (Difference-in-Difference):**
- `Maturation_G32A_specific`: (D65_G32A - D35_G32A) - (D65_Ctrl - D35_Ctrl)
- `Maturation_R403C_specific`: (D65_R403C - D35_R403C) - (D65_Ctrl - D35_Ctrl)

#### 3. Gene Set Enrichment Analysis (GSEA)

**Software:** fgsea (v1.26+) via clusterProfiler wrapper

**Method:**
- Ranking metric: t-statistic from limma (preserves direction and significance)
- Permutations: 10,000
- Gene ID conversion: HGNC symbols to Entrez IDs via org.Hs.eg.db
- Multiple testing: Benjamini-Hochberg FDR

**Databases (12 total):**

| Database | Source | Pathways | Description |
|----------|--------|----------|-------------|
| Hallmark | MSigDB | 50 | Curated biological states |
| KEGG | MSigDB | 186 | Metabolic & signaling pathways |
| Reactome | MSigDB | 1,615 | Curated pathway reactions |
| GO:BP | MSigDB | 7,658 | Biological processes |
| GO:CC | MSigDB | 1,006 | Cellular components |
| GO:MF | MSigDB | 1,738 | Molecular functions |
| WikiPathways | MSigDB | 664 | Community-curated pathways |
| Canonical | MSigDB | 2,922 | Combined canonical pathways |
| CGP | MSigDB | 3,358 | Chemical/genetic perturbations |
| TF | MSigDB | 1,137 | Transcription factor targets |
| SynGO | SynGO v1.1 | ~300 | Synaptic gene ontology (CC) |
| MitoCarta | MitoCarta 3.0 | 149 | Mitochondrial pathways |

#### 4. Gene Set Variation Analysis (GSVA)

**Software:** GSVA R package (v1.48+)

**Method:**
- Single-sample enrichment scoring (Gaussian kernel)
- Sample-level pathway activity scores
- Enables trajectory analysis across development
- Pathway size filter: 10-500 genes

---

## Trajectory Framework

Contrasts are mapped to a developmental trajectory framework for pattern analysis:

| Original Contrast | Trajectory Name | Stage | Description |
|-------------------|-----------------|-------|-------------|
| G32A_vs_Ctrl_D35 | Early_G32A | Early | Initial mutation effect |
| R403C_vs_Ctrl_D35 | Early_R403C | Early | Initial mutation effect |
| G32A_vs_Ctrl_D65 | Late_G32A | Late | Mature neuron state |
| R403C_vs_Ctrl_D65 | Late_R403C | Late | Mature neuron state |
| Maturation_G32A_specific | TrajDev_G32A | TrajDev | Trajectory deviation |
| Maturation_R403C_specific | TrajDev_R403C | TrajDev | Trajectory deviation |

### Pattern Classifications

Pathways are classified into 8 mutually exclusive patterns based on Early → TrajDev → Late dynamics. The framework distinguishes **active** patterns (requiring significant trajectory deviation) from **passive** patterns (developmental buffering without active compensation).

**Full 8-Pattern Taxonomy:**

| Pattern | Active? | Criteria | Interpretation |
|---------|---------|----------|----------------|
| **Compensation** | Active | Early defect + TrajDev opposes + Late improved | Active adaptive response; system compensates |
| **Sign_reversal** | Active | Early defect + TrajDev opposes + Late opposite sign | Trajectory reversal; defect direction completely flipped |
| **Progressive** | Active | Early defect + TrajDev amplifies + Late worsened | Cumulative damage; active maladaptive response |
| **Natural_worsening** | Passive | Early defect + no TrajDev + Late worsened | Passive deterioration; no adaptive capacity |
| **Natural_improvement** | Passive | Early defect + no TrajDev + Late improved | Passive recovery; developmental buffering |
| **Late_onset** | - | No early defect + Late defect emerges | Maturation-dependent dysfunction |
| **Transient** | - | Early defect + resolved by Late | Developmental delay that recovers |
| **Complex** | - | Non-linear or multiphasic | Requires individual inspection |

**Super-Categories (Simplified for Main Text):**

| Super-Category | Includes | Use Case |
|----------------|----------|----------|
| **Active_Compensation** | Compensation | Main narrative - active adaptive responses |
| **Active_Reversal** | Sign_reversal | Main narrative - trajectory reversal (sign flip) |
| **Active_Progression** | Progressive | Main narrative - active worsening |
| **Passive** | Natural_improvement, Natural_worsening | Main narrative - developmental buffering |
| **Late_onset** | Late_onset | Maturation-dependent effects |
| **Other** | Transient, Complex | Requires individual inspection |

**Canonical reference:** See `docs/PATTERN_CLASSIFICATION.md` for full specifications and `01_Scripts/Python/pattern_definitions.py` for implementation.

---

## Key Findings

### Core Discovery: The Translation Paradox

DRP1 mutations cause a paradoxical energetic-translational crisis:

1. **Ribosome biogenesis INCREASED** (NES +2.25, FDR < 10⁻¹²)
2. **Cytoplasmic translation DECREASED** (NES -2.07, FDR < 10⁻¹²)
3. **Synaptic translation FAILED** (NES -2.9 to -3.0, FDR < 10⁻¹²)

**Mechanism:** Mitochondrial dysfunction causes local ATP depletion at synapses, preventing ribosome function despite compensatory biogenesis.

### Cross-Database Validation

The translation crisis finding replicates across 6 independent databases (50+ significant pathways), confirming robustness of the core observation.

### Mutation-Specific Effects

| Finding | G32A | R403C |
|---------|------|-------|
| Ribosome biogenesis | Stronger (+2.25) | Moderate (+1.59) |
| Synaptic compensation | Strong | "Empty synapses" pattern |
| OXPHOS compensation | Not significant | Strong (+1.60) |

---

## Master Data Tables

Three comprehensive master tables provide complete analysis results:

### Summary Comparison

| Table | Method | Scope | Rows | Size | Use Case |
|-------|--------|-------|------|------|----------|
| **master_gsea_table.csv** | GSEA | All pathways × contrasts | 109,990 | 43 MB | Differential enrichment discovery |
| **master_gsva_focused_table.csv** | GSVA | 7 key modules | 42 | 12 KB | Quick trajectory analysis |
| **master_gsva_all_table.csv** | GSVA | All pathways × groups | 87,001 | 25 MB | Comprehensive trajectory exploration |

### 1. Master GSEA Table

**File:** `master_gsea_table.csv`

**Overview:**
- **109,990 rows**: pathway × contrast combinations
- **12,221 unique pathways** across 12 databases
- **25 columns**: IDs, statistics, pattern classifications
- **43% significant**: At least one contrast shows FDR < 0.05

**Key Columns:**

| Column Group | Key Columns | Description |
|--------------|-------------|-------------|
| **Identification** | pathway_id, database, Description | Pathway metadata |
| **GSEA Statistics** | NES, pvalue, p.adjust, setSize | Enrichment results |
| **Pattern Classification** | Pattern_G32A, Pattern_R403C | One of 8 temporal patterns per mutation |
| **Pattern Confidence** | Confidence_G32A, Confidence_R403C | "High" (p<0.05) or "Medium" (p<0.10) |
| **Super-Categories** | Super_Category_G32A, Super_Category_R403C | Simplified 6-category grouping |
| **Trajectory NES** | NES_Early/TrajDev/Late_[Mutation] | Stage-specific enrichment scores |
| **Pattern Consistency** | change_consistency | Whether G32A and R403C show same pattern |
| **Significance Flags** | ever_significant, ever_significant_trajectory | Any contrast significant |

**Generation Script:** `02_Analysis/4.1.create_master_pathway_table.py`

**Usage Examples:**

```python
import pandas as pd

# Load table
df = pd.read_csv('03_Results/02_Analysis/master_gsea_table.csv')

# Find MitoCarta compensation pathways (high confidence)
mito_comp = df[
    (df['database'] == 'MitoCarta') &
    (df['Pattern_G32A'] == 'Compensation') &
    (df['Confidence_G32A'] == 'High')
]
print(f"Found {len(mito_comp)} MitoCarta compensation pathways")

# Find pathways showing Sign_reversal in either mutation
sign_rev = df[
    ((df['Pattern_G32A'] == 'Sign_reversal') |
     (df['Pattern_R403C'] == 'Sign_reversal')) &
    (df['category'] == 'Early')
].drop_duplicates('pathway_id')
print(f"Found {len(sign_rev)} pathways with trajectory reversal")

# Compare Active vs Passive patterns
active = df[df['Super_Category_G32A'].isin(['Active_Compensation', 'Active_Reversal', 'Active_Progression'])]
passive = df[df['Super_Category_G32A'] == 'Passive']
print(f"Active patterns: {len(active.drop_duplicates('pathway_id'))}")
print(f"Passive patterns: {len(passive.drop_duplicates('pathway_id'))}")
```

### 2. Master GSVA Table (Focused - 7 Key Modules)

**Files:**
- `master_gsva_focused_table.csv` (long format, 42 rows) - Per-group GSVA scores and statistics
- `gsva_pattern_summary.csv` (wide format, 7 rows) - Pattern classifications for each module
- `gsva_statistics_summary.txt` (summary statistics)

**Overview:**
- **7 modules**: Carefully curated biological modules
- **42 rows** (master table): 7 modules × 3 genotypes × 2 timepoints
- **7 rows** (pattern summary): 1 row per module with trajectory-wide pattern classification
- **Small & fast**: Ideal for focused hypothesis testing

**The 7 Trajectory Modules:**

| Module | Source | N Genes | Panel | Biological Function |
|--------|--------|---------|-------|---------------------|
| **Ribosome Biogenesis** | GO:BP | 158 | A | Nuclear ribosome assembly |
| **Cytoplasmic Translation** | GO:BP | 76 | B | Cytoplasmic protein synthesis |
| **Synaptic Ribosomes** | SynGO | 65 | C | Synapse-localized translation |
| **Mitochondrial Ribosome** | MitoCarta | 77 | D | Mitoribosome structure |
| **Mito Ribosome Assembly** | MitoCarta | 24 | E | Mitoribosome biogenesis |
| **mtDNA Maintenance** | MitoCarta | 29 | F | Mitochondrial genome integrity |
| **OXPHOS** | MitoCarta | 139 | G | Oxidative phosphorylation |

**Key Columns (master_gsva_focused_table.csv):**

| Column Group | Key Columns | Description |
|--------------|-------------|-------------|
| **Identification** | Module, Display_Name, Source_Database | Module metadata |
| **GSVA Scores** | Mean_GSVA, SD_GSVA, SE_GSVA, N | Enrichment statistics |
| **Trajectory** | Expression_vs_CtrlD35 | Deviation from baseline |
| **Divergence** | Divergence_vs_Ctrl | Difference from matched control |
| **Statistics** | t_statistic, p_value, p_adjusted, significant | T-test results |

**Key Columns (gsva_pattern_summary.csv):**

| Column Group | Key Columns | Description |
|--------------|-------------|-------------|
| **Identification** | Module, Display_Name, Source_Database | Module metadata |
| **Trajectory Metrics** | Expression_vs_CtrlD35_[Mutation]_[Day] | Deviation from baseline at each stage |
| **Divergence Metrics** | Divergence_vs_Ctrl_[Mutation]_[Day] | Instantaneous difference from control |
| **Significance** | p_adjusted_[Mutation]_[Day], significant_* | Statistical test results |
| **Trajectory Deviation** | TrajDev_G32A, TrajDev_R403C | Calculated as (D65_mut - D35_mut) - (D65_ctrl - D35_ctrl) |
| **Pattern Classification** | Pattern_G32A, Pattern_R403C | One of 8 temporal patterns per mutation |
| **Pattern Confidence** | Confidence_G32A, Confidence_R403C | "High" or "Medium" (or NA for Complex) |
| **Super-Categories** | Super_Category_G32A, Super_Category_R403C | Simplified 6-category grouping |
| **Pattern Consistency** | Change_Consistency | Whether both mutations show same pattern |

**Generation Script:** `02_Analysis/1.7.create_master_gsva_table.R`

**Usage Examples:**

```r
library(dplyr)

# Load long format table (GSVA scores)
df <- read.csv('03_Results/02_Analysis/master_gsva_focused_table.csv')

# View OXPHOS trajectory
df %>%
  filter(Module == 'OXPHOS') %>%
  select(Genotype, Day, Mean_GSVA, Expression_vs_CtrlD35, Divergence_vs_Ctrl)

# Load pattern summary (classifications)
patterns <- read.csv('03_Results/02_Analysis/gsva_pattern_summary.csv')

# Find modules with compensation in G32A
patterns %>%
  filter(Pattern_G32A == 'Compensation') %>%
  select(Display_Name, Pattern_G32A, Confidence_G32A, Super_Category_G32A)

# Compare pattern consistency across mutations
patterns %>%
  select(Display_Name, Pattern_G32A, Pattern_R403C, Change_Consistency)
```

### 3. Master GSVA Table (Comprehensive - All Pathways)

**File:** `master_gsva_all_table.csv`

**Overview:**
- **87,001 rows**: pathway × genotype × timepoint combinations
- **~14,500 unique pathways** (after 10-500 gene size filter)
- **All 12 databases**: Complete coverage for exploration
- **Powers interactive explorer**: Used by `DRP1_Pathway_Explorer.html`

**Key Columns:**

| Column Group | Key Columns | Description |
|--------------|-------------|-------------|
| **Identification** | pathway_id, database, pathway_name, n_genes | Pathway metadata |
| **GSVA Scores** | Mean_GSVA, SD_GSVA, SE_GSVA, N | Enrichment statistics |
| **Trajectory** | Expression_vs_CtrlD35 | Cumulative change from baseline |
| **Divergence** | Divergence_vs_Ctrl | Instantaneous difference from control |
| **Statistics** | t_statistic, p_value, p_adjusted, significant | T-test results (FDR corrected) |

**Generation Script:** `02_Analysis/1.6.gsva_analysis.R`

**Computation Time:** ~15-40 minutes (cached in checkpoint)

**Usage Example:**

```r
# Load table
df <- read.csv('03_Results/02_Analysis/master_gsva_all_table.csv')

# Find G32A pathways with largest D65 divergence
top_divergent <- df %>%
  filter(Genotype == 'G32A', Day == 65, !is.na(Divergence_vs_Ctrl)) %>%
  arrange(desc(abs(Divergence_vs_Ctrl))) %>%
  head(20)
```

### How Pattern Classifications Are Stored

Pattern classifications appear in two places with slightly different structures:

**1. master_gsea_table.csv (GSEA results):**
- One row per pathway × contrast combination (109,990 rows)
- Pattern columns apply to the entire trajectory (Early→TrajDev→Late)
- Each pathway has Pattern_G32A, Pattern_R403C, Super_Category_G32A, Super_Category_R403C
- Pattern is the same across all contrast rows for a given pathway (classification uses all trajectory stages)

**2. gsva_pattern_summary.csv (GSVA focused results):**
- One row per module (7 rows total)
- Contains trajectory-wide metrics: TrajDev_G32A, TrajDev_R403C
- Pattern columns: Pattern_G32A, Pattern_R403C, Confidence_*, Super_Category_*
- Companion to master_gsva_focused_table.csv which has per-group GSVA scores

**3. master_gsva_all_table.csv (GSVA comprehensive):**
- Does NOT contain pattern classifications (87,001 rows would be unwieldy)
- Contains per-group GSVA scores and statistics only
- Use for exploratory trajectory analysis, then apply patterns post-hoc if needed

### GSEA vs GSVA: When to Use Which

| Feature | GSEA (master_gsea_table) | GSVA (focused) | GSVA (all) |
|---------|--------------------------|----------------|------------|
| **Method** | Pathway enrichment (rank-based) | Single-sample scoring | Single-sample scoring |
| **Granularity** | Contrast-level | Sample-level | Sample-level |
| **Scope** | 12,221 pathways | 7 modules | 14,500 pathways |
| **Question** | "Is pathway differentially enriched?" | "How do key modules change?" | "Which pathways show trajectories?" |
| **Advantage** | Statistical power for DE | Fast, focused | Comprehensive exploration |
| **Statistical test** | Permutation test | T-test | T-test |
| **Pattern classifications** | YES (all pathways) | YES (7 modules) | NO (too large) |
| **Best for** | Discovering enriched pathways | Hypothesis testing | Hypothesis generation |

---

## Data Files

### Differential Expression Results

**Location:** `DE_results/`

- **Format:** CSV with gene-level statistics
- **Columns:** gene, logFC, AveExpr, t, P.Value, adj.P.Val, B
- **Files:** 9 contrasts × 2 versions = 18 files

### Checkpoint Files

**Location:** `checkpoints/`

- **Format:** RDS (R serialized objects)
- **Purpose:** Cache expensive computations for reproducibility
- **Total size:** ~140 MB
- **Key files:**
  - `qc_variables.rds` - Expression matrix and metadata
  - `all_gsea_results.rds` - MSigDB GSEA results
  - `mitocarta_gsea_results.rds` - MitoCarta GSEA
  - `syngo_gsea_results.rds` - SynGO GSEA
  - `gsva_module_scores.rds` - GSVA scores (7 modules)
  - `gsva_all_pathways.rds` - GSVA scores (all pathways)

---

## Generating Scripts

| Script | Purpose | Runtime | Output |
|--------|---------|---------|--------|
| `1.1.main_pipeline.R` | Core DE + GSEA | 30-60 min | DE_results/, checkpoints/, Plots/GSEA/, Plots/Volcano/ |
| `1.2.generate_contrast_tables.R` | Extract DE tables | <1 min | DE_results/*.csv |
| `1.3.add_mitocarta.R` | MitoCarta GSEA | 2-5 min | checkpoints/mitocarta_*.rds, Plots/GSEA/*/MitoCarta/ |
| `1.4.export_gsea_for_python.R` | Export for Python | <1 min | Python_exports/*.csv |
| `4.1.create_master_pathway_table.py` | Master GSEA table | <1 min | master_gsea_table.csv |
| `1.7.create_master_gsva_table.R` | Focused GSVA table | <1 min | master_gsva_focused_table.csv |
| `1.6.gsva_analysis.R` | Comprehensive GSVA | 15-40 min | master_gsva_all_table.csv |
| `4.4.prepare_explorer_data.py` | Interactive explorer | 2-5 min | Explorer/DRP1_Pathway_Explorer.html |
| `viz_*.R` | Specialized visualizations | Varies | Plots/* subdirectories |
| `3.1.publication_figures.py` | Publication figures | <5 min | Plots/Publication_Figures/ |

---

## Reproducibility

### Checkpoint Caching

The pipeline uses checkpoint caching to avoid recomputation:

```r
result <- load_or_compute(
  checkpoint_file = "path/to/checkpoint.rds",
  compute_fn = function() { expensive_computation() },
  force_recompute = config$force_recompute,
  description = "Computation name"
)
```

**To force recomputation:** Set `config$force_recompute = TRUE` in `1.1.main_pipeline.R`

### Re-running the Analysis

```bash
# Full pipeline (30-60 min, uses caching)
Rscript 02_Analysis/1.1.main_pipeline.R

# Generate master tables
python3 02_Analysis/4.1.create_master_pathway_table.py          # GSEA master table
Rscript 02_Analysis/1.7.create_master_gsva_table.R             # Focused GSVA table
Rscript 02_Analysis/1.6.gsva_analysis.R          # Comprehensive GSVA table (slow)

# Generate interactive explorer
python3 02_Analysis/4.4.prepare_explorer_data.py

# Visualizations (R)
Rscript 02_Analysis/2.1.viz_ribosome_paradox.R
Rscript 02_Analysis/2.2.viz_mito_translation_cascade.R
Rscript 02_Analysis/2.4.viz_critical_period_trajectories_gsva.R

# Publication figures (Python)
python3 02_Analysis/3.1.publication_figures.py
python3 02_Analysis/3.4.pattern_summary_normalized.py
```

---

## Software Versions

| Package | Version | Purpose |
|---------|---------|---------|
| R | 4.3+ | Statistical computing |
| limma | 3.56+ | Linear modeling |
| edgeR | 3.42+ | Normalization |
| clusterProfiler | 4.8+ | GSEA framework |
| fgsea | 1.26+ | Fast GSEA algorithm |
| GSVA | 1.48+ | Single-sample enrichment |
| msigdbr | 7.5+ | MSigDB gene sets |
| Python | 3.10+ | Visualization |
| pandas | 2.0+ | Data manipulation |
| matplotlib | 3.7+ | Plotting |
| seaborn | 0.12+ | Statistical graphics |

---

## Citation

When using these results, cite:

1. The analysis pipeline and pattern classification framework
2. Relevant pathway databases (MSigDB, SynGO, MitoCarta)
3. Statistical methods (limma, fgsea, GSVA)

**Key References:**
- **limma**: Ritchie et al., Nucleic Acids Research, 2015
- **fgsea**: Korotkevich et al., bioRxiv, 2016
- **GSVA**: Hänzelmann et al., BMC Bioinformatics, 2013
- **MSigDB**: Liberzon et al., Cell Systems, 2015
- **MitoCarta 3.0**: Rath et al., Nucleic Acids Research, 2021
- **SynGO**: Koopmans et al., Neuron, 2019

---

## Related Documentation

- `CLAUDE.md` - Project overview and analysis instructions
- `Plots/README.md` - Comprehensive figure documentation
- `DE_results/README.md` - Differential expression tables
- `checkpoints/README.md` - Checkpoint file documentation
- `Python_exports/README.md` - Python export files
- `CHANGELOG.md` - Analysis change tracking
- `ISSUES.md` - Known issues and limitations
