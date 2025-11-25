# Differential Expression Analysis Results

**Project:** DRP1 Mutation Bulk RNA-seq Analysis
**Last Updated:** 2025-11-25
**Pipeline Version:** Post-migration validated (v1.0+)

---

## Overview

This directory contains the complete differential expression and pathway enrichment analysis of iPSC-derived cortical neurons carrying DRP1 mutations (G32A and R403C) compared to isogenic controls at two developmental timepoints (Day 35 and Day 65).

### Experimental Design

| Factor | Levels | Description |
|--------|--------|-------------|
| **Genotype** | Control, G32A, R403C | DRP1 mutation status |
| **Timepoint** | D35, D65 | Days of neuronal maturation |
| **Replicates** | n=3 per condition | Biological replicates |
| **Total samples** | 18 | 3 genotypes x 2 timepoints x 3 replicates |

### Mutation Biology

- **G32A**: GTPase domain mutation affecting DRP1 enzymatic activity
- **R403C**: Stalk domain mutation affecting oligomerization and membrane binding

---

## Directory Structure

```
03_Results/02_Analysis/
├── README.md                    # This file
├── MASTER_TABLE_README.md       # Documentation for master_pathway_table.csv
│
├── DE_results/                  # Differential expression statistics (9 contrasts)
│   └── README.md
│
├── checkpoints/                 # Cached R objects for reproducibility
│   └── README.md
│
├── Plots/                       # All visualizations
│   ├── README.md               # Comprehensive figure overview
│   ├── General/                # QC plots (MDS, correlation, UpSet)
│   ├── Volcano/                # Volcano plots (8 variants)
│   ├── GSEA/                   # Pathway enrichment (12 databases x 9 contrasts)
│   ├── Publication_Figures/    # Final manuscript figures (Python)
│   ├── Ribosome_paradox/       # Core finding visualization
│   ├── Mito_translation_cascade/
│   ├── Synaptic_ribosomes/
│   ├── Critical_period_trajectories/
│   ├── Cross_database_validation/
│   └── Pattern_Summary_Normalized/
│
├── Calcium_genes/               # Calcium signaling gene analysis
│   └── README.md
│
├── Python_exports/              # Data exports for Python visualization
│   └── README.md
│
├── Verification_reports/        # QC reports and gene lists
│   └── README.md
│
├── master_pathway_table.csv     # Comprehensive GSEA results (110K rows)
└── Summary/                     # DE gene count summaries
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
| Parameter | Value | Rationale |
|-----------|-------|-----------|
| p_cutoff | 0.05 | Standard significance threshold |
| fc_cutoff | 2 (|log2FC| > 1) | Biologically meaningful fold-change |
| FDR threshold | 0.05 | Benjamini-Hochberg adjusted |

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

Pathways are classified based on Early → TrajDev → Late dynamics:

| Pattern | Criteria | Interpretation |
|---------|----------|----------------|
| **Compensation** | Early defect, TrajDev opposes, Late improved | Active adaptive response |
| **Progressive** | Early defect, TrajDev amplifies, Late worsened | Cumulative damage |
| **Natural_worsening** | Early defect, no TrajDev, Late worsened | Passive deterioration |
| **Natural_improvement** | Early defect, no TrajDev, Late improved | Passive recovery |
| **Late_onset** | No early defect, Late defect emerges | Maturation-dependent |
| **Transient** | Early defect, resolved by Late | Developmental delay |
| **Persistent** | Stable defect, no trajectory change | Unchanging dysfunction |
| **Complex** | Non-linear pattern | Multi-phase dynamics |

---

## Key Findings

### Core Discovery: The Translation Paradox

DRP1 mutations cause a paradoxical energetic-translational crisis:

1. **Ribosome biogenesis INCREASED** (NES +2.25, FDR < 10^-12)
2. **Cytoplasmic translation DECREASED** (NES -2.07, FDR < 10^-12)
3. **Synaptic translation FAILED** (NES -2.9 to -3.0, FDR < 10^-12)

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

## Data Files

### Master Data Table

**File:** `master_pathway_table.csv` (44 MB)

- **Rows:** 109,989 (pathway x contrast combinations)
- **Unique pathways:** 12,221
- **Columns:** 25 (IDs, statistics, classifications)
- **Documentation:** See `MASTER_TABLE_README.md`

### DE Result Tables

**Location:** `DE_results/`

- **Format:** CSV with gene-level statistics
- **Columns:** gene, logFC, AveExpr, t, P.Value, adj.P.Val, B
- **Files:** 9 contrasts x 2 versions = 18 files

### Checkpoint Files

**Location:** `checkpoints/`

- **Format:** RDS (R serialized objects)
- **Purpose:** Cache expensive computations for reproducibility
- **Total size:** ~140 MB

---

## Generating Scripts

| Script | Purpose | Output |
|--------|---------|--------|
| `02_Analysis/1a.Main_pipeline.R` | Core DE + GSEA | DE_results/, checkpoints/, Plots/GSEA/, Plots/Volcano/ |
| `02_Analysis/1b.generate_contrast_tables.R` | Extract DE tables | DE_results/*.csv |
| `02_Analysis/2.add_MitoCarta.R` | MitoCarta GSEA | checkpoints/mitocarta_*.rds, Plots/GSEA/*/MitoCarta/ |
| `02_Analysis/3.export_gsea_for_python.R` | Export for Python | Python_exports/*.csv |
| `02_Analysis/viz_*.R` | Specialized visualizations | Plots/* subdirectories |
| `02_Analysis/6.publication_figures.py` | Publication figures | Plots/Publication_Figures/ |
| `02_Analysis/9.create_master_pathway_table.py` | Master table | master_pathway_table.csv |

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

**To force recomputation:** Set `config$force_recompute = TRUE` in `1a.Main_pipeline.R`

### Re-running the Analysis

```bash
# Full pipeline (30-60 min, uses caching)
Rscript 02_Analysis/1a.Main_pipeline.R

# Supplementary modules
Rscript 02_Analysis/1b.generate_contrast_tables.R
Rscript 02_Analysis/2.add_MitoCarta.R
Rscript 02_Analysis/3.export_gsea_for_python.R

# Visualizations (R)
Rscript 02_Analysis/viz_ribosome_paradox.R
Rscript 02_Analysis/viz_mito_translation_cascade.R
# ... other viz_*.R scripts

# Publication figures (Python)
python3 02_Analysis/6.publication_figures.py
python3 02_Analysis/8.pattern_summary_normalized.py
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
3. Statistical methods (limma, fgsea)

---

## Related Documentation

- `CLAUDE.md` - Project overview and analysis instructions
- `Plots/README.md` - Comprehensive figure documentation
- `MASTER_TABLE_README.md` - Master pathway table schema
- Individual folder READMEs for specific outputs
