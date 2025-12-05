# Analysis Scripts Inventory

This document catalogs all analysis scripts, their purposes, and their status (active vs deprecated).

**Last Updated**: 2025-12-04

---

## Active Scripts

### Main Pipeline (Run in Order)

| Script | Purpose | Runtime | Output Location | Dependencies |
|--------|---------|---------|-----------------|--------------|
| **1.1.main_pipeline.R** | Core DE analysis, GSEA, checkpointing | 30-60 min (first run), 5-10 min (cached) | `03_Results/02_Analysis/` | Count matrix, metadata |
| **1.2.generate_contrast_tables.R** | Generate contrast-specific CSV tables | 2-5 min | `03_Results/02_Analysis/DE_results/` | Checkpoints from 1.1 |
| **1.3.add_mitocarta.R** | Annotate results with MitoCarta pathways | 5-10 min | Adds MitoCarta columns to DE tables | MitoCarta 3.0 database |
| **1.4.export_gsea_for_python.R** | Export GSEA results to CSV for Python | 2-3 min | `03_Results/02_Analysis/Python_exports/` | GSEA checkpoints |

**Usage**:
```bash
# Full pipeline
Rscript 02_Analysis/1.1.main_pipeline.R
Rscript 02_Analysis/1.2.generate_contrast_tables.R
Rscript 02_Analysis/1.3.add_mitocarta.R
Rscript 02_Analysis/1.4.export_gsea_for_python.R
```

---

### Visualization Scripts (R)

Run independently after main pipeline completes. Each generates plots in `03_Results/02_Analysis/Plots/[folder]/`.

| Script | Output Folder | Purpose | Key Plots |
|--------|---------------|---------|-----------|
| **2.1.viz_ribosome_paradox.R** | `Ribosome_paradox/` | Translation paradox visualization | Three-pool trajectory plot, temporal dynamics |
| **2.2.viz_mito_translation_cascade.R** | `Mito_translation_cascade/` | Energy→translation→synapse cascade | Mechanistic cascade heatmap, module summary |
| **2.3.viz_synaptic_ribosomes.R** | `Synaptic_ribosomes/` | Pre/postsynaptic ribosome analysis | Expression heatmap, compartment comparison |
| **2.4.viz_critical_period_trajectories_gsva.R** | `Critical_period_trajectories/gsva/` | GSVA temporal trajectories | Pathway trajectory heatmaps across D35-D65 |
| **viz_developmental_framework.R** | `Developmental_Framework/` | Developmental pattern analysis | Maturation-specific pathway patterns |
| **2.5.viz_complex_v_analysis.R** | `Complex_V_analysis/` | ATP synthase (Complex V) deep-dive | Complex V subunit analysis, assembly factors |
| **2.6.viz_calcium_genes.R** | `Calcium_genes/` | Calcium signaling gene analysis | Volcano plots with calcium genes highlighted |
| **viz_pooled_dotplots.R** | `GSEA/Cross_database_pooled/` | Cross-database GSEA dotplots | Pooled dotplots showing top pathways |
| **3.6.viz_alluvial_ggalluvial.R** | `Trajectory_Flow/` | Classical alluvial diagram (ggalluvial) | Flow diagram: Early Defect → Mechanism → Outcome |

**Usage**:
```bash
# Run individual visualizations
Rscript 02_Analysis/2.1.viz_ribosome_paradox.R
Rscript 02_Analysis/2.2.viz_mito_translation_cascade.R
Rscript 02_Analysis/2.3.viz_synaptic_ribosomes.R
Rscript 02_Analysis/3.6.viz_alluvial_ggalluvial.R
# ... (run any script independently)
```

---

### Python Visualization Scripts

Run after R exports complete (`1.4.export_gsea_for_python.R`).

| Script | Output Folder | Purpose | Key Figures |
|--------|---------------|---------|-------------|
| **3.1.publication_figures.py** | `Publication_Figures/` | Main publication-ready figures | Fig1-4 (ribosome paradox, MitoCarta, patterns, semantic overview) |
| **3.2.publication_figures_dotplot.py** | `Publication_Figures_Dotplot/` | Dotplot version of publication figures | Alternative visualization style for same data |
| **3.3.ribosome_upset_plot.py** | `Ribosome_paradox/` | Ribosome gene set overlap UpSet plot | UpSet plot showing pre/post/mito ribosome overlaps |
| **3.4.pattern_summary_normalized.py** | `Pattern_Summary_Normalized/` + `Cross_database_validation/` | Temporal pattern classification and visualization | Trajectory comparative plots, pattern summaries |
| **3.5.viz_trajectory_flow.py** | `Trajectory_Flow/` | Trajectory flow visualization (alluvial diagrams) | Sankey/alluvial diagrams, trajectory heatmaps, violin plots |
| **3.7.viz_bump_chart.py** | `Trajectory_Flow/bump/` | Static bump chart generator (matplotlib) | Pathway trajectory bump charts with 5 variants (uniform, weighted, highlighted, curved, combined) |
| **3.8.viz_interactive_bump.py** | `Trajectory_Flow/bump/` | Interactive bump charts (Plotly) | Standalone HTML interactive bump charts with tooltips and pattern filtering |
| **3.8.viz_interactive_bump_dashboard.py** | `Trajectory_Flow/` | Unified interactive dashboard (Plotly) | Single-page HTML dashboard for pathway trajectory exploration |

**Usage**:
```bash
# Run Python visualizations (after R export)
python3 02_Analysis/3.1.publication_figures.py
python3 02_Analysis/3.2.publication_figures_dotplot.py
python3 02_Analysis/3.3.ribosome_upset_plot.py
python3 02_Analysis/3.4.pattern_summary_normalized.py
python3 02_Analysis/3.5.viz_trajectory_flow.py

# Bump chart visualizations (static and interactive)
python3 02_Analysis/3.7.viz_bump_chart.py          # Static matplotlib bump charts
python3 02_Analysis/3.8.viz_interactive_bump.py    # Interactive Plotly bump charts (multiple HTML files)
python3 02_Analysis/3.8.viz_interactive_bump_dashboard.py  # Unified dashboard (single HTML)
```

#### Bump Chart Variants

The bump chart scripts generate multiple visualization variants to explore pathway trajectories:

**Static Bump Charts (`3.7.viz_bump_chart.py`):**
- **Uniform**: Baseline visualization with equal line weights
- **Weighted**: Line thickness reflects pattern frequency/importance
- **Highlighted**: Top significant pathways labeled with text annotations
- **Curved**: TrajDev visualization using Bezier curves to show trajectory deviation
- **Combined (FINAL)**: Paper-ready version combining weighted + curved + highlighted

**Interactive Bump Charts (`3.8.viz_interactive_bump.py`):**
- All static variants above but with Plotly interactivity
- Tooltips showing pathway name, database, NES values, p-values, pattern classification
- Legend toggling for pattern categories
- Outputs multiple standalone HTML files (one per variant)

**Dashboard (`3.8.viz_interactive_bump_dashboard.py`):**
- Single unified HTML file with all data
- Interactive filtering by mutation, pattern, database
- Toggle between NES and rank y-axes
- Weighted/curved visualization options
- Comprehensive tooltips with statistical details

**Scopes:** All bump charts support two data scopes:
- **Focused**: 7 key trajectory modules (ribosome, mitochondrial, synaptic pathways)
- **Significant**: All pathways meeting statistical thresholds (p.adjust < 0.05)

---

### Data Export & Master Tables

Run after main pipeline to create comprehensive summary tables.

| Script | Purpose | Runtime | Output | Dependencies |
|--------|---------|---------|--------|--------------|
| **1.5.create_master_pathway_table.py** | Create comprehensive GSEA pathway table | 2-5 min | `master_gsea_table.csv` (109K rows) | Python exports, pattern classifications |
| **1.6.gsva_analysis.R** | Run GSVA on all pathways (all databases) | 15-40 min | `checkpoints/gsva_all_pathways.rds` (87K rows) | Expression data, all databases |
| **1.7.create_master_gsva_table.R** | Create comprehensive GSVA master tables | 1-2 min | Multiple GSVA output files | GSVA checkpoint from 1.6, QC variables |

**Usage**:
```bash
# Create master pathway table (GSEA results, all databases)
python3 02_Analysis/1.5.create_master_pathway_table.py

# Run comprehensive GSVA analysis (all pathways, all databases)
Rscript 02_Analysis/1.6.gsva_analysis.R

# Create GSVA master tables (focused + comprehensive)
Rscript 02_Analysis/1.7.create_master_gsva_table.R
```

**Outputs**:
- `master_gsea_table.csv`: All GSEA pathways across all contrasts with pattern classifications (109K rows)
- `master_gsva_focused_table.csv`: GSVA scores for 7 key modules in long format with statistics (42 rows)
- `master_gsva_all_table.csv`: GSVA scores for all pathways across all databases (87K rows)
- `gsva_pattern_summary.csv`: GSVA pattern classifications in wide format (7 modules)
- `checkpoints/gsva_all_pathways.rds`: Cached GSVA scores for all pathways (reusable)
- `gsva_statistics_summary.txt`: Summary statistics and usage examples

**Note:** The script numbering was recently updated:
- `4.1.create_master_pathway_table.py` → `1.5.create_master_pathway_table.py`
- `4.3.comprehensive_gsva_analysis.R` → `1.6.gsva_analysis.R`
- `4.2.create_master_gsva_table.R` → `1.7.create_master_gsva_table.R`

See `03_Results/02_Analysis/README.md` for detailed documentation.

---

### Utility Scripts (Run On-Demand)

| Script | Purpose | When to Use | Output |
|--------|---------|-------------|--------|
| **0.1.runtime_installs.R** | Install missing R packages | Initial setup, package troubleshooting | N/A (package installation) |
| **1.8.extract_syngo_ribosome_genes.R** | Extract SynGO ribosome gene lists | Manual gene set curation | Console output, gene lists |
| **3.9.viz_pooled_dotplots.R** | Generate cross-database pooled dotplots | Visualize top pathways across databases | `03_Results/02_Analysis/Plots/GSEA/Cross_database_pooled/` |
| **regenerate_gsea_plots.R** | Regenerate GSEA plots from checkpoints | Updating plot styles without recomputing GSEA | Replaces existing GSEA plots |
| **Supp1.verify_enrichments.R** | Validate GSEA enrichments, report stats | Debugging, validation checks | `03_Results/02_Analysis/Verification_reports/` |
| **Supp2.diagnose_calcium_genes.R** | Check calcium gene presence in data | Debugging calcium gene analyses | Console output, diagnostic stats |
| **Supp3.analyze_de_thresholds.R** | Analyze effects of different DE thresholds | Exploring threshold sensitivity | Console output, threshold comparison |
| **Supp4.sensitivity_analysis.py** | Pattern classification sensitivity analysis | Validate pattern thresholds | Console output, sensitivity stats |
| **Supp5.prepare_explorer_data.py** | Prepare data for interactive explorer | Generate data files for web explorer | Explorer data files |
| **Supp6.app_bump_chart_explorer.py** | Streamlit app for bump chart exploration | Interactive exploration (run locally) | Web app (localhost) |

**Usage**:
```bash
# Example: Verify enrichments
Rscript 02_Analysis/Supp1.verify_enrichments.R

# Example: Regenerate GSEA plots
Rscript 02_Analysis/regenerate_gsea_plots.R

# Example: Run sensitivity analysis
python3 02_Analysis/Supp4.sensitivity_analysis.py

# Example: Launch interactive explorer (requires Streamlit)
streamlit run 02_Analysis/Supp6.app_bump_chart_explorer.py
```

---

## Pattern Classification Workflow

### Overview

The pattern classification system distinguishes between **active** adaptive responses (requiring significant trajectory deviation) and **passive** developmental changes. This framework is central to interpreting pathway dynamics across neuronal maturation.

### Canonical Reference

**Single source of truth:** `/workspaces/Gama_Vivian_DRP1_bulkRNAseq/01_Scripts/Python/pattern_definitions.py`

This module defines:
- Pattern thresholds (significance, effect size, improvement/worsening ratios)
- Classification logic for 8 trajectory patterns
- Super-category mappings for simplified reporting
- Helper functions for pattern assignment

**Documentation:** See `/workspaces/Gama_Vivian_DRP1_bulkRNAseq/docs/PATTERN_CLASSIFICATION.md` for complete framework description.

### The 8-Pattern Classification System

| Pattern | Active? | Criteria | Biological Interpretation |
|---------|---------|----------|---------------------------|
| **Compensation** | Active | Early defect + TrajDev opposes + Late improved | Adaptive plasticity compensates for mutation |
| **Sign_reversal** | Active | Early defect + TrajDev opposes + Late OPPOSITE SIGN | Trajectory completely reversed defect direction |
| **Progressive** | Active | Early defect + TrajDev amplifies + Late worsened | Cumulative damage via maladaptive response |
| **Natural_worsening** | Passive | Early defect + TrajDev NS + Late worsened | Passive deterioration without compensation |
| **Natural_improvement** | Passive | Early defect + TrajDev NS + Late improved | Developmental buffering without active response |
| **Late_onset** | - | No Early defect + Late defect emerges | Maturation-dependent dysfunction |
| **Transient** | - | Strong Early defect + Late resolved | Developmental delay that recovers |
| **Complex** | - | Inconsistent or multiphasic patterns | Requires individual inspection |

**Key distinction:** Active patterns require **significant TrajDev** (p.adjust < 0.05, |NES| > 0.5), indicating transcriptional plasticity beyond normal development.

### Super-Categories (Simplified Reporting)

For main figures and text, patterns are grouped into 6 super-categories:

| Super-Category | Includes | Use Case |
|----------------|----------|----------|
| **Active_Compensation** | Compensation | Active adaptive responses |
| **Active_Reversal** | Sign_reversal | Trajectory sign reversal |
| **Active_Progression** | Progressive | Active worsening (rare) |
| **Passive** | Natural_improvement, Natural_worsening | Developmental buffering |
| **Late_onset** | Late_onset | Maturation-dependent effects |
| **Other** | Transient, Complex | Individual inspection required |

### Pattern Classification in Scripts

**GSEA (NES-based):**
- Applied in: `1.5.create_master_pathway_table.py`
- Uses: `pattern_definitions.classify_pattern()`
- Output: `master_gsea_table.csv` with pattern columns

**GSVA (enrichment score-based):**
- Applied in: `1.7.create_master_gsva_table.R`
- Uses: R implementation of `classify_pattern()` logic
- Output: `master_gsva_all_table.csv`, `gsva_pattern_summary.csv`

**Synchronization:** Pattern classification logic must stay synchronized between Python and R implementations. See `CLAUDE.md` "Pattern System Synchronization Checklist" for update protocol.

### Trajectory Components

| Stage | GSEA Contrast | Biological Meaning |
|-------|---------------|-------------------|
| **Early** | `{Mutation}_vs_Ctrl_D35` | Initial mutation effect (D35, immature) |
| **TrajDev** | `Maturation_{Mutation}_specific` | Mutation-specific deviation from normal maturation |
| **Late** | `{Mutation}_vs_Ctrl_D65` | Mature state effect (D65) |

**TrajDev formula:**
```
TrajDev = (D65_mut - D35_mut) - (D65_ctrl - D35_ctrl)
```

Significant TrajDev indicates the mutant's maturation actively differs from control maturation.

### Key Thresholds

```python
# Significance
PADJ_SIGNIFICANT = 0.05      # High confidence
PADJ_TRENDING = 0.10         # Medium confidence

# Effect size (GSEA NES)
NES_EFFECT = 0.5             # Minimum detectable effect
NES_STRONG = 1.0             # Strong defect threshold

# Magnitude changes
IMPROVEMENT_RATIO = 0.7      # ≥30% reduction = improved
WORSENING_RATIO = 1.3        # ≥30% increase = worsened

# GSVA-specific (scaled for -1 to 1 range)
GSVA_EFFECT = 0.15           # Equivalent to NES 0.5
GSVA_STRONG = 0.30           # Equivalent to NES 1.0
```

### Accessing Pattern Functions

**Python:**
```python
from Python.pattern_definitions import (
    classify_pattern,           # Single pathway classification
    add_pattern_classification, # Batch classification for dataframe
    add_super_category_columns, # Add super-category columns
    get_pattern_colors,         # Get color mappings
    MEANINGFUL_PATTERNS,        # List of meaningful patterns
    ACTIVE_PATTERNS,            # List of active patterns
    PATTERN_DEFINITIONS         # Full pattern metadata
)
```

**R:**
```r
# Pattern classification implemented directly in:
# 02_Analysis/1.7.create_master_gsva_table.R
# Function: classify_gsva_pattern()
```

### Visualization Support

Pattern colors are consistent across all visualization scripts:

**Static visualizations:** Use `get_pattern_colors()` from `pattern_definitions.py`
**Interactive visualizations:** Colors embedded in HTML dashboards
**Cross-script consistency:** All scripts import from central definition

---

### Helper Scripts (01_Scripts/R_scripts/)

These are sourced by other scripts, not run directly.

| Script | Purpose | Sourced By |
|--------|---------|------------|
| **run_syngo_gsea.R** | SynGO GSEA wrapper function | `1.1.main_pipeline.R` |
| **run_syngo_only.R** | Standalone SynGO GSEA (recovery) | Manual recovery if pipeline fails |
| **parse_mitocarta_gmx.R** | Parse MitoCarta .gmx files | Visualization scripts |
| **run_mitocarta_gsea.R** | MitoCarta GSEA wrapper | `1.1.main_pipeline.R` |
| **read_count_matrix.R** | Load and validate count matrix | `1.1.main_pipeline.R` |
| **generate_vertical_volcanos.R** | Generate vertical volcano plots | `1.1.main_pipeline.R` |
| **syngo_running_sum_plot.R** | SynGO-specific running sum plots | Visualization scripts |

---

## Deprecated Scripts

Moved to `02_Analysis/.deprecated/` - no longer used in active pipeline.

| Script | Deprecated Date | Reason | Replaced By |
|--------|-----------------|--------|-------------|
| **4.visualize_trajectory_patterns.py** | 2025-11-24 | Sparse heatmaps, superseded | `3.1.publication_figures.py` |
| **5.focused_pathway_analysis.py** | 2025-11-24 | Intermediate analysis, superseded | `3.1.publication_figures.py` |
| **viz_critical_period_trajectories.R** | 2025-11-20 | Old trajectory method | `2.4.viz_critical_period_trajectories_gsva.R` (GSVA-based) |
| **viz_critical_period_trajectories_with_mitocarta.R** | 2025-11-20 | MitoCarta variant, merged | `2.4.viz_critical_period_trajectories_gsva.R` |
| **viz_temporal_trajectory.R** | 2025-11-20 | Redundant trajectory viz | `3.1.publication_figures.py` |
| **regenerate_volcanoes.R** | 2025-11-20 | Manual volcano regeneration | Integrated into `1.1.main_pipeline.R` |
| **Analysis_pipeline.R** | 2025-11-19 | Old pipeline version | `1.1.main_pipeline.R` (renamed, restructured) |
| **experimental.R** | Unknown | Experimental code, not used | N/A |

**Note**: Deprecated scripts are kept for reference but should not be used. See `.deprecated/` folder for historical code.

---

## Script Dependencies

### Dependency Graph

```
1.1.main_pipeline.R
    ├── Generates checkpoints/*.rds
    └── Depends on:
        ├── 01_Scripts/R_scripts/run_syngo_gsea.R
        ├── 01_Scripts/R_scripts/run_mitocarta_gsea.R
        ├── 01_Scripts/R_scripts/read_count_matrix.R
        └── 01_Scripts/RNAseq-toolkit/* (submodule)

1.2.generate_contrast_tables.R
    └── Depends on: checkpoints from 1.1

1.3.add_mitocarta.R
    └── Depends on: DE results from 1.2

1.4.export_gsea_for_python.R
    └── Depends on: GSEA checkpoints from 1.1

1.5.create_master_pathway_table.py
    ├── Depends on: Python exports from 1.4
    └── Uses: 01_Scripts/Python/pattern_definitions.py

1.6.gsva_analysis.R
    └── Depends on: Expression data, GSVA checkpoint from 1.1

1.7.create_master_gsva_table.R
    └── Depends on: GSVA checkpoint from 1.6, QC variables checkpoint from 1.1

1.8.extract_syngo_ribosome_genes.R
    └── Depends on: SynGO database, checkpoints from 1.1

All 2.x viz_*.R scripts
    └── Depend on: checkpoints from 1.1

All 3.x Python visualization scripts
    ├── Depend on: Python exports from 1.4
    └── Use: 01_Scripts/Python/pattern_definitions.py

3.6.viz_alluvial_ggalluvial.R
    └── Depends on: GSVA checkpoint from 1.1 or 1.6

3.7.viz_bump_chart.py
    ├── Depends on: master_gsea_table.csv from 1.5
    └── Uses: 01_Scripts/Python/viz_bump_charts.py

3.8.viz_interactive_bump*.py (both variants)
    ├── Depends on: master_gsea_table.csv from 1.5
    └── Uses: 01_Scripts/Python/viz_bump_charts.py

3.9.viz_pooled_dotplots.R
    └── Depends on: GSEA checkpoints from 1.1
```

### Checkpoint Files

Generated by `1.1.main_pipeline.R`, used by downstream scripts:

| Checkpoint | Purpose | Used By |
|------------|---------|---------|
| `fit_object.rds` | Limma fit object | Most viz scripts |
| `de_results.rds` | DE test results | 1.2, viz scripts |
| `all_gsea_results.rds` | Complete GSEA results | All viz scripts, Python export |
| `syngo_gsea_results.rds` | SynGO-specific GSEA | SynGO viz scripts |
| `mitocarta_gsea_results.rds` | MitoCarta GSEA | MitoCarta viz scripts |
| `voom_object.rds` | Voom transformation | QC, diagnostics |
| `DGE_object.rds` | DGEList object | Normalization checks |
| `contrasts_matrix.rds` | Contrast definitions | 1.2, reference |
| `gene_intersections.rds` | Gene overlap data | Ribosome overlap analyses |
| `gsva_module_scores.rds` | GSVA scores | Trajectory viz, 4.2 master table |
| `qc_variables.rds` | QC variables and sample metadata | GSVA master table generation |

---

## Script Naming Conventions

### Active Scripts

- **Main pipeline**: `[major].[minor].name_with_underscores.R/py`
  - Example: `1.1.main_pipeline.R`, `1.4.export_gsea_for_python.R`
- **Visualization**: `[major].[minor].viz_descriptive_name.R`
  - Example: `2.1.viz_ribosome_paradox.R`
- **Python**: `[major].[minor].descriptive_name.py`
  - Example: `3.1.publication_figures.py`
- **Helpers**: `function_name.R` (verb-noun format)
  - Example: `run_syngo_gsea.R`, `parse_mitocarta_gmx.R`

### Deprecated Scripts

- Original naming preserved
- Moved to `.deprecated/` folder
- May have `.bak` extension removed for cleaner organization

---

## Adding New Scripts

When adding a new analysis script:

1. **Follow naming convention**: Use appropriate prefix (`viz_`, number, etc.)
2. **Add script header**: Include purpose, output, runtime, dependencies
3. **Document here**: Add to appropriate section in this file
4. **Update CLAUDE.md**: Add to relevant sections if it's a major pipeline script
5. **Create output README**: If script generates plots, create `README.md` in output folder

### Script Header Template

```r
###############################################################################
## Script: viz_[name].R
## Purpose: [One-line description]
## Output: 03_Results/02_Analysis/Plots/[folder]/
## Runtime: ~[X] minutes
## Dependencies: Checkpoints from 1.1.main_pipeline.R
## Author: [Name]
## Date: [YYYY-MM-DD]
## See: README in output folder for plot descriptions
###############################################################################
```

---

## Quick Reference

### Complete Pipeline Run Order

```bash
# ============================================================================
# PHASE 1: Core Analysis Pipeline
# ============================================================================
Rscript 02_Analysis/1.1.main_pipeline.R                   # 30-60 min first run
Rscript 02_Analysis/1.2.generate_contrast_tables.R        # 2-5 min
Rscript 02_Analysis/1.3.add_mitocarta.R                   # 5-10 min
Rscript 02_Analysis/1.4.export_gsea_for_python.R          # 2-3 min

# ============================================================================
# PHASE 2: Master Tables (Optional but recommended)
# ============================================================================
python3 02_Analysis/1.5.create_master_pathway_table.py    # 2-5 min (GSEA master table)
Rscript 02_Analysis/1.6.gsva_analysis.R                   # 15-40 min (GSVA all pathways)
Rscript 02_Analysis/1.7.create_master_gsva_table.R        # 1-2 min (GSVA master tables)

# ============================================================================
# PHASE 3: R Visualizations (run in parallel or sequentially)
# ============================================================================
Rscript 02_Analysis/2.1.viz_ribosome_paradox.R
Rscript 02_Analysis/2.2.viz_mito_translation_cascade.R
Rscript 02_Analysis/2.3.viz_synaptic_ribosomes.R
Rscript 02_Analysis/2.4.viz_critical_period_trajectories_gsva.R
Rscript 02_Analysis/2.5.viz_complex_v_analysis.R
Rscript 02_Analysis/2.6.viz_calcium_genes.R

# ============================================================================
# PHASE 4: Python Visualizations (after R exports complete)
# ============================================================================
python3 02_Analysis/3.1.publication_figures.py
python3 02_Analysis/3.2.publication_figures_dotplot.py
python3 02_Analysis/3.3.ribosome_upset_plot.py
python3 02_Analysis/3.4.pattern_summary_normalized.py
python3 02_Analysis/3.5.viz_trajectory_flow.py

# ============================================================================
# PHASE 5: Advanced Visualizations (after master tables generated)
# ============================================================================
Rscript 02_Analysis/3.6.viz_alluvial_ggalluvial.R         # Classical alluvial
python3 02_Analysis/3.7.viz_bump_chart.py                 # Static bump charts
python3 02_Analysis/3.8.viz_interactive_bump.py           # Interactive bump charts
python3 02_Analysis/3.8.viz_interactive_bump_dashboard.py # Unified dashboard

# Optional utilities
Rscript 02_Analysis/3.9.viz_pooled_dotplots.R             # Cross-database dotplots
```

**Total runtime**:
- Core pipeline (Phase 1): ~40-75 min first run, ~10-20 min cached
- Master tables (Phase 2): ~20-50 min first run, ~3-7 min cached
- All visualizations (Phases 3-5): ~15-30 min
- **Full pipeline**: ~1.5-2.5 hours (first run), ~30-60 min (subsequent runs)

---

**For additional documentation**:
- Pipeline overview: `README.md` (root)
- Claude Code instructions: `CLAUDE.md`
- Visualization guide: `02_Analysis/README_visualization_scripts.md`
- Plot-specific docs: `03_Results/02_Analysis/Plots/[folder]/README.md`
