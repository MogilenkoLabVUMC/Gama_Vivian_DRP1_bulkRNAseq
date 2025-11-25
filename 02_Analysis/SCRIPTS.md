# Analysis Scripts Inventory

This document catalogs all analysis scripts, their purposes, and their status (active vs deprecated).

**Last Updated**: 2025-11-25

---

## Active Scripts

### Main Pipeline (Run in Order)

| Script | Purpose | Runtime | Output Location | Dependencies |
|--------|---------|---------|-----------------|--------------|
| **1a.Main_pipeline.R** | Core DE analysis, GSEA, checkpointing | 30-60 min (first run), 5-10 min (cached) | `03_Results/02_Analysis/` | Count matrix, metadata |
| **1b.generate_contrast_tables.R** | Generate contrast-specific CSV tables | 2-5 min | `03_Results/02_Analysis/DE_results/` | Checkpoints from 1a |
| **2.add_MitoCarta.R** | Annotate results with MitoCarta pathways | 5-10 min | Adds MitoCarta columns to DE tables | MitoCarta 3.0 database |
| **3.export_gsea_for_python.R** | Export GSEA results to CSV for Python | 2-3 min | `03_Results/02_Analysis/Python_exports/` | GSEA checkpoints |

**Usage**:
```bash
# Full pipeline
Rscript 02_Analysis/1a.Main_pipeline.R
Rscript 02_Analysis/1b.generate_contrast_tables.R
Rscript 02_Analysis/2.add_MitoCarta.R
Rscript 02_Analysis/3.export_gsea_for_python.R
```

---

### Visualization Scripts (R)

Run independently after main pipeline completes. Each generates plots in `03_Results/02_Analysis/Plots/[folder]/`.

| Script | Output Folder | Purpose | Key Plots |
|--------|---------------|---------|-----------|
| **viz_ribosome_paradox.R** | `Ribosome_paradox/` | Translation paradox visualization | Three-pool trajectory plot, temporal dynamics |
| **viz_mito_translation_cascade.R** | `Mito_translation_cascade/` | Energy→translation→synapse cascade | Mechanistic cascade heatmap, module summary |
| **viz_synaptic_ribosomes.R** | `Synaptic_ribosomes/` | Pre/postsynaptic ribosome analysis | Expression heatmap, compartment comparison |
| **viz_critical_period_trajectories_gsva.R** | `Critical_period_trajectories/gsva/` | GSVA temporal trajectories | Pathway trajectory heatmaps across D35-D65 |
| **viz_developmental_framework.R** | `Developmental_Framework/` | Developmental pattern analysis | Maturation-specific pathway patterns |
| **viz_complex_v_analysis.R** | `Complex_V_analysis/` | ATP synthase (Complex V) deep-dive | Complex V subunit analysis, assembly factors |
| **viz_calcium_genes.R** | `Calcium_genes/` | Calcium signaling gene analysis | Volcano plots with calcium genes highlighted |
| **viz_pooled_dotplots.R** | `GSEA/Cross_database_pooled/` | Cross-database GSEA dotplots | Pooled dotplots showing top pathways |

**Usage**:
```bash
# Run individual visualizations
Rscript 02_Analysis/viz_ribosome_paradox.R
Rscript 02_Analysis/viz_mito_translation_cascade.R
Rscript 02_Analysis/viz_synaptic_ribosomes.R
# ... (run any script independently)
```

---

### Python Visualization Scripts

Run after R exports complete (`3.export_gsea_for_python.R`).

| Script | Output Folder | Purpose | Key Figures |
|--------|---------------|---------|-------------|
| **6.publication_figures.py** | `Publication_Figures/` | Main publication-ready figures | Fig1-4 (ribosome paradox, MitoCarta, patterns, semantic overview) |
| **6b.publication_figures_dotplot.py** | `Publication_Figures_Dotplot/` | Dotplot version of publication figures | Alternative visualization style for same data |
| **7.ribosome_upset_plot.py** | `Ribosome_paradox/` | Ribosome gene set overlap UpSet plot | UpSet plot showing pre/post/mito ribosome overlaps |
| **8.pattern_summary_normalized.py** | `Pattern_Summary_Normalized/` + `Cross_database_validation/` | Temporal pattern classification and visualization | Trajectory comparative plots, pattern summaries |

**Usage**:
```bash
# Run Python visualizations (after R export)
python3 02_Analysis/6.publication_figures.py
python3 02_Analysis/6b.publication_figures_dotplot.py
python3 02_Analysis/7.ribosome_upset_plot.py
python3 02_Analysis/8.pattern_summary_normalized.py
```

---

### Utility Scripts (Run On-Demand)

| Script | Purpose | When to Use | Output |
|--------|---------|-------------|--------|
| **verify_enrichments.R** | Validate GSEA enrichments, report stats | Debugging, validation checks | `03_Results/02_Analysis/Verification_reports/` |
| **diagnose_calcium_genes.R** | Check calcium gene presence in data | Debugging calcium gene analyses | Console output, diagnostic stats |
| **analyze_DE_thresholds.R** | Analyze effects of different DE thresholds | Exploring threshold sensitivity | Console output, threshold comparison |
| **extract_syngo_ribosome_genes.R** | Extract SynGO ribosome gene lists | Manual gene set curation | Console output, gene lists |
| **runtime_installs.R** | Install missing R packages | Initial setup, package troubleshooting | N/A (package installation) |

**Usage**:
```bash
# Example: Verify enrichments
Rscript 02_Analysis/verify_enrichments.R

# Example: Diagnose calcium genes
Rscript 02_Analysis/diagnose_calcium_genes.R
```

---

### Helper Scripts (01_Scripts/R_scripts/)

These are sourced by other scripts, not run directly.

| Script | Purpose | Sourced By |
|--------|---------|------------|
| **run_syngo_gsea.R** | SynGO GSEA wrapper function | `1a.Main_pipeline.R` |
| **run_syngo_only.R** | Standalone SynGO GSEA (recovery) | Manual recovery if pipeline fails |
| **parse_mitocarta_gmx.R** | Parse MitoCarta .gmx files | Visualization scripts |
| **run_mitocarta_gsea.R** | MitoCarta GSEA wrapper | `1a.Main_pipeline.R` |
| **read_count_matrix.R** | Load and validate count matrix | `1a.Main_pipeline.R` |
| **generate_vertical_volcanos.R** | Generate vertical volcano plots | `1a.Main_pipeline.R` |
| **syngo_running_sum_plot.R** | SynGO-specific running sum plots | Visualization scripts |

---

## Deprecated Scripts

Moved to `02_Analysis/.deprecated/` - no longer used in active pipeline.

| Script | Deprecated Date | Reason | Replaced By |
|--------|-----------------|--------|-------------|
| **4.visualize_trajectory_patterns.py** | 2025-11-24 | Sparse heatmaps, superseded | `6.publication_figures.py` |
| **5.focused_pathway_analysis.py** | 2025-11-24 | Intermediate analysis, superseded | `6.publication_figures.py` |
| **viz_critical_period_trajectories.R** | 2025-11-20 | Old trajectory method | `viz_critical_period_trajectories_gsva.R` (GSVA-based) |
| **viz_critical_period_trajectories_with_mitocarta.R** | 2025-11-20 | MitoCarta variant, merged | `viz_critical_period_trajectories_gsva.R` |
| **viz_temporal_trajectory.R** | 2025-11-20 | Redundant trajectory viz | `6.publication_figures.py` |
| **regenerate_volcanoes.R** | 2025-11-20 | Manual volcano regeneration | Integrated into `1a.Main_pipeline.R` |
| **Analysis_pipeline.R** | 2025-11-19 | Old pipeline version | `1a.Main_pipeline.R` (renamed, restructured) |
| **experimental.R** | Unknown | Experimental code, not used | N/A |

**Note**: Deprecated scripts are kept for reference but should not be used. See `.deprecated/` folder for historical code.

---

## Script Dependencies

### Dependency Graph

```
1a.Main_pipeline.R
    ├── Generates checkpoints/*.rds
    └── Depends on:
        ├── 01_Scripts/R_scripts/run_syngo_gsea.R
        ├── 01_Scripts/R_scripts/run_mitocarta_gsea.R
        ├── 01_Scripts/R_scripts/read_count_matrix.R
        └── 01_Scripts/RNAseq-toolkit/* (submodule)

1b.generate_contrast_tables.R
    └── Depends on: checkpoints from 1a

2.add_MitoCarta.R
    └── Depends on: DE results from 1b

3.export_gsea_for_python.R
    └── Depends on: GSEA checkpoints from 1a

All viz_*.R scripts
    └── Depend on: checkpoints from 1a

All Python scripts (6.py, 6b.py, 7.py, 8.py)
    └── Depend on: Python exports from 3.export_gsea_for_python.R
```

### Checkpoint Files

Generated by `1a.Main_pipeline.R`, used by downstream scripts:

| Checkpoint | Purpose | Used By |
|------------|---------|---------|
| `fit_object.rds` | Limma fit object | Most viz scripts |
| `de_results.rds` | DE test results | 1b, viz scripts |
| `all_gsea_results.rds` | Complete GSEA results | All viz scripts, Python export |
| `syngo_gsea_results.rds` | SynGO-specific GSEA | SynGO viz scripts |
| `mitocarta_gsea_results.rds` | MitoCarta GSEA | MitoCarta viz scripts |
| `voom_object.rds` | Voom transformation | QC, diagnostics |
| `DGE_object.rds` | DGEList object | Normalization checks |
| `contrasts_matrix.rds` | Contrast definitions | 1b, reference |
| `gene_intersections.rds` | Gene overlap data | Ribosome overlap analyses |
| `gsva_module_scores.rds` | GSVA scores | Trajectory viz |

---

## Script Naming Conventions

### Active Scripts

- **Main pipeline**: `[number][letter].Name_with_underscores.R/py`
  - Example: `1a.Main_pipeline.R`, `3.export_gsea_for_python.R`
- **Visualization**: `viz_descriptive_name.R`
  - Example: `viz_ribosome_paradox.R`
- **Python**: `[number][letter].descriptive_name.py`
  - Example: `6.publication_figures.py`
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
## Dependencies: Checkpoints from 1a.Main_pipeline.R
## Author: [Name]
## Date: [YYYY-MM-DD]
## See: README in output folder for plot descriptions
###############################################################################
```

---

## Quick Reference

### Complete Pipeline Run Order

```bash
# 1. Main analysis pipeline
Rscript 02_Analysis/1a.Main_pipeline.R          # 30-60 min first run
Rscript 02_Analysis/1b.generate_contrast_tables.R  # 2-5 min
Rscript 02_Analysis/2.add_MitoCarta.R           # 5-10 min
Rscript 02_Analysis/3.export_gsea_for_python.R  # 2-3 min

# 2. R visualizations (run in parallel or sequentially)
Rscript 02_Analysis/viz_ribosome_paradox.R
Rscript 02_Analysis/viz_mito_translation_cascade.R
Rscript 02_Analysis/viz_synaptic_ribosomes.R
Rscript 02_Analysis/viz_critical_period_trajectories_gsva.R
Rscript 02_Analysis/viz_developmental_framework.R
Rscript 02_Analysis/viz_complex_v_analysis.R
Rscript 02_Analysis/viz_calcium_genes.R
Rscript 02_Analysis/viz_pooled_dotplots.R

# 3. Python visualizations (after R exports)
python3 02_Analysis/6.publication_figures.py
python3 02_Analysis/6b.publication_figures_dotplot.py
python3 02_Analysis/7.ribosome_upset_plot.py
python3 02_Analysis/8.pattern_summary_normalized.py
```

**Total runtime**: ~1-2 hours (first run with caching), ~20-30 minutes (subsequent runs)

---

**For additional documentation**:
- Pipeline overview: `README.md` (root)
- Claude Code instructions: `CLAUDE.md`
- Visualization guide: `02_Analysis/README_visualization_scripts.md`
- Plot-specific docs: `03_Results/02_Analysis/Plots/[folder]/README.md`
