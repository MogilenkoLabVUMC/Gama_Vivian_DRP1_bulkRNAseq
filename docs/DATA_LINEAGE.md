# Data Lineage: Table Production and Consumption

**Last Updated**: 2025-12-04
**Purpose**: Documents which scripts produce which data tables, how they're consumed downstream, and the overall data architecture.

## Executive Summary

The DRP1 bulk RNA-seq analysis pipeline follows a **two-tier master table architecture**:

1. **GSEA Master Table**: `master_gsea_table.csv` (produced by `1.5.create_master_pathway_table.py`)
   - Primary source for all downstream Python visualizations
   - Contains 110K rows with pattern classifications

2. **GSVA Master Tables**: `master_gsva_*.csv` (produced by `1.7.create_master_gsva_table.R`)
   - Focused table: 7 key modules (42 rows)
   - Comprehensive table: All pathways (87K rows)

**Architecture Status**: Clean "one script produces one master table" design with minimal intermediate files.

---

## Data Flow Architecture

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                            DATA PRODUCTION LAYER                            │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│  ┌─────────────────────┐                                                    │
│  │ 1.1.main_pipeline.R │ ◄── Counts matrix + Metadata                       │
│  │                     │                                                    │
│  │ PRODUCES:           │                                                    │
│  │  • checkpoints/*.rds│ ─────────────────────────────────────────┐         │
│  │  • DE_results/*.csv │                                          │         │
│  └──────────┬──────────┘                                          │         │
│             │                                                     │         │
│             ▼                                                     │         │
│  ┌─────────────────────┐         ┌─────────────────────┐          │         │
│  │ 1.3.add_mitocarta.R │         │ 1.6.gsva_analysis.R │ ◄────────┘         │
│  │                     │         │                     │                    │
│  │ PRODUCES:           │         │ PRODUCES:           │                    │
│  │  • mitocarta_*.rds  │         │  • gsva_*.rds       │                    │
│  └──────────┬──────────┘         └──────────┬──────────┘                    │
│             │                               │                               │
│             ▼                               ▼                               │
│  ┌─────────────────────────────────────────────────────────────┐            │
│  │ 1.4.export_gsea_for_python.R                                │            │
│  │                                                             │            │
│  │ PRODUCES (Python_exports/):                                 │            │
│  │  • gsea_results_long.csv ────────────────┐                  │            │
│  │  • gsea_results_wide.csv ─────────────┐  │                  │            │
│  └───────────────────────────────────────┼──┼──────────────────┘            │
│                                          │  │                               │
└──────────────────────────────────────────┼──┼───────────────────────────────┘
                                           │  │
┌──────────────────────────────────────────┼──┼───────────────────────────────┐
│                        MASTER TABLE LAYER│  │                               │
├──────────────────────────────────────────┼──┼───────────────────────────────┤
│                                          │  │                               │
│  ┌───────────────────────────────────────┘  │                               │
│  │                                          │                               │
│  ▼                                          ▼                               │
│  ┌──────────────────────────────┐    ┌────────────────────────────────┐     │
│  │ 1.5.create_master_pathway_   │    │ 1.7.create_master_gsva_table.R │     │
│  │     table.py                 │    │                                │     │
│  │                              │    │ PRODUCES:                      │     │
│  │ PRODUCES:                    │    │  • master_gsva_focused_*.csv   │     │
│  │  • master_gsea_table.csv ◄───┼────┤  • master_gsva_all_table.csv   │     │
│  │    (50MB, 110K rows)         │    │  • gsva_pattern_summary.csv    │     │
│  └──────────────┬───────────────┘    └────────────────────────────────┘     │
│                 │                                                           │
└─────────────────┼───────────────────────────────────────────────────────────┘
                  │
┌─────────────────┼───────────────────────────────────────────────────────────┐
│                 │              VISUALIZATION LAYER                          │
├─────────────────┼───────────────────────────────────────────────────────────┤
│                 │                                                           │
│                 ├──► 3.1.publication_figures.py (via data_loader.py)        │
│                 ├──► 3.2.publication_figures_dotplot.py (via data_loader)   │
│                 ├──► 3.4.pattern_summary_normalized.py (direct load)        │
│                 ├──► 3.5.viz_trajectory_flow.py (direct + wide CSV)         │
│                 ├──► 3.6.viz_alluvial_ggalluvial.R (direct load)            │
│                 ├──► 3.7.viz_bump_chart.py (direct load)                    │
│                 ├──► 3.8.viz_interactive_bump_dashboard.py                  │
│                 └──► Supp5.prepare_explorer_data.py                         │
│                                                                             │
│  [R viz scripts use RDS checkpoints directly, not CSVs]                     │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

## File Production and Consumption Matrix

### Python Exports (`03_Results/02_Analysis/Python_exports/`)

Produced by: `1.4.export_gsea_for_python.R`

| File | Size | Status | Consumers |
|------|------|--------|-----------|
| `gsea_results_long.csv` | 31MB | **ACTIVE** | `1.5.create_master_pathway_table.py` |
| `gsea_results_wide.csv` | 5MB | **ACTIVE** | `data_loader.py`, `3.5.viz_trajectory_flow.py`, `Supp4.sensitivity_analysis.py` |

### Master Tables (`03_Results/02_Analysis/`)

| File | Size | Producer | Consumers | Purpose |
|------|------|----------|-----------|---------|
| `master_gsea_table.csv` | 50MB | `1.5.create_master_pathway_table.py` | 3.1, 3.2, 3.4, 3.5, 3.6, 3.7, 3.8, data_loader, Supp5 | **Primary downstream source** |
| `master_gsva_focused_table.csv` | 11KB | `1.7.create_master_gsva_table.R` | External reference, explorers | 7 key modules × 6 groups |
| `master_gsva_all_table.csv` | 26MB | `1.7.create_master_gsva_table.R` | `Supp5.prepare_explorer_data.py` | Comprehensive GSVA |
| `gsva_pattern_summary.csv` | 3KB | `1.7.create_master_gsva_table.R` | Reference | Pattern classifications |

### RDS Checkpoints (`03_Results/02_Analysis/checkpoints/`)

| File | Producer | Consumers | Purpose |
|------|----------|-----------|---------|
| `model_objects.rds` | `1.1.main_pipeline.R` | `1.2`, all R viz | Core model objects |
| `contrast_tables.rds` | `1.1.main_pipeline.R` | `1.3`, `1.4`, `1.6`, R viz | DE statistics |
| `all_gsea_results.rds` | `1.1.main_pipeline.R` | `1.4`, `2.x` R viz | GSEA results (MSigDB) |
| `syngo_gsea_results.rds` | `1.1.main_pipeline.R` | `1.4`, `2.x` R viz | SynGO enrichments |
| `mitocarta_gsea_results.rds` | `1.3.add_mitocarta.R` | `1.4`, `2.x` R viz | MitoCarta enrichments |
| `gsva_module_scores.rds` | `1.6.gsva_analysis.R` | `1.7`, `2.4` | GSVA scores |
| `qc_variables.rds` | `1.1.main_pipeline.R` | `1.7`, R viz | Sample metadata |

---

## Critical Dependency Chains

### Chain 1: Main Pipeline
```
1.1.main_pipeline.R
  └─► checkpoints/*.rds
        ├─► 1.3.add_mitocarta.R → mitocarta_gsea_results.rds
        ├─► 1.4.export_gsea_for_python.R → Python_exports/*.csv
        ├─► 1.6.gsva_analysis.R → gsva_module_scores.rds
        └─► 2.x R visualization scripts
```

### Chain 2: GSEA Master Table
```
1.4.export_gsea_for_python.R
  └─► gsea_results_long.csv
        └─► 1.5.create_master_pathway_table.py
              └─► master_gsea_table.csv ─────►┬─► 3.x Python viz
                                              ├─► data_loader.py
                                              └─► Explorers
```

### Chain 3: GSVA Master Tables
```
1.6.gsva_analysis.R
  └─► gsva_module_scores.rds
        └─► 1.7.create_master_gsva_table.R
              ├─► master_gsva_focused_table.csv
              ├─► master_gsva_all_table.csv
              └─► gsva_pattern_summary.csv
```

---

## Data Loader Module Pattern

The `01_Scripts/Python/data_loader.py` module provides centralized data access:

```python
from Python.data_loader import load_classified_pathways

# Loads from master_gsea_table.csv (primary)
# Merges p.adjust columns from gsea_results_wide.csv (secondary)
df = load_classified_pathways()
```

**Why two files are still needed:**
- `master_gsea_table.csv`: Contains NES, patterns, classifications
- `gsea_results_wide.csv`: Contains p.adjust values in trajectory format (p.adjust_Early_G32A, etc.)

---

## Verification Checklist

Before paper submission, verify:

- [ ] `master_gsea_table.csv` is up-to-date (run `1.5.create_master_pathway_table.py`)
- [ ] `master_gsva_*.csv` tables are up-to-date (run `1.7.create_master_gsva_table.R`)
- [ ] Pattern classifications are consistent between GSEA and GSVA tables
- [ ] All visualization outputs use the master tables as source

---

## Appendix: Script Numbering Convention

```
1.x - Core Pipeline (Data Production)
    1.1 - Main DE pipeline
    1.2 - Quick contrast regeneration
    1.3 - MitoCarta enrichment
    1.4 - Export GSEA to CSV (2 files)
    1.5 - Create master GSEA table
    1.6 - GSVA analysis
    1.7 - Create master GSVA tables
    1.8 - SynGO gene extraction

2.x - R Visualizations (Consume RDS checkpoints)
    2.1 - Ribosome paradox
    2.2 - Mito translation cascade
    2.3 - Synaptic ribosomes
    2.4 - Critical period trajectories (GSVA)
    2.5 - Complex V analysis
    2.6 - Calcium genes

3.x - Python/Publication Visualizations (Consume master tables)
    3.1 - Publication figures (heatmaps)
    3.2 - Publication figures (dotplots)
    3.3 - Ribosome UpSet plot
    3.4 - Pattern summary normalized
    3.5 - Trajectory flow (Sankey/alluvial)
    3.6 - Alluvial (ggalluvial)
    3.7 - Bump chart
    3.8 - Interactive dashboard
    3.9 - Pooled dotplots

Supp - Supplementary/Exploratory
    Supp1 - Verify enrichments
    Supp2 - Diagnose calcium genes
    Supp3 - Analyze DE thresholds
    Supp4 - Sensitivity analysis
    Supp5 - Prepare explorer data
    Supp6 - Bump chart explorer app
```

---

## Change Log

| Date | Change |
|------|--------|
| 2025-12-04 | Removed unused intermediate files (gsea_significant, gsea_trajectory*, contrast_mapping) |
| 2025-12-04 | Comprehensive audit of data flow; renamed from SCRIPT_DEPENDENCY_ANALYSIS.md |
| 2025-12-01 | Initial investigation of 1.4 vs 1.5 relationship |
