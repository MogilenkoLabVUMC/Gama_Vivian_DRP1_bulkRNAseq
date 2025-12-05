# Visualization Overview: DRP1 Mutation Analysis

**Last Updated:** 2025-12-05
**Analysis:** Bulk RNA-seq of DRP1 mutant cortical neurons (G32A, R403C) at D35 and D65

---

## Overview

This directory contains all visualization outputs from the DRP1 mutation RNA-seq analysis pipeline. The visualizations are organized into subdirectories by analysis type, with each folder containing its own detailed README.

**Core Finding:** DRP1 mutations cause a translation paradox where ribosome biogenesis increases while synaptic translation fails, representing a fundamental energy-translation coupling failure during neuronal maturation.

---

## Directory Structure

```
03_Results/02_Analysis/Plots/
├── README.md (this file)
│
├── Chord_Diagrams/              # NEW: Gene-pathway chord diagrams (GOChord style)
├── Complex_V_analysis/          # ATP synthase (Complex V) deep-dive
├── Critical_period_trajectories/ # Temporal trajectory patterns (D35→D65)
├── General/                     # General exploratory plots
├── GSEA/                        # Gene Set Enrichment Analysis results
├── Mito_translation_cascade/    # Mitochondrial-translation cascade heatmaps
├── Pattern_Summary_Normalized/  # Pattern classification visualizations
├── Publication_Figures/         # Main manuscript figures (heatmap versions)
├── Publication_Figures_Dotplot/ # Dotplot versions of publication figures
├── Ribosome_paradox/           # Core ribosome paradox finding
├── Synaptic_ribosomes/         # Pre/postsynaptic translation analysis
├── Trajectory_Flow/            # Alluvial and bump chart trajectory visualizations
└── Volcano/                    # Volcano plots for differential expression
```

---

## Key Subdirectories

### Chord_Diagrams/ [NEW]
**Generator Script:** `02_Analysis/3.7.viz_chord_diagrams.py`
**Purpose:** GOChord-style chord diagrams showing gene-pathway leading edge membership for synaptic compartments

**Contents:**
- Gene-pathway relationship visualizations for G32A and R403C mutations
- Stacked rectangles show D35/D65 fold changes for each gene
- Colored ribbons indicate leading-edge genes for each pathway
- Includes presynaptic/postsynaptic ribosomes and cytosol pathways from SynGO

See [Chord_Diagrams/README.md](Chord_Diagrams/README.md) for detailed documentation.

---

### Complex_V_analysis/
**Purpose:** Deep-dive analysis of ATP synthase (Complex V) pathway components

See [Complex_V_analysis/README.md](Complex_V_analysis/README.md) for details.

---

### Critical_period_trajectories/
**Generator Script:** `02_Analysis/2.4.viz_critical_period_trajectories_gsva.R`
**Purpose:** GSVA-based temporal trajectory analysis showing module-level changes during D35→D65 maturation

**Key Finding:** Mutants fail to support synaptic maturation during critical period. Includes sample-level GSVA points for data transparency.

See [Critical_period_trajectories/README.md](Critical_period_trajectories/README.md) for full documentation.

---

### General/
**Purpose:** General exploratory plots and preliminary visualizations

See [General/README.md](General/README.md) for contents.

---

### GSEA/
**Generator Scripts:**
- `02_Analysis/1.1.main_pipeline.R` - Main GSEA pipeline
- `02_Analysis/1.3.add_mitocarta.R` - MitoCarta GSEA
- `02_Analysis/3.9.viz_pooled_dotplots.R` - Cross-database pooled visualizations

**Purpose:** Comprehensive GSEA results across 12 databases and 9 contrasts

**Structure:**
- Per-contrast subdirectories (e.g., `G32A_vs_Ctrl_D35/`)
- Per-database subfolders (hallmark, kegg, reactome, gobp, gocc, gomf, wiki, canon, cgp, tf, SynGO, MitoCarta)
- `Cross_database_pooled/` - Combined visualizations across databases

See [GSEA/README.md](GSEA/README.md) for complete documentation.

---

### Mito_translation_cascade/
**Generator Script:** `02_Analysis/2.2.viz_mito_translation_cascade.R`
**Purpose:** Mechanistic cascade heatmaps showing energy crisis → ribosome biogenesis → translation failure → synaptic dysfunction

**Key Visualizations:**
- Module-level summary heatmap (5 functional modules)
- Mechanistic cascade heatmap (gene-level expression)
- Module gene lists with documented sources

See [Mito_translation_cascade/README.md](Mito_translation_cascade/README.md) for details.

---

### Pattern_Summary_Normalized/
**Generator Script:** `02_Analysis/3.4.pattern_summary_normalized.py`
**Purpose:** Pattern classification visualizations showing temporal trajectory patterns (Compensation, Progressive, Sign_reversal, etc.)

**Contents:**
- Normalized pattern comparison plots for G32A and R403C
- Dual-panel visualizations showing Early→Late transitions with TrajDev curvature
- Pattern distribution summaries

See [Pattern_Summary_Normalized/README.md](Pattern_Summary_Normalized/README.md) for details.

---

### Publication_Figures/
**Generator Scripts:**
- `02_Analysis/3.1.publication_figures.py` - Main publication figures
- `02_Analysis/3.3.ribosome_upset_plot.py` - Ribosome gene overlap UpSet plot

**Purpose:** Publication-ready integrated figures combining multiple analyses (heatmap versions)

**Key Figures:**
- Fig1_Ribosome_Paradox.pdf - Core ribosome paradox finding
- Fig2_MitoCarta_Trajectory_Patterns.pdf - Mitochondrial pathway trajectory heatmap
- Fig3_Pattern_Classification_Summary.pdf - Cross-database pattern distribution
- Fig4_Semantic_Pathway_Overview.pdf - Semantic category overview heatmap

See [Publication_Figures/README.md](Publication_Figures/README.md) for complete figure inventory.

---

### Publication_Figures_Dotplot/
**Generator Script:** `02_Analysis/3.2.publication_figures_dotplot.py`
**Purpose:** Dotplot versions of publication figures (dot size = significance, color = NES)

**Advantages:**
- Dot size explicitly encodes statistical significance
- Easier to identify significant pathways compared to heatmaps
- Standard format for pathway enrichment publications

See [Publication_Figures_Dotplot/README.md](Publication_Figures_Dotplot/README.md) for details.

---

### Ribosome_paradox/
**Generator Script:** `02_Analysis/2.1.viz_ribosome_paradox.R`
**Purpose:** Core ribosome paradox finding showing opposing trajectories of ribosome biogenesis vs synaptic translation

**Key Finding:**
- Cytoplasmic ribosome biogenesis UP (NES +2.25)
- Presynaptic ribosomes DOWN (NES -2.9)
- Postsynaptic ribosomes DOWN (NES -3.0)
- Represents failed compensation due to ATP/translation bottleneck at synapses

See [Ribosome_paradox/README.md](Ribosome_paradox/README.md) for detailed interpretation.

---

### Synaptic_ribosomes/
**Generator Script:** `02_Analysis/2.3.viz_synaptic_ribosomes.R`
**Purpose:** Pre/postsynaptic translation analysis showing compartment-specific expression patterns

**Key Finding:** Translation machinery failure is pan-synaptic (affects both pre and postsynaptic compartments)

See [Synaptic_ribosomes/README.md](Synaptic_ribosomes/README.md) for details.

---

### Trajectory_Flow/
**Generator Scripts:**
- `02_Analysis/3.5.viz_trajectory_flow.py` - Alluvial diagrams (Plotly)
- `02_Analysis/3.6.viz_alluvial_ggalluvial.R` - Classical ggalluvial diagrams
- `02_Analysis/3.7.viz_bump_chart.py` - Static bump charts
- `02_Analysis/3.8.viz_interactive_bump.py` - Interactive bump charts
- `02_Analysis/3.8.viz_interactive_bump_dashboard.py` - Main interactive dashboard

**Purpose:** Pathway trajectory dynamics visualizations (alluvial diagrams, bump charts, flow visualizations)

**Key Visualizations:**
- `interactive_bump_dashboard.html` - Comprehensive interactive explorer
- `bump_focused_curved_nes.pdf` - Focused modules with trajectory curves
- `alluvial/` - All alluvial/Sankey diagrams

**Key Finding:** 99.9% of early defects improve during maturation (60% via active compensation, 40% via passive buffering)

See [Trajectory_Flow/README.md](Trajectory_Flow/README.md) for complete documentation.

---

### Volcano/
**Purpose:** Volcano plots for differential expression results across all contrasts

See [Volcano/README.md](Volcano/README.md) for details.

---

## Analysis Scripts Reference

### R Scripts

| Script | Purpose | Output Directory |
|--------|---------|------------------|
| `1.1.main_pipeline.R` | Main DE + GSEA pipeline | Multiple (checkpoints, DE_results, GSEA/) |
| `1.3.add_mitocarta.R` | MitoCarta GSEA | GSEA/*/MitoCarta/ |
| `2.1.viz_ribosome_paradox.R` | Core ribosome paradox | Ribosome_paradox/ |
| `2.2.viz_mito_translation_cascade.R` | Mechanistic cascade heatmaps | Mito_translation_cascade/ |
| `2.3.viz_synaptic_ribosomes.R` | Synaptic compartment analysis | Synaptic_ribosomes/ |
| `2.4.viz_critical_period_trajectories_gsva.R` | GSVA temporal trajectories | Critical_period_trajectories/ |
| `2.5.viz_complex_v_analysis.R` | ATP synthase deep-dive | Complex_V_analysis/ |
| `2.6.viz_calcium_genes.R` | Calcium signaling genes | General/ or Calcium/ |
| `3.6.viz_alluvial_ggalluvial.R` | Classical alluvial diagrams | Trajectory_Flow/alluvial/ |
| `3.9.viz_pooled_dotplots.R` | Cross-database pooled dotplots | GSEA/Cross_database_pooled/ |

### Python Scripts

| Script | Purpose | Output Directory |
|--------|---------|------------------|
| `1.5.create_master_pathway_table.py` | Master GSEA table | master_gsea_table.csv |
| `3.1.publication_figures.py` | Main publication figures | Publication_Figures/ |
| `3.2.publication_figures_dotplot.py` | Dotplot versions | Publication_Figures_Dotplot/ |
| `3.3.ribosome_upset_plot.py` | Ribosome gene overlap UpSet | Publication_Figures/ |
| `3.4.pattern_summary_normalized.py` | Pattern visualizations | Pattern_Summary_Normalized/ |
| `3.5.viz_trajectory_flow.py` | Alluvial diagrams (Plotly) | Trajectory_Flow/alluvial/ |
| `3.7.viz_chord_diagrams.py` | GOChord-style diagrams | Chord_Diagrams/ |
| `3.7.viz_bump_chart.py` | Static bump charts | Trajectory_Flow/ |
| `3.8.viz_interactive_bump.py` | Interactive bump variants | Trajectory_Flow/bump/ |
| `3.8.viz_interactive_bump_dashboard.py` | Main interactive dashboard | Trajectory_Flow/ |

---

## Data Sources

### Enrichment Databases (12 total)

1. **Hallmark** - MSigDB curated hallmark gene sets
2. **KEGG** - Kyoto Encyclopedia of Genes and Genomes
3. **Reactome** - Pathway database
4. **GO BP** - Gene Ontology Biological Process
5. **GO CC** - Gene Ontology Cellular Component
6. **GO MF** - Gene Ontology Molecular Function
7. **WikiPathways** - Community-curated pathways
8. **Canonical Pathways** - MSigDB canonical pathways
9. **CGP** - MSigDB Chemical & Genetic Perturbations
10. **TF** - Transcription Factor Targets
11. **SynGO** - Synapse-specific Gene Ontology
12. **MitoCarta 3.0** - Mitochondrial pathways

### Experimental Contrasts (9 total)

**Baseline mutation effects (4):**
- G32A_vs_Ctrl_D35
- R403C_vs_Ctrl_D35
- G32A_vs_Ctrl_D65
- R403C_vs_Ctrl_D65

**Time-course maturation effects (3):**
- Time_Ctrl (Control D65 vs D35)
- Time_G32A (G32A D65 vs D35)
- Time_R403C (R403C D65 vs D35)

**Mutation-specific maturation (2):**
- Maturation_G32A_specific (interaction effect)
- Maturation_R403C_specific (interaction effect)

### Expression Data

- **Genes analyzed:** ~15,000 protein-coding genes
- **Samples:** 18 total (3 genotypes × 2 timepoints × 3 replicates)
- **Checkpoints:** `03_Results/02_Analysis/checkpoints/`
- **Master tables:**
  - `master_gsea_table.csv` - 109,989 pathway-contrast combinations
  - `master_gsva_focused_table.csv` - 42 rows (7 modules × 6 groups)
  - `gsva_pattern_summary.csv` - Pattern classifications

---

## Key Findings Summary

### Core Finding: The Translation Paradox

**Question:** How do DRP1 mutations cause epilepsy?

**Answer:** Mitochondrial dysfunction → Local ATP depletion → Synaptic translation failure → Failed neuronal maturation → Seizures

**Evidence:**
- Cytoplasmic ribosome biogenesis UP (NES +2.25, FDR < 10⁻¹²)
- Presynaptic ribosomes DOWN (NES -2.9, FDR < 10⁻¹²)
- Postsynaptic ribosomes DOWN (NES -3.0, FDR < 10⁻¹²)

### Pattern Classification

Pathways are classified into 8 temporal trajectory patterns:

| Pattern | Description | G32A | R403C |
|---------|-------------|------|-------|
| **Compensation** | Active trajectory opposing early defect | ~1,400 pathways | ~1,200 pathways |
| **Sign_reversal** | Sign flipped between Early and Late | ~100 pathways | ~80 pathways |
| **Progressive** | Active trajectory amplifying early defect | ~20 pathways | ~15 pathways |
| **Natural_improvement** | Passive improvement (TrajDev NS) | ~800 pathways | ~700 pathways |
| **Natural_worsening** | Passive worsening (TrajDev NS) | ~50 pathways | ~40 pathways |
| **Late_onset** | New defect emerging at D65 | ~60 pathways | ~45 pathways |
| **Transient** | Early defect resolved by D65 | ~40 pathways | ~30 pathways |
| **Complex** | Inconsistent/multiphasic patterns | ~500 pathways | ~400 pathways |

See `docs/PATTERN_CLASSIFICATION.md` for complete framework documentation.

---

## Trajectory Framework

All temporal analyses follow the Early → TrajDev → Late framework:

- **Early (D35):** Initial mutation effects at early neuronal stage
- **TrajDev (Trajectory Deviation):** Mutation-specific developmental trajectory changes during D35→D65 maturation
- **Late (D65):** Mature neuron state at D65

**Key Insight:** ~60% of pathways show active compensation (significant TrajDev opposing early defect), while ~40% show passive improvement through normal developmental processes.

---

## How to Regenerate Visualizations

### Full Pipeline

```bash
# 1. Run main analysis pipeline (generates GSEA results and checkpoints)
Rscript 02_Analysis/1.1.main_pipeline.R

# 2. Add MitoCarta GSEA
Rscript 02_Analysis/1.3.add_mitocarta.R

# 3. Create master tables (optional, for comprehensive queries)
python3 02_Analysis/1.5.create_master_pathway_table.py
Rscript 02_Analysis/1.7.create_master_gsva_table.R

# 4. Generate R visualizations
Rscript 02_Analysis/2.1.viz_ribosome_paradox.R
Rscript 02_Analysis/2.2.viz_mito_translation_cascade.R
Rscript 02_Analysis/2.3.viz_synaptic_ribosomes.R
Rscript 02_Analysis/2.4.viz_critical_period_trajectories_gsva.R
Rscript 02_Analysis/3.6.viz_alluvial_ggalluvial.R
Rscript 02_Analysis/3.9.viz_pooled_dotplots.R

# 5. Generate Python visualizations
python3 02_Analysis/3.1.publication_figures.py
python3 02_Analysis/3.2.publication_figures_dotplot.py
python3 02_Analysis/3.3.ribosome_upset_plot.py
python3 02_Analysis/3.4.pattern_summary_normalized.py
python3 02_Analysis/3.5.viz_trajectory_flow.py
python3 02_Analysis/3.7.viz_chord_diagrams.py
python3 02_Analysis/3.7.viz_bump_chart.py
python3 02_Analysis/3.8.viz_interactive_bump.py
python3 02_Analysis/3.8.viz_interactive_bump_dashboard.py
```

### Individual Subdirectories

Each subdirectory README contains specific regeneration instructions for that analysis.

---

## File Formats

**Vector graphics (publication):**
- PDF - Editable vector format for journal submission
- Suitable for infinite scaling without quality loss

**Raster graphics (presentations):**
- PNG - High-resolution raster format (300 DPI)
- Suitable for presentations and web viewing

**Interactive HTML:**
- HTML - Interactive visualizations with hover tooltips and filtering
- Requires web browser, best for exploratory analysis

**Data exports:**
- CSV - Raw data tables for downstream analysis

---

## Related Documentation

- `CLAUDE.md` - Project overview and pipeline documentation
- `docs/PATTERN_CLASSIFICATION.md` - Pattern classification framework
- `01_Scripts/Python/pattern_definitions.py` - Canonical pattern classification code
- `03_Results/02_Analysis/README.md` - Results directory overview
- Individual subdirectory READMEs for detailed documentation

---

**Last Updated:** 2025-12-05
**Status:** Active analysis - all subdirectories documented

For questions about specific visualizations, see the README in each subdirectory.
