# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [2.0.1] - 2025-12-04

### Summary

Data lineage cleanup and documentation update. Removed unused intermediate files to streamline the data pipeline.

---

### Removed

- **Unused intermediate files from Python_exports/**:
  - `gsea_significant.csv` (14MB) - superseded by master_gsea_table.csv
  - `gsea_trajectory.csv` (7MB) - never consumed
  - `gsea_trajectory_all.csv` (22MB) - never consumed
  - `contrast_mapping.csv` (490B) - documentation only

- **Unused function from data_loader.py**:
  - `load_gsea_trajectory_data()` - never called by any script

### Changed

- **1.4.export_gsea_for_python.R**: Now only produces 2 files (was 6):
  - `gsea_results_long.csv` - consumed by 1.5
  - `gsea_results_wide.csv` - consumed by data_loader.py

- **Documentation renamed**:
  - `docs/SCRIPT_DEPENDENCY_ANALYSIS.md` → `docs/DATA_LINEAGE.md`

### Documentation

- Updated `docs/DATA_LINEAGE.md` with clean data flow architecture
- Updated `03_Results/02_Analysis/Python_exports/README.md`
- Updated `01_Scripts/Python/README.md`
- Updated `CLAUDE.md` script numbering (4.x → 1.x references)

---

## [2.0.0] - 2025-11-25

### Summary

Major update for manuscript resubmission with comprehensive visualization suite, cross-database validation framework, publication-ready figures, and extensive documentation. This release pushes the analysis from being a preliminary exploration to a complete, reproducible analytical pipeline with biological interpretation.

**Key Achievement**: Identification and validation of the **Ribosome Paradox** - a translation crisis where neurons increase ribosome biogenesis (NES +2.25, FDR < 1e-11) while synaptic translation programs collapse (NES -2.9 to -3.0, FDR < 1e-12), validated across 10 independent pathway databases.

---

### 🎯 Major Features Added

#### 1. Core Scientific Discovery: The Ribosome Paradox

**New Analysis Framework**: Three-pool ribosome model revealing paradoxical compensation

- **Pool 1 (Cytoplasmic Biogenesis)**: Compensatory upregulation during maturation (NES +2.25)
- **Pool 2 (Synaptic Ribosomes)**: Progressive collapse at synapses (NES -2.9 to -3.0)
- **Pool 3 (Mitochondrial Ribosomes)**: Compensatory upregulation (NES +1.9)

**Statistical Support**: FDR < 1e-11 for key effects, replicated across 10 databases

**Biological Interpretation**: Failed compensation due to unresolved ATP delivery failure at synapses

**Files**:
- `02_Analysis/2.1.viz_ribosome_paradox.R` - Three-pool trajectory visualization
- `03_Results/02_Analysis/Plots/Ribosome_paradox/` - Core finding plots and data

#### 2. Cross-Database Validation Framework

**Purpose**: Validate findings across 10 independent pathway databases to distinguish robust biological signals from database-specific artifacts

**Implementation**: Pattern classification system categorizing pathways into 7 temporal patterns:
- **Compensation** (dominant): Early disruption → recovery by D65
- **Progressive**: Gradual worsening D35 → D65
- **Transient**: D35 only, resolved by D65
- **Late_onset**: Normal D35, disrupted D65
- **Natural_improvement**: Enhanced at both timepoints
- **Natural_worsening**: Disrupted at both timepoints
- **Complex**: Non-linear/inconsistent patterns

**Databases Validated**:
1. Gene Ontology (BP, CC, MF) - 3000+ pathways
2. MSigDB Hallmark - 50 pathways
3. KEGG - 186 pathways
4. Reactome - 674 pathways
5. WikiPathways - 454 pathways
6. MSigDB Canonical - 1329 pathways
7. MSigDB TF Targets - 610 pathways
8. Chemical & Genetic Perturbations - 3396 pathways
9. **SynGO** (synaptic ontology) - 127 pathways
10. **MitoCarta** (mitochondrial) - 142 pathways

**Key Finding**: Compensation is dominant pattern (~60% of significant pathways), with G32A showing stronger compensation than R403C

**Files**:
- `02_Analysis/3.4.pattern_summary_normalized.py` - Pattern classification and trajectory plots
- `03_Results/02_Analysis/Plots/Cross_database_validation/` - 10 trajectory plots + pattern summary
- `03_Results/02_Analysis/Plots/Cross_database_validation/README.md` - Comprehensive documentation

#### 3. Publication Figure Suite (Python Pipeline)

**New Python Scripts** for manuscript-ready figures:

- `3.1.publication_figures.py` - Main publication figures (Fig 1-4)
  - Fig 1: Ribosome paradox three-pool visualization
  - Fig 2: MitoCarta trajectory patterns
  - Fig 3: Pattern classification summary across databases
  - Fig 4: Semantic pathway overview heatmap

- `3.2.publication_figures_dotplot.py` - Alternative dotplot visualizations
  - Same data, dotplot presentation style
  - Optimized for colorblind-safe palettes

- `3.3.ribosome_upset_plot.py` - Ribosome gene overlap UpSet plots
  - Pre/postsynaptic/mitochondrial ribosome overlap
  - Gene set intersection analysis

- `3.4.pattern_summary_normalized.py` - Pattern classification suite
  - 10 database-specific trajectory plots
  - Pattern distribution summary
  - Cross-database validation figures

**Output**: `03_Results/02_Analysis/Plots/Publication_Figures/` (publication-ready PDFs and PNGs)

#### 4. Critical Period Trajectory Analysis (GSVA-Based)

**New Script**: `2.4.viz_critical_period_trajectories_gsva.R`

**Enhancements over previous version**:
- GSVA (Gene Set Variation Analysis) for robust pathway scores
- Individual sample-level visualization (n=3 per condition)
- Mitochondria/ATP module added (ATP synthase, OXPHOS, dynamics)
- Mechanistic arrows annotating energy → translation → synapse cascade
- Divergence plots showing mutation-specific maturation failures

**Key Modules**:
1. Ribosome Biogenesis - Control: declines normally; Mutants: excessive increase
2. Cytoplasmic Translation - Control: stable; Mutants: progressive decline
3. Synaptic Ribosomes - Control: increase (synaptogenesis); Mutants: failed increase
4. Mitochondrial/ATP - Energy crisis module

**Output**: `03_Results/02_Analysis/Plots/Critical_period_trajectories/gsva/`

**Documentation**: Comprehensive 250+ line README explaining:
- Two different y-axis scales (absolute vs relative)
- How to read trajectory plots
- Sample variability interpretation
- Module-specific biological interpretation

#### 5. Mechanistic Cascade Visualization

**New Script**: `2.2.viz_mito_translation_cascade.R`

**Purpose**: Shows complete pathway from mitochondria → translation → synapse → calcium dysregulation

**Five Modules** (366 genes total):
1. **Energy Crisis** (35 genes): ATP synthase, OXPHOS, mitochondrial dynamics/transport
2. **Ribosome Biogenesis** (158 genes): UPREGULATED - compensatory response
3. **Cytoplasmic Translation** (76 genes): DOWNREGULATED - functional failure
4. **Synaptic Function** (65 genes): DOWNREGULATED - synaptic crisis
5. **Calcium Dysregulation** (32 genes): Mixed patterns

**Visualization Features**:
- Gene-level heatmap with hierarchical clustering
- Module-level flow diagram
- Increased font sizes for readability
- Color-coded modules

**Output**: `03_Results/02_Analysis/Plots/Mito_translation_cascade/`

#### 6. Synaptic Ribosome Compartment Analysis

**New Script**: `2.3.viz_synaptic_ribosomes.R`

**Purpose**: Pre/postsynaptic ribosome-specific analysis using SynGO database

**Key Analyses**:
- **52 presynaptic ribosome genes** (SynGO)
- **70 postsynaptic ribosome genes** (SynGO)
- **18 postsynaptic-specific genes** (not in presynaptic)
- Expression heatmaps by compartment
- Lollipop plots showing enrichment by compartment

**Key Finding**: Both compartments show severe downregulation during maturation (NES ~ -3.0), suggesting global synaptic translation failure rather than compartment-specific effects

**Output**: `03_Results/02_Analysis/Plots/Synaptic_ribosomes/`

#### 7. Additional Visualization Scripts

- `viz_developmental_framework.R` - Developmental pattern analysis
- `2.5.viz_complex_v_analysis.R` - ATP synthase (Complex V) deep-dive
- `2.6.viz_calcium_genes.R` - Calcium signaling gene analysis (13 validated genes)
- `viz_pooled_dotplots.R` - Cross-database GSEA dotplots

---

### 🔧 Major Pipeline Improvements

#### Volcano Plot Improvements (Session 2)

**Problem**: Y-axis ambiguity (p-value vs FDR mode) and calcium gene visibility issues

**Solution**:
- Added explicit subtitles distinguishing p-value vs FDR modes
- Created dedicated ggrepel layer for calcium genes (`max.overlaps = Inf`, bold black, `force = 5`)
- Corrected calcium gene list to 13 validated symbols (NNAT, CACNG3, CACNA1S, ATP2A1, RYR1, MYLK3, VDR, STIM1, STIM2, ORAI1_1, CALB1, CALR, PNPO)
- Standardized FDR-mode plots to use `label_method = "top"` to avoid empty plots
- Created diagnostic script `5.3.diagnose_calcium_genes.R` for validation

**Files Modified**:
- `01_Scripts/RNAseq-toolkit/scripts/DE/plot_standard_volcano.R`
- `02_Analysis/1.1.main_pipeline.R` (formerly Analysis_pipeline_fin.R)
- `02_Analysis/2.6.viz_calcium_genes.R`

**Testing**: Integration tests added (`tests/test_subtitle_volcano.R`, `tests/test_integrated_volcano.R`)

#### Pipeline Refactoring and Checkpointing

**Major Refactor**: Applied DRY (Don't Repeat Yourself) principles throughout main pipeline

**Performance Improvements**:
- **Eliminated 8+ duplicate `topTable()` calls** - Contrasts now processed once, reused everywhere
- **Consolidated volcano plot generation** - 60+ lines → 20 lines (single parameterized function)
- **Refactored CSV export** - 14 lines → 2 lines using utility functions
- **Improved GSEA/SynGO loops** - Now use pre-computed contrast tables

**New Utility Functions** (`01_Scripts/RNAseq-toolkit/scripts/utils_plotting.R`):
- `ensure_dir()` - Safe directory creation
- `save_plot()` - PDF save wrapper with error handling
- `load_checkpoint()` - Checkpoint loading utility
- `process_all_contrasts()` - Single-pass contrast processing
- `save_contrast_csvs()` - Batch CSV export
- `get_plot_params()` - Database-specific plot dimensions
- `close_all_devices()` - Graphics device cleanup
- `log_message()` - Timestamped logging
- `format_contrast_name()` - Display formatting

**Code Reduction**: ~150 lines of duplication eliminated

**Progress Tracking**: Added clear progress messages with emojis for better pipeline tracking

#### Checkpoint System Enhancement

**New Checkpoints** (10 total in `03_Results/02_Analysis/checkpoints/`):
1. `fit_object.rds` - Limma fit object (DEG modeling)
2. `de_results.rds` - decideTests results
3. `voom_object.rds` - Voom-transformed expression
4. `DGE_object.rds` - DGEList with metadata
5. `contrasts_matrix.rds` - Contrast definitions
6. `all_gsea_results.rds` - **Complete MSigDB GSEA results** (computationally expensive, 50-100 MB)
7. `syngo_gsea_results.rds` - SynGO GSEA results (2.4 MB)
8. `syngo_lists.rds` - SynGO term2gene/term2name (13 KB)
9. `qc_variables.rds` - QC metrics (samples, annotation, colors, logCPM)
10. `gene_intersections.rds` - Gene overlap data (baseline9, mat38)
11. `mitocarta_gsea_results.rds` - MitoCarta GSEA results
12. `gsva_module_scores.rds` - GSVA scores for trajectory analysis

**Benefits**:
- First run: 30-60 minutes
- Subsequent runs: 5-10 minutes (uses cached results)
- Set `config$force_recompute = TRUE` to bypass cache

**Critical Fix**: Missing SynGO checkpoints
- Root cause: SynGO analysis section hadn't been executed in pipeline
- Missing helper `format_pathway_name` added to sourcing list
- Generated syngo_gsea_results.rds (2.4 MB) and syngo_lists.rds (13 KB)
- All 9 contrasts processed successfully (32 SynGO pathways per contrast)
- Fixed viz_pooled_dotplots.R and viz_ribosome_pathways.R dependency issues

#### Pipeline Structure Reorganization

**Script Renaming**:
- `Analysis_pipeline.R` → `1.1.main_pipeline.R` (clearer ordering)
- Deprecated scripts moved to `02_Analysis/.deprecated/`

**New Pipeline Structure**:
```
1.1.main_pipeline.R                 # Core DE + GSEA (checkpointed)
1.2.generate_contrast_tables.R      # Generate contrast CSVs
1.3.add_mitocarta.R                  # Annotate with MitoCarta pathways
1.4.export_gsea_for_python.R         # Export for Python pipeline
2.x.viz_*.R                          # R visualizations (6 scripts)
3.x.*.py                             # Python publication figures (4 scripts)
4.x.*.py/R                           # Master table creation (4 scripts)
5.x.*.R                              # Sensitivity and verification (5 scripts)
```

**Total Runtime**: ~1-2 hours (first run), ~20-30 minutes (subsequent runs with caching)

---

### 📚 Documentation Overhaul

#### New Documentation Files

1. **`02_Analysis/SCRIPTS.md`** - Complete script inventory
   - Active vs deprecated scripts
   - Dependencies and run order
   - Runtime estimates
   - Checkpoint usage
   - Quick reference commands

2. **`docs/SESSION_HISTORY.md`** - Development log
   - Session 2: Volcano plot improvements
   - Session 3: Biological context & mechanistic model
   - Session 4: Translation crisis verification
   - Session 5: Visualization improvements
   - Milestone tracking

3. **`docs/NEXT_STEPS.md`** - Priorities and future directions
   - Immediate priorities (now completed in v2.0)
   - Medium-term goals
   - Future directions (post-manuscript)
   - Success criteria

4. **`docs/bio_notes.md`** - Consolidated biological context (replaces `docs/biological_context.md`)
   - Clinical context (DRP1 mutations → epilepsy)
   - Synaptic ribosome biology
   - Mitochondria-translation coupling
   - Critical period window (D35-D65)
   - Mechanistic model

5. **`03_Results/02_Analysis/Plots/README.md`** - Visualization overview
   - Figure organization
   - Core findings summary
   - Plot folder descriptions
   - Reading guide

#### Plot Folder READMEs (7 comprehensive documents)

Each major plot folder now contains extensive documentation:

1. **`Cross_database_validation/README.md`** (200+ lines)
   - Pattern classification framework
   - 10 database descriptions
   - Plot structure and interpretation
   - Key findings with statistics
   - Data file descriptions

2. **`Ribosome_paradox/README.md`** (250+ lines)
   - Three-pool ribosome model
   - Core finding explanation
   - Plot-by-plot breakdown
   - Statistical annotations guide
   - Biological interpretation

3. **`Publication_Figures/README.md`** (150+ lines)
   - Figure-by-figure descriptions
   - Manuscript figure mapping
   - Generation workflow
   - File format descriptions

4. **`Mito_translation_cascade/README.md`** (180+ lines)
   - Five-module cascade explanation
   - Module gene lists (366 genes total)
   - Visualization methods
   - Interpretation guide

5. **`Synaptic_ribosomes/README.md`** (200+ lines)
   - Compartment-specific analysis
   - Gene counts (52 pre, 70 post)
   - SynGO integration details
   - Key findings

6. **`Critical_period_trajectories/README.md`** (300+ lines)
   - Two y-axis scale explanations (absolute vs relative)
   - Visual elements guide (points, lines, trajectories)
   - Reading strategy step-by-step
   - GSVA methodology
   - Sample variability interpretation

7. **`Pattern_Summary_Normalized/README.md`** (150+ lines)
   - Pattern classification details
   - Normalization methods
   - Database-specific insights

**README Template Standardization**: All plot READMEs follow consistent structure:
- Purpose
- Generating script(s)
- Runtime and dependencies
- Plots generated (with detailed descriptions)
- Data files
- Methods (statistical tests, thresholds)
- Interpretation guide
- Key findings with statistics
- Related analyses

#### Enhanced Core Documentation

**Updated `CLAUDE.md`**:
- Architecture section expanded
- RNAseq-toolkit submodule documentation
- Checkpoint system guide
- SynGO integration patterns
- Common tasks with code examples
- Troubleshooting section

**Updated `README.md`**:
- Quick start guide
- Key findings summary with statistics
- Repository structure diagram
- Script organization
- Documentation roadmap
- What to explore first (for different audiences)

---

### 🧪 Analysis Enhancements

#### Contrast Framework Refinement

**9 Statistical Contrasts** organized into three categories:

**Mutation Effects** (4 contrasts):
- `G32A_vs_Ctrl_D35` - G32A mutation effect at early timepoint
- `R403C_vs_Ctrl_D35` - R403C mutation effect at early timepoint
- `G32A_vs_Ctrl_D65` - G32A mutation effect at late timepoint
- `R403C_vs_Ctrl_D65` - R403C mutation effect at late timepoint

**Maturation Effects** (3 contrasts):
- `Time_Ctrl` - Normal maturation trajectory (Control D65 vs D35)
- `Time_G32A` - G32A-specific maturation (G32A D65 vs D35)
- `Time_R403C` - R403C-specific maturation (R403C D65 vs D35)

**Interaction Effects** (2 contrasts):
- `Maturation_G32A_specific` - G32A-specific maturation changes
- `Maturation_R403C_specific` - R403C-specific maturation changes

**Trajectory Framework Mapping**:
- **Early**: D35 mutation effects
- **TrajDev**: Mutation-specific maturation changes (the critical period crisis)
- **Late**: D65 mutation effects

#### Gene Set Enrichment Analysis (GSEA) Expansion

**11 Pathway Databases** (up from original 6):
1. Hallmark (50 pathways)
2. KEGG (186 pathways)
3. Reactome (674 pathways)
4. GO:BP (7658 pathways)
5. GO:CC (996 pathways)
6. GO:MF (1738 pathways)
7. WikiPathways (454 pathways)
8. MSigDB Canonical (1329 pathways)
9. Chemical & Genetic Perturbations (3396 pathways)
10. Transcription Factor Targets (610 pathways)
11. **SynGO** (synaptic ontology) - 127 pathways
12. **MitoCarta 3.0** (mitochondrial) - 142 pathways

**Total**: ~17,000 pathway × 9 contrast = ~153,000 GSEA tests

**Custom Database Integration**:
- **SynGO** (Synapse Gene Ontology): Pre/postsynaptic compartment-specific enrichment
- **MitoCarta 3.0**: Mitochondrial pathway-specific enrichment
- Custom helper functions: `run_syngo_gsea.R`, `run_mitocarta_gsea.R`, `parse_mitocarta_gmx.R`

#### Calcium Gene Analysis Refinement

**Corrected Gene List** (13 validated genes):
- Removed: CACNA1C, CASR (not in dataset or not reliable)
- Added: ORAI1_1 (correct symbol in dataset)
- **Final list**: NNAT, CACNG3, CACNA1S, ATP2A1, RYR1, MYLK3, VDR, STIM1, STIM2, ORAI1_1, CALB1, CALR, PNPO

**Validation**: `5.3.diagnose_calcium_genes.R` confirms 13/13 genes present in dataset

**Visualization**: Dedicated calcium gene highlighting in volcano plots with dedicated ggrepel layer

#### Statistical Thresholds Standardization

**Differential Expression**:
- Significance: FDR < 0.05
- Fold-change: |log2FC| > 1 (i.e., fold-change > 2)
- Method: limma-voom with empirical Bayes moderation

**GSEA**:
- Significance: FDR < 0.05
- Permutations: 10,000
- Minimum gene set size: 15 genes
- Maximum gene set size: 500 genes

**Visualization**:
- Volcano plots: p-value mode (p < 0.05) or FDR mode (FDR < 0.1)
- GSEA plots: FDR < 0.05 threshold
- Pattern classification: FDR < 0.05 for "significant" designation

---

### 🐛 Bug Fixes

1. **Calcium gene visibility in volcano plots**
   - **Issue**: Calcium genes disappeared due to ggrepel `max.overlaps` default
   - **Fix**: Dedicated geom_text_repel layer with `max.overlaps = Inf`

2. **Empty FDR-mode volcano plots**
   - **Issue**: Maturation contrasts had no significant genes at FDR < 0.05
   - **Fix**: Use `label_method = "top"` to show top-ranked genes regardless of significance

3. **ORAI1 gene symbol mismatch**
   - **Issue**: Dataset contains `ORAI1_1` not `ORAI1`
   - **Fix**: Updated calcium gene list to use correct symbol

4. **Missing SynGO checkpoints**
   - **Issue**: SynGO analysis section wasn't executing, causing downstream script failures
   - **Fix**: Added missing helper function sourcing, generated checkpoints
   - **Result**: All 9 contrasts × 32 SynGO pathways processed successfully

5. **Path portability issues**
   - **Issue**: Hardcoded `/workspaces/GVDRP1/` paths broke on migration
   - **Fix**: Converted all paths to `here::here()` relative paths

6. **Duplicate computation inefficiency**
   - **Issue**: 8+ duplicate `topTable()` calls, redundant volcano generation
   - **Fix**: Process contrasts once, reuse results
   - **Impact**: Reduced runtime and improved code maintainability

---

### 🔬 Utility Scripts Added

1. **`Supp1.verify_enrichments.R`** - Systematic validation and documentation
   - Calcium pathway enrichment verification across all databases
   - Gene list intersection verification (baseline9, mat38)
   - Ribosome pathway statistics for manuscript
   - Calcium gene differential expression summary
   - Markdown summary report with interpretations

2. **`Supp2.diagnose_calcium_genes.R`** - Calcium gene presence validation
   - Checks all 13 calcium genes in dataset
   - Reports missing genes, synonyms, alternatives
   - Reproducible diagnostic for future datasets

3. **`1.8.extract_syngo_ribosome_genes.R`** - SynGO gene list extraction
   - Extracts pre/postsynaptic ribosome genes from SynGO
   - Gene counts and overlap analysis
   - Used for manual gene set curation

4. **`run_syngo_only.R`** - Standalone SynGO GSEA runner
   - Loads checkpoints from main pipeline
   - Runs SynGO analysis independently
   - Useful for checkpoint recovery or re-running just SynGO

5. **`Supp3.analyze_de_thresholds.R`** - Threshold sensitivity analysis
   - Explores effects of different FDR and FC thresholds
   - Helps justify threshold choices for manuscript

---

### 🗂️ Data Export for Python Pipeline

**New Script**: `1.4.export_gsea_for_python.R`

**Purpose**: Export R GSEA results to CSV format for Python visualization pipeline

**Exports**:
1. `gsea_trajectory.csv` - Complete GSEA results with trajectory metadata
2. `gsea_by_contrast.csv` - Contrast-specific GSEA results
3. `pathway_annotations.csv` - Pathway metadata and classifications
4. `gene_set_leadingEdge.csv` - Leading edge gene lists

**Output Directory**: `03_Results/02_Analysis/Python_exports/`

**Benefits**:
- Enables Python-based visualization (matplotlib, seaborn)
- Cross-language reproducibility
- Easier integration with external tools

---

### 📦 Repository Organization

#### Deprecated Scripts Archive

**Moved to** `02_Analysis/.deprecated/`:
- `Analysis_pipeline.R` (old version)
- `viz_critical_period_trajectories.R` (old method)
- `viz_critical_period_trajectories_with_mitocarta.R` (merged into new version)
- `viz_temporal_trajectory.R` (redundant)
- `regenerate_volcanoes.R` (integrated into pipeline)
- `4.visualize_trajectory_patterns.py` (superseded by v6)
- `5.focused_pathway_analysis.py` (superseded by v6)
- `experimental.R` (unused code)

**Rationale**: Keep for reference but prevent accidental use

#### Git Submodule Update

**RNAseq-toolkit** submodule updated to include:
- New utility functions (`utils_plotting.R`)
- Enhanced GSEA plotting functions
- Improved volcano plot functions with subtitle support
- Custom theme improvements

#### Container Configuration

**Docker Compose** setup with:
- Image: `scdock-r-dev:v0.5.1`
- Environment variable support (`.env` file)
- Volume mounts for data persistence
- UID/GID mapping for file permissions

---

### 🧪 Testing

#### Integration Tests Added

**Test Suite**:
1. `tests/test_subtitle_volcano.R` - Volcano plot subtitle functionality
2. `tests/test_integrated_volcano.R` - Calcium gene visibility and label methods
3. `tests/test_checkpoint_loading.R` - Checkpoint system validation

**Coverage**:
- Volcano plot generation (both p-value and FDR modes)
- Calcium gene highlighting
- Label method behavior (top vs significant)
- Checkpoint loading and validation

#### Diagnostic Scripts

**Validation Tools**:
- `Supp2.diagnose_calcium_genes.R` - Gene presence validation
- `Supp1.verify_enrichments.R` - GSEA enrichment verification
- `Supp3.analyze_de_thresholds.R` - Threshold sensitivity

**Output**: Console reports and CSV summaries for manuscript

---

### 📊 Key Findings Documented

#### The Ribosome Paradox

**Quantified**:
- Pool 1 (Cytoplasmic Biogenesis): Early NES -2.6 → TrajDev NES +2.25 → Late NES ~0
- Pool 2 (Synaptic Ribosomes): Early NES +2.3 → TrajDev NES -3.0 → Late NES -2.5
- Pool 3 (Mitochondrial Ribosomes): Early NES -2.2 → TrajDev NES +1.9 → Late NES +0.5

**Validated**: FDR < 1e-11 for biogenesis, FDR < 1e-12 for synaptic collapse

#### Mitochondrial Compensation

**4 Shared Pathways** compensating in both mutations:
1. Mitochondrial ribosome
2. Mitochondrial ribosome assembly
3. mtDNA maintenance
4. OXPHOS

**Mutation Specificity**:
- G32A: 45 MitoCarta compensation pathways
- R403C: 29 MitoCarta compensation pathways
- Pattern: G32A shows stronger synaptic pathway compensation
- R403C shows broader compensation profile (100 vs 74 pathways total)

#### Pattern Distribution Across Databases

**Compensation Dominance**:
- ~60% of significant pathways show compensation pattern
- Validated across all 10 databases
- MitoCarta and GO:BP show strongest compensation signals

**Database Consistency**:
- SynGO: 19-22 compensation pathways per mutation
- MitoCarta: 45 (G32A) vs 29 (R403C) compensation pathways
- GO:BP: 500+ compensation pathways (filtered to top 50 for visualization)

---

### 🔄 Breaking Changes

1. **Script Naming**:
   - `Analysis_pipeline.R` → `1.1.main_pipeline.R`
   - Old scripts moved to `.deprecated/`
   - Update any external scripts referencing old names

2. **Checkpoint Structure**:
   - Added new checkpoints: `gsva_module_scores.rds`, `mitocarta_gsea_results.rds`
   - Checkpoint directory: `03_Results/02_Analysis/checkpoints/`
   - May need to regenerate checkpoints if upgrading from v1.x

3. **Calcium Gene List**:
   - Changed from 15 genes to 13 genes (removed CACNA1C, CASR; corrected ORAI1 to ORAI1_1)
   - Scripts using old list need updating

4. **Python Dependencies**:
   - New Python scripts require: matplotlib, seaborn, pandas, numpy, scipy, upsetplot
   - Install: `pip install matplotlib seaborn pandas numpy scipy upsetplot`

---

### 🚀 Migration from v1.0

**If upgrading from v1.0** (manuscript submission version):

1. **Pull latest code**:
   ```bash
   git checkout main
   git pull
   git submodule update --init --recursive
   ```

2. **Delete old checkpoints** (optional, will regenerate):
   ```bash
   rm 03_Results/02_Analysis/checkpoints/*.rds
   ```

3. **Run updated pipeline**:
   ```bash
   Rscript 02_Analysis/1.1.main_pipeline.R
   Rscript 02_Analysis/1.2.generate_contrast_tables.R
   Rscript 02_Analysis/1.3.add_mitocarta.R
   Rscript 02_Analysis/1.4.export_gsea_for_python.R
   ```

4. **Generate new visualizations**:
   ```bash
   # R visualizations
   Rscript 02_Analysis/2.1.viz_ribosome_paradox.R
   Rscript 02_Analysis/2.4.viz_critical_period_trajectories_gsva.R
   # ... (run others as needed)

   # Python visualizations
   python3 02_Analysis/3.1.publication_figures.py
   python3 02_Analysis/3.4.pattern_summary_normalized.py
   ```

**Total runtime**: ~1-2 hours (first run with cache generation)

---

### 📖 Documentation Changes

#### New Documentation

- `02_Analysis/SCRIPTS.md` - Complete script inventory
- `docs/PATTERN_CLASSIFICATION.md` - Pattern classification system documentation
- `docs/bio_notes.md` - Curated biological notes (consolidated from prior literature files)
- `docs/biological_context.md` - Additional biological context (now consolidated into `docs/bio_notes.md`)
- 7 comprehensive plot folder READMEs
- `03_Results/02_Analysis/Plots/README.md` - Visualization overview

#### Updated Documentation

- `README.md` - Enhanced with key findings, script organization, exploration guide
- `CLAUDE.md` - Expanded architecture, checkpoint system, SynGO integration, troubleshooting
- `MIGRATION.md` - Container setup and migration notes
- `CHANGELOG.md` - This comprehensive v2.0 changelog

#### Documentation Statistics

- **New markdown files**: 15
- **Total documentation**: ~5,000 lines
- **Plot folder READMEs**: 7 folders × ~200 lines average = ~1,400 lines
- **Code comments**: Enhanced throughout all scripts

---

### 🎓 Reproducibility Improvements

1. **Checkpoint System**: 30-60 min → 5-10 min re-runs
2. **Docker Container**: Standardized environment across workstations
3. **Comprehensive READMEs**: Every major analysis documented
4. **Script Inventory**: `SCRIPTS.md` documents all active scripts
5. **Path Portability**: `here::here()` throughout
6. **Git Tags**: v2.0 tagged at this commit for manuscript resubmission
7. **Session Info**: Saved with all major checkpoint files

---

### 🔮 Future Directions

**Medium-Term Goals** (1-2 weeks):
- Protein-level validation (if proteomics available)
- Single-cell RNA-seq for cell-type-specific effects
- Network analysis (co-expression modules)
- Comparative analyses with other DRP1 mutations

**Long-Term Goals** (Post-manuscript):
- Functional assays (ATP measurements, translation rates)
- Drug screening (rescue compounds)
- Patient-derived tissue comparison
- Related mitochondrial disorder comparisons

---

### 👥 Contributors

**Development Sessions**:
- Session 2 (2025-11-20): Volcano plot improvements
- Session 3 (2025-11-20): Biological context research
- Session 4 (2025-11-20): Translation crisis verification
- Session 5 (2025-11-20): Visualization improvements
- Session 6-8 (2025-11-21 to 2025-11-25): Publication figures, cross-database validation, documentation

**Code Development**: Claude Code (Anthropic) + Human collaboration

**Scientific Context**: Dr. Vivian Gama (PI), literature review, biological interpretation

---

### 📝 Notes for Reviewers

This version represents a major analytical and visualization overhaul for manuscript resubmission. Key improvements for reviewers:

1. **Robustness**: Findings validated across 10 independent databases
2. **Transparency**: Comprehensive documentation of methods and interpretation
3. **Reproducibility**: Complete checkpoint system, containerized environment
4. **Clarity**: Publication-ready figures with detailed captions
5. **Depth**: Multiple analysis layers (gene, pathway, pattern, temporal)

**For Questions**: See individual plot folder READMEs for detailed methods and interpretation guidance.

---

## [1.0.0] - 2025-06-18

### Summary

Frozen version corresponding to initial manuscript submission. Tagged at commit `d6ec164`.

### Included Features

- Complete bulk RNA-seq analysis pipeline for DRP1 mutations (G32A, R403C)
- Differential expression analysis using edgeR/limma voom
- Gene set enrichment analysis (GSEA) with multiple databases
- SynGO synaptic ontology enrichment
- Calcium signaling gene analysis
- Comprehensive visualization suite:
  - Volcano plots (standard and vertical comparison)
  - MA plots
  - Heatmaps (correlation, gene expression)
  - MDS plots
  - GSEA enrichment plots

### Analysis Details

- **Samples**: 18 (3 genotypes × 2 timepoints × 3 biological replicates)
- **Contrasts**: 9 statistical comparisons
  - 4 mutation effects (G32A vs Ctrl, R403C vs Ctrl at D35 and D65)
  - 3 maturation effects (Time_Ctrl, Time_G32A, Time_R403C)
  - 2 interaction effects (Maturation_G32A_specific, Maturation_R403C_specific)
- **Normalization**: TMM (Trimmed Mean of M-values)
- **Statistical method**: limma-voom with empirical Bayes moderation
- **Thresholds**: FDR < 0.05, |log2FC| > 1

### Data Files

- Counts matrix: `03_Results/01_Preprocessing/04_FeatureCounts/count_matrices_fc/sorted_counts_matrix.txt`
- Metadata: `03_Results/01_Preprocessing/04_FeatureCounts/count_matrices_fc/metadata.csv`
- SynGO database: `00_Data/SynGO_bulk_20231201/`

### Known Limitations (Addressed in v2.0)

- Limited cross-database validation
- No systematic pattern classification
- Calcium gene visibility issues in volcano plots
- Limited documentation of biological interpretation
- No checkpoint system (slow re-runs)
- Hardcoded file paths (portability issues)

---

## Version Numbering

- **Major version (X.0.0)**: Substantial changes to analysis pipeline or methodology
- **Minor version (1.X.0)**: New analyses, figure updates, or reviewer-requested changes
- **Patch version (1.0.X)**: Bug fixes, documentation updates, minor corrections

---

## Commit Message Convention

When making changes for reviewer comments, use descriptive commit messages:

```
[Reviewer 2] Add mitochondrial gene enrichment analysis

- Perform GSEA on mitochondrial GO terms
- Create supplementary figure S5
- Update methods section accordingly

Addresses Reviewer #2, Comment #3
```

---

## Notes for Maintainers

### Tracking Reviewer Comments

For each reviewer comment being addressed:
1. Add entry to CHANGELOG.md under appropriate section
2. Reference the specific reviewer and comment number
3. Describe what was changed and why
4. Note any new files or figures generated
5. Commit with descriptive message linking to reviewer comment

### Reproducibility Notes

- Always document the exact version of code used for each submission
- Tag git commits corresponding to manuscript revisions (e.g., `v2.0-resubmission`)
- Include R session info in analysis outputs when reporting new results
- Document any changes to statistical thresholds or parameters

---

**For pattern classification details, see**: `docs/PATTERN_CLASSIFICATION.md`
**For comprehensive literature notes, see**: `docs/bio_notes.md`
**For additional biological context, see**: `docs/bio_notes.md`
