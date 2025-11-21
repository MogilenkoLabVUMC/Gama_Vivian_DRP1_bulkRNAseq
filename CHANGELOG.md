# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added - 2025-11-20 Volcano Plot Improvements
- **NEW:** `create_standard_volcano()` accepts a `subtitle` argument so the pipeline can explicitly describe p-value vs FDR modes and avoids ambiguity for downstream viewers.
- **NEW:** `diagnose_calcium_genes.R` provides a reproducible check that every curated calcium symbol exists in the dataset (13 corrected genes) before plotting.
- **NEW:** Integration tests `tests/test_subtitle_volcano.R` and `tests/test_integrated_volcano.R` cover the new subtitles, calcium gene highlighting, and label-method behavior across disease/maturation contrasts.

### Changed - 2025-11-20 Volcano Plot Improvements
- **Updated:** `01_Scripts/RNAseq-toolkit/scripts/DE/plot_standard_volcano.R` draws calcium genes via a dedicated `geom_text_repel` layer (`max.overlaps = Inf`, bold black, `force = 5`) and surfaces the provided subtitle.
- **Updated:** `02_Analysis/Analysis_pipeline.R` and `02_Analysis/viz_calcium_genes.R` use the corrected 13-gene calcium list (`ORAI1_1` replaces `ORAI1`, drop `CACNA1C`/`CASR`), forward mode-specific subtitles, and standardize the FDR plotting block to `label_method = "top"` for consistency.

### Fixed - 2025-11-20 Volcano Plot Improvements
- Calcium genes no longer disappear: every curated gene now survives ggrepel thanks to the dedicated layer and force settings.
- FDR-mode volcano plots for maturation contrasts (and any weak effect) now show top-ranked genes instead of empty canvases because `label_method = "top"`.
- The calcium gene list matches the actual dataset so the diagnostic script reports 13/13 present genes (including `ORAI1_1`).

### Testing - 2025-11-20
- Added focused integration tests to validate subtitles, calcium visibility, and label behavior instead of relying solely on visual inspection.
- Diagnostic tooling (`diagnose_calcium_genes.R`) confirms that calcium symbols and synonyms remain reproducible with future datasets.

### Documentation - 2025-11-20
- Captured the volcano improvements, testing notes, and next-session plan in `VOLCANO_IMPROVEMENTS_SUMMARY.md` plus the session/next-session handoff notes so future contributors can follow the remaining visualization priorities without re-running the investigative work.

### Added - 2025-11-20 Refactoring Session
- **NEW:** `utils_plotting.R` in RNAseq-toolkit - Reusable DRY helper functions for RNA-seq plotting
  - `ensure_dir()` - Safe directory creation
  - `save_plot()` - PDF save wrapper with error handling
  - `load_checkpoint()` - Checkpoint loading utility
  - `process_all_contrasts()` - Single-pass contrast processing (eliminates duplicate topTable calls)
  - `save_contrast_csvs()` - Batch CSV export
  - `get_plot_params()` - Database-specific plot dimensions
  - `close_all_devices()` - Graphics device cleanup
  - `log_message()` - Timestamped logging
  - `format_contrast_name()` - Display formatting
- **NEW:** `run_syngo_only.R` - Standalone SynGO GSEA analysis runner with comprehensive error handling
  - Loads checkpoints from main pipeline
  - Runs SynGO analysis independently
  - Generates missing syngo_gsea_results.rds and syngo_lists.rds
  - Useful for checkpoint recovery or re-running just SynGO analysis

### Changed - 2025-11-20 Refactoring Session
- **MAJOR REFACTOR:** `Analysis_pipeline.R` - Applied DRY principles throughout
  - **Eliminated 8+ duplicate `topTable()` calls** - Contrasts now processed once, reused everywhere
  - **Consolidated volcano plot generation** - Single parameterized function replaces 60+ lines of duplication
  - **Refactored CSV export** - 14 lines → 2 lines using utility functions
  - **Refactored volcano plots** - 60 lines → 20 lines (p/fdr modes)
  - **Refactored MD plots** - 35 lines → 21 lines
  - **Refactored FC-B plots** - 20 lines → 17 lines
  - **Updated GSEA loops** - Now use pre-computed contrast_tables
  - **Updated SynGO loops** - Now use pre-computed contrast_tables
  - Added clear progress messages with emojis for better tracking
  - Improved code organization with consistent patterns
- **IMPROVED:** Volcano plot y-axis titles now include explicit threshold information
  - P-value mode: "(Threshold: p ≤ 0.05)"
  - FDR mode: "(Threshold: FDR ≤ 0.1)"
- Sources utils_plotting.R from RNAseq-toolkit for reusability

### Fixed - 2025-11-20 Refactoring Session
- **CRITICAL FIX:** Missing SynGO checkpoints issue
  - Root cause: SynGO analysis section hadn't been executed
  - Missing helper function `format_pathway_name` added to sourcing list
  - Generated syngo_gsea_results.rds (2.4 MB) and syngo_lists.rds (13 KB)
  - All 9 contrasts processed successfully (32 SynGO pathways per contrast)
  - viz_pooled_dotplots.R and viz_ribosome_pathways.R now work correctly

### Performance Improvements - 2025-11-20
- **Reduced redundant computation:** Contrast results computed once instead of 8+ times
- **Improved code maintainability:** Eliminated ~150 lines of code duplication
- **Better error handling:** Added tryCatch blocks and clear error messages
- **Cleaner output:** Progress messages clearly indicate which step is running

---

### Added - Previous Sessions
- MIGRATION.md documentation describing workstation migration process
- CHANGELOG.md for tracking changes (this file)
- .env file support for container configuration
- Docker Compose configuration for dev container

### Changed - Previous Sessions
- **BREAKING:** Updated all hardcoded `/workspaces/GVDRP1/` paths to relative paths using `here::here()`
  - Affected file: `02_Analysis/Analysis_pipeline_fin.R` (7 path references)
- Updated RNAseq-toolkit git submodule path to `/data1/users/antonz/pipeline/RNAseq-toolkit/`
- Migrated from direct Docker image to Docker Compose-based dev container
- Updated docker-compose.yml to remove duplicate volume mounts

### Removed
- External mount dependencies for raw sequencing data (no longer needed for analysis)
  - Raw FASTQ files mount
  - Reference genome mount
  - Bash preprocessing scripts mount

### Fixed
- Corrected RNAseq-toolkit submodule initialization
- Standardized workspace paths for portability across workstations

### Documentation
- Updated README.md with migration notes and new setup instructions
- Added troubleshooting section to MIGRATION.md

### Summary
Updated the main analysis pipleine script with checkpoints creation for faster re-runs. 

## 📝 Files Created/Modified

### 1. Updated Main Pipeline
**File:** `Analysis_pipeline_fin.R`

**Changes Made:**
1. ✅ Added PNPO to calcium genes (line 29) – now 15 genes total
2. ✅ Created checkpoint directory (lines 35-37)
3. ✅ Added checkpoint saving after DEG modeling (lines 149-158)
4. ✅ Added checkpoint saving after GSEA (lines 477-482)
5. ✅ Added checkpoint saving after SynGO GSEA (lines 538-544)
6. ✅ Added checkpoint saving for QC variables (lines 226-233)
7. ✅ Fixed baseline9/mat38 loading with verification (lines 632-679)

**Checkpoints Saved (8 files):**
- `fit_object.rds` – limma fit object
- `de_results.rds` – decideTests results
- `voom_object.rds` – voom-transformed expression
- `DGE_object.rds` – DGEList with metadata
- `contrasts_matrix.rds` – contrast matrix
- `all_gsea_results.rds` – **ALL MSigDB GSEA results (EXPENSIVE!)**
- `syngo_gsea_results.rds` – SynGO GSEA results
- `syngo_lists.rds` – SynGO term2gene/term2name
- `qc_variables.rds` – ordered_samples, annot, colors, logCPM
- `gene_intersections.rds` – baseline9 & mat38 (loaded + computed)

---

### 2. New Visualization Scripts (3 files)

#### A. `viz_ribosome_pathways.R` (382 lines)
**Purpose:** THE KEY FINDING – Ribosome pathway dysregulation

**Generates:**
- Running sum plots for Presynaptic/Postsynaptic Ribosome pathways
- Ribosomal gene expression heatmaps (logCPM)
- Average expression heatmap by condition
- Pathway context dotplots (ribosome + mitochondria + synapse)
- Summary text report

**Output Directory:** `03_Results/02_Analysis/Plots/Ribosome_analysis/`

**Runtime:** ~15 minutes (uses checkpoints, no GSEA re-run)

---

#### B. `viz_pooled_dotplots.R` (330 lines)
**Purpose:** Cross-database integrated pathway view

**Generates TWO versions for each priority contrast:**
1. **Comprehensive** – All 11 databases (Hallmark, GO-BP/CC/MF, KEGG, Reactome, Wiki, Canonical, CGP, TF, SynGO)
2. **Focused** – 4 key databases (Hallmark, KEGG, Reactome, SynGO)

**Priority Contrasts:**
- Maturation_G32A_specific
- Maturation_R403C_specific
- Time_Ctrl
- G32A_vs_Ctrl_D65
- R403C_vs_Ctrl_D65

**Output Directories:**
- `03_Results/02_Analysis/Plots/GSEA/Cross_database_pooled/comprehensive/`
- `03_Results/02_Analysis/Plots/GSEA/Cross_database_pooled/focused/`

**Runtime:** ~20 minutes

---

#### C. `verify_enrichments.R` (470 lines)
**Purpose:** Systematic verification and documentation

**Task 1: Calcium Pathway Enrichment Verification**
- Searches ALL databases for calcium-related terms
- Outputs: `calcium_pathways_all.csv`, `calcium_pathways_significant.csv`
- **Answers:** Can we claim pathway-level calcium dysregulation?

**Task 2: Gene List Intersection Verification**
- Compares loaded vs computed baseline9 & mat38
- Outputs: Comparison tables, text files
- **Answers:** Are gene lists consistent with previous analysis?

**Task 3: Ribosome Pathway Statistics**
- Documents all ribosome pathway enrichments
- Outputs: `ribosome_pathway_statistics.csv`, `ribosome_pathway_significant.csv`
- **Provides:** Exact stats for manuscript (NES, FDR, gene counts)

**Task 4: Calcium Gene Analysis**
- Which calcium genes are DE across contrasts
- Outputs: `calcium_genes_DE_summary.csv`

**Task 5: Summary Report**
- **`verification_summary_report.md`** ⭐ **READ THIS FIRST**
- Markdown report with interpretations and recommendations

**Output Directory:** `03_Results/02_Analysis/Verification_reports/`

**Runtime:** ~10 minutes

## [1.0] - 2025-06-18

### Summary
Frozen version corresponding to manuscript submission/publication. Tagged at commit `d6ec164`.

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
- **Samples:** 26 (3 genotypes × 2 timepoints × biological replicates)
- **Contrasts:** 9 statistical comparisons
  - 4 pairwise (mutation vs control at each timepoint)
  - 3 maturation (time effects within genotype)
  - 2 interaction (mutation-specific maturation changes)
- **Normalization:** TMM (Trimmed Mean of M-values)
- **Statistical method:** limma-voom with empirical Bayes moderation

### Data Files
- Counts matrix: `03_Results/01_Preprocessing/04_FeatureCounts/count_matrices_fc/sorted_counts_matrix.txt`
- Metadata: `03_Results/01_Preprocessing/04_FeatureCounts/count_matrices_fc/metadata.csv`
- SynGO database: `00_Data/SynGO_bulk_20231201/`

---

## Future Sections (Template for Reviewer Comments)

Use the sections below to track changes made in response to reviewer feedback:

## [1.1] - YYYY-MM-DD (In Progress)

### Added - Reviewer Response
<!-- Example:
- New analysis of mitochondrial genes per Reviewer #2 comment
- Additional GSEA database (WikiPathways) per Reviewer #3 suggestion
-->

### Changed - Reviewer Response
<!-- Example:
- Updated FDR threshold from 0.05 to 0.1 per Reviewer #1
- Re-ran GSEA with updated gene sets
- Modified volcano plot color scheme for better visibility
-->

### Figures - Reviewer Response
<!-- Example:
- Regenerated Figure 3 with higher resolution
- Added supplementary figure showing cell type markers
- Updated Figure 5 to include NNAT (Neuronatin) as requested
-->

---

## Notes for Maintainers

### Version Numbering
- **Major version (X.0):** Substantial changes to analysis pipeline or methodology
- **Minor version (1.X):** New analyses, figure updates, or reviewer-requested changes
- **Patch version (1.0.X):** Bug fixes, documentation updates, minor corrections

### Commit Message Convention
When making changes for reviewer comments, use descriptive commit messages:
```
[Reviewer 2] Add mitochondrial gene enrichment analysis

- Perform GSEA on mitochondrial GO terms
- Create supplementary figure S5
- Update methods section accordingly

Addresses Reviewer #2, Comment #3
```

### Tracking Reviewer Comments
For each reviewer comment being addressed:
1. Add entry to CHANGELOG.md under appropriate section
2. Reference the specific reviewer and comment number
3. Describe what was changed and why
4. Note any new files or figures generated
5. Commit with descriptive message linking to reviewer comment

### Reproducibility Notes
- Always document the exact version of code used for each submission
- Tag git commits corresponding to manuscript revisions (e.g., `v1.1-revision1`)
- Include R session info in analysis outputs when reporting new results
- Document any changes to statistical thresholds or parameters
