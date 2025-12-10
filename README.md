# DRP1 Bulk RNA-seq Analysis

# Bioinformatics analyst and github repo maintainer 
[Anton Zhelonkin, MD](https://github.com/tony-zhelonkin)
# Project supervisor 
[Denis Mogilenko, PhD](https://github.com/MogilenkoLab)
# Principal Investigator
[Vivian Gama](https://medschool.vanderbilt.edu/cdb/person/vivian-gama-ph-d/)


# Intro

Transcriptional analysis of DRP1 mutations (G32A, R403C) in iPSC-derived cortical neurons across neuronal maturation (Day 35 → Day 65). This repository contains a complete downstream analysis pipeline starting from count matrices to hopefully insightfull visualisations.

## 📂 Repository Structure

```
├── 00_Data/                       # Reference databases (SynGO, MitoCarta)
├── 01_Scripts/
│   ├── RNAseq-toolkit/           # Git submodule with GSEA & DE helper functions
│   ├── R_scripts/                # Project-specific helpers (SynGO, MitoCarta integration)
│   └── Python/                   # Python modules (pattern_definitions, viz_bump_charts)
├── 02_Analysis/                   # Analysis pipeline scripts (see SCRIPTS.md)
│   ├── 1.1-1.7.*.R/py            # Main pipeline: DE, GSEA, GSVA, master tables
│   ├── 2.1-2.6.viz_*.R           # R visualization scripts (6 scripts)
│   ├── 3.1-3.9.*.py/R            # Python/R visualizations including bump charts
│   ├── Supp*.R/py                # Supplementary/utility scripts
│   ├── SCRIPTS.md                # Complete script inventory and dependencies
│   └── .deprecated/              # Archived old scripts
├── 03_Results/
│   ├── 01_Preprocessing/         # Count matrices and metadata (input data)
│   └── 02_Analysis/
│       ├── checkpoints/          # Cached computation results (.rds files)
│       ├── DE_results/           # Differential expression tables (9 contrasts)
│       ├── Python_exports/       # GSEA data exported for Python (CSV files)
│       ├── master_gsea_table.csv # Comprehensive GSEA results (109K rows)
│       ├── master_gsva_*.csv     # GSVA enrichment scores (focused + all pathways)
│       └── Plots/                # All visualizations (see folder READMEs)
│           ├── Trajectory_Flow/           # Bump charts & alluvial diagrams
│           ├── Cross_database_validation/ # Pattern validation across 10 databases
│           ├── Ribosome_paradox/          # Translation crisis core finding
│           ├── Publication_Figures/       # Manuscript-ready figures
│           └── ...                        # Additional visualization folders
├── docs/                          # Extended documentation
│   ├── bio_notes.md              # Curated biological notes and mechanistic model
│   ├── PATTERN_CLASSIFICATION.md # Pattern classification system documentation
│   ├── CHANGELOG.md              # Analysis evolution and version history
│   └── MIGRATION.md              # Container setup and migration notes
└── .deprecated/                   # Archived materials (old docs, backups)
```

## 🚀 Quick Start

### Running the Analysis Pipeline

**Prerequisites**: to run the analysis in the same computational enviroment it was tested in you will need to build and install the [scbio-dock v0.5.1](https://github.com/tony-zhelonkin/scbio-docker/tree/v0.5.1)

```bash
# Phase 1: Core analysis pipeline (run in order)
Rscript 02_Analysis/1.1.main_pipeline.R          
Rscript 02_Analysis/1.2.generate_contrast_tables.R
Rscript 02_Analysis/1.3.add_mitocarta.R
Rscript 02_Analysis/1.4.export_gsea_for_python.R

# Phase 2: Master tables (injested downstream)
python3 02_Analysis/1.5.create_master_pathway_table.py   # Master GSEA table
Rscript 02_Analysis/1.6.gsva_analysis.R                  # GSVA all pathways
Rscript 02_Analysis/1.7.create_master_gsva_table.R       # GSVA master tables

# Phase 3: R visualizations (run any/all, independent)
Rscript 02_Analysis/2.1.viz_ribosome_paradox.R
Rscript 02_Analysis/2.4.viz_critical_period_trajectories_gsva.R
# ... (see 02_Analysis/SCRIPTS.md for complete list)

# Phase 4: Python visualizations (after R exports)
python3 02_Analysis/3.1.publication_figures.py
python3 02_Analysis/3.4.pattern_summary_normalized.py

# Phase 5: Trajectory visualizations
python3 02_Analysis/3.7.viz_bump_chart.py                # Static bump charts
python3 02_Analysis/3.8.viz_interactive_bump_dashboard.py # Interactive explorer
```


### Script Organization

- **Main pipeline** (1.1-1.7): Core analysis, master tables - run in order
- **R visualizations** (2.1-2.6): Independent plots - run after Phase 1
- **Python visualizations** (3.1-3.9): Publication figures, bump charts - run after exports
- **Utilities** (Supp*.R/py): Verification, sensitivity analysis, explorers
- **See**: `02_Analysis/SCRIPTS.md` for complete inventory and dependencies

## 🧬 Experimental Design

| Component | Details |
|-----------|---------|
| **Cell type** | iPSC-derived cortical excitatory neurons |
| **Genotypes** | Ctrl, G32A (GTPase domain), R403C (stalk domain) - heterozygous DRP1 mutations |
| **Timepoints** | Day 35 (early maturation) and Day 65 (mature neurons) |
| **Replicates** | N=3 biological replicates per group (18 samples total) |
| **Contrasts** | 9 contrasts: mutation effects (6), maturation effects (3), interactions (2) |
| **Question** | How do domain-specific DRP1 mutations alter maturation trajectories and mitochondrial-synaptic coupling? |

**Contrast Framework**:
- **Early/Late effects**: G32A vs Ctrl (D35, D65), R403C vs Ctrl (D35, D65)
- **Maturation trajectories**: Time_Ctrl, Time_G32A, Time_R403C (D35 → D65 changes)
- **Mutation-specific maturation**: Maturation_G32A_specific, Maturation_R403C_specific (interactions)

## 📊 Key Findings

### Discoveries

DRP1 mutations seem to trigger a **translation crisis** where neurons attempt a compensation by increasing ribosome production, increasing mitochondrial function, but fail to maintain functional synaptic translation, probably because of diminished mitochondrial network at the synapse?:

| Pool | Early (D35) | Maturation (TrajDev) | Late (D65) | Pattern |
|------|-------------|----------------------|------------|---------|
| **Cytoplasmic ribosome biogenesis** | ↓ (NES -2.6) | **↑ (NES +2.25)** | Normal | Compensation |
| **Synaptic ribosomes (pre/post)** | ↑ (NES +2.3-2.5) | **↓ (NES -2.9 to -3.0)** | ↓ (NES -2.5) | Collapse |
| **Mitochondrial ribosomes** | ↓ (NES -2.2) | **↑ (NES +1.9)** | Normal | Compensation |

**Statistical support**: FDR < 0.05 for key effects, validated across 10 independent pathway databases

### Potential Mechanistic Model

```
DRP1 Mutation → Mitochondrial Positioning Failure
    ↓
Synaptic ATP Depletion
    ↓
┌─────────────────────────────────┬──────────────────────────────────┐
│ COMPENSATION (Pools 1 & 3)      │ FAILURE (Pool 2)                 │
├─────────────────────────────────┼──────────────────────────────────┤
│ ↑ Ribosome biogenesis           │                                  │
│ ↑ Mitochondrial translation     |                                  │
│ ↑ Mitochondrial function        │ ↓ Synaptic translation           │
│ ↑ Mitochondrial ribosomes       │ ↓ Pre/postsynaptic programs      │
│ (Adaptive response)             │ (Energy bottleneck)              │
└─────────────────────────────────┴──────────────────────────────────┘
    ↓
Translation Crisis at Synapses → Synaptic Dysfunction → Epilepsy
```

**Key insight**: Compensation fails presumably because the root cause (ATP delivery) remains unresolved—making more ribosomes doesn't help when they lack energy to function.

### Cross-Database Validation

Findings replicated across 10 independent databases:
- **SynGO** (synaptic ontology): 19-22 compensation pathways per mutation
- **MitoCarta** (mitochondrial): 45 (G32A) and 29 (R403C) compensation pathways
- **GO:BP/CC/MF** (Gene Ontology): 3000+ pathways show consistent patterns
- **MSigDB** (Hallmark, KEGG, Reactome, WikiPathways, TF targets): Corroborate findings

**Pattern distribution**: Compensation is the dominant response (~60% of significant pathways), with G32A showing stronger compensation than R403C.

## 🔄 Pattern Classification System

Pathway trajectories are classified into 8 mutually exclusive patterns based on three stages:
- **Early**: Mutation effect at D35 (immature neurons)
- **TrajDev**: Mutation-specific deviation from normal maturation
- **Late**: Mutation effect at D65 (mature neurons)

### Pattern Types

| Pattern | Active? | Description |
|---------|---------|-------------|
| **Compensation** | Active | Early defect + TrajDev opposes + Late improved |
| **Sign_reversal** | Active | Early defect + TrajDev opposes + sign flipped |
| **Progressive** | Active | Early defect + TrajDev amplifies + Late worsened |
| **Natural_improvement** | Passive | Early defect + no TrajDev + Late improved |
| **Natural_worsening** | Passive | Early defect + no TrajDev + Late worsened |
| **Late_onset** | - | No Early defect + Late defect emerges |
| **Transient** | - | Strong Early defect + Late resolved |
| **Complex** | - | Inconsistent or multiphasic |

**Active vs Passive**: Active patterns require significant TrajDev (p < 0.05, |NES| > 0.5), indicating transcriptional plasticity beyond normal development.

**Canonical source**: `01_Scripts/Python/pattern_definitions.py`
**Full documentation**: `docs/PATTERN_CLASSIFICATION.md`

## 📈 Trajectory Visualizations

### Bump Charts

Bump charts show pathway NES trajectories from Early to Late, with curvature representing TrajDev magnitude:

| Visualization | Script | Description |
|--------------|--------|-------------|
| **Static bump charts** | `3.7.viz_bump_chart.py` | PDF/PNG with weighted lines, curves, labels |
| **Interactive dashboard** | `3.8.viz_interactive_bump_dashboard.py` | HTML explorer with filtering, tooltips |

**How to read curves**:
- **Upward bulge**: Positive TrajDev (pathway upregulated during maturation)
- **Downward bulge**: Negative TrajDev (pathway downregulated during maturation)
- **Straight line**: No significant TrajDev

**Key outputs** (in `03_Results/02_Analysis/Plots/Trajectory_Flow/`):
- `bump_focused_FINAL_paper_combined.pdf` - Publication-ready figure
- `interactive_bump_dashboard.html` - Interactive explorer

**See**: `03_Results/02_Analysis/Plots/Trajectory_Flow/README.md` for visualization details.

## 📖 Key References

### Core Documentation

| File | Purpose |
|------|---------|
| [02_Analysis/SCRIPTS.md](02_Analysis/SCRIPTS.md) | Script inventory: active vs deprecated, dependencies, usage |
| [docs/PATTERN_CLASSIFICATION.md](docs/PATTERN_CLASSIFICATION.md) | Pattern classification framework (canonical reference) |
| [docs/CHANGELOG.md](docs/CHANGELOG.md) | Version history, reviewer tracking, major changes |

### Scientific Context

| File | Purpose |
|------|---------|
| [docs/bio_notes.md](docs/bio_notes.md) | Curated biological notes and mechanistic context |
| [docs/MIGRATION.md](docs/MIGRATION.md) | Container setup, environment configuration |

### Plot-Specific Documentation

Each major plot folder contains a comprehensive `README.md` with generating scripts, plot descriptions, methods, and interpretation guides.

**Key documented folders**:
- `Trajectory_Flow/` - Bump charts, alluvial diagrams, interactive explorers
- `Cross_database_validation/` - Pattern validation across 10 databases
- `Ribosome_paradox/` - Translation crisis core finding
- `Publication_Figures/` - Manuscript-ready integrated figures
- `Mito_translation_cascade/` - Energy → translation → synapse cascade

## 🛠️ Setup

### Prerequisites

- Docker and VS Code with Dev Containers extension
- Git with submodule support
- ~10 GB disk space for container and results

### Installation

```bash
# Clone repository
git clone <repo-url>
cd Gama_Vivian_DRP1_bulkRNAseq

# Initialize submodules (RNAseq-toolkit)
git submodule update --init --recursive

# Configure environment (optional, has defaults)
cp .env.example .devcontainer/.env
# Edit .devcontainer/.env to set UID/GID if needed

# Launch container
# Open folder in VS Code → Command Palette → "Dev Containers: Reopen in Container"
```

**Container**: [`scdock-r-dev:v0.5.1`](https://github.com/tony-zhelonkin/scbio-docker/tree/v0.5.1) (R 4.3+, Python 3.x, Bioconductor packages)

## 🔍 What to Explore First

### For Biological Interpretation

1. **Interactive Trajectory Explorer** (`03_Results/02_Analysis/Plots/Trajectory_Flow/`)
   - `interactive_bump_dashboard.html` - Explore all pathway trajectories interactively
   - `bump_focused_FINAL_paper_combined.pdf` - Key publication figure

2. **Publication Figures** (`03_Results/02_Analysis/Plots/Publication_Figures/`)
   - Fig1: Ribosome paradox visualization
   - Fig2: MitoCarta trajectory patterns
   - Fig3: Pattern classification summary across databases

3. **Cross-Database Validation** (`03_Results/02_Analysis/Plots/Cross_database_validation/`)
   - Pattern validation across 10 pathway databases
   - Compensation dominance (~60% of significant pathways)

### For Methods and Reproducibility

1. **Main pipeline script**: `02_Analysis/1.1.main_pipeline.R`
   - Checkpoint-based analysis (caches results for fast re-runs)
   - Complete DE and GSEA workflow

2. **Script inventory**: `02_Analysis/SCRIPTS.md`
   - Active vs deprecated scripts, dependencies, run order
   - Pattern classification workflow documentation

3. **Pattern framework**: `docs/PATTERN_CLASSIFICATION.md`
   - 8-pattern classification system
   - Active vs Passive distinction, thresholds

### For Data Access

- **Master tables** (comprehensive analysis outputs):
  - `03_Results/02_Analysis/master_gsea_table.csv` (109K rows, all pathways + patterns)
  - `03_Results/02_Analysis/master_gsva_*.csv` (GSVA scores, focused + all)
- **DE results**: `03_Results/02_Analysis/DE_results/*.csv` (9 contrast tables)
- **GSEA results**: `03_Results/02_Analysis/checkpoints/*.rds` (R) or `Python_exports/*.csv` (flat files)

## 🔬 Methodology

### Statistical Analysis Pipeline

#### 1. Data Preprocessing

**Upstream Processing (Alignment & Quantification):**

Raw FASTQ files were processed using custom bulk RNA-seq pipeline ([bulkRNAseq_pipeline_scripts](https://github.com/tony-zhelonkin/bulkRNAseq_pipeline_scripts)) with tools version-locked at [scbio-docker v0.2.0](https://github.com/tony-zhelonkin/scbio-docker/tree/v0.2.0).

**Reference genome:**
- **Assembly:** GRCh38.p14 (hg38, GCF_000001405.40)
- **Release date:** February 3, 2022
- **Annotation:** GCF_000001405.40_GRCh38.p14_genomic.gff.gz

**Reference preparation:**
```bash
# 1. Convert GTF to refFlat format
./GTFtoRefFlat.sh -i GCF_000001405.40_GRCh38.p14_genomic.gtf.gz \
                  -o GCF_000001405.40_GRCh38.p14_genomic.refflat

# 2. Extract ribosomal intervals for QC
./getRibosomalIntervals_from_gtf.sh \
    -i GCF_000001405.40_GRCh38.p14_genomic.gtf.gz \
    -r GCF_000001405.40_GRCh38.p14_genomic.fna.gz \
    -o GCF_000001405.40_GRCh38.p14_genomic.ribosomal_intervals

# 3. Convert GTF to BED12 format
./GTFtoBED12.sh -i GCF_000001405.40_GRCh38.p14_genomic.gtf.gz \
                -o GCF_000001405.40_GRCh38.p14_genomic.bed12

# 4. Create ReSeQC-compatible BED
./BEDtoRefSeqBED_human.sh -i GCF_000001405.40_GRCh38.p14_genomic.bed12 \
                          -o GCF_000001405.40_GRCh38.p14_genomic.reseqc.bed12
```

**Alignment:**
- **Software:** STAR (two-pass mode)
- **Index generation:**
```bash
STAR --runMode genomeGenerate \
     --runThreadN 8 \
     --genomeDir ../ref_index100 \
     --genomeFastaFiles ./GCF_000001405.40_GRCh38.p14_genomic.fna \
     --sjdbGTFfile ./GCF_000001405.40_GRCh38.p14_genomic.gtf \
     --sjdbOverhang 100
```

**Post-alignment QC:**
- **Strand inference:** RSeQC `infer_experiment.py` (determined non-stranded library)
- **Alignment metrics:** Picard `CollectRnaSeqMetrics`, `MarkDuplicates`
- **Read statistics:** samtools `flagstat`
- **Comprehensive reporting:** MultiQC aggregation

**Feature quantification:**
- **Software:** featureCounts (Subread package)
- **Parameters:**
  - Feature type: exon
  - Grouping attribute: gene_id
  - Strandedness: 0 (unstranded)
  - Mode: union (paired-end)
  - Threads: 8
- **Output:** Gene-level count matrix (`03_Results/01_Preprocessing/04_FeatureCounts/`)

**Normalization:**
- **Method:** TMM (Trimmed Mean of M-values) via edgeR
- **Filtering:** filterByExpr() via edgeR
- **Transformation:** edgeR::voomLmFit with sample.weights = FALSE
-
- voom log2-CPM with observation  weights

**Quality Control:**
- Sample correlation heatmaps (Pearson r > 0.9 within replicates)
- MDS plots for outlier detection
- Expression distribution assessment

#### 2. Differential Expression Analysis

**Software:** limma-voom (v3.56+)

**Model:**
```
Design matrix: ~0 + group
Where group = genotype × timepoint (6 levels)
```

**Statistical framework:**
1. Linear model fitting: `voomLmFit()` without precision weights
2. Contrast estimation: `contrasts.fit()` for 9 biological comparisons
3. Empirical Bayes moderation: `eBayes(robust = TRUE)`
4. Multiple testing: Benjamini-Hochberg FDR correction

**Significance thresholds:**
| Parameter | Value | Rationale |
|-----------|-------|-----------|
| p-value | 0.05 | Standard significance level |
| FDR | 0.05 | Controls false discovery rate |
| |log2FC| | > 1 | 2-fold change minimum for biological relevance |

#### 3. Gene Set Enrichment Analysis (GSEA)

**Software:** fgsea (v1.26+) via clusterProfiler wrapper

**Algorithm:**
- Ranking metric: Moderated t-statistic (preserves sign and precision)
- Permutations: 10,000 (gene set randomization)
- Gene set sizes: 10-500 genes
- Normalization: Enrichment Score normalized by set size

**Databases analyzed (12 total):**

| Category | Database | Gene Sets | Source |
|----------|----------|-----------|--------|
| Curated | Hallmark | 50 | MSigDB |
| Pathway | KEGG | 186 | MSigDB |
| Pathway | Reactome | 1,615 | MSigDB |
| Ontology | GO:BP | 7,658 | MSigDB |
| Ontology | GO:CC | 1,006 | MSigDB |
| Ontology | GO:MF | 1,738 | MSigDB |
| Pathway | WikiPathways | 664 | MSigDB |
| Curated | Canonical | 2,922 | MSigDB |
| Perturbation | CGP | 3,358 | MSigDB |
| Regulatory | TF targets | 1,137 | MSigDB |
| Synaptic | SynGO | ~300 | SynGO v1.1 |
| Mitochondrial | MitoCarta | 149 | MitoCarta 3.0 |

**Interpretation:**
- NES > 0: Pathway genes enriched in upregulated direction
- NES < 0: Pathway genes enriched in downregulated direction
- |NES| > 1.5: Moderate enrichment
- |NES| > 2.0: Strong enrichment

#### 4. Trajectory Pattern Classification

**Framework:** Early → TrajDev → Late

Contrasts mapped to developmental stages:
| Stage | G32A Contrast | R403C Contrast | Biological Question |
|-------|---------------|----------------|---------------------|
| Early | G32A_vs_Ctrl_D35 | R403C_vs_Ctrl_D35 | Initial mutation effect |
| TrajDev | Maturation_G32A_specific | Maturation_R403C_specific | Trajectory deviation |
| Late | G32A_vs_Ctrl_D65 | R403C_vs_Ctrl_D65 | Mature state effect |

**Pattern system:** 8 mutually exclusive patterns distinguishing active (significant TrajDev) vs passive (normal developmental buffering) responses. See [Pattern Classification System](#-pattern-classification-system) section above for pattern definitions.

**Key thresholds:**
- Significance: p.adjust < 0.05 (High), < 0.10 (Medium)
- Effect size: |NES| > 0.5 (minimum), > 1.0 (strong)
- Improvement: |Late|/|Early| < 0.7 (≥30% reduction)
- Worsening: |Late|/|Early| > 1.3 (≥30% increase)

**Full specification:** `docs/PATTERN_CLASSIFICATION.md`

### Reproducibility Features

**Checkpoint caching:**
- Expensive computations saved as RDS files
- Automatic loading on re-runs (10x speedup)
- Force recompute: `config$force_recompute = TRUE`

**Version control:**
- Analysis scripts tracked in git
- Helper functions via git submodule (`RNAseq-toolkit`)
- Container environment locked at [scbio-dock v0.5.1](https://github.com/tony-zhelonkin/scbio-docker/tree/v0.5.1)
- Runtime R packages documented in `R_session_info.txt` and `R_packages.txt`
- Python dependencies: `requirements.txt` (core packages) and `python_requirements_freeze.txt` (full freeze) 

### Software Versions

**R environment:**
- **Base:** R 4.3+ (via scbio-dock v0.5.1 container)
- **Package manifest:** `R_session_info.txt` (complete sessionInfo output)
- **Package list:** `R_packages.txt` (simple version listing)
- **Runtime installs:** See `02_Analysis/0.1.runtime_installs.R` for additional packages

**Python environment:**
- **Base:** Python 3.9+ (via scbio-dock v0.5.1 container)
- **Core packages:** `requirements.txt` (install with `pip install -r requirements.txt`)
- **Full freeze:** `python_requirements_freeze.txt` (exact versions for reproducibility)

**Key Python packages:**
- Data manipulation: pandas (1.3+), numpy (1.21+)
- Visualization: matplotlib (3.4+), seaborn (0.11+)
- Specialized plots: upsetplot (0.6+)

**To reproduce exact environment:**
```bash
# R packages - manually verify against R_session_info.txt
Rscript 02_Analysis/0.1.runtime_installs.R

# Python packages - install exact versions
pip install -r python_requirements_freeze.txt
```

## Supporting tools

### [scbio-dock v0.5.1](https://github.com/tony-zhelonkin/scbio-docker/tree/v0.5.1) 
RNAseq docker container for reproducible computational environments across lab projects 

### [RNAseq-toolkit](https://github.com/tony-zhelonkin/RNAseq-toolkit/tree/dev-GVDRP1) 
Custom RNAseq script automating DE, GSEA analysis and visualisations 

### [Preprocessing scripts](https://github.com/tony-zhelonkin/bulkRNAseq_pipeline_scripts)
Raw fastq files pre-processed with the custom scripts and bioinformatics tools version locked at [scbio-docker v0.2.0](https://github.com/tony-zhelonkin/scbio-docker/tree/v0.2.0)

### AI-assistants

This project utilized AI tools to accelerate analysis pipeline development, biological interpretation, and documentation:

**[Claude Code](https://github.com/anthropics/claude-code) (Anthropic):**
- **Primary tool** for refactoring initial scripts into a streamlined pipeline, debugging, and visualization implementation
- **Custom agents** developed for specialized tasks:
  - `bio-research-visualizer`: Deep web research for biological mechanism interpretation, literature synthesis, and visualization strategy recommendations. Used for cross-database validation framework development and ribosome paradox biological interpretation.
  - `rnaseq-insight-explorer`: RNAseq results exploration, pattern discovery, critical evaluation of biological claims, and data-driven visualization recommendations. Used for rough cross-database GSEA results querry and publication figure design.
  - `readme-auditor`: for quick alignment of the output files, code with a directory specific README

**[Gemini CLI](https://github.com/google-gemini/gemini-cli) (Google):**
- Primarily used for the implementation of the interactive pathway dashboard `3.8.viz_interactive_bump_dashboard.py`
- Supporting tool for specific documentation restructure and code review, web

**[Codex](https://github.com/openai/codex) (OpenAI):**
- Supporting tool for git management

**Scope:** AI tools assisted with code implementation, documentation generation, and preliminary biological context research at the stage of the RNAseq analysis. All scientific interpretations were validated against peer-reviewed literature, and all statistical analyses follow standard bioinformatics best practices. Human oversight guided all analytical decisions and biological conclusions.

**Transparency:** Analysis scripts, custom agent definitions (`.claude/agents/`), and project instructions (`CLAUDE.md`) are version-controlled, and committed to the repo for reproducibility .

## Documentation

### Quick Start
- **[SETUP.md](SETUP.md)** - Quick setup guide (post-migration)
- **[INSTALL.md](INSTALL.md)** - Detailed installation instructions

### Analysis Reference
- **[CLAUDE.md](CLAUDE.md)** - Claude Code instructions, architecture, troubleshooting
- **[02_Analysis/SCRIPTS.md](02_Analysis/SCRIPTS.md)** - Complete script inventory, dependencies, pattern workflow
- **[docs/PATTERN_CLASSIFICATION.md](docs/PATTERN_CLASSIFICATION.md)** - 8-pattern trajectory classification framework

### Scientific Documentation
See **[docs/](docs/)** directory for:
- **[Biological Notes](docs/bio_notes.md)** - Curated literature notes and mechanistic model
- **[CHANGELOG](docs/CHANGELOG.md)** - Analysis improvements since v1.0
- **[Migration Notes](docs/MIGRATION.md)** - Technical migration details

### Results Documentation
- **[Results README](03_Results/02_Analysis/README.md)** - Master tables documentation
- **[Trajectory Flow README](03_Results/02_Analysis/Plots/Trajectory_Flow/README.md)** - Bump charts and interactive explorers
- **[Plot READMEs](03_Results/02_Analysis/Plots/)** - Individual visualization guides

### Data Availability

**Master tables** (comprehensive summary outputs):
- `master_gsea_table.csv` - All GSEA pathways with pattern classifications (109K rows)
- `master_gsva_focused_table.csv` - GSVA scores for 7 key modules
- `master_gsva_all_table.csv` - GSVA scores for all pathways (87K rows)

**Raw outputs:**
- DE tables: `03_Results/02_Analysis/DE_results/*.csv`
- GSEA results: `03_Results/02_Analysis/checkpoints/*.rds`
- Visualizations: `03_Results/02_Analysis/Plots/`

**See:** Individual folder READMEs for detailed file documentation
