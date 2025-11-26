# DRP1 Bulk RNA-seq Analysis

# Repositry analyst and maintainer 
[Anton Zhelonkin, MD](https://github.com/tony-zhelonkin)

# Intro

Transcriptional analysis of DRP1 mutations (G32A, R403C) in iPSC-derived cortical neurons across neuronal maturation (Day 35 → Day 65). This repository contains a complete downstream analysis pipeline starting from count matrices, revealing a translation paradox where ribosome biogenesis increases while synaptic translation fails.

## 📂 Repository Structure

```
├── 00_Data/                       # Reference databases (SynGO, MitoCarta)
├── 01_Scripts/
│   ├── RNAseq-toolkit/           # Git submodule with GSEA & DE helper functions
│   └── R_scripts/                # Project-specific helpers (SynGO, MitoCarta integration)
├── 02_Analysis/                   # Analysis pipeline scripts (see SCRIPTS.md)
│   ├── 1a.Main_pipeline.R        # Core DE analysis + GSEA with checkpointing
│   ├── 1b-3.*.R                  # Contrast tables, MitoCarta, Python exports
│   ├── viz_*.R                   # R visualization scripts (8 scripts)
│   ├── *.py                      # Python visualization scripts (4 scripts)
│   ├── SCRIPTS.md                # Complete script inventory and dependencies
│   └── .deprecated/              # Archived old scripts
├── 03_Results/
│   ├── 01_Preprocessing/         # Count matrices and metadata (input data)
│   └── 02_Analysis/
│       ├── checkpoints/          # Cached computation results (.rds files)
│       ├── DE_results/           # Differential expression tables (9 contrasts)
│       ├── Python_exports/       # GSEA data exported for Python (CSV files)
│       └── Plots/                # All visualizations (see folder READMEs)
│           ├── Cross_database_validation/  # Pattern validation across 10 databases
│           ├── Ribosome_paradox/          # Translation crisis core finding
│           ├── Publication_Figures/        # Manuscript-ready figures
│           ├── Mito_translation_cascade/   # Energy→translation→synapse cascade
│           └── Synaptic_ribosomes/         # Pre/postsynaptic translation analysis
├── docs/                          # Extended documentation
│   ├── SESSION_HISTORY.md        # Session-by-session development log
│   ├── NEXT_STEPS.md             # Current priorities and future directions
│   ├── MIGRATION.md              # Container setup and migration notes
│   └── biological_context.md     # Literature review and mechanistic model
└── .deprecated/                   # Archived materials (old docs, backups)
```

## 🚀 Quick Start

### Running the Analysis Pipeline

**Prerequisites**: to run the analysis in the same computational enviroment it was tested in you will need to build and install the [scbio-dock v0.5.1](https://github.com/tony-zhelonkin/scbio-docker/tree/v0.5.1)

```bash
# 1. Main differential expression and GSEA (uses checkpointing)
Rscript 02_Analysis/1a.Main_pipeline.R          # ~30-60 min first run, 5-10 min cached

# 2. Generate contrast tables and annotations
Rscript 02_Analysis/1b.generate_contrast_tables.R
Rscript 02_Analysis/2.add_MitoCarta.R
Rscript 02_Analysis/3.export_gsea_for_python.R

# 3. Generate R visualizations (run any/all, independent)
Rscript 02_Analysis/viz_ribosome_paradox.R
Rscript 02_Analysis/viz_mito_translation_cascade.R
Rscript 02_Analysis/viz_synaptic_ribosomes.R
# ... (see 02_Analysis/SCRIPTS.md for complete list)

# 4. Generate Python publication figures
python3 02_Analysis/6.publication_figures.py
python3 02_Analysis/8.pattern_summary_normalized.py
```

**Total runtime**: ~1-2 hours (first run), ~20-30 minutes (subsequent runs with caching)

### Script Organization

- **Main pipeline** (1a-3.R): Core analysis, run in order
- **Visualization scripts** (viz_*.R): Independent R plots, run in any order
- **Python scripts** (6-8.py): Publication figures, run after R exports
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

DRP1 mutations seem to trigger a **translation crisis** where neurons attempt a compensation by increasing ribosome production, but fail to maintain functional synaptic translation:

| Pool | Early (D35) | Maturation (TrajDev) | Late (D65) | Pattern |
|------|-------------|----------------------|------------|---------|
| **Cytoplasmic ribosome biogenesis** | ↓ (NES -2.6) | **↑ (NES +2.25)** | Normal | Compensation |
| **Synaptic ribosomes (pre/post)** | ↑ (NES +2.3-2.5) | **↓ (NES -2.9 to -3.0)** | ↓ (NES -2.5) | Collapse |
| **Mitochondrial ribosomes** | ↓ (NES -2.2) | **↑ (NES +1.9)** | Normal | Compensation |

**Statistical support**: FDR < 1e-11 for key effects, validated across 10 independent pathway databases

### Potential Mechanistic Model

```
DRP1 Mutation → Mitochondrial Positioning Failure
    ↓
Synaptic ATP Depletion
    ↓
┌─────────────────────────────────┬──────────────────────────────────┐
│ COMPENSATION (Pools 1 & 3)     │ FAILURE (Pool 2)                 │
├─────────────────────────────────┼──────────────────────────────────┤
│ ↑ Ribosome biogenesis          │ ↓ Synaptic translation           │
│ ↑ Mitochondrial ribosomes      │ ↓ Pre/postsynaptic programs      │
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

## 📖 Documentation

### Core Documentation

| File | Purpose |
|------|---------|
| [02_Analysis/SCRIPTS.md](02_Analysis/SCRIPTS.md) | Script inventory: active vs deprecated, dependencies, usage |
| [docs/CHANGELOG.md](CHANGELOG.md) | Version history, reviewer tracking, major changes |

### Scientific Context

| File | Purpose |
|------|---------|
| [docs/biological_context.md](docs/biological_context.md) | Literature review, mechanistic model, clinical relevance |
| [docs/MIGRATION.md](docs/MIGRATION.md) | Container setup, environment configuration |
| Plot folder READMEs | Scientific documentation for each analysis (see below) |

### Plot-Specific Documentation

Each major plot folder contains a comprehensive `README.md` with:
- Generating script and runtime
- Plot descriptions (axes, colors, statistics)
- Methods (GSEA parameters, statistical tests, sample sizes)
- Interpretation guide
- Key findings and biological context

**Documented folders**:
- `Cross_database_validation/` - Pattern validation across databases
- `Ribosome_paradox/` - Translation crisis core finding
- `Publication_Figures/` - Manuscript-ready integrated figures
- `Mito_translation_cascade/` - Energy → translation → synapse cascade
- `Synaptic_ribosomes/` - Compartment-specific ribosome analysis
- *(More READMEs available in other plot folders)*

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

1. **Publication Figures** (`03_Results/02_Analysis/Plots/Publication_Figures/`)
   - Fig1: Ribosome paradox visualization
   - Fig2: MitoCarta trajectory patterns
   - Fig3: Pattern classification summary across databases

2. **Ribosome Paradox** (`03_Results/02_Analysis/Plots/Ribosome_paradox/`)
   - Three-pool trajectory plot showing divergent compensation vs collapse
   - Data tables with statistics for all pathways

3. **Cross-Database Validation** (`03_Results/02_Analysis/Plots/Cross_database_validation/`)
   - 10 trajectory comparative plots (one per database)
   - Pattern summary showing compensation dominance

### For Methods and Reproducibility

1. **Main pipeline script**: `02_Analysis/1a.Main_pipeline.R`
   - Checkpoint-based analysis (caches results for fast re-runs)
   - Complete DE and GSEA workflow

2. **Script inventory**: `02_Analysis/SCRIPTS.md`
   - Active vs deprecated scripts
   - Dependencies and run order
   - Quick reference commands

### For Data Access

- **DE results**: `03_Results/02_Analysis/DE_results/*.csv` (9 contrast tables)
- **GSEA results**: `03_Results/02_Analysis/checkpoints/*.rds` (R objects) or `Python_exports/*.csv` (flat files)
- **Checkpoints**: `03_Results/02_Analysis/checkpoints/` (cached computation for fast re-analysis)

## 🔬 Methodology

### Statistical Analysis Pipeline

#### 1. Data Preprocessing

**Upstream Processing (Alignment & Quantification):**

Raw FASTQ files were processed using a standardized bulk RNA-seq pipeline ([bulkRNAseq_pipeline_scripts](https://github.com/tony-zhelonkin/bulkRNAseq_pipeline_scripts)) with tools version-locked at [scbio-docker v0.2.0](https://github.com/tony-zhelonkin/scbio-docker/tree/v0.2.0).

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
- **Strand inference:** RSeQC `infer_experiment.py` (determined: non-stranded library)
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
1. Linear model fitting: `voomLmFit()` with precision weights
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
| Late | G32A_vs_Ctrl_D65 | R403C_vs_Ctrl_D65 | Mature state effect |
| TrajDev | Maturation_G32A_specific | Maturation_R403C_specific | Trajectory deviation |

**Pattern definitions:**
| Pattern | Criteria | Interpretation |
|---------|----------|----------------|
| Compensation | Early defect + opposing TrajDev + Late improvement | Active adaptive response |
| Progressive | Early defect + amplifying TrajDev + Late worsening | Cumulative damage |
| Persistent | Stable defect, no TrajDev change | Static dysfunction |
| Late_onset | No Early defect, Late defect emerges | Maturation-dependent |
| Transient | Early defect, resolved by Late | Developmental delay |

**Thresholds:**
- Defect: |NES| > 0.5 minimum, > 1.0 for clear effect
- Change: |ΔNES| > 0.5 between stages

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
- **Runtime installs:** See `02_Analysis/0.runtime_installs.R` for additional packages

**Key R packages:**
- Differential expression: edgeR (3.42+), limma (3.56+), DESeq2 (1.40+)
- GSEA: clusterProfiler (4.8+), fgsea (1.26+), msigdbr (7.5+)
- Annotation: org.Hs.eg.db (3.17+)
- Co-expression: WGCNA (1.72+)
- Visualization: ggplot2 (3.4+), pheatmap (1.0+), patchwork (1.1+)

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
Rscript 02_Analysis/0.runtime_installs.R

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
- **Primary tool** for pipeline development, refactoring, and visualization implementation
- **Custom agents** developed for specialized tasks:
  - `bio-research-visualizer` (November 2025): Deep web research for biological mechanism interpretation, literature synthesis, and visualization strategy recommendations. Used for cross-database validation framework development and ribosome paradox biological interpretation.
  - `rnaseq-insight-explorer` (November 2025): RNAseq results exploration, pattern discovery, critical evaluation of biological claims, and data-driven visualization recommendations. Used for rough cross-database GSEA results querry and publication figure design.

**[Gemini CLI](https://github.com/google-gemini/gemini-cli) (Google):**
- Supporting tool for specific documentation restructure and code review 

**[Codex](https://github.com/openai/codex) (OpenAI):**
- Supporting tool for git management

**Scope:** AI tools assisted with code implementation, documentation generation, and preliminary biological context research at the stage of the RNAseq analysis. All scientific interpretations were validated against peer-reviewed literature, and all statistical analyses follow standard bioinformatics best practices. Human oversight guided all analytical decisions and biological conclusions.

**Transparency:** Analysis scripts, custom agent definitions (`.claude/agents/`), and project instructions (`CLAUDE.md`) are version-controlled, and committed to the repo for reproducibility .

### Data Availability

**Input data:**
- Count matrices: `03_Results/01_Preprocessing/04_FeatureCounts/count_matrices_fc/`
- Sample metadata: Same location

**Output data:**
- DE tables: `03_Results/02_Analysis/DE_results/*.csv`
- GSEA results: `03_Results/02_Analysis/checkpoints/*.rds`
- Visualizations: `03_Results/02_Analysis/Plots/`

**See:** Individual folder READMEs for detailed file documentation