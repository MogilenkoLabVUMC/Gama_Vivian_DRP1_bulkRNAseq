# Vivian Gama: DRP1 BulkRNAseq Experiment

## Overview
This repository contains the computational analysis pipeline for investigating the transcriptional effects of DRP1 mutations on neuronal maturation using bulk RNA-seq data.

## Migration Note (2025-11-19)

This repository has been migrated to a new workstation environment. The frozen version corresponding to the published/submitted analysis is tagged as **v1.0** (commit `d6ec164`).

**Key Changes:**
- Updated to [scdock-r-dev:v0.5.1](https://github.com/tony-zhelonkin/scbio-docker/tree/v0.5.1) container with Docker Compose configuration
- All paths now use relative references for cross-workstation portability
- RNAseq-toolkit integrated as a git submodule (no external mounts required)
- Raw sequencing data dependencies removed (preprocessing complete, counts matrices available)

**Container Details:**
- **Image:** `scdock-r-dev:v0.5.1` from [scbio-docker](https://github.com/tony-zhelonkin/scbio-docker/tree/v0.5.1)
- **Features:** R 4.x with Bioconductor, Python 3.x with scientific stack, Jupyter, RStudio Server
- **Configuration:** Docker Compose with environment variable support

**For setup instructions, see:** `MIGRATION.md`

**For change tracking:** `CHANGELOG.md`

## Background
Mitochondrial fission, mediated by DRP1, is vital for the rapid metabolic shifts that occur during cortical development and plays a critical role in neuronal development and function. De-novo DNM1L loss-of-function mutations (G32A in the GTP-ase domain, R403C in the stalk domain) are emerging in paediatric neuro-developmental cohorts but their domain-specific impact on neuronal maturation is unclear. This study examines two domain-specific DRP1 mutations: G32A (GTPase domain) and R403C (stalk domain), both identified in pediatric neurodevelopmental disorders. The analysis compares transcriptional profiles between mutant and control iPSC-derived cortical neurons at two developmental timepoints (35 and 65 days in vitro).

## Experimental Design
- **Cell lines**: iPSC-derived cortical neurons with heterozygous G32A or R403C mutations vs. controls
- **Timepoints**: Day 35 and Day 65 of differentiation
- **Sample size**: 26 samples total across conditions and timepoints
- **Sequencing**: Bulk RNA-seq, paired-end reads
- **Reference genome**: GRCh38.p14 (hg38)


## Thoughts on analysis approaches

**Synaptic Gene Ontologies.** I found this SynGO database, one of the reviewers mentioned. The ontologies are downloadable. I can integrate this database into our pipeline and test DEG lists against this databases. That should highlight presynaptic vs postsynaptic signatures more cleanly than generic GO.

**Transcription factor programs.** I'm thinking if any transcription factors are involved. I hope I could try and run some inference tool (there are couple like decoupleR) that predicts Transcription factor activity from bulk RNAseq expression data. This might give a mechanistic bridge from elongated mitochondria -> TF programmes -> synaptic genes.

**Specific genes of interest.** Have a specific focus on a list of genes that are biologically relevant to keep an eye on.

**Gene Networks.** We could run weighted-gene correlation network (WGCNA) analysis on all the samples. This might uncover modules of co-regulated genes in your data. We could further test enrichment of these modules in SynGO or any other database of interest. We could also correlate the modules with some phenotype data, that you might have, say mitochondrial length, Ca²⁺ peak, synapse density or whatever, or make a correlation heatmap matrix of "co-regulated network module" X "phenotypic trait". 

Also I don't think you exploit your data enough in terms of statistical modelling. We could build a linear model to interrogate it with our biological questions, that would span beyond just the standard pair-wise comparisons. 

### Statistical Design and Contrasts

For example, we may of course ask ourselves, what's the difference between each kind of mutant vs control at each time-point.

Mathematically, this would encode into these 4 separate model questions:

```r
G32A_vs_Ctrl_D35  =  D35_G32A  - D35_Control,
R403C_vs_Ctrl_D35 =  D35_R403C - D35_Control,
G32A_vs_Ctrl_D65  =  D65_G32A  - D65_Control,
R403C_vs_Ctrl_D65 =  D65_R403C - D65_Control,
```

Each such contrast (G32A_vs_Ctrl_D35 etc.) gives its own list of DEGs, that we could test in a pathway analysis of any kind (GO, GSEA doesn't matter).

And of course we might ask, what's the maturation effect inside each genotype:

```r
Time_Ctrl   =  D65_Control - D35_Control, # What genes change during normal maturation?
Time_G32A   =  D65_G32A    - D35_G32A,    # What genes change during G32A mutant maturation?
Time_R403C  =  D65_R403C   - D35_R403C,   # What genes change during R403C mutant maturation?
```

Again, each question gives us its separate DEG list, which we could dig deeper into pathways further.

But we also might go for more complex **interaction questions** - "Difference-in-difference" questions:

```r
Int_G32A = (D65_G32A - D65_Control) - (D35_G32A - D35_Control),   # Does G32A's maturation trajectory differ from control trajectory?
Int_R403C = (D65_R403C - D65_Control) - (D35_R403C - D35_Control), # Does R403C's maturation trajectory differ from control?
```

Positive DEGs: are the ones differentially up-regulated over time only in the mutant.

These type of questions can be modelled too. They isolate genes whose time-course is selectively altered by the mutation, not just genes that differ at one time-point.

Combination of such Q&A would give us some insight. For example:

| Result pattern | Time_Ctrl | Time_G32A | Int_G32A | Interpretation |
|----------------|-----------|-----------|----------|----------------|
| gene doubles in Ctrl & doubles in G32A | +1 | +1 | 0 | same maturation → not highlighted |
| gene doubles in Ctrl but stays flat in G32A | +1 | 0 | –1 | G32A fails to up-regulate this gene |
| gene flat in Ctrl but doubles in G32A | 0 | +1 | +1 | gene activation is mutant-specific |

We could go and model other questions like this. For example, we might be interested in common DEGs mis-expressed similarly in both G32A and R403C (all DRP1 mutants) vs control. We also might include phenotypic traits in the model.


## Hypothesis 
Domain-specific DRP1 mutations differentially disturb mitochondrial dynamics, which in turn impairs synaptic development, calcium handling and the transcriptional programme of maturing cortical neurons.

## Experimental design
Generate iPSC lines harbouring heterozygous G32A or R403C mutations → dual-SMAD differentiation to glutamatergic cortical cultures → Phenotype mitochondria (SIM, live spinning-disk), synapses (SYP/PSD-95 SIM), and Ca²⁺ dynamics (Fluo-4 imaging) at 35, 65, 100 DIV

## Experimental findings
Both mutants retain the expected cortical-neuron marker profile but exhibit persistently hyper-elongated axonal mitochondria. Size-dependent, mutation-specific motility phenotypes (e.g. faster small mitochondria in R403C axons). Functional read-outs confirm diminished pre/post-synaptic marker volume and markedly exaggerated glutamate-evoked Ca²⁺ transients in both mutants

### Novelty
**Domain-resolved design** – analysing two distinct DRP1 mutations enables structure–function inference and partially explains heterogeneous patient presentations.  

**Multi-modal validation** – transcriptomic, imaging and functional assays converge on impaired synapse/Ca²⁺ regulation, lending biological credibility to RNA-seq findings.  


### Keep an eye out on 
**Calcium signaling genes:** NEURONATIN (NNAT), CACNG3, CACNA1C, CACNA1S, ATP2A1, RYR1, MYLK3, CASR, VDR, STIM1, STIM2, ORAI1, CALB1 and CALR

### Collaborator questions
1. The key question is: what are the effects of the G32A and R403C mutations in the transcriptional landscape of neurons?
2. We have 2 time points (35 DIV and 65 DIV), thus we were also wondering whether the mutations had any effects on maturation (compared to controls)


### General Approach Notes
The study asks **how two domain-specific DRP1 mutations (G32A = GTPase, R403C = stalk) disturb cortical-neuron maturation**.  

**Phenotypes observed:**

- hyper-elongated, sluggishly trafficked mitochondria
- widespread DEGs at 35 DIV (ion-channel, synapse, ROS) that largely normalise by 65 DIV
- smaller/rarer synapses and exaggerated glutamate-evoked Ca²⁺ transients
    

The reviewers want a **mechanistic bridge** from _mitochondrial shape → nuclear transcriptome → synaptic outcome_.




## Repository Structure

Below is a high-level overview of the repository:

```
├── 00_Data/                    # Reference data (SynGO database)
├── 01_Scripts/                 # Analysis scripts and custom functions
│   ├── RNAseq-toolkit/        # Git submodule: analysis toolkit and GSEA module
│   └── R_scripts/             # Custom R functions for this project
├── 02_Analysis/               # Main analysis pipeline
├── 03_Results/                # Output files and results
│   ├── 01_Preprocessing/      # Data preprocessing results (includes counts matrices)
│   └── 02_Analysis/           # Analysis results and plots
├── .devcontainer/             # Docker container configuration
├── MIGRATION.md               # Migration documentation
├── CHANGELOG.md               # Change tracking
└── README.md                  # This file
```

**Note:** The RNAseq-toolkit is a git submodule. Initialize it after cloning:
```bash
git submodule update --init --recursive
```
## Data Processing Workflow

### 1. Data Preprocessing



#### Preparing index

GRCh38.p14 reference genome (hg38, GCF_000001405.40, release date Feb 3, 2022) for read mapping. 'GCF_000001405.40_GRCh38.p14_genomic.gff.gz' was used as the annotation. To run our custom QC scripts we pre-processed the annotation in the following manner. 

**1. Prep refFlat from gtf**
```bash
./GTFtoRefFlat.sh -i GCF_000001405.40_GRCh38.p14_genomic.gtf.gz -o GCF_000001405.40_GRCh38.p14_genomic.refflat
```

**2. Prep ribosomal intervals**
```bash
./getRibosomalIntervals_from_gtf.sh -i GCF_000001405.40_GRCh38.p14_genomic.gtf.gz -r GCF_000001405.40_GRCh38.p14_genomic.fna.gz -o GCF_000001405.40_GRCh38.p14_genomic.ribosomal_intervals
```

**3. Get BED12 from GTF**
```bash
./GTFtoBED12.sh -i GCF_000001405.40_GRCh38.p14_genomic.gtf.gz -o GCF_000001405.40_GRCh38.p14_genomic.bed12
```

**4. (Optional) Change to ReSeQC compatible bed**
```bash
./BEDtoRefSeqBED_human.sh -i GCF_000001405.40_GRCh38.p14_genomic.bed12 -o GCF_000001405.40_GRCh38.p14_genomic.reseqc.bed12
```

**5. Alignment**

STAR two-pass alignment with genome index generation:
```bash
STAR --runMode genomeGenerate \
     --runThreadN 8 \
     --genomeDir ../ref_index100 \
     --genomeFastaFiles ./GCF_000001405.40_GRCh38.p14_genomic.fna \
     --sjdbGTFfile ./GCF_000001405.40_GRCh38.p14_genomic.gtf \
     --sjdbOverhang 100  
```

**6. Post-Alignment QC**
Post-alignment QC using multiple tools:
- **Strand inference**: RSeQC infer_experiment.py (determined non-stranded)
- **Alignment metrics**: Picard CollectRnaSeqMetrics, MarkDuplicates
- **Read statistics**: samtools flagstat
- **Comprehensive reporting**: MultiQC aggregation

**Infer experiment**. 
Experiment non-stranded
```bash
for file in 03_Results/01_Preprocessing/02_Alignment/aligned_bam/*.bam; do \
	echo "Strandedness of"
	infer_experiment.py \
		-i "$file" \
		-r /data/human_ref/ref_genome/GCF_000001405.40_GRCh38.p14_genomic.bed12 \
		-s 2000000 \
		-q 30
done
```

**Alignment metrics capture**
Running post-alignment QC scripts
```bash
docker run --rm -it \
  --name post_alignQC \
  -v "$HOME/projects/GVDRP1/03_Results/01_Preprocessing/02_Alignment/aligned_bam":/in:ro \
  -v "$HOME/projects/GVDRP1/03_Results/01_Preprocessing/03_Post_AlignmentQC":/out \
  -v "$HOME/projects/GVDRP1/03_Results/01_Preprocessing/03_Post_AlignmentQC/tmp":/out/tmp \
  -e TMPDIR=/out/tmp \
  -v "$HOME/data/GRCh38_hs_genome/ref_genome":/genome:ro \
  -v "$HOME/pipeline/bulkRNAseq_scripts":/pipeline:ro \
  scdock-r-dev:v0.2 \
  bash -c "/pipeline/scripts_bash/getPostAlignmentQC.sh \
             -i /in \
             -o /out \
             -j 2 \
             -t 4 \
             -r /genome/GCF_000001405.40_GRCh38.p14_genomic.fna \
             -e /genome/GCF_000001405.40_GRCh38.p14_genomic.ribosomal_intervals \
             -f /genome/GCF_000001405.40_GRCh38.p14_genomic.refflat \
             -b /genome/GCF_000001405.40_GRCh38.p14_genomic.reseqc.bed12 \
             --strand NONE \
             --copy-sorted"
```

**7. Feature quantification**

Gene-level quantification using featureCounts:

```bash
docker run --rm -it \
  --name feature_counts_run \
  -v "$HOME/projects/GVDRP1/03_Results/01_Preprocessing/02_Alignment/aligned_bam":/in:ro \
  -v "$HOME/projects/GVDRP1/03_Results/01_Preprocessing/04_FeatureCounts":/out \
  -v "$HOME/data/GRCh38_hs_genome/ref_genome":/genome:ro \
  -v "$HOME/pipeline/bulkRNAseq_scripts":/pipeline:ro \
  scdock-r-dev:v0.2 \
  bash -c "/pipeline/scripts_bash/runPostAlignmentQC.sh \
             -i /in \
             -o /out \
             -a /genome/GCF_000001405.40_GRCh38.p14_genomic.gtf \
             -s 0 \
             -t 8 \
             -f "exon" \
             -g "gene_id" \
             -p yes"
```

## 2. Analysis Pipeline

The main analysis is performed by `02_Analysis/Analysis_pipeline_fin.R`, which implements a comprehensive differential expression and pathway analysis workflow.

### Analysis Overview

The pipeline performs the following key analyses:

1. **Data preprocessing and quality control**
2. **Differential expression analysis** using limma-voom
3. **Multiple visualization approaches** (volcano plots, heatmaps, PCA)
4. **Pathway enrichment analysis** using MSigDB databases
5. **Specialized SynGO analysis** for synaptic gene ontologies
6. **Calcium signaling gene focus analysis**

### Statistical Design

#### Experimental Design Matrix
Factorial design 'genotype × timepoint' with the following contrasts:
**Pairwise Comparisons (Mutation vs Control)**
- `G32A_vs_Ctrl_D35`: G32A mutation effect at day 35
- `R403C_vs_Ctrl_D35`: R403C mutation effect at day 35  
- `G32A_vs_Ctrl_D65`: G32A mutation effect at day 65
- `R403C_vs_Ctrl_D65`: R403C mutation effect at day 65

**Maturation Effects (Time Course Within Genotype)**
- `Time_Ctrl`: Normal maturation trajectory (D65 vs D35 in controls)
- `Time_G32A`: Maturation in G32A mutants
- `Time_R403C`: Maturation in R403C mutants

**Interaction Effects (Mutation-Specific Maturation Changes)**
- `Maturation_G32A_specific`: G32A-specific alterations in maturation trajectory
- `Maturation_R403C_specific`: R403C-specific alterations in maturation trajectory


#### Statistical Methods
- **Normalization**: TMM (Trimmed Mean of M-values)
- **Differential expression**: limma-voom with empirical Bayes moderation
- **Multiple testing correction**: Benjamini-Hochberg FDR
- **Significance thresholds**: p < 0.05, |log2FC| > 1

### Custom Functions and Scripts

The analysis pipeline uses several custom functions located in `01_Scripts/R_scripts/`:


#### `read_count_matrix.R` - Data Processing Function
**File**: `01_Scripts/R_scripts/read_count_matrix.R`
**Function:** `process_rnaseq_data()`

**Purpose**: Integrates count matrix with sample metadata and creates analysis-ready DGEList object.

**Key Features:**
- Handles count matrix and metadata integration with flexible sample name matching
- Performs low-expression gene filtering using `filterByExpr()`
- Applies TMM normalization for library size differences
- Creates properly structured DGEList object with sample annotations

**Input**: Count matrix file, metadata CSV file
**Output**: Normalized DGEList object with sample annotations

**Usage:**
```r
DGE <- process_rnaseq_data(config$counts_file, config$metadata_file, annotate = FALSE)
```

#### `generate_vertical_volcanos.R` - Specialized Volcano Plot Generation
**File**: `01_Scripts/R_scripts/generate_vertical_volcanos.R`
**Function:** `generate_vertical_volcano_sets()`

**Purpose**: Creates comprehensive volcano plot collections for systematic comparison across contrasts.

**Key Features**:
- Predefined contrast groupings for biological interpretation
- Dual significance modes (p-value and FDR-based)
- Optional gene highlighting (e.g., calcium signaling genes)
- Multi-panel layout generation
- Automated file organization and naming

**Contrast Groupings**:
- D35 comparisons (early timepoint effects)
- D65 comparisons (late timepoint effects)
- Combined disease vs control comparisons
- Temporal effects within genotypes
- Interaction effects (mutation-specific maturation changes)

**Usage:**
```r
generate_vertical_volcano_sets(contrast_tables, config, highlight_calcium = TRUE)
```

#### `run_syngo_gsea.R` - Synaptic Gene Ontology Analysis

**File**: `01_Scripts/R_scripts/run_syngo_gsea.R`
**Function:** `run_syngo_gsea()`

**Purpose**: Performs GSEA using the SynGO database for synapse-specific pathway analysis.

**Key Operations**:
- SynGO cellular component ontology integration
- clusterProfiler GSEA execution with custom gene sets
- Comprehensive visualization generation (dotplots, barplots, running sum plots)
- Directional analysis (up/down regulated pathways)
- Results serialization and structured output

**Visualization Outputs**:
- Combined and directional dotplots
- NES (Normalized Enrichment Score) barplots
- Faceted up/down regulation plots
- Running sum plots for top pathways

**Helper Function:** `syngo_gmt()`
- Processes SynGO Excel files into TERM2GENE and TERM2NAME mappings
- Handles GO domain filtering (CC = Cellular Component)
- Cleans pathway descriptions for better visualization

**Usage:**
```r
syngo_gsea_results[[contrast]] <- run_syngo_gsea(
  tbl, contrast,
  T2G = syngo_lists$T2G,
  T2N = syngo_lists$T2N,
  sample_annotation = annot
)
```

#### `syngo_running_sum_plot.R` - Specialized GSEA Visualization
**File**: `01_Scripts/R_scripts/syngo_running_sum_plot.R`
**Function:** `syngo_running_sum_plot()`

**Purpose**: Creates publication-ready running sum plots specifically optimized for SynGO GSEA results.

**Key Features:**
- Handles both SYNGO: and GO: prefixed pathway identifiers
- Generates multi-panel running sum plots with:
  - Enrichment score curves
  - Gene hit locations
  - Ranked gene list visualization
- Applies consistent styling and color schemes
- Truncates long pathway names for readability
- Uses patchwork for professional multi-panel layout

**Usage:**
```r
syngo_running_sum_plot(
  gsea_obj = res, 
  gene_set_ids = top5,
  base_size = plot_par$font_size
)
```

## Analysis Outputs

The pipeline generates comprehensive outputs organized in `03_Results/02_Analysis/`:

#### Differential Expression Results
- **DE_results/**: CSV files with complete differential expression statistics for each contrast
- **Plots/Volcano/**: Multiple volcano plot variants
  - Standard p-value and FDR-based plots
  - Vertical comparison panels
  - Mean-difference (MD) plots
  - Fold-change vs B-statistic plots
- **Plots/Volcano/MD/**: Mean-difference plots for each contrast
- **Plots/Volcano/FC-B/**: Fold-change vs B-statistic plots

#### Quality Control and Overview
- **Plots/General/**: Sample correlation heatmaps, MDS plots, DEG count summaries
- **Summary/**: Analysis summary tables and optional HTML reports

#### Pathway Analysis
- **Plots/GSEA/**: MSigDB-based pathway analysis results for each contrast
  - Hallmark pathways, GO Biological Process, GO Cellular Component
  - Reactome and KEGG pathway databases
  - Individual contrast subdirectories with comprehensive plot sets

#### SynGO Analysis
- **Plots/GSEA/[contrast]/SynGO/**: Synapse-specific pathway analysis
  - Dotplots showing enriched synaptic processes
  - Running sum plots for top pathways
  - Directional analysis (up vs down regulated)

#### Calcium Gene Analysis
- **Calcium_genes/**: Focused analysis of calcium signaling genes
  - Expression heatmaps across conditions
  - Individual gene boxplots showing condition effects
  - Differential expression results for calcium genes
  - Volcano plots highlighting calcium genes

## Technical Implementation

#### Reproducible Workflow
- Centralized configuration management
- Robust path handling with here::here()
- Automatic package installation and loading
- Comprehensive error handling and logging

## Reproducibility

All analysis parameters are centralized in the configuration section of the main pipeline script. The workflow uses:
- Consistent file path handling with `here::here()`
- Automatic package installation and version checking
- Comprehensive logging and error reporting
- Structured output organization for easy navigation

## Dependencies

**Core Analysis**: edgeR, limma, dplyr, ggplot2, pheatmap, RColorBrewer, viridis
**Pathway Analysis**: clusterProfiler, msigdbr, fgsea, org.Hs.eg.db
**Visualization**: patchwork, UpSetR, VennDiagram, enrichplot
**Utilities**: here, reshape2, readxl

## Results

The analysis pipeline generates multiple complementary views of the data:

### Differential Expression Summary
- **Baseline effects**: Direct comparison of mutants vs controls at each timepoint
- **Maturation effects**: Time-course analysis within each genotype  
- **Interaction effects**: Mutation-specific alterations in maturation trajectory

### Pathway Analysis Insights
- **MSigDB analysis**: Broad pathway coverage including Hallmark, GO, Reactome, KEGG
- **SynGO analysis**: Synapse-specific processes and cellular components
- **Custom gene sets**: Focused analysis of calcium signaling pathways

### Visualization Outputs
- **Quality control plots**: Sample relationships and batch effect assessment
- **Volcano plots**: Multiple styles for different biological questions
- **Pathway plots**: Comprehensive GSEA visualization suite
- **Gene-specific analysis**: Calcium signaling gene expression patterns

## Future Directions

Potential extensions to the current analysis:

1. **Transcription Factor Analysis**: Integration with decoupleR for TF activity inference
2. **Network Analysis**: WGCNA for co-expression module identification
3. **Integration Analysis**: Correlation with phenotypic measurements (mitochondrial length, Ca²⁺ dynamics)
4. **Time Series Analysis**: More sophisticated temporal modeling approaches
5. **Single Cell Integration**: If scRNA-seq data becomes available
