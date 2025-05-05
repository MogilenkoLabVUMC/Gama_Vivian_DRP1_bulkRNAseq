# Yasmine-RetroT: Temperature-Induced Transposable Element Dynamics in T Helper Cells

## Overview
This repository contains the code and analysis pipeline used to investigate the effects of temperature exposure on transposable element (TE) expression in Th1 and Th17 helper T cell populations. We examine how fever-range temperature affects TE expression dynamics across multiple time points and correlate these changes with pathway enrichment.

I\`ve already explored the gene set pathways behaviour along the fever-range temperature exposure on this data [here](https://github.com/MogilenkoLabVUMC/T_cells_temperature-Yasmine), but we never tested the transposable elements hypothesis

## Background
Transposable elements (TEs) constitute a significant portion of mammalian genomes and are increasingly recognized for their roles in gene regulation and immune responses. This study explores how physiological temperature shifts (37°C to 39°C) influence TE expression patterns in distinct T helper cell subsets, which may contribute to understanding fever-associated immune modulation.

## Repository Structure

Below is a high-level overview of the repository:

* `1_Scripts/`
  * `runTrim.sh`, `runSTARalign_TEt.sh` — Preprocessing scripts used in the pipeline. Some archive scripts that were not directly used in this project.
  * `R_scripts/` — R scripts for TE parsing, enrichment, GSEA, correlation, etc.

* `2_Analysis/`
  * `1b.TE_Enrich_analysis_all_samples.R` — Initial TE analysis on all samples.
  * `1a.TE_Enrich_analysis_cleaned.R` — Analysis excluding A22, A14.
  * `2.TE_GSEA_clean.R` — GSEA pooling & TE–pathway correlation.

* `3_Results/`
  * `Plots/`
  * `GSEA/`
  * `TE/`
  * (… other result directories.)

* `README.md` (this file)

## Data Processing Workflow

### 1. Data Preprocessing

1. Read Trimming  
`1_Scripts/runTrim.sh`

* Raw sequencing reads were preprocessed using `Trimmomatic` with the following parameters:
```
LEADING=10
TRAILING=10
SLIDINGWINDOW="4:20"
MINLEN=36
ILLUMINACLIP_SETTINGS="2:30:10"
ADAPTERS="/usr/share/Trimmomatic-0.39/adapters/TruSeq3-PE.fa"
```
  * Leading/trailing quality trimming (threshold: 10)
  * Sliding window quality filtering (4:20)
  * Minimum read length: 36bp
  * Adapter removal using TruSeq3-PE adapters
  * Quality control was performed using FastQC and MultiQC both before and after trimming.

1. Alignment  
`1_Scripts/runSTARalign_TEt.sh`
```
"STAR --runThreadN 10 \
      --genomeDir \"$GENOME_DIR\" \
      --readFilesIn \"$R1_FILE\" \"$R2_FILE\" \
      --readFilesCommand zcat \
      --sjdbOverhang 100 \
      --limitBAMsortRAM 50000000000 \
      --twopassMode Basic \
      --outSAMattributes Standard \
      --outSAMattrRGline ID:$ID LB:$ID SM:$SM PL:$PL PU:$PU \
      --outFileNamePrefix \"$OUTPUT_DIR/${BASE_NAME}_\" \
      --outSAMunmapped None \
      --outSAMtype BAM Unsorted \
      --alignEndsType EndToEnd \
      --outMultimapperOrder Random \
      --runRNGseed 777 \
      --outFilterMultimapNmax 100 \
      --winAnchorMultimapNmax 200 \
      --outFilterMismatchNmax 10 \
      --outSAMprimaryFlag AllBestScore \
      --outFilterType BySJout \
      --outFilterScoreMinOverLread 0.4 \
      --outFilterMatchNminOverLread 0.4"
```
* Trimmed reads were aligned to the mouse reference genome (GRCm39) using STAR with parameters optimized for transposable element detection:
* Two-pass mode for improved splice junction detection
* outSAMtype BAM Unsorted for further TEtranscripts reads quantification
* Maximum multimapper count: 100
* Random assignment of multi-mapped reads with seed for reproducibility
* Primary alignment flagging: AllBestScore


3. Quantification
TE and gene expression was quantified using `TEtranscripts`, which handles the ambiguous assignment of reads mapping to repetitive elements. 
```
TEtranscripts --format BAM \
    --mode multi \
    --stranded no \
    -t A21_Aligned.out.bam A22_Aligned.out.bam A23_Aligned.out.bam A24_Aligned.out.bam \
    -c A37_Aligned.out.bam A38_Aligned.out.bam A39_Aligned.out.bam A40_Aligned.out.bam \
    --GTF /workspaces/Yasmine-retroT/0_Data/ref/ref_genome/GCF_000001635.27_GRCm39_genomic.gtf \
    --TE /workspaces/Yasmine-retroT/0_Data/ref/ref_genome/GRCm39_RefSeqAcc_rmsk_TE.gtf \
    --padj 0.05 \
    --minread 1 \
    --iteration 100 \
    --project Th17_24h_vs_0.5h
```
* “-t” references the treatment group (e.g., 0.5 vs 24h).
* “-c” references the control group.
* “--GTF” uses the curated gene GTF.
* “--TE” uses a TE annotation GTF.

Multiple contrasts were quantified separately and then merged into a unified count matrix for downstream analysis.  
The exact parameters have been written down in the `1_Scripts/TETranscripts.md`  
The merge parameters and checks are available in the `2_Analysis/1b.TE_Enrich_analysis_all_samples.R`  


### 2. Analysis Pipeline

1. Initial Data Analysis  
`2_Analysis/1b.TE_Enrich_analysis_all_samples.R`

* Complete dataset analysis including all samples.

2. Refined Data Analysis (Excluding Outliers)  
`2_Analysis/1a.TE_Enrich_analysis_cleaned.R`

* After quality assessment, PCA and heatmap visualisations – samples A22 and A14 were identified as outliers and excluded from subsequent analyses. 

These samples represented Th17 samples exposed to 39 degrees fro 24h, and 96h respectively, and seemed to deviate from the average Th17 pattern in terms of their overall gene and TE expression. 

1. Pathway and Correlation Analysis  
`2_Analysis/2.TE_GSEA_clean.R`

Gene Set Enrichment Analysis (GSEA) and correlation with TE expression patterns.

### Key Custom Functions

We used a couple of custom scripts to parse TE information, test for TE enrichment and visualize TE enrichment and mean expression. 

#### TE Annotation and Parsing 
`1_Scripts/R_scripts/parse_te_info.R`
**Function:** parse_te_info(name)
**Purpose:** Parses TE identifiers from annotations into structured information. Splits row names (e.g., “chr:subfamily:family”) to identify if a feature is a TE vs. a gene. Returns subfamily and family if it’s a TE

* TE status (True/False)
* Unique identifier
* Subfamily classification
* Family classification
* Example format: TE_ID:subfamily:family

#### TE Enrichment Analysis
`1_Scripts/R_scripts/te_enrichment_multi_contrast.R`
**Function:** `create_te_enrichment_plot(...)`
**Purpose:** Tests for the enrichment of TEs in different contrasts. Implements limma’s “geneSetTest” approach on TE expression data, grouping TEs by subfamily or family with multiple testing correction functionality.  The function calculates a p-value for each TE group by checking whether its log-fold changes deviate systematically from zero

* Uses gene set testing to detect collective expression changes of TE groups
* Calculates enrichment significance and direction (up/down-regulation)
* Multiple testing correction across contrasts or globally
* Filters by replication mode (cytoplasmic vs. nuclear)
* Visualizes results in bubble plots with size indicating significance and color showing fold change


#### TE Expression Visualization
`1_Scripts/R_scripts/create_te_heatmap.R` 
**Function:** `create_te_heatmap(...)`
**Purpose:** Creates heatmaps of TE expression with flexible options:

* Aggregation (mean expression) by TE subfamily or family
* Optional filtering for significant TEs
* Z-score normalization or raw expression values
* Custom sample annotations (Time, CellType, Temperature) and colorblind-friendly palettes
* Custom clustering options

#### Gene Set Enrichment Analysis
We leverage the [clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html) and [msigdbr](https://igordot.github.io/msigdbr/) packages:

`1_Scripts/R_scripts/runGSEA.R`  
`clusterProfiler::GSEA` wrapper function
**Function:** `runGSEA(...)`
**Purpose**: run `clusterProfiler::GSEA` against `MSigDB` gene sets

`1_Scripts/R_scripts/runGSEA_pool.R`  
**Function:** `run_pooled_gsea(...)`
**Purpose:** Automates GSEA over multiple contrasts, returning a combined result list.
• Summarizes significant pathways across all contrasts, merges them, and calculates pathway scores.

`1_Scripts/R_scripts/runGSEA_pool_helper.R`  
Functions provided:  
`get_significant_pathways(...)`  
`get_pathway_genes(...)` & `get_pathway_genes_all(...)`  
`calculate_pathway_scores(...)`  
**Purpose:** Collect significant pathways, extract the genes involved, and compute average expression for each pathway and each sample.

Implements a comprehensive GSEA pipeline for multi-contrast experiments:

* Uses preranked gene lists based on differential expression statistics
* Supports multiple pathway databases (Hallmark, KEGG, GO, Reactome, etc.)
* Pools significant pathways across contrasts for integrated analysis
* Calculates pathway scores for samples using core enrichment genes
* Handles errors robustly with detailed logging

**The mathematical approach:**
* Genes are ranked by statistical significance (t-statistic)
* GSEA iteratively evaluates whether the genes in a given pathway are concentrated at the top (or bottom) of the ranked list.
Permutations or sample-label shuffling estimates significance

`1_Scripts/R_scripts/plot_pathway_heatmap.R` is used to visualize the average expression scores of the genes involved in the top pathways identified on a heatmap. 

#### TE-Pathway Correlation Analysis
`1_Scripts/R_scripts/get_bulk_te_gsea_expr_correlation.R`

**Purpose:** Correlates TE expression with pathway activity:

* Calculates pairwise correlations between TE subfamily expression (mean expression score of the genes involved in the TE subfamily) and pathway scores (mean expression score of the core genes involved in the particular pathway)
* Uses Spearman correlation to capture monotonic relationships
* Generates heatmaps showing correlation patterns
* Supports multiple visualization options with colorblind-friendly palettes

**Steps:  **
1. Extract logCPM expression for TEs from the DGE object.
2. Calculate average pathway expression for each sample.
3. Compute correlation, visualized as a heatmap (using pheatmap).

**Values visualized:**
* Correlation coefficients between each TE subfamily/family and pathway
* Color intensity indicates correlation strength and direction
* Clustered to reveal patterns across TE-pathway relationships

## Results
The analysis reveals:  

* Distinct TE expression patterns between Th1 and Th17 cells – Th1 seem to derepress their TEs with cytosolic replication life cycle early in the fever-range temperature exposure 
* Derepressed TE subfamilies correlate with chromatin remodelling pathways and immune-related pathways such as cGAS. 

## Citation
If you use this code or analysis in your research, please cite this GitHub repository (and the forthcoming paper, eventually)

## Contact
For questions or issues, please open an issue in this repository or contact the repository maintainers.


# Future directions

1. Update the `1_Scripts/R_scripts/runGSEA_pool.R` to include a list of the top of the top pathways across all databases tested, to later visualise the top key pathways under one heatmap. 

2. I was thinking about doing kinda of a multifactorial correspondence analysis to visualise the TE subfamilies along with the GSEA pathways on a single latent space 2D plot to visually identify patterns of correspondence. But that might require some digging and data wrangling. (One could the develop the idea into visualizing a sort of a causal graph or DAG, connecting TEs and their relevant pathways).

3. Test `TEtranscripts` against the `featureCounts`

I planned to test `featureCounts` vs `TEtranscripts` downstream analysis implications. For that I created a unified gtf annotation, including the standard GRCm39 mouse genome annotation along with the transposable elements annotation provided with the `TEtranscripts`.

`cat GRCm39_RefSeqAcc_rmsk_TE.gtf | head -10`

```
NC_000067.7	mm39_rmsk	exon	8387807	8388657	3777	+	.	gene_id "Lx2B2"; transcript_id "Lx2B2"; family_id "L1"; class_id "LINE";
NC_000067.7	mm39_rmsk	exon	41942995	41943142	595	+	.	gene_id "B3"; transcript_id "B3"; family_id "B2"; class_id "SINE";
NC_000067.7	mm39_rmsk	exon	50331619	50332377	1796	-	.	gene_id "Lx7"; transcript_id "Lx7"; family_id "L1"; class_id "LINE";
NC_000067.7	mm39_rmsk	exon	58720078	58721182	5180	+	.	gene_id "L1MdV_III"; transcript_id "L1MdV_III"; family_id "L1"; class_id "LINE";
NC_000067.7	mm39_rmsk	exon	100663165	100663479	1316	+	.	gene_id "MLTR14"; transcript_id "MLTR14"; family_id "ERV1"; class_id "LTR";
NC_000067.7	mm39_rmsk	exon	109051878	109052234	2314	-	.	gene_id "Lx4B"; transcript_id "Lx4B"; family_id "L1"; class_id "LINE";
NC_000067.7	mm39_rmsk	exon	117440147	117440529	2330	-	.	gene_id "Lx2B"; transcript_id "Lx2B"; family_id "L1"; class_id "LINE";
NC_000067.7	mm39_rmsk	exon	142606236	142606423	328	+	.	gene_id "X9_LINE"; transcript_id "X9_LINE"; family_id "L1"; class_id "LINE";
NC_000067.7	mm39_rmsk	exon	176160666	176160805	560	+	.	gene_id "B1F1"; transcript_id "B1F1"; family_id "Alu"; class_id "SINE";
NC_000067.7	mm39_rmsk	exon	3145577	3147149	5715	+	.	gene_id "Lx9"; transcript_id "Lx9"; family_id "L1"; class_id "LINE";
```

`cat GCF_000001635.27_GRCm39_genomic.gtf | head -10`

```
#gtf-version 2.2
#!genome-build GRCm39
#!genome-build-accession NCBI_Assembly:GCF_000001635.27
#!annotation-date 02/01/2024
#!annotation-source NCBI RefSeq GCF_000001635.27-RS_2024_02
NC_000067.7	cmsearch	gene	3172239	3172348	.	+	.	gene_id "Gm26206"; transcript_id ""; db_xref "GeneID:115487594"; db_xref "MGI:MGI:5455983"; description "predicted gene, 26206"; gbkey "Gene"; gene "Gm26206"; gene_biotype "snRNA";
NC_000067.7	cmsearch	transcript	3172239	3172348	.	+	.	gene_id "Gm26206"; transcript_id "XR_004936710.1"; db_xref "GeneID:115487594"; db_xref "RFAM:RF00026"; gbkey "ncRNA"; gene "Gm26206"; inference "COORDINATES: nucleotide motif:Rfam:12.0:RF00026"; inference "COORDINATES: profile:INFERNAL:1.1.1"; product "U6 spliceosomal RNA"; transcript_biotype "snRNA";
NC_000067.7	cmsearch	exon	3172239	3172348	.	+	.	gene_id "Gm26206"; transcript_id "XR_004936710.1"; db_xref "GeneID:115487594"; db_xref "RFAM:RF00026"; gene "Gm26206"; inference "COORDINATES: nucleotide motif:Rfam:12.0:RF00026"; inference "COORDINATES: profile:INFERNAL:1.1.1"; product "U6 spliceosomal RNA"; transcript_biotype "snRNA"; exon_number "1";
NC_000067.7	BestRefSeq%2CGnomon	gene	3269956	3741733	.	-	.	gene_id "Xkr4"; transcript_id ""; db_xref "GeneID:497097"; db_xref "MGI:MGI:3528744"; description "X-linked Kx blood group related 4"; gbkey "Gene"; gene "Xkr4"; gene_biotype "protein_coding"; gene_synonym "Gm210"; gene_synonym "mKIAA1889"; gene_synonym "XRG4";
NC_000067.7	Gnomon	transcript	3269956	3741733	.	-	.	gene_id "Xkr4"; transcript_id "XM_006495550.5"; db_xref "GeneID:497097"; gbkey "mRNA"; gene "Xkr4"; model_evidence "Supporting evidence includes similarity to: 3 mRNAs, 3 ESTs, 4 Proteins, 766 long SRA reads, and 100% coverage of the annotated genomic feature by RNAseq alignments, including 17 samples with support for all annotated introns"; product "X-linked Kx blood group related 4, transcript variant X1"; transcript_biotype "mRNA";
```

Get the gene GTF without headers (skip lines starting with #)
`grep -v '^#' GCF_000001635.27_GRCm39_genomic.gtf > GRCm39_GCF_000001635.27_TE_RefSeq_combined.gtf`

Append the TE GTF
`cat GRCm39_RefSeqAcc_rmsk_TE.gtf >> GRCm39_GCF_000001635.27_TE_RefSeq_combined.gtf`

For now I ve only tested `TEtranscripts`. # Gama_Vivian_DRP1_bulkRNAseq
