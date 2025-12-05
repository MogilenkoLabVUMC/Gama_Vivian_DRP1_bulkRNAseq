# Methods

## Sample Collection and RNA Sequencing

Human induced pluripotent stem cells (iPSCs) from a healthy donor (PGP1 cell line) were differentiated into cortical neurons following established protocols. We generated isogenic lines carrying two DRP1 mutations using CRISPR-Cas9 gene editing: G32A (affecting the GTPase domain) and R403C (affecting the stalk domain). Cortical neurons were harvested at two developmental timepoints: Day 35 (early maturation stage) and Day 65 (late maturation stage). The experimental design consisted of three genotypes (Control, G32A, R403C), two timepoints (D35, D65), with 3 biological replicates per condition, yielding 18 samples total.

Total RNA was extracted using standard protocols and sequenced on an Illumina platform with paired-end 150 bp reads. Sequencing reads were aligned to the human reference genome GRCh38.p14 using STAR aligner. Gene-level quantification was performed using featureCounts from the Subread package, producing raw count matrices for downstream analysis.

## Differential Expression Analysis

### Normalization and Transformation

Raw count matrices were processed using the edgeR and limma R packages (v3.42+ and v3.56+, respectively). We applied Trimmed Mean of M-values (TMM) normalization via edgeR to account for sequencing depth and compositional biases across samples. Normalized counts were transformed to log2 counts-per-million (CPM) using the voom method from limma, which estimates precision weights for each observation based on the mean-variance relationship.

### Statistical Modeling

We used a factorial design to model the effects of genotype (Control, G32A, R403C) and timepoint (D35, D65) simultaneously:

```
design ~ 0 + genotype:timepoint
```

This design matrix allowed us to estimate effects for each genotype-timepoint combination directly. Linear models were fitted using `voomLmFit()`, which combines voom transformation and linear modeling in a single step for improved stability. Empirical Bayes moderation was applied using `eBayes()` with the robust setting to guard against outlier samples.

### Contrast Definitions

We defined nine contrasts organized into three biological categories to address distinct research questions:

**1. Mutation Effects (Baseline Comparisons)**

These contrasts quantify the direct effect of each mutation compared to control at each timepoint:

- `G32A_vs_Ctrl_D35`: G32A mutation effect at Day 35
- `R403C_vs_Ctrl_D35`: R403C mutation effect at Day 35
- `G32A_vs_Ctrl_D65`: G32A mutation effect at Day 65
- `R403C_vs_Ctrl_D65`: R403C mutation effect at Day 65

**2. Maturation Effects (Within-Genotype Time Comparisons)**

These contrasts capture developmental changes within each genotype:

- `Time_Ctrl`: Control maturation (D65_Control - D35_Control)
- `Time_G32A`: G32A maturation (D65_G32A - D35_G32A)
- `Time_R403C`: R403C maturation (D65_R403C - D35_R403C)

**3. Mutation-Specific Maturation Effects (Interaction Terms)**

These difference-in-differences contrasts identify how mutations alter the normal maturation trajectory:

- `Maturation_G32A_specific`: (D65_G32A - D35_G32A) - (D65_Control - D35_Control)
- `Maturation_R403C_specific`: (D65_R403C - D35_R403C) - (D65_Control - D35_Control)

These interaction contrasts represent the deviation of mutant maturation from control maturation. Algebraically, the Late mutation effect equals the sum of the Early mutation effect plus the mutation-specific maturation trajectory: Late = Early + TrajDev. This relationship is central to our trajectory-based pattern classification framework (see Pattern Classification System below).

### Multiple Testing Correction

We applied Benjamini-Hochberg false discovery rate (FDR) correction to control for multiple testing. Genes were classified as differentially expressed using limma's `decideTests()` function with FDR < 0.05. For visualization purposes, volcano plots additionally display fold-change thresholds (|log2 fold-change| ≥ 2, equivalent to 4-fold change) to highlight genes with both statistical significance and biological effect size.

## Gene Set Enrichment Analysis

### GSEA Implementation

We performed gene set enrichment analysis (GSEA) using the fgsea algorithm (v1.26+) via the clusterProfiler R package (v4.8+) wrapper. For each contrast, genes were ranked by the t-statistic from limma, which preserves both the direction and statistical significance of differential expression. We used 100,000 permutations to estimate enrichment p-values. Gene identifiers were converted from HGNC symbols to Entrez IDs using the org.Hs.eg.db annotation package, with multi-mapping genes assigned to their primary Entrez ID.

### Pathway Databases

We tested enrichment across 12 pathway databases to ensure comprehensive coverage and cross-validation of findings:

**MSigDB Collections (v7.5+, Homo sapiens):**
- Hallmark (50 pathways): Curated gene sets representing well-defined biological states
- KEGG (186 pathways): Metabolic and signaling pathways
- Reactome (1,615 pathways): Curated biological pathway reactions
- GO Biological Process (7,658 pathways): Gene Ontology biological processes
- GO Cellular Component (1,006 pathways): Gene Ontology cellular localization
- GO Molecular Function (1,738 pathways): Gene Ontology molecular activities
- WikiPathways (664 pathways): Community-curated pathway database
- Canonical Pathways (2,922 pathways): Combined curated pathway collection
- Chemical and Genetic Perturbations (3,358 pathways): Experimental perturbation signatures
- Transcription Factor Targets (1,137 pathways): Predicted regulatory targets

**Specialized Databases:**
- SynGO v1.1 (approximately 300 pathways, Cellular Component namespace): Synaptic gene ontology providing synapse-specific annotations
- MitoCarta 3.0 (149 pathways): Mitochondrial pathway annotations from the comprehensive mitochondrial proteome catalog

GSEA results were considered statistically significant at FDR < 0.05. The Normalized Enrichment Score (NES) was used as the primary metric for effect size, with |NES| > 0.5 indicating a biologically meaningful effect.

## Developmental Trajectory Framework

To characterize how pathway dysregulation evolves across neuronal maturation, we developed a three-stage trajectory framework. This descriptive framework facilitates the identification of biologically interesting temporal patterns but does not constitute formal statistical hypothesis testing.

### Trajectory Stages

Each pathway's enrichment trajectory is characterized by three stages:

1. **Early**: Initial mutation effect at the immature neuronal stage (Day 35), measured by the baseline mutation contrasts (e.g., `G32A_vs_Ctrl_D35`, `R403C_vs_Ctrl_D35`)

2. **TrajDev (Trajectory Deviation)**: How the mutation-specific maturation trajectory deviates from normal control maturation, measured by the interaction contrasts (e.g., `Maturation_G32A_specific`, `Maturation_R403C_specific`). This represents the mutation-specific component of developmental change: TrajDev = (D65_mut - D35_mut) - (D65_ctrl - D35_ctrl)

3. **Late**: Net mutation effect at the mature neuronal stage (Day 65), measured by the late-stage mutation contrasts (e.g., `G32A_vs_Ctrl_D65`, `R403C_vs_Ctrl_D65`)

### Biological Interpretation of TrajDev

A statistically significant TrajDev indicates that the mutant's maturation actively differs from control maturation, representing active transcriptional plasticity rather than passive developmental changes. This distinction is critical for identifying compensatory mechanisms: pathways with significant TrajDev that opposes the Early defect direction suggest active adaptive responses, whereas pathways that improve or worsen without significant TrajDev indicate passive developmental buffering.

## Pattern Classification System

We classified pathway enrichment trajectories into eight mutually exclusive patterns to facilitate discovery of biologically interesting dynamics. This classification is **descriptive** - it summarizes trajectory characteristics using predefined criteria but does not represent formal statistical hypothesis testing of pattern membership. Claims about specific pathway trajectories require individual validation and visual inspection of underlying data.

### Classification Criteria

All patterns require both statistical significance and biological effect size. We used two confidence levels:

- **High confidence**: FDR < 0.05 and |NES| > 0.5
- **Medium confidence**: FDR < 0.10 and |NES| > 0.5 (used for sensitivity analysis)

Primary results report high-confidence classifications unless otherwise noted.

### The Eight-Pattern Taxonomy

**Active Patterns** (require significant TrajDev opposing or amplifying the Early defect):

1. **Compensation**: Early shows a significant defect (FDR < 0.05, |NES| > 0.5), TrajDev is significant and opposes the Early defect direction (opposite sign, FDR < 0.05, |NES| > 0.5), and Late shows improvement (|Late NES| / |Early NES| < 0.7, representing at least 30% reduction in effect size, or |Late NES| < 0.5). This pattern indicates active adaptive plasticity where the system compensates for the initial mutation effect.

2. **Sign_reversal**: Early shows a significant defect, TrajDev is significant and opposes the Early defect direction, and Late shows the opposite sign from Early with substantial effect (|Late NES| > 0.5). This pattern captures complete trajectory reversals where the defect direction flips during development. Biological interpretation depends on pathway context: for synaptic ribosomes, the transition from Early upregulation to Late downregulation may represent failed compensation leading to disease phenotype; for mitochondrial Complex I, the transition from Early deficiency to Late upregulation may represent successful recovery.

3. **Progressive**: Early shows a significant defect, TrajDev is significant and amplifies the Early defect direction (same sign, FDR < 0.05, |NES| > 0.5), and Late shows worsening (|Late NES| / |Early NES| > 1.3, representing at least 30% increase). This pattern indicates cumulative damage with active maladaptive responses.

**Passive Patterns** (no significant TrajDev; changes follow normal developmental buffering):

4. **Natural_improvement**: Early shows a significant defect, TrajDev is not significant (FDR >= 0.05 or |NES| <= 0.5), and Late shows improvement (|Late NES| / |Early NES| < 0.7 or |Late NES| < 0.5). The defect resolves through normal developmental processes without requiring active compensatory transcriptional programs.

5. **Natural_worsening**: Early shows a significant defect, TrajDev is not significant, and Late shows worsening (|Late NES| / |Early NES| > 1.3). The defect worsens passively during development, indicating lack of adaptive capacity.

**Special Developmental Patterns**:

6. **Late_onset**: No significant Early defect (FDR >= 0.10 or |NES| <= 0.5), but a strong Late defect emerges (FDR < 0.05, |NES| > 1.0). These pathways exhibit maturation-dependent dysfunction, appearing normal at early stages but becoming significantly dysregulated in mature neurons.

7. **Transient**: Strong Early defect (FDR < 0.05, |NES| > 1.0) that fully resolves by Late stage (|Late NES| < 0.5). These pathways show developmental delays that recover by maturation, indicating timing-specific vulnerability.

8. **Complex**: Pathways that do not fit the above patterns, typically showing multiphasic or non-linear dynamics. These require individual inspection and may involve multiple regulatory mechanisms. For detailed analysis, Complex patterns are further subtyped into Stagnant (significant TrajDev but no outcome change), Weak_early (sign flip but subthreshold Early effect), Multiphasic (inconsistent trajectory dynamics), and Weak_signal (insufficient early effect size).

### Super-Category Grouping

For main text interpretation, we grouped the eight patterns into six super-categories to simplify the narrative:

- **Active_Compensation**: Compensation pattern
- **Active_Reversal**: Sign_reversal pattern
- **Active_Progression**: Progressive pattern
- **Passive**: Natural_improvement and Natural_worsening patterns
- **Late_onset**: Late_onset pattern
- **Other**: Transient and Complex patterns

Both the full eight-pattern taxonomy and super-category groupings are reported in Results and Supplementary Materials.

### Classification Algorithm Order

Pattern classification follows a hierarchical decision tree where evaluation order matters:

1. First, assess Late_onset (pathways with no Early defect but significant Late dysfunction)
2. For pathways with Early defects:
   - Check Active patterns first (significant TrajDev): Compensation → Sign_reversal → Progressive
   - Then check Transient (strong Early defect that fully resolves, evaluated after active patterns because significant opposing TrajDev takes precedence)
   - Then check Passive patterns (no significant TrajDev): Natural_improvement → Natural_worsening
3. Remaining pathways are classified as Complex

This order ensures that active transcriptional responses (indicated by significant TrajDev) are prioritized over passive developmental changes.

## Gene Set Variation Analysis

To complement the differential enrichment approach of GSEA, we performed Gene Set Variation Analysis (GSVA) using the GSVA R package (v1.48+). GSVA transforms gene expression data from a gene-by-sample matrix into a pathway-by-sample matrix of enrichment scores, enabling sample-level pathway activity assessment.

We applied GSVA with the Gaussian kernel method using two analysis modes. For comprehensive pathway-level analysis across all databases, we filtered pathways to those with 10-500 genes (size filter ensures statistical robustness while excluding very large gene sets that may dilute signal). For focused module analysis of seven pre-selected biological pathways, no upper size limit was applied since modules were chosen for biological relevance to the core findings. GSVA scores represent single-sample enrichment and range approximately from -1 to 1. For GSVA-based pattern classification, we used scaled thresholds: |GSVA score| > 0.15 (equivalent to |NES| > 0.5) for minimum effect size, and |GSVA score| > 0.30 (equivalent to |NES| > 1.0) for strong effects.

The seven focused modules represent the core biological findings: Ribosome Biogenesis (GO:BP, 158 genes), Cytoplasmic Translation (GO:BP, 76 genes), Synaptic Ribosomes (SynGO, 65 genes), Mitochondrial Ribosome (MitoCarta, 77 genes), Mitochondrial Ribosome Assembly (MitoCarta, 24 genes), mtDNA Maintenance (MitoCarta, 29 genes), and OXPHOS (MitoCarta, 139 genes).

Statistical comparisons of GSVA scores between groups were performed using two-sample t-tests with Benjamini-Hochberg FDR correction. Pattern classifications for GSVA modules used the same eight-pattern framework as GSEA, with one key difference: GSVA TrajDev represents a calculated difference-of-differences rather than a formally tested contrast, so the TrajDev significance criterion uses magnitude (|TrajDev| > 0.15) rather than p-values.

## Visualization and Statistical Software

All statistical analyses were performed in R version 4.3 or later. Key R packages included: edgeR (v3.42+), limma (v3.56+), fgsea (v1.26+), GSVA (v1.48+), clusterProfiler (v4.8+), msigdbr (v7.5+), org.Hs.eg.db (v3.17+), ggplot2 (v3.4+), pheatmap (v1.0+), and patchwork (v1.1+) for composite figures.

Publication figures were generated using Python 3.10+ with pandas (v2.0+), matplotlib (v3.7+), seaborn (v0.12+), and custom visualization functions. Pattern classification was implemented in Python using numpy (v1.24+) with the canonical pattern definitions defined in `01_Scripts/Python/pattern_definitions.py`. All visualization scripts applied colorblind-safe palettes from the Okabe-Ito and ColorBrewer sets.

The complete analysis pipeline, including all source code, pattern classification algorithms, and detailed documentation, is available in the project repository. Master summary tables containing all GSEA results with pattern classifications (`master_gsea_table.csv`, 109,990 pathway-contrast combinations), focused GSVA trajectory modules (`master_gsva_focused_table.csv`, 7 modules), and comprehensive GSVA pathway scores (`master_gsva_all_table.csv`, 87,001 pathway-group combinations) are provided as Supplementary Data files.

## Data Availability

The raw count matrices, metadata, differential expression results, and complete GSEA/GSVA outputs are available in the project repository. Full reproducibility is enabled through checkpoint caching of intermediate results and detailed documentation of all analysis parameters. The analysis pipeline version used for manuscript submission is permanently tagged as v1.0 (commit d6ec164) in the git repository.
