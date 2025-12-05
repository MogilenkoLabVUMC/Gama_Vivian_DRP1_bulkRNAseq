# Figure Captions: DRP1 Mutation RNA-seq Analysis

**Manuscript:** DRP1 mutations in iPSC-derived cortical neurons
**Last Updated:** 2025-12-05

---

## Figure X1: Sample Quality and Differential Expression Overview

**Full Caption:**

Transcriptomic analysis of iPSC-derived cortical neurons harboring DRP1 mutations. **(A)** Sample-to-sample Pearson correlation heatmap based on log2 counts per million (CPM) expression values. Samples (n=18) are ordered by timepoint (Day 35, Day 65) and genotype (Control, G32A, R403C). Color intensity indicates correlation coefficient (yellow-brown scale; r range 0.85-1.0). Annotation bars indicate genotype (green=Control, orange=G32A, purple=R403C) and timepoint (pink=D35, green=D65). **(B)** Principal component analysis showing sample clustering. Points represent individual samples colored by genotype and shaped by timepoint (circle=D35, triangle=D65). PC1 (35.1% variance) separates samples primarily by timepoint, while PC2 (15.3% variance) captures genotype-associated variation. **(C)** Volcano plots displaying differential expression results for nine contrasts organized by experimental comparison type. Top row: D35 mutation effects (G32A vs. Control, R403C vs. Control at Day 35). Second row: D65 mutation effects at Day 65. Third row: within-genotype maturation effects (D65 vs. D35). Bottom row: mutation-specific maturation effects (interaction contrasts). Calcium signaling genes (highlighted in orange/red) were prioritized based on known DRP1-calcium relationships. Each plot shows log2 fold-change (x-axis) versus -log10(raw p-value) (y-axis). Significance decisions use FDR (FDR≤0.1 shown); horizontal dashed line indicates the raw p-value corresponding to the FDR threshold; vertical dashed lines indicate |log2FC|=2 (4-fold change). Genes passing both FDR and FC thresholds shown in color. **(D)** Bar chart showing counts of differentially expressed genes (DEGs; FDR<0.05, no FC cutoff) per contrast. Orange=upregulated, blue=downregulated. Notably, mutation-specific maturation contrasts (Maturation_G32A_specific, Maturation_R403C_specific) yield zero DEGs at this stringency, motivating the use of gene set enrichment analysis to detect pathway-level patterns. **(E)** UpSet plot displaying intersection of DEG sets across contrasts. Horizontal bars show total DEGs per contrast; vertical bars show intersection sizes. Connected dots indicate which contrasts share DEGs. Time_G32A and Time_R403C show minimal overlap with Time_Ctrl, indicating distinct maturation trajectories in mutant neurons.

**Generating Scripts:**
- Panel A: `02_Analysis/1.1.main_pipeline.R` (sample correlation heatmap)
- Panel B: `02_Analysis/1.1.main_pipeline.R` (PCA/MDS plot)
- Panel C: `02_Analysis/1.1.main_pipeline.R` (volcano plots)
- Panel D: `02_Analysis/1.1.main_pipeline.R` (DE gene counts)
- Panel E: `02_Analysis/1.1.main_pipeline.R` (UpSet plot)

**Source Files:**
- A: `03_Results/02_Analysis/Plots/General/sample_correlation_heatmap_ordered.pdf`
- B: `03_Results/02_Analysis/Plots/General/MDS_plot.pdf`
- C: `03_Results/02_Analysis/Plots/Volcano/fdr/*.pdf`
- D: `03_Results/02_Analysis/Plots/General/DE_gene_counts.pdf`
- E: `03_Results/02_Analysis/Plots/General/UpSet_plot_all_contrasts.pdf`

---

## Figure X2: Pathway Trajectory Analysis and Synaptic Ribosome Deficits

**Full Caption:**

Gene set enrichment analysis reveals developmental trajectory patterns and synaptic ribosome dysfunction. **(A)** Normalized pattern distribution across 12 pathway databases for pathways with classifiable trajectory dynamics (excluding Complex patterns). Each database scaled to 100% to enable comparison of pattern proportions across databases of different sizes. Bars show proportion of pathways exhibiting each pattern; absolute counts (N) and percentages shown inside bars; total classifiable pathways per database shown as "n=X" at right. Pattern classification based on trajectory dynamics across Early (D35 mutation effect), TrajDev (mutation-specific maturation deviation), and Late (D65 mutation effect) stages. Compensation (green) represents active trajectory deviation opposing early defects; Sign_reversal (purple) captures complete direction reversal; Late_onset (pink) indicates maturation-dependent deficits emerging without early problems; Natural_improvement (cyan) represents passive developmental buffering without significant trajectory deviation. **(B)** Bump chart showing pathway trajectory dynamics for all pathways with at least one significant enrichment (p<0.05 at any stage; n=4,142 pathways per mutation). Lines connect Early (D35) to Late (D65) normalized enrichment scores (NES); line curvature reflects TrajDev magnitude and direction (upward bulge = positive trajectory deviation during maturation). Colors indicate classified patterns. Legend shows pattern counts per mutation. An interactive version of this visualization is available at [URL/supplementary materials]. **(C)** Focused bump chart showing trajectory dynamics for MitoCarta (mitochondrial) and SynGO (synaptic) pathways specifically (n=104 pathways per mutation). Key pathways labeled. Highlights include: presynaptic/postsynaptic ribosome pathways showing strong downward trajectories (green, Compensation pattern) transitioning from early upregulation to late downregulation; mitochondrial central dogma and ribosome assembly pathways showing upward compensation toward baseline. **(D)** GSEA running enrichment score plots for top 5 SynGO pathways at D35 for G32A (left) and R403C (right). Running enrichment score (y-axis) versus gene rank (x-axis; ranked by t-statistic from differential expression). Positive enrichment indicates coordinated upregulation; negative indicates coordinated downregulation. **(E)** GSEA running enrichment score plots for top 5 SynGO pathways at D65, showing the transition to negative enrichment for synaptic ribosome pathways in mature neurons. **(F)** Chord diagrams showing gene-pathway relationships for synaptic compartment pathways (G32A left, R403C right). Left arcs: individual genes with stacked rectangles showing log2 fold-change at D35 (inner) and D65 (outer); color scale from blue (downregulated) to orange (upregulated); gray indicates FDR>0.25. Right arcs: pathway categories (Presynaptic Ribosome, Postsynaptic Ribosome, Presynaptic Cytosol, Postsynaptic Cytosol, Postsynaptic Density, Postsynaptic Membrane). Ribbons connect leading-edge genes to their enriched pathways. **(X)** Expression heatmap of synaptic ribosome genes across developmental trajectory. Rows: ribosomal protein genes grouped by synaptic compartment localization (postsynaptic-only at top, both pre- and postsynaptic below). Columns: trajectory stages (Early, TrajDev, Late) for G32A and R403C mutations. Color indicates log2 fold-change (blue=downregulated, orange=upregulated; scale -0.6 to +0.6). Hierarchical clustering (Euclidean distance, complete linkage) reveals coordinated regulation of ribosomal subunits.

**Pattern Classification Framework:**
Patterns were classified based on significance (p.adjust<0.05) and effect size (|NES|>0.5) thresholds. Active patterns (Compensation, Sign_reversal, Progressive) require significant TrajDev indicating mutation-specific transcriptional adaptation. Passive patterns (Natural_improvement, Natural_worsening) show improvement or worsening without significant trajectory deviation. See Methods for complete pattern definitions.

**Generating Scripts:**
- Panel A: `02_Analysis/3.4.pattern_summary_normalized.py`
- Panels B-C: `02_Analysis/3.7.viz_bump_chart.py`, `02_Analysis/3.8.viz_interactive_bump_dashboard.py`
- Panels D-E: `02_Analysis/1.1.main_pipeline.R` (GSEA running sum plots)
- Panel F: `02_Analysis/3.7.viz_chord_diagrams.py`
- Panel X: `02_Analysis/2.3.viz_synaptic_ribosomes.R`

**Source Files:**
- A: `03_Results/02_Analysis/Plots/Pattern_Summary_Normalized/pattern_summary_normalized.pdf`
- B: `03_Results/02_Analysis/Plots/Trajectory_Flow/bump_curved_nes_significant.pdf`
- C: `03_Results/02_Analysis/Plots/Trajectory_Flow/bump_focused_FINAL_paper_combined.pdf`
- D: `03_Results/02_Analysis/Plots/GSEA/*/SynGO/*_running_sum.pdf` (D35 contrasts)
- E: `03_Results/02_Analysis/Plots/GSEA/*/SynGO/*_running_sum.pdf` (D65 contrasts)
- F: `03_Results/02_Analysis/Plots/Chord_Diagrams/chord_diagram_*.pdf`
- X: `03_Results/02_Analysis/Plots/Synaptic_ribosomes/Panel_C_Expression_Heatmap.pdf`

---

## Figure X3: Mechanistic Pathway Overview and Calcium Gene Expression

**Full Caption:**

Comprehensive pathway analysis reveals coordinated mitochondrial compensation and calcium signaling dysregulation. **(A)** Dotplot showing GSEA results for biologically-curated semantic pathway categories across all trajectory stages. Pathways organized by functional category (Synapse, Neuronal Development, ATP Synthase, Mitochondrial Metabolism, Mitochondrial Function, Mitochondrial Ribosome, Mitochondrial Translation, Ribosome Biogenesis, Cytoplasmic Ribosome, Cytoplasmic Translation, Calcium Signaling, Other). Each dot represents a pathway-contrast combination. Dot position: rows=pathways, columns=trajectory stage (Early, TrajDev, Late) for G32A and R403C. Dot color: Normalized Enrichment Score (NES; blue=downregulated, orange=upregulated; colorblind-safe Blue-White-Orange diverging palette). Dot size: statistical significance (-log10(FDR); larger dots=more significant). Black outline: FDR<0.05 (highly significant); gray outline: FDR 0.05-0.10. Filtered to show pathways with at least one significant result. Note coordinated patterns within functional modules: Synaptic pathways show early upregulation transitioning to late downregulation; Mitochondrial modules show compensation patterns with positive TrajDev values. **(B)** Gene expression heatmap for three mechanistic modules: ATP Synthase (Complex V; 20 genes from KEGG oxidative phosphorylation pathway), Calcium Signaling (12 curated genes including NNAT, CACNG3, STIM1/2), and Mitochondrial Central Dogma (25 top genes from MitoCarta pathway by effect size). Rows: individual genes with hierarchical clustering within each module. Columns: trajectory stages for G32A and R403C. Color: log2 fold-change (blue=downregulated, white=unchanged, orange=upregulated; scale -0.6 to +0.6). ATP Synthase and Mitochondrial Central Dogma modules display compensation patterns (early deficit followed by positive TrajDev); Calcium Signaling shows progressive downregulation that worsens from Early to Late, particularly for NNAT (neuronatin) and PNPO. **(C)** Sample-level expression heatmap of calcium signaling genes. Rows: calcium-related genes with hierarchical clustering. Columns: individual samples grouped by genotype and timepoint. Color: z-score normalized expression (blue=low, red=high). Clustering reveals coordinated expression modules: NNAT, CACNG3, ORAI1, and PNPO cluster together, showing progressive downregulation in both mutations.

**Module Curation:**
- ATP Synthase genes: F1 catalytic subunits (ATP5F1A-E), FO proton channel subunits (ATP5MC1-3, ATP5PB/D/F/O, ATP5ME/F/G/J/K), and assembly factors (ATPAF1/2, TMEM70). Curated from KEGG hsa00190 to exclude chromatin-remodeling ATPases.
- Calcium genes: Study-prioritized list based on known DRP1-calcium relationships and prior literature.
- Mitochondrial Central Dogma: Top 25 genes from MitoCarta 3.0 pathway ranked by mean |logFC| across TrajDev contrasts.

**Generating Scripts:**
- Panel A: `02_Analysis/3.2.publication_figures_dotplot.py`
- Panel B: `02_Analysis/2.2.viz_mito_translation_cascade.R`
- Panel C: `02_Analysis/2.6.viz_calcium_genes.R`

**Source Files:**
- A: `03_Results/02_Analysis/Plots/Publication_Figures_Dotplot/Fig4_Semantic_Pathway_Overview_dotplot.pdf`
- B: `03_Results/02_Analysis/Plots/Mito_translation_cascade/Mechanistic_Cascade_Heatmap.pdf`
- C: `03_Results/02_Analysis/Calcium_genes/calcium_genes_expression_heatmap.pdf`

---

## Supplementary Figure X1: Cross-Database GSEA Dotplots

**Full Caption:**

Cross-database pooled GSEA results for all contrasts. Dotplots display top 10 significantly enriched pathways (FDR<0.05) per database pooled across four curated pathway databases: KEGG (metabolic/signaling pathways), Reactome (biological processes), SynGO (synaptic ontology), and MitoCarta 3.0 (mitochondrial pathways). Database selection prioritizes curated resources most relevant to iPSC-derived cortical neuron biology. Non-neuronal tissue-specific pathways (pancreatic, cardiac, kidney, intestinal, lung, liver) were filtered out.

**Panel Organization:**
- **Top row:** D35 mutation effects (G32A_vs_Ctrl_D35, R403C_vs_Ctrl_D35) showing early pathway dysregulation
- **Second row:** D65 mutation effects (G32A_vs_Ctrl_D65, R403C_vs_Ctrl_D65) showing mature state
- **Bottom row:** Within-genotype maturation effects (Time_G32A, Time_R403C) showing D65 vs. D35 changes in mutant neurons

**Visual Encoding:**
- X-axis: Gene Ratio (leading edge genes / set size); higher values indicate larger fraction of pathway genes contributing to enrichment
- Y-axis: Pathway descriptions sorted by Gene Ratio
- Dot color: NES (blue=downregulated, orange=upregulated; colorblind-safe diverging palette)
- Dot size: -log10(FDR); larger dots indicate stronger statistical significance
- Black outline: FDR<0.05 (highly significant)

**Key Observations:**
- D35 mutation effects: Both mutations show upregulation of SynGO presynaptic/postsynaptic ribosome pathways (early stress response), cell cycle/DNA replication pathways, and mitochondrial pathways
- D65 mutation effects: Transition to downregulation of SynGO synaptic ribosome pathways (translation failure), continued mitochondrial pathway activation
- Maturation effects: Strong synaptic pathway involvement (SynGO terms) and DNA replication/cell cycle pathway activation during maturation in mutants

**Generating Script:** `02_Analysis/3.9.viz_pooled_dotplots.R`

**Source Files:**
- `03_Results/02_Analysis/Plots/GSEA/Cross_database_pooled/focused/G32A_vs_Ctrl_D35_pooled_focused.pdf`
- `03_Results/02_Analysis/Plots/GSEA/Cross_database_pooled/focused/R403C_vs_Ctrl_D35_pooled_focused.pdf`
- `03_Results/02_Analysis/Plots/GSEA/Cross_database_pooled/focused/G32A_vs_Ctrl_D65_pooled_focused.pdf`
- `03_Results/02_Analysis/Plots/GSEA/Cross_database_pooled/focused/R403C_vs_Ctrl_D65_pooled_focused.pdf`
- `03_Results/02_Analysis/Plots/GSEA/Cross_database_pooled/focused/Time_G32A_pooled_focused.pdf`
- `03_Results/02_Analysis/Plots/GSEA/Cross_database_pooled/focused/Time_R403C_pooled_focused.pdf`

---

## Technical Notes

### Statistical Methods

**Differential Expression:**
- TMM normalization (edgeR) followed by voom transformation (limma)
- Linear model with factorial design: `~ 0 + Group` where Group = Genotype_Timepoint
- Nine contrasts tested: 4 mutation effects (D35/D65 x G32A/R403C), 3 maturation effects (Time_Ctrl/G32A/R403C), 2 interaction contrasts (Maturation_G32A/R403C_specific)
- DEG classification (for counts/UpSet): FDR<0.05 via limma::decideTests with BH adjustment (no FC cutoff applied)
- Volcano plot thresholds: FDR≤0.1 (fdr mode) or p≤0.05 (p mode); |log2FC|≥2 (4-fold) for vertical dashed lines

**Gene Set Enrichment Analysis:**
- fgsea package with 100,000 permutations
- Genes ranked by moderated t-statistic from limma
- 12 pathway databases: MSigDB (Hallmark, KEGG, Reactome, GO:BP/CC/MF, WikiPathways, Canonical, CGP, TF), SynGO (synaptic ontology), MitoCarta 3.0 (mitochondrial)
- Significance threshold: FDR<0.05

**Pattern Classification:**
- Eight mutually exclusive patterns based on trajectory dynamics
- Significance requirements: p.adjust<0.05 (High confidence) or p.adjust<0.10 (Medium confidence) with |NES|>0.5
- Active patterns (Compensation, Sign_reversal, Progressive): Require significant TrajDev
- Passive patterns (Natural_improvement, Natural_worsening): TrajDev not significant
- Full specification: `docs/PATTERN_CLASSIFICATION.md`

### Reproducibility

All figures can be regenerated from the analysis pipeline:
```bash
# Main pipeline (generates GSEA checkpoints)
Rscript 02_Analysis/1.1.main_pipeline.R

# Create master tables
python3 02_Analysis/1.5.create_master_pathway_table.py    # Master GSEA table
Rscript 02_Analysis/1.6.gsva_analysis.R                   # Comprehensive GSVA
Rscript 02_Analysis/1.7.create_master_gsva_table.R        # Master GSVA tables

# Generate visualizations
Rscript 02_Analysis/2.1.viz_ribosome_paradox.R
Rscript 02_Analysis/2.2.viz_mito_translation_cascade.R
Rscript 02_Analysis/2.3.viz_synaptic_ribosomes.R
Rscript 02_Analysis/2.4.viz_critical_period_trajectories_gsva.R
python3 02_Analysis/3.1.publication_figures.py
python3 02_Analysis/3.2.publication_figures_dotplot.py
python3 02_Analysis/3.4.pattern_summary_normalized.py
python3 02_Analysis/3.7.viz_chord_diagrams.py
python3 02_Analysis/3.7.viz_bump_chart.py
Rscript 02_Analysis/3.9.viz_pooled_dotplots.R
```

---

**Document Version:** 1.0
**Generated:** 2025-12-05
