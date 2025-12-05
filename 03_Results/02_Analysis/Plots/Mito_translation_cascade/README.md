# Mitochondrial & Calcium Gene Expression Trajectories

## Purpose

Focused visualization of gene expression trajectories for three key functional modules relevant to DRP1 mutation pathology:

1. **ATP Synthase (Complex V)** - Mitochondrial ATP production machinery
2. **Calcium Signaling** - Curated genes of interest including NNAT and PNPO
3. **Mitochondrial Central Dogma** - mtDNA maintenance, transcription, and replication

This figure complements the synaptic ribosome analysis (Panel_C_Expression_Heatmap) by showing mitochondrial and signaling gene expression patterns across the developmental trajectory.

---

## Experimental Context

### Study Design
- **Cell type**: iPSC-derived cortical neurons
- **Mutations**: DRP1 G32A (GTPase domain) and R403C (stalk domain)
- **Timepoints**: Day 35 (early maturation) and Day 65 (late maturation)
- **Replicates**: n=3 biological replicates per condition
- **Sequencing**: Bulk RNA-seq, paired-end 150bp

### Contrasts Displayed
| Column Label | Full Contrast Name | Biological Meaning |
|--------------|-------------------|-------------------|
| G32A Early | G32A_vs_Ctrl_D35 | G32A mutation effect at Day 35 |
| G32A TrajDev | Maturation_G32A_specific | G32A-specific change during D35→D65 maturation (interaction term) |
| G32A Late | G32A_vs_Ctrl_D65 | G32A mutation effect at Day 65 |
| R403C Early | R403C_vs_Ctrl_D35 | R403C mutation effect at Day 35 |
| R403C TrajDev | Maturation_R403C_specific | R403C-specific change during D35→D65 maturation (interaction term) |
| R403C Late | R403C_vs_Ctrl_D65 | R403C mutation effect at Day 65 |

### TrajDev (Trajectory Development) Explained
The TrajDev column represents the **interaction term** from the linear model:
```
TrajDev = (Mutant_D65 - Mutant_D35) - (Control_D65 - Control_D35)
```
This captures mutation-specific developmental changes beyond normal maturation. A positive TrajDev indicates the mutation causes greater upregulation during maturation than controls; negative indicates greater downregulation.

---

## Generating Script

- **Script**: `02_Analysis/2.2.viz_mito_translation_cascade.R`
- **Runtime**: ~30 seconds
- **Dependencies**:
  - `03_Results/02_Analysis/checkpoints/mitocarta_gsea_results.rds` - MitoCarta GSEA results
  - `03_Results/02_Analysis/checkpoints/fit_object.rds` - limma model fit with logFC coefficients

---

## Figure Files

### Mechanistic_Cascade_Heatmap.pdf

**Main figure for publication**

#### How to Read This Figure

**Structure:**
- **Rows**: Individual genes (n=57 total) grouped into three functional modules
- **Columns**: Six trajectory conditions organized by mutation (G32A, R403C) and developmental stage (Early, TrajDev, Late)
- **Color**: Log2 fold-change (logFC) relative to control at each stage
  - **Blue** (negative logFC): Downregulated compared to control
  - **White** (logFC ≈ 0): No change compared to control
  - **Orange** (positive logFC): Upregulated compared to control

**Key visual elements:**
- **Dendrograms** (left of each module): Hierarchical clustering showing co-regulated gene groups within each module
- **Module color bar** (right): Color-coded identification of functional module membership
- **Column splits**: Visual separation between G32A and R403C mutations for direct comparison
- **Module labels** (left): Bold text identifying each functional module

#### Interpretation Guide

**Reading the trajectory (left to right within each mutation):**
- **Early**: Initial mutation effect at Day 35 (before full maturation)
- **TrajDev**: How the mutation specifically affects the maturation process (D35→D65)
- **Late**: Final mutation effect at Day 65 (mature neurons)

**Pattern recognition:**
| Pattern | Visual Signature | Biological Interpretation |
|---------|-----------------|--------------------------|
| Compensation | Blue → Orange → White | Early deficit rescued by developmental upregulation |
| Progressive | Blue → Blue → Darker Blue | Deficit worsens over development |
| Late onset | White → Any → Blue | Deficit emerges only in mature neurons |
| Persistent | Blue → White → Blue | Consistent deficit, no developmental rescue |

#### Module-Specific Findings

**ATP Synthase (Complex V) - 20 genes:**
- Shows **compensation pattern** in G32A: Early downregulation (blue) → TrajDev upregulation (orange) → Late near-baseline
- ATP5 subunit genes (F1 and FO components) cluster together showing coordinated regulation
- Assembly factors (ATPAF1/2, TMEM70) cluster with their targets
- R403C shows similar but weaker pattern

**Calcium Signaling - 12 genes:**
- Shows **progressive downregulation** pattern in both mutations
- **NNAT** (neuronatin): Strong progressive downregulation (increasingly blue from Early→Late)
- **PNPO**: Similar progressive deficit pattern
- Store-operated calcium entry genes (STIM1, STIM2, CALR) form a distinct cluster
- Indicates calcium homeostasis dysfunction that worsens with neuronal maturation

**Mt Central Dogma - 25 genes:**
- Shows **strong compensation pattern** in G32A
- Early deficit in mtDNA maintenance/transcription genes (POLQ, DNA2, RECQL4)
- Robust TrajDev upregulation (orange) suggesting adaptive transcriptional response
- Mitochondrial ribosomal proteins (MRPL/MRPS genes) show coordinated upregulation
- Pattern suggests compensatory mitochondrial biogenesis response

### Module_Summary_Heatmap.pdf

**Supplementary figure showing module-level averages**

- Shows mean logFC across all genes within each module
- Numeric values displayed in cells for quantitative reference
- Useful for statistical comparisons between modules and mutations

---

## Methods

### Differential Expression Analysis

**Pipeline:**
1. Raw counts from featureCounts (GRCh38.p14 reference)
2. TMM normalization (edgeR)
3. Voom transformation with sample quality weights (limma)
4. Linear model fit with factorial design: `~ 0 + Group` where Group = Genotype_Timepoint
5. Contrasts extracted for mutation effects and interaction terms

**Statistical model:**
```r
design <- model.matrix(~ 0 + Group, data = metadata)
contrasts <- makeContrasts(
  G32A_vs_Ctrl_D35 = G32A_D35 - Ctrl_D35,
  G32A_vs_Ctrl_D65 = G32A_D65 - Ctrl_D65,
  Maturation_G32A_specific = (G32A_D65 - G32A_D35) - (Ctrl_D65 - Ctrl_D35),
  # ... similar for R403C
  levels = design
)
```

### Gene Selection

**ATP Synthase (Complex V) - 20 genes:**

Curated from KEGG Oxidative Phosphorylation pathway (hsa00190), specifically Complex V (ATP synthase) subunits:

| Subcomplex | Genes | Function |
|------------|-------|----------|
| F1 (catalytic) | ATP5F1A, ATP5F1B, ATP5F1C, ATP5F1D, ATP5F1E | ATP synthesis catalysis |
| FO (proton channel) | ATP5PB, ATP5MC1, ATP5MC2, ATP5MC3, ATP5PD, ATP5PF, ATP5MF, ATP5MG, ATP5PO, ATP5ME | Proton translocation |
| Peripheral stalk | ATP5MJ, ATP5MK | Structural connection |
| Assembly factors | ATPAF1, ATPAF2, TMEM70 | Complex V assembly |

**Important note:** This gene set is distinct from GO:CC "ATPase complex" (GO:0043234) which includes SWI/SNF chromatin remodeling ATPases and would be inappropriate for mitochondrial analysis.

**Calcium Signaling - 12 genes:**

Curated list of calcium-related genes prioritized for this study:

| Gene | Full Name | Calcium Role |
|------|-----------|--------------|
| **NNAT** | Neuronatin | Imprinted gene; regulates ER calcium; key study target |
| **PNPO** | Pyridoxamine 5'-phosphate oxidase | Vitamin B6 metabolism; affects calcium via neurotransmitter synthesis |
| CACNG3 | Voltage-dependent calcium channel gamma-3 | AMPA receptor trafficking |
| CACNA1S | Voltage-dependent L-type calcium channel alpha-1S | Excitation-contraction coupling |
| ATP2A1 | SERCA1 | ER calcium reuptake |
| RYR1 | Ryanodine receptor 1 | ER calcium release |
| MYLK3 | Myosin light chain kinase 3 | Calcium-dependent kinase |
| VDR | Vitamin D receptor | Calcium homeostasis transcription factor |
| STIM1 | Stromal interaction molecule 1 | ER calcium sensor (SOCE) |
| STIM2 | Stromal interaction molecule 2 | ER calcium sensor (SOCE) |
| CALB1 | Calbindin 1 | Calcium buffering protein |
| CALR | Calreticulin | ER calcium buffering |

**Note:** ORAI1 was in the original config but not detected in expression data.

**Mt Central Dogma - 25 genes:**

Selected from MitoCarta 3.0 "Mitochondrial_central_dogma" pathway:
- **Source**: GSEA core enrichment genes from Maturation_G32A_specific contrast
- **GSEA statistics**: NES = 2.41, p.adjust = 2.14e-13, pathway size = 216 genes
- **Selection**: Top 25 genes ranked by mean |logFC| across TrajDev contrasts
- **Rationale**: Limits visual complexity while retaining highest-effect genes

Functional categories represented:
| Category | Genes |
|----------|-------|
| mtDNA maintenance | POLQ, DNA2, RECQL4, PIF1, UNG |
| mtDNA replication | PRIMPOL, EXOG, ENDOG |
| Mitochondrial transcription | MTERF3, MTG1 |
| mt-tRNA modification | PUS1, MRM2, TRMT61B, DARS2 |
| Mitochondrial ribosome | MRPL1, MRPL11, MRPL15, MRPS11, MRPS17, MRPS18C |
| Other | DDX28, GUF1, MPV17L2, PDF, RBFA |

### Visualization Parameters

**Clustering:**
- **Method**: Hierarchical clustering (hclust)
- **Distance metric**: Euclidean distance
- **Linkage**: Complete linkage
- **Scope**: Applied independently WITHIN each module (not across modules)
- **Module order**: Fixed (ATP Synthase → Calcium → Mt Central Dogma), not affected by clustering

**Color scale:**
- **Type**: Symmetric diverging gradient
- **Range**: -0.6 to +0.6 logFC
- **Colors**: Blue (#2166AC) → White (#F7F7F7) → Orange (#B35806)
- **Rationale**: Narrower than typical ±1.5 scale to reveal subtle expression changes; values outside range are clipped to endpoints

**Software:**
- R version 4.x
- ComplexHeatmap v2.24.1
- circlize v0.4.16

---

## Paper Caption Suggestions

### Main Text Figure Caption (Short)

> **Figure X. Mitochondrial and calcium gene expression trajectories in DRP1 mutant neurons.**
> Heatmap showing log2 fold-change (logFC) of 57 genes across three functional modules in G32A and R403C DRP1 mutations. Columns represent developmental trajectory stages: Early (D35 mutation effect), TrajDev (mutation-specific maturation change), and Late (D65 mutation effect). Color indicates logFC vs. control (blue: downregulated; orange: upregulated; scale: ±0.6). Dendrograms show hierarchical clustering within modules. ATP Synthase and Mt Central Dogma show compensation patterns; Calcium Signaling shows progressive downregulation including NNAT.

### Main Text Figure Caption (Detailed)

> **Figure X. Gene expression trajectories reveal compensatory and progressive patterns in mitochondrial and calcium signaling pathways.**
> (A) Heatmap displaying log2 fold-change (logFC) values for genes in three functional modules across the developmental trajectory in DRP1 G32A and R403C mutant iPSC-derived cortical neurons. Rows represent individual genes (n=57 total): ATP Synthase/Complex V (n=20, curated KEGG hsa00190 subunits), Calcium Signaling (n=12, study-prioritized genes including NNAT and PNPO), and Mitochondrial Central Dogma (n=25, top genes from MitoCarta pathway by effect size). Columns represent trajectory stages: Early (D35 mutation vs. control), TrajDev (mutation-specific maturation effect, interaction term), and Late (D65 mutation vs. control). Color scale indicates logFC (blue: downregulated; white: unchanged; orange: upregulated; range ±0.6). Dendrograms (left) show hierarchical clustering of genes within each module, revealing co-regulated subgroups. ATP Synthase and Mt Central Dogma modules display compensation patterns characterized by Early deficits (blue) followed by TrajDev upregulation (orange). Calcium Signaling genes, including NNAT (neuronatin), show progressive downregulation that worsens from Early to Late stages, suggesting impaired calcium homeostasis in mature mutant neurons.

### Supplementary Figure Caption

> **Supplementary Figure X. Module-level summary of expression trajectories.**
> Mean log2 fold-change (logFC) values averaged across all genes within each functional module. Values displayed in cells. ATP Synthase (Complex V) and Mt Central Dogma modules show positive TrajDev values (G32A: +0.14 and +0.69 respectively) indicating compensatory upregulation during maturation. Calcium Signaling module shows consistently negative values across all stages (Early: -0.39, TrajDev: -0.52, Late: -0.91 for G32A), indicating progressive deficit without compensation.

---

## Peer Reviewer Notes

### Methodological Transparency

#### 1. Gene Selection Rationale

**Why curated gene sets rather than pathway databases?**

| Module | Selection Method | Justification |
|--------|-----------------|---------------|
| ATP Synthase | Manual curation from KEGG | GO:CC "ATPase complex" includes chromatin remodelers (SWI/SNF); KEGG hsa00190 Complex V is specific to mitochondrial ATP synthase |
| Calcium Signaling | Study-specific priority list | Captures genes of biological interest (NNAT, PNPO) not in standard calcium pathways; hypothesis-driven |
| Mt Central Dogma | GSEA core enrichment + effect size filter | Data-driven selection from significant pathway; limited to top 25 for visual clarity |

#### 2. Statistical Considerations

**Expression values:**
- logFC values are coefficients from limma-voom linear model
- Represent log2(Mutant/Control) at each condition
- TrajDev is the interaction term, not a direct measurement

**Significance:**
- Individual gene significance not shown (no asterisks)
- Module-level significance from GSEA: Mt Central Dogma p.adj = 2.14e-13
- Calcium and ATP Synthase modules are curated (not GSEA-derived)

**Multiple testing:**
- No correction applied at gene level for this visualization
- Genes selected by biological criteria, not statistical thresholds

#### 3. Clustering Methodology

- Clustering is applied WITHIN modules only
- Does NOT affect module order or cross-module interpretation
- Purpose: reveal co-regulated gene subgroups for biological insight
- Default ComplexHeatmap parameters (Euclidean distance, complete linkage)

#### 4. Color Scale Choice

- Range ±0.6 chosen to match synaptic ribosome heatmap (Panel_C) for visual consistency
- Typical range (±1.5 or ±2) would compress most values near white
- Extreme values (|logFC| > 0.6) are clipped to endpoints
- This is a **visualization choice** and does not affect underlying data

### Reproducibility

**To regenerate this figure:**
```bash
cd /workspaces/Gama_Vivian_DRP1_bulkRNAseq
Rscript 02_Analysis/2.2.viz_mito_translation_cascade.R
```

**Required checkpoints** (generated by `02_Analysis/1.1.main_pipeline.R`):
- `03_Results/02_Analysis/checkpoints/fit_object.rds` - limma model fit with logFC coefficients
- `03_Results/02_Analysis/checkpoints/mitocarta_gsea_results.rds` - MitoCarta GSEA results

**Gene list definitions in script:**
- Calcium genes: lines 50-67 (CALCIUM_GENES_CONFIG)
- ATP Synthase genes: lines 143-151 (atp_f1, atp_fo, atp_stalk, atp_assembly)
- Mt Central Dogma: extracted from GSEA core enrichment, limited to top 25 by |logFC|

### Potential Reviewer Questions

**Q: Why not show statistical significance (asterisks) on individual genes?**
A: This figure emphasizes pattern visualization across the trajectory. Gene-level statistics are available in DE results tables. The Mt Central Dogma module as a whole is highly significant (GSEA p.adj = 2.14e-13).

**Q: Why limit Mt Central Dogma to 25 genes when GSEA found 132 core enrichment genes?**
A: Visual clarity. Showing 132 genes would make patterns difficult to discern. Top 25 by effect size captures the strongest signals. Full gene list available in GSEA results.

**Q: Why is the calcium gene list different from standard calcium signaling pathways?**
A: The list is hypothesis-driven based on prior literature linking DRP1 to calcium dysregulation. NNAT and PNPO were specifically prioritized. This is stated transparently in Methods.

**Q: Could the compensation pattern be an artifact of the color scale?**
A: No. The pattern (negative Early → positive TrajDev → near-zero Late) is evident in the numeric logFC values regardless of visualization. Module Summary Heatmap shows actual values.

### Limitations

1. **Calcium gene selection is hypothesis-driven** - not an unbiased pathway analysis
2. **Mt Central Dogma genes filtered by effect size** - may miss biologically important low-effect genes
3. **Color scale clips extreme values** - |logFC| > 0.6 displayed as maximum color
4. **No individual gene statistics shown** - refer to DE results for p-values
5. **Clustering is descriptive** - does not imply regulatory relationships

---

## Related Visualizations

| Figure | Location | Relationship |
|--------|----------|-------------|
| Synaptic Ribosome Heatmap | `../Synaptic_ribosomes/Panel_C_Expression_Heatmap.pdf` | Same style; shows synaptic translation genes |
| Complex V Deep-Dive | `../Complex_V_analysis/` | Detailed ATP synthase subunit analysis |
| GSVA Trajectories | `../Critical_period_trajectories/` | Pathway-level GSVA scores over time |
| Cross-Database Validation | `../Cross_database_validation/` | Pattern validation across GSEA databases |

---

## Complete Gene Lists

### ATP Synthase (Complex V) - 20 genes
```
ATP5F1A, ATP5F1B, ATP5F1C, ATP5F1D, ATP5F1E, ATP5MC1, ATP5MC2, ATP5MC3,
ATP5ME, ATP5MF, ATP5MG, ATP5MJ, ATP5MK, ATP5PB, ATP5PD, ATP5PF, ATP5PO,
ATPAF1, ATPAF2, TMEM70
```

### Calcium Signaling - 12 genes
```
ATP2A1, CACNA1S, CACNG3, CALB1, CALR, MYLK3, NNAT, PNPO, RYR1,
STIM1, STIM2, VDR
```

### Mt Central Dogma - 25 genes (top by |logFC|)
```
DARS2, DDX28, DNA2, ENDOG, EXOG, GUF1, MPV17L2, MRM2, MRPL1, MRPL11,
MRPL15, MRPS11, MRPS17, MRPS18C, MTERF3, MTG1, PDF, PIF1, POLQ,
PRIMPOL, PUS1, RBFA, RECQL4, TRMT61B, UNG
```

---

**Last Updated**: 2025-12-05
**Analysis Version**: v2.0 - Focused 3-module heatmap with Panel_C clustering style
**Script**: `02_Analysis/2.2.viz_mito_translation_cascade.R`

**Key Changes from v1.0**:
- Reduced from 6 modules to 3 focused modules
- Replaced incorrect GOCC_ATPASE_COMPLEX with curated ATP5* genes
- Added hierarchical clustering within modules (dendrograms)
- Narrower logFC scale (-0.6 to 0.6) for pattern visibility
- Included NNAT and PNPO in calcium signaling module
- Comprehensive documentation for peer review
