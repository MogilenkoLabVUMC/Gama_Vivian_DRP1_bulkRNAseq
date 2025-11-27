# Critical Period Trajectories: D35 to D65 Maturation Dynamics

## Overview

This analysis visualizes how ribosome-related gene expression changes during the critical developmental period (D35 to D65) in control vs. DRP1 mutant neurons. The plots reveal a paradoxical compensatory failure: mutants increase ribosome production but fail to maintain functional translation capacity.

---

## IMPORTANT: Two Different Y-Axis Scales Used

This analysis uses **two complementary visualization approaches**:

### **Panels A-C: Expression Trajectories (Absolute Scale)**
- **Y-axis**: "Mean Expression (logFC from Ctrl D35)"
- **Shows**: How expression changes over time for ALL genotypes
- **Baseline**: Control D35 = 0 (reference timepoint)
- **Control trajectory**: Shows actual developmental changes (CAN move away from zero)

### **Panel D: Divergence Trajectories (Relative Scale)**
- **Y-axis**: "Divergence from Control (Module Score)"
- **Shows**: How mutants differ from control AFTER accounting for normal development
- **Baseline**: Control mean = 0 at EACH timepoint (by definition)
- **Control trajectory**: Not shown (would always be zero)

---

## How to Read These Plots

### Panels A-C: Visual Elements Guide

Each plot contains multiple layers showing absolute expression relative to Ctrl D35:

#### **1. Y-Axis: "Mean Expression (logFC from Ctrl D35)"**
- **What it means**: Expression level relative to Control D35 baseline
- **Zero line (dashed gray)**: Represents Ctrl D35 reference timepoint (NOT control at all times)
  - Points/lines **at zero** = same expression as Ctrl D35
  - Points/lines **above zero** = higher expression than Ctrl D35
  - Points/lines **below zero** = lower expression than Ctrl D35
- **Scale**: Log2 fold-change from baseline
  - Example: +0.5 means ~1.4× higher expression than Ctrl D35
  - Example: -0.5 means ~1.4× lower expression than Ctrl D35
- **KEY**: Control can (and does) move away from zero as it develops from D35→D65

#### **2. Individual Sample Points (Small, Semi-Transparent Dots)**
- **What they are**: Each individual biological replicate (n=3 per group/timepoint)
- **Why semi-transparent**: To show overlapping values without obscuring data
- **Why jittered**: X-positions are slightly randomized (±0.5 days) to prevent exact overlap
- **Color coding**:
  - Gray = Control samples
  - Orange = G32A mutant samples
  - Green = R403C mutant samples
- **How to interpret**:
  - **Tight clustering** = low biological variability (consistent response)
  - **Wide spread** = high biological variability (variable response)
  - **Distance from zero** = how different this sample is from control mean

#### **3. Mean Trajectory Lines (Thick, Prominent Lines)**
- **What they represent**: Average expression across n=3 replicates at each timepoint
- **How to read**:
  - **Line slope** = direction of change over time
  - **Upward slope** = expression increasing from D35 to D65
  - **Downward slope** = expression decreasing from D35 to D65
  - **Flat line** = stable expression across development
- **Why important**: Shows the overall trend independent of individual variability

#### **4. Mean Trajectory Points (Large Dots on Lines)**
- **What they are**: Mean values at each timepoint (D35 and D65)
- **Why larger/bolder**: Emphasize the actual data points that define the trajectory
- **Color coding**: Same as individual points (Gray/Orange/Green)

#### **5. Reference Line (Dashed Gray Horizontal)**
- **Position**: Always at y = 0
- **Meaning**: Represents the control group mean expression
- **Why important**: Provides visual reference for interpreting divergence magnitude

### Reading Strategy: What to Look For

**Step 1: Check the trajectory slope**
- Is the line going up, down, or staying flat from D35 to D65?
- This tells you the **direction** of developmental change

**Step 2: Check Control trajectory first**
- Where does the gray (Control) line go from D35 to D65?
- This shows the **normal developmental pattern** for this module
- Example: In Panel A, Control goes from 0 to -0.46 (ribosome biogenesis decreases normally)

**Step 3: Compare mutant trajectories to control**
- Do mutants (orange/green) follow the same developmental pattern as control (gray)?
- Or do they show a different trajectory?
- This reveals genotype-specific developmental failures

**Step 4: Check individual point spread**
- Are the small dots tightly clustered or widely spread?
- This tells you how **consistent** the response is across biological replicates

**Step 5: Compare across panels A-C**
- Do all three gene modules show coordinated changes, or opposite patterns?
- This reveals whether biological processes are coupled or uncoupled during development

### Example Interpretation

**Actual Panel A (Ribosome Biogenesis) reading**:
> "In Panel A, I see:
> - Control (gray) starts at zero (D35 reference) and drops to -0.46 at D65 → ribosome biogenesis normally **decreases** during maturation
> - G32A (orange) starts at -0.49 (below control) but only drops to -0.45 → shows less decrease than control
> - R403C (green) starts at -0.28 (below control) and drops to -0.38 → also shows less decrease than control
> - Individual dots are tightly clustered → consistent response across all 3 replicates
> - Interpretation: Control normally downregulates ribosome production during maturation, but mutants show **compensatory resistance** to this downregulation (attempting to maintain higher ribosome production)"

---

## Plot Descriptions

### Panel_A_Biogenesis_Trajectory.pdf
**Ribosome Biogenesis: Compensatory Response**

- **What it shows**:
  - Mean trajectory: Expression of ribosome biogenesis genes (n=158) relative to Ctrl D35 baseline
  - Individual points: All n=3 biological replicates per group/timepoint
  - Y-axis: Mean Expression (logFC from Ctrl D35)
- **What to observe**:
  - Control (gray): 0 → -0.46 → Normal developmental **downregulation** of ribosome production
  - G32A (orange): -0.49 → -0.45 → Shows less downregulation than control
  - R403C (green): -0.28 → -0.38 → Also shows less downregulation than control
  - Individual points are tightly clustered → consistent response across replicates
- **Key finding**: Control normally downregulates ribosome biogenesis during D35→D65; mutants resist this downregulation
- **Interpretation**: Mutants show **compensatory maintenance** of ribosome production (failing to execute normal developmental downregulation)

### Panel_B_Translation_Trajectory.pdf
**Cytoplasmic Translation: Functional Failure**

- **What it shows**:
  - Mean trajectory: Expression of cytoplasmic translation genes (n=76) relative to Ctrl D35 baseline
  - Individual points: All n=3 biological replicates per group/timepoint
  - Y-axis: Mean Expression (logFC from Ctrl D35)
- **What to observe**:
  - Control (gray): 0 → +0.26 → Normal developmental **increase** in translation capacity
  - G32A (orange): +0.22 → +0.15 → Starts similar to control but **fails to increase**
  - R403C (green): +0.17 → +0.13 → Also starts similar but **fails to increase**
  - Mutants show flat/declining trajectories while control increases
- **Key finding**: Control normally upregulates translation during maturation; mutants show **failure to upregulate**
- **Interpretation**: Despite maintaining ribosome production (Panel A), mutants cannot increase functional translation capacity during development

### Panel_C_Synaptic_Trajectory.pdf
**Synaptic Ribosomes: Critical Period Crisis**

- **What it shows**:
  - Mean trajectory: Expression of synaptic ribosome genes (n=70) relative to Ctrl D35 baseline
  - Individual points: All n=3 biological replicates per group/timepoint
  - Y-axis: Mean Expression (logFC from Ctrl D35)
- **What to observe**:
  - Control (gray): 0 → +0.18 → Normal developmental **increase** for synaptogenesis
  - G32A (orange): +0.15 → +0.07 → Starts elevated but **declines** during maturation
  - R403C (green): +0.11 → +0.04 → Also starts elevated but **declines** during maturation
  - Control trajectory is upward; mutant trajectories are downward → **opposite directions**
- **Key finding**: Control increases synaptic translation during critical period; mutants show **progressive decline**
- **Interpretation**: Mutants cannot maintain synaptic translation capacity during the D35-D65 synaptogenesis window (critical period failure)

### Panel_D_Divergence_All_Modules.pdf
**Divergence from Control During Critical Period**

- **What it shows**:
  - All three gene modules (Panels A-C) overlaid on one plot
  - Individual sample points (semi-transparent, jittered) show biological variability for all mutant samples
  - Lines/shapes show mean trajectories for G32A and R403C mutants
  - Different colors represent different gene modules (not genotypes in this panel)
- **Color coding** (by module):
  - Green = Ribosome biogenesis (n=158 genes)
  - Orange = Cytoplasmic translation (n=76 genes)
  - Blue = Synaptic ribosomes (n=70 genes)
- **Shape/line coding** (by genotype):
  - Circles, solid lines = G32A mutant
  - Triangles, dashed lines = R403C mutant
- **What to observe**:
  - Three modules show **opposite trajectories** in mutants:
    - Green (biogenesis): Upward trajectory → compensatory increase
    - Orange (translation): Downward trajectory → functional decline
    - Blue (synaptic): Downward trajectory → critical period failure
  - G32A and R403C show similar patterns → mutation-independent DRP1 phenotype
  - Individual points show tight clustering → reproducible biological response
- **Key finding**: Biological processes are uncoupled in DRP1 mutants (production vs. function)
- **Interpretation**: Compensatory upregulation of ribosome production is insufficient to prevent translation failure

### Combined_Trajectories_3panel.pdf
**All Trajectories Together**

- **What it shows**: All three panels (A, B, C) stacked vertically for direct comparison
- **Layout**:
  - Top: Panel A (Ribosome Biogenesis)
  - Middle: Panel B (Cytoplasmic Translation)
  - Bottom: Panel C (Synaptic Ribosomes)
- **Use case**: Publication-ready figure showing the complete paradoxical response pattern
- **Advantage**: See all three biological processes simultaneously with shared legend and axis scales

---

## Biological Interpretation

### The Critical Period (D35 to D65)
The D35→D65 window is when cortical neurons undergo rapid synaptogenesis and circuit maturation. This requires:
1. Increased ribosome production (biogenesis)
2. Enhanced translation capacity (cytoplasmic translation)
3. Localized synaptic translation (synaptic ribosomes)

### Control Response (Normal Development)
- **D35 baseline**: Moderate ribosome biogenesis and translation
- **D35→D65 trajectory**: Coordinated increase in ribosome production AND translation capacity
- **Result**: Successful synaptic maturation

### DRP1 Mutant Response (Paradoxical Failure)
- **D35 baseline**: Already showing ribosome stress (higher biogenesis attempt)
- **D35→D65 trajectory**:
  - Ribosome biogenesis increases MORE than control (compensation attempt)
  - Cytoplasmic translation DECLINES (functional failure)
  - Synaptic ribosomes FAIL to increase (critical period crisis)
- **Result**: Energy-translation mismatch prevents synaptic maturation → seizures

### Key Mechanistic Insight
DRP1 mutants show a **dissociation between ribosome production and translation function**. They attempt to compensate by making more ribosomes, but the underlying mitochondrial dysfunction prevents these ribosomes from functioning effectively. This creates an energy crisis during the critical period when synaptic demands are highest.

---

## Methodological Rationale: Why Two Different Y-Axes?

### The Problem: Developmental vs. Genotype Effects

When studying **developmental time-series with genetic perturbations**, expression changes have two components:
1. **Developmental component**: Normal changes during maturation (happens in control)
2. **Genotype component**: Mutation-specific effects (different from control)

### Two Complementary Approaches

#### **Approach 1: Absolute Expression (Panels A-C)**
**Formula**: Expression = sample_mean - Ctrl_D35_mean

**What it shows**:
- Total expression change from a fixed baseline (Ctrl D35)
- Includes both developmental AND genotype effects
- Control trajectory shows normal development

**When to use**:
- To see how each genotype changes over time
- To compare developmental trajectories across genotypes
- To identify if control itself changes during development

**Standard in**: Developmental transcriptomics, time-series RNA-seq

#### **Approach 2: Divergence from Control (Panel D)**
**Formula**: Divergence = sample_mean - control_mean_at_same_timepoint

**What it shows**:
- Mutation-specific effects AFTER removing developmental baseline
- Isolates genotype effects independent of normal development
- Control would always be zero (excluded from plot)

**When to use**:
- To compare mutation effects across different developmental stages
- To identify genotype-specific dysregulation
- To ask "how different is the mutant from control at each age?"

**Standard in**: Treatment/genotype comparison studies, clinical transcriptomics

### Why Use Both?

**Panels A-C answer**: "Does control change during development? Do mutants follow this pattern?"
- Reveals: Control ribosome biogenesis **decreases** -0.46 logFC
- Reveals: Mutants resist this decrease (compensatory response)

**Panel D answers**: "How do mutants diverge from control at each stage?"
- Reveals: At D65, mutants show +0.02 to +0.09 divergence in biogenesis (relative upregulation)
- Cleaner comparison of mutation effects across timepoints

### Is This Standard Practice?

**YES** - This hybrid approach is common in developmental perturbation studies:

**Literature examples**:
1. **Developmental time-series + perturbation**: Use absolute expression to show developmental dynamics
2. **Cross-stage genotype comparison**: Use divergence/difference scores to isolate treatment effects
3. **Clinical studies**: Compare patient samples to age-matched controls (divergence approach)

**Key reference approaches**:
- DESeq2 tutorial recommends both "time effects" and "genotype:time interaction" contrasts
- Bioconductor workflows show both absolute changes and treatment-specific deviations
- Standard practice in longitudinal RNA-seq with interventions

### Mathematical Formulas

**Panels A-C (Expression)**:
```
Expression_ij = mean(logCPM_ij) - mean(logCPM_Ctrl_D35)

where:
  i = genotype (Ctrl, G32A, R403C)
  j = timepoint (D35, D65)
  Ctrl D35 is reference = 0
```

**Panel D (Divergence)**:
```
Divergence_ij = mean(logCPM_mutant_ij) - mean(logCPM_Ctrl_j)

where:
  i = mutant genotype (G32A or R403C only)
  j = timepoint (D35 or D65)
  Control mean at same timepoint j is subtracted
  Control divergence would be 0 by definition (not shown)
```

---

## Methods Note

### Data Processing
- **Module scores**: Mean log2 counts per million (logCPM) across all genes in each functional module
- **Individual sample points**: Each biological replicate shown as semi-transparent dot (n=3 per group/timepoint)
- **X-axis jitter**: Sample points randomized by ±0.5 days to prevent overplotting
- **Trajectory lines**: Connect mean values at D35 and D65 for each genotype
- **Gene sets**:
  - Ribosome biogenesis (n=158 genes): From GO:BP enrichment in Maturation_G32A_specific comparison
  - Cytoplasmic translation (n=76 genes): From GO:BP enrichment in Maturation_G32A_specific comparison
  - Synaptic ribosomes (n=70 genes): From SynGO synaptic annotations

### Statistical Design
- **Normalization**: edgeR TMM normalization + log2 transformation
- **Modeling**: limma-voom linear modeling with genotype × timepoint design
- **Contrasts**: Time-course effects (D65 vs D35) within each genotype
- **Module scores**: Simple mean of normalized expression (NOT GSVA, NOT enrichment scores)
  - Justified because: Gene modules are functionally coherent (from GO/SynGO enrichment)
  - Mean logCPM is biologically interpretable as "average pathway activity"

### Quality Metrics
- **Biological replicates**: n=3 per group/timepoint (18 samples total)
- **Data transparency**: All individual replicates shown on plots
- **Reproducibility**: Tight clustering of replicates indicates consistent biological response

## Related Visualizations

See also:
- [../Publication_Figures/](../Publication_Figures/README.md) - Main manuscript figures
- [../Ribosome_paradox/](../Ribosome_paradox/README.md) - Core ribosome downregulation finding
- [../Mito_translation_cascade/](../Mito_translation_cascade/README.md) - Mechanistic cascade visualization
- [../Synaptic_ribosomes/](../Synaptic_ribosomes/README.md) - Synaptic translation deep-dive
- [../Cross_database_validation/](../Cross_database_validation/README.md) - Pattern validation framework
- [../Pattern_Summary_Normalized/](../Pattern_Summary_Normalized/README.md) - Pattern classification summary

---

## Generated by
`02_Analysis/viz_critical_period_trajectories.R`

**Dependencies**: Requires checkpoint files from main differential expression analysis
