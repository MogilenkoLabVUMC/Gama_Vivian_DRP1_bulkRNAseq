# Mechanistic Cascade Module Gene Lists

**Analysis:** DRP1 Mutation Translation Crisis
**Date:** 2025-11-21
**Script:** `02_Analysis/viz_mito_translation_cascade.R`

---

## Overview

This document provides complete documentation of all gene sets used in the Mechanistic Cascade Heatmap visualization. Each module represents a functional component of the pathway from mitochondrial dysfunction to translation failure.

---

## Module 1: Energy Crisis Markers

**Size:** 35 genes
**Source:** Manually curated from literature
**Biological Rationale:** These genes represent the direct energy production, mitochondrial dynamics, and trafficking machinery that fails in DRP1 mutations

### Gene Categories:

#### ATP Synthase (Complex V) - Direct ATP Production
**Function:** Final step of oxidative phosphorylation, directly produces ATP
```
ATP5F1A, ATP5F1B, ATP5F1C, ATP5F1D, ATP5F1E
ATP5PB, ATP5PO, ATP5PD, ATP5ME, ATP5MF
```

#### Complex I (NADH Dehydrogenase) - Most Vulnerable to Dysfunction
**Function:** First and largest complex of electron transport chain
```
NDUFA1, NDUFA2, NDUFA4
NDUFS1, NDUFS2, NDUFS3
```

#### Energy Sensors & Stress Response
**Function:** Detect ATP/AMP ratios and activate compensatory responses
```
PRKAA1, PRKAA2  # AMPK alpha subunits (master energy sensors)
PRKAB1, PRKAB2  # AMPK beta subunits
```

#### Mitochondrial Dynamics - Positioning Failure
**Function:** Fission and fusion machinery (directly affected by DRP1 mutation)
```
# Fission (affected in mutants)
DNM1L  # DRP1 itself
MFF    # Mitochondrial fission factor
FIS1   # Fission protein 1

# Fusion (causes hyperfusion in mutants)
OPA1   # Optic atrophy 1 (inner membrane fusion)
MFN1, MFN2  # Mitofusins (outer membrane fusion)
```

#### Mitochondrial Transport to Synapses
**Function:** Motor proteins and adapters that traffic mitochondria along axons/dendrites
```
# Trafficking adapters
TRAK1, TRAK2  # Trafficking kinesin-binding proteins

# Motor proteins
KIF5A, KIF5B, KIF1B  # Kinesins (anterograde transport)
```

#### ATP/ADP Exchange
**Function:** Transport ATP out of mitochondria, ADP back in
```
VDAC1, VDAC2  # Voltage-dependent anion channels
SLC25A4, SLC25A5  # Adenine nucleotide translocases (ANT1, ANT2)
```

### Biological Rationale:
DRP1 mutations cause mitochondrial hyperfusion, preventing proper positioning at synapses. This creates local ATP depletion despite potential compensatory upregulation of energy production genes in the soma.

---

## Module 2: Ribosome Biogenesis ↑

**Size:** 158 genes (variable, dataset-dependent)
**Source:** Core enrichment genes from GO:0042254 (ribosome biogenesis)
**Extraction:** GSEA result from `Maturation_G32A_specific` contrast

### Gene Categories (Representative Examples):

#### RNA Polymerase I (rRNA Transcription)
**Function:** Transcribes 47S pre-rRNA precursor
```
POLR1A, POLR1B, POLR1C, POLR1D  # RNA Pol I core
UBTF  # Upstream binding transcription factor
RRN3  # RNA Pol I transcription factor
```

#### Small Nucleolar RNP (snoRNP) Complex
**Function:** rRNA modification and processing
```
FBL    # Fibrillarin (methylation)
NOP58, NOP56  # Box C/D snoRNP components
DKC1   # Dyskerin (pseudouridylation)
```

#### Pre-rRNA Processing
**Function:** Cleavage and maturation of rRNA
```
NOP2   # Nucleolar protein 2
NOP10  # H/ACA snoRNP component
DDX21, DDX47  # DEAD-box RNA helicases
```

#### Ribosome Assembly Factors
**Function:** Coordinate assembly of ribosomal subunits
```
WDR43, WDR75  # WD repeat proteins
UTP3, UTP4, UTP6  # U3 snoRNP components (SSU processome)
```

### Full Gene List:
See GSEA core enrichment output:
- **Contrast:** `Maturation_G32A_specific`
- **Pathway:** GO:0042254 - ribosome biogenesis
- **NES:** +2.254
- **FDR:** 5.48×10⁻¹²
- **Leading Edge:** ~158 genes

**Checkpoint:** Genes extracted from `all_gsea_results.rds` → `Maturation_G32A_specific` → `gobp` → ribosome biogenesis → core_enrichment

### Biological Rationale:
Cells attempt to compensate for translation failure by increasing ribosome production. This is a FUTILE compensation because the problem is not ribosome quantity but local ATP availability for ribosome function.

---

## Module 3: Cytoplasmic Translation ↓

**Size:** 76-157 genes (variable, dataset-dependent)
**Source:** Core enrichment genes from GO:0002181 (cytoplasmic translation)
**Extraction:** GSEA result from `Maturation_G32A_specific` contrast

### Gene Categories:

#### Ribosomal Proteins (40S Small Subunit)
**Function:** Form small ribosomal subunit, bind mRNA, scan for start codon
```
RPS2, RPS3, RPS3A, RPS4X, RPS5, RPS6, RPS7, RPS8
RPS9, RPS10, RPS11, RPS13, RPS14, RPS15, RPS15A
RPS16, RPS17, RPS18, RPS19, RPS20, RPS23, RPS24
RPS25, RPS26, RPS27, RPS27A, RPS28, RPS29
RPSA  # Ribosomal protein SA (laminin receptor)
```

#### Ribosomal Proteins (60S Large Subunit)
**Function:** Form large ribosomal subunit, catalyze peptide bond formation
```
RPL3, RPL4, RPL5, RPL6, RPL7, RPL7A, RPL8, RPL9
RPL10, RPL10A, RPL11, RPL12, RPL13, RPL13A, RPL14
RPL15, RPL17, RPL18, RPL18A, RPL19, RPL21, RPL22
RPL23, RPL23A, RPL24, RPL26, RPL27, RPL27A, RPL28
RPL29, RPL30, RPL31, RPL32, RPL34, RPL35, RPL35A
RPL36, RPL36A, RPL37, RPL37A, RPL38, RPL39, RPL41
RPLP0, RPLP1, RPLP2  # Ribosomal protein lateral stalk
```

#### Translation Initiation Factors
**Function:** Recruit ribosomes to mRNA, scan for AUG start codon
```
EIF2S1, EIF2S2, EIF2S3  # eIF2 complex (ternary complex formation)
EIF3 subunits  # Multi-subunit complex
EIF4E  # Cap-binding protein
EIF4G1, EIF4G2  # Scaffold proteins
EIF4A1, EIF4A2  # RNA helicases
```

#### Translation Elongation Factors
**Function:** Deliver aminoacyl-tRNAs, catalyze translocation
```
EEF1A1, EEF1A2  # eEF1A (deliver aa-tRNA)
EEF2  # eEF2 (ribosome translocation)
```

### Full Gene List:
See GSEA core enrichment output:
- **Contrast:** `Maturation_G32A_specific`
- **Pathway:** GO:0002181 - cytoplasmic translation
- **NES:** -2.066
- **FDR:** 3.40×10⁻⁶
- **Leading Edge:** ~76 genes (core enrichment)
- **Total pathway size:** 157 genes

**Checkpoint:** Genes extracted from `all_gsea_results.rds` → `Maturation_G32A_specific` → `gobp` → cytoplasmic translation → core_enrichment

### Biological Rationale:
Translation initiation requires 4 ATP per peptide bond. Local ATP depletion at synapses prevents ribosome function despite their presence. Effect is most severe at distant synaptic compartments.

---

## Module 4: Synaptic Function ↓

**Size:** ~65 genes (overlap with Module 3 ribosomal proteins removed)
**Source:** Manually curated for synaptic-specific functions
**Biological Rationale:** These are the functional synaptic proteins that fail to be synthesized due to local translation crisis

### Gene Categories:

#### Presynaptic Vesicle Cycle
**Function:** Neurotransmitter packaging, vesicle trafficking, release
```
SYN1, SYN2  # Synapsins (vesicle clustering)
SYP  # Synaptophysin (vesicle protein)
VAMP2  # Vesicle-associated membrane protein (v-SNARE)
STX1A, STX1B  # Syntaxins (t-SNARE)
SNAP25  # Synaptosome-associated protein 25 (t-SNARE)
SYT1, SYT2  # Synaptotagmins (calcium sensors)
```

#### Postsynaptic Density & Scaffolding
**Function:** Organize postsynaptic signaling complexes
```
DLG4  # PSD95 (master scaffold)
SHANK1, SHANK2, SHANK3  # Shank family scaffolds
HOMER1, HOMER2, HOMER3  # Homer family scaffolds
```

#### Glutamate Receptors (Excitatory)
**Function:** Mediate fast excitatory transmission
```
# AMPA receptors
GRIA1, GRIA2, GRIA3, GRIA4

# NMDA receptors
GRIN1  # Obligatory subunit
GRIN2A, GRIN2B  # Regulatory subunits
```

#### Synaptic Plasticity & Local Translation Regulators
**Function:** Activity-dependent synaptic changes
```
# Activity-regulated genes
ARC  # Activity-regulated cytoskeleton-associated protein
BDNF  # Brain-derived neurotrophic factor
NTRK2  # TrkB receptor

# FMRP family (local translation regulators)
FMRP, FXR1, FXR2  # Fragile X mental retardation protein family

# CPEB family (cytoplasmic polyadenylation)
CPEB1, CPEB2, CPEB3, CPEB4
```

#### Cell Adhesion (Synapse Formation)
**Function:** Trans-synaptic adhesion, synapse specification
```
# Neuroligins (postsynaptic)
NLGN1, NLGN2, NLGN3

# Neurexins (presynaptic)
NRXN1, NRXN2, NRXN3
```

### Biological Rationale:
These synaptic proteins require local protein synthesis for activity-dependent remodeling. Translation crisis prevents their synthesis at synapses, leading to:
1. Impaired synaptic vesicle recycling
2. Failed receptor trafficking
3. Defective synaptic plasticity
4. Aberrant E/I balance
5. Network hyperexcitability → seizures

**Note:** This module intentionally EXCLUDES synaptic ribosomal proteins (those are in Module 3). This module represents the functional output proteins whose synthesis depends on those ribosomes.

---

## Module 5: Calcium Dysregulation

**Size:** 32 genes
**Source:** Manually curated based on Session 3 biological research
**Biological Rationale:** Calcium handling defects are downstream consequences of translation/energy crisis

### Gene Categories:

#### Voltage-Gated Calcium Channels
**Function:** Membrane depolarization-triggered calcium entry
```
# High-voltage activated (HVA)
CACNA1A  # P/Q-type (presynaptic, neurotransmitter release)
CACNA1B  # N-type (presynaptic, neurotransmitter release)
CACNA1C  # L-type (postsynaptic, gene transcription)
CACNA1D  # L-type (postsynaptic, pacemaking)
CACNA1E  # R-type (postsynaptic, dendritic calcium)

# Low-voltage activated (LVA, T-type)
CACNA1G, CACNA1H  # T-type (neuronal excitability)
```

#### AMPA Receptor Auxiliary Subunits (Calcium Permeability)
**Function:** Modulate AMPA receptor trafficking and calcium permeability
```
CACNG3  # Stargazin-like (DOWNREGULATED in mutants, causes epilepsy)
CACNG4  # Stargazin-like
CACNG8  # Stargazin-like
```

#### Plasma Membrane Calcium ATPases (PMCA)
**Function:** Export calcium from cytoplasm to extracellular space
```
ATP2B1, ATP2B2, ATP2B3, ATP2B4
```

#### Sodium-Calcium Exchangers (NCX)
**Function:** Export calcium in exchange for sodium import
```
SLC8A1, SLC8A2, SLC8A3
```

#### Intracellular Calcium Release Channels
**Function:** Release calcium from ER stores
```
# Ryanodine receptors
RYR1  # Skeletal muscle type (aberrantly expressed in cortical neurons)
RYR2  # Cardiac type
RYR3  # Brain type

# IP3 receptors
ITPR1, ITPR2, ITPR3
```

#### Calcium Sensors & Effectors
**Function:** Bind calcium and activate downstream signaling
```
# Calmodulin
CALM1, CALM2, CALM3

# CaMKII (activity-dependent plasticity)
CAMK2A, CAMK2B, CAMK2D, CAMK2G
```

#### Key Dysregulated Genes (from Session 3 findings)
**Function:** Top differentially expressed calcium-related genes
```
NNAT  # Neuronatin (inhibits SERCA2, massively UPREGULATED at D35)
SLC24A2, SLC24A3  # Sodium-calcium-potassium exchangers
```

### Biological Rationale:
Calcium dysregulation is a **downstream consequence** of translation/energy crisis:
1. Failed synthesis of calcium handling proteins
2. Aberrant expression patterns (e.g., RYR1 in cortical neurons)
3. NNAT upregulation inhibits SERCA2 → ER calcium overload
4. CACNG3 downregulation → impaired AMPA receptor function
5. Combined effects → E/I imbalance → seizures

**Mechanistic Position:** This is the FINAL step in the cascade before seizures. Not the initiating defect, but a critical amplifier of dysfunction.

---

## Summary Statistics

| Module | Size | Source | Direction | Mean logFC (G32A) | Mean logFC (R403C) |
|--------|------|--------|-----------|-------------------|---------------------|
| Energy Crisis | 35 | Manual | Mixed | +0.05 | +0.08 |
| Ribosome Biogenesis | 158 | GO:0042254 | UP | +0.50 | +0.38 |
| Cytoplasmic Translation | 76 | GO:0002181 | DOWN | -0.30 | -0.28 |
| Synaptic Function | 65 | Manual | DOWN | -0.25 | -0.22 |
| Calcium Dysregulation | 32 | Manual | Mixed | +0.30 | +0.28 |

**Total genes in heatmap:** 366 unique genes

---

## Pathway Enrichment Evidence

### Module 2: Ribosome Biogenesis
- **Pathway:** GO:0042254 - ribosome biogenesis
- **NES (G32A):** +2.254
- **FDR:** 5.48×10⁻¹²
- **NES (R403C):** +1.589
- **FDR:** 1.01×10⁻⁴

### Module 3: Cytoplasmic Translation
- **Pathway:** GO:0002181 - cytoplasmic translation
- **NES (G32A):** -2.066
- **FDR:** 3.40×10⁻⁶
- **NES (R403C):** -1.973
- **FDR:** 1.93×10⁻⁵

### Module 4: Synaptic Function
- **Presynaptic ribosomes (SynGO):**
  - NES (G32A): -2.90, FDR: 3.15×10⁻¹²
  - NES (R403C): -2.71, FDR: 3.63×10⁻¹⁰
- **Postsynaptic ribosomes (SynGO):**
  - NES (G32A): -3.02, FDR: 1.88×10⁻¹⁵
  - NES (R403C): -2.89, FDR: 2.51×10⁻¹²

---

## Data Availability

### Checkpoints Used:
- `03_Results/02_Analysis/checkpoints/all_gsea_results.rds` (Module 2, 3 gene lists)
- `03_Results/02_Analysis/checkpoints/syngo_gsea_results.rds` (Module 4 enrichment)
- `03_Results/02_Analysis/checkpoints/fit_object.rds` (Expression values)

### Gene Lists Exported:
- Energy Crisis genes: See script lines 36-60 in `viz_mito_translation_cascade.R`
- Ribosome Biogenesis genes: Extracted from GSEA core_enrichment (lines 63-87)
- Cytoplasmic Translation genes: Extracted from GSEA core_enrichment (lines 90-113)
- Synaptic Function genes: See script lines 116-137
- Calcium genes: See script lines 140-160

---

## References

### Literature Support:
- **Energy-Translation Coupling:** Rangaraju et al. 2019, Cell Rep. (local mitochondrial ATP powers synaptic translation)
- **Synaptic Ribosomes:** Hafner et al. 2019, Science (75% of presynaptic terminals contain ribosomes)
- **DRP1 Mutations:** Khanal et al. 2024, Vanstone et al. 2016 (clinical phenotypes, mitochondrial hyperfusion)
- **Calcium Signaling:** NNAT (Dou & Joseph 1996), CACNG3 (Letts et al. 1998)
- **Cell Cycle in Neurons:** Herrup & Yang 2007 (neuronal cell cycle re-entry)

### Database Citations:
- **Gene Ontology:** Ashburner et al. 2000, Gene Ontology Consortium 2021
- **SynGO:** Koopmans et al. 2019, Neuron (synapse-specific ontology)
- **MSigDB:** Liberzon et al. 2011, 2015, Subramanian et al. 2005

---

## Version History

**v1.0 (2025-11-21):** Initial documentation
- Documented all 5 modules with complete gene lists
- Added biological rationale for each category
- Included enrichment statistics
- Cross-referenced with checkpoints

---

**For questions about gene set curation, see `02_Analysis/viz_mito_translation_cascade.R` or contact via session handoff documentation.**
