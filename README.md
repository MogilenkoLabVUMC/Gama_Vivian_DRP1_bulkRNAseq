# Vivian Gama: DRP1 BulkRNAseq Experiment

## Overview
This repository contains the code and analysis pipeline used to investigate the effects of DRP1 mutations on neuronal maturation.

## My comments 
I read the paper and reviewer`s comments. I feel like most of them we could address with our standard bulkRNAseq pipeline with standard DE and GSEA against a couple of databases in MSigDbr, such as Hallmark pathway, GO:BP, GO:CC, Reactome, KEGG etc. -> this mostly covers the reviewer’s request for “mechanistic insights”

But I think we can squeeze much more biology out of the data than the original analysis did

Synaptic Gene Ontologies. I found this SynGO database, one of the reviewers mentioned. The ontologies are downloadable. I canintegrate this database into our pipeline and  test DEG lists against this databases. That should highlight presynaptic vs postsynaptic signatures more cleanly than generic GO.

Transcription factor programs. I`m thinking if any transcription factors are involved. I hope I could try and run some inference tool (there are couple like decoupleR) that predicts Transcription factor activity from bulk RNAseq expression data. This might give a mechanistic bridge from elongated mitochondria -> TF programmes -> synaptic genes.

Specific genes of interest. Concerning AMPARs and any other genes of interest specifically, I could just make a per-sample heatmap of LogFC or Relative expression for the main AMPAR subunit genes and its auxiliary genes. If you could give me a list of genes that you always want to have an eye on, I will keep them on all of the volcanos explicitly labeled. 

Gene Networks. You have 26 samples. While splitting it into subgroups per-mutant per-day would yield too little samples, but we could run weighted-gene correlation network (WGCNA) analysis on all the samples. This might uncover modules of co-regulated genes in your data. We could further test enrichment of these modules in SynGO or any other database of interest. We could also correlate the modules with some phenotype data, that you might have, say mitocondrial length, Ca²⁺ peak, synapse density or whatever, or make a correlation heatmap matrix of "co-regulated network module" X "phenotypic trait". 

Also I don`t think you exploit your data enough in terms of statistical modelling. We could build a linear model to interrogate it with our biological questions, that would span beyond just the standard pair-wise comparisons. 

For example, 
We may of course ask ourselves, what`s the difference between each kind of mutant vs control at each time-point.
Mathematically, this would encode into this 4 separate model questions:
  G32A_vs_Ctrl_D35  =  D35_G32A  - D35_Control,
  R403C_vs_Ctrl_D35 =  D35_R403C - D35_Control,
  G32A_vs_Ctrl_D65  =  D65_G32A  - D65_Control,
  R403C_vs_Ctrl_D65 =  D65_R403C - D65_Control,
Each such contrast (G32A_vs_Ctrl_D35 etc.) gives its own list of DEGs, that we could test in a pathway analysis of any kind (GO, GSEA doesn`t matter).

And of course we might ask, what`s the maturation effect inside each genotype
  Time_Ctrl   =  D65_Control - D35_Control, i.e. What genes change during normal maturation?
  Time_G32A   =  D65_G32A    - D35_G32A, i.e. What genes change during G32A mutant maturation?
  Time_R403C  =  D65_R403C   - D35_R403C, i.e. What genes change during R403C mutant maturation?
Again, each question gives us its separate DEG list, which we could dig deeper into pathways further 

But we also might go for more complex interaction questions
"Difference-in-difference" questions
Int_G32A = (D65_G32A - D65_Control) - (D35_G32A - D35_Control), i.e. Does G32A’s maturation trajectory differ from control trajectory?
Int_R403C = (D65_R403C - D65_Control) - (D35_R403C - D35_Control), i.e. Does R403C`s maturation trajectory differ from control?
Positive DEGs: are the ones differentially up-regulated over time only in the mutant 

These type of questions can be modelled too. They isolate genes whose time-course is selectively altered by the mutation, not just genes that differ at one time-point

Combination of such Q&A would give us some insight. For example, 

result pattern
Time_Ctrl
Time_G32A
Int_G32A
interpretation
gene doubles in Ctrl & doubles in G32A
+1
+1
0
same maturation → not highlighted
gene doubles in Ctrl but stays flat in G32A
+1
0
–1
G32A fails to up-regulate this gene
gene flat in Ctrl but doubles in G32A
0
+1
+1
gene activation is mutant-specific


We could go and model other questions like this. For example, we might be interested in common DEGs mis-expressed similarly in both G32A and R403C (all DRP1 mutants) vs control. We also might include phenotypic traits in the model. 

## Background
Mitochondrial fission, mediated by DRP1, is vital for the rapid metabolic shifts that occur during cortical development. De-novo DNM1L loss-of-function mutations (G32A in the GTP-ase domain, R403C in the stalk domain) are emerging in paediatric neuro-developmental cohorts but their domain-specific impact on neuronal maturation is unclear.

## Hypothesis 
Domain-specific DRP1 mutations differentially disturb mitochondrial dynamics, which in turn impairs synaptic development, calcium handling and the transcriptional programme of maturing cortical neurons.

## Experimental design
Generate iPSC lines harbouring heterozygous G32A or R403C mutations → dual-SMAD differentiation to glutamatergic cortical cultures.</li><li>Phenotype mitochondria (SIM, live spinning-disk), synapses (SYP/PSD-95 SIM), and Ca²⁺ dynamics (Fluo-4 imaging) at 35, 65, 100 DIV

## Experimental findings
Both mutants retain the expected cortical-neuron marker profile but exhibit persistently hyper-elongated axonal mitochondria. Size-dependent, mutation-specific motility phenotypes (e.g. faster small mitochondria in R403C axons). Functional read-outs confirm diminished pre/post-synaptic marker volume and markedly exaggerated glutamate-evoked Ca²⁺ transients in both mutants

### Novelty
**Domain-resolved design** – analysing two distinct DRP1 mutations enables structure–function inference and partially explains heterogeneous patient presentations.  

**Multi-modal validation** – transcriptomic, imaging and functional assays converge on impaired synapse/Ca²⁺ regulation, lending biological credibility to RNA-seq findings.  


### Keep an eye out on 
NEURONATIN, CACNG3, CACNA1C, CACNA1S, ATP2A1, RYR1, MYLK3, CASR, VDR, STIM1, STIM2, ORAI1, CALB1 and CALR

Genes associated with calcium signaling

### Collaborator questions
1. The key question is: what are the effects of the G32A and R403C mutations in the transcriptional landscape of neurons?
2. We have 2 time points (35 DIV and 65 DIV), thus we were also wondering whether the mutations had any effects on maturation (compared to controls)


### General Approach Notes
The study asks **how two domain-specific DRP1 mutations (G32A = GTPase, R403C = stalk) disturb cortical-neuron maturation**.  
Phenotypes observed:

- hyper-elongated, sluggishly trafficked mitochondria
- widespread DEGs at 35 DIV (ion-channel, synapse, ROS) that largely normalise by 65 DIV
- smaller/rarer synapses and exaggerated glutamate-evoked Ca²⁺ transients
    

The reviewers want a **mechanistic bridge** from _mitochondrial shape → nuclear transcriptome → synaptic outcome_.




## Repository Structure

Below is a high-level overview of the repository:

* `00_Data/`

* `01_Scripts/`

* `02_Analysis/`

* `03_Results/`

* `README.md` (this file)



## Data Processing Workflow

### 1. Data Preprocessing



#### Preparing index

GRCh38.p14 reference genome (hg38, GCF_000001405.40, release date Feb 3, 2022) for read mapping. 'GCF_000001405.40_GRCh38.p14_genomic.gff.gz' was used as the annotation. To run our custom QC scripts we pre-processed the annotation in the following manner. 

**1. Prep refFlat from gtf**
```
./GTFtoRefFlat.sh -i GCF_000001405.40_GRCh38.p14_genomic.gtf.gz -o GCF_000001405.40_GRCh38.p14_genomic.refflat
```

**2. Prep ribosomal intervals**
```
./getRibosomalIntervals_from_gtf.sh -i GCF_000001405.40_GRCh38.p14_genomic.gtf.gz -r GCF_000001405.40_GRCh38.p14_genomic.fna.gz -o GCF_000001405.40_GRCh38.p14_genomic.ribosomal_intervals
```

**3. Get BED12 from GTF**
```
./GTFtoBED12.sh -i GCF_000001405.40_GRCh38.p14_genomic.gtf.gz -o GCF_000001405.40_GRCh38.p14_genomic.bed12
```

**4. (Optional) Change to ReSeQC compatible bed**
```
./BEDtoRefSeqBED_human.sh -i GCF_000001405.40_GRCh38.p14_genomic.bed12 -o GCF_000001405.40_GRCh38.p14_genomic.reseqc.bed12
```


```
STAR --runMode genomeGenerate \
     --runThreadN 8 \
     --genomeDir ../ref_index100 \
     --genomeFastaFiles ./GCF_000001405.40_GRCh38.p14_genomic.fna \
     --sjdbGTFfile ./GCF_000001405.40_GRCh38.p14_genomic.gtf \
     --sjdbOverhang 100  
```


3. Quantification


### 2. Analysis Pipeline

## Results


## Citation


## Contact



# Future directions
