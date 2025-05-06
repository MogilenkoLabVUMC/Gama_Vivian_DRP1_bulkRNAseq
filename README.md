# Vivian Gama: DRP1 BulkRNAseq Experiment

## Overview
This repository contains the code and analysis pipeline used to investigate the effects of DRP1 mutations on neuronal maturation.

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
