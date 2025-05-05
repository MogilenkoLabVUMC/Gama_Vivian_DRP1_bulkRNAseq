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


3. Quantification


### 2. Analysis Pipeline

## Results


## Citation


## Contact



# Future directions
