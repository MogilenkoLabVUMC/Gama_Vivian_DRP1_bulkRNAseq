# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- MIGRATION.md documentation describing workstation migration process
- CHANGELOG.md for tracking changes (this file)
- .env file support for container configuration
- Docker Compose configuration for dev container

### Changed
- **BREAKING:** Updated all hardcoded `/workspaces/GVDRP1/` paths to relative paths using `here::here()`
  - Affected file: `02_Analysis/Analysis_pipeline_fin.R` (7 path references)
- Updated RNAseq-toolkit git submodule path to `/data1/users/antonz/pipeline/RNAseq-toolkit/`
- Migrated from direct Docker image to Docker Compose-based dev container
- Updated docker-compose.yml to remove duplicate volume mounts

### Removed
- External mount dependencies for raw sequencing data (no longer needed for analysis)
  - Raw FASTQ files mount
  - Reference genome mount
  - Bash preprocessing scripts mount

### Fixed
- Corrected RNAseq-toolkit submodule initialization
- Standardized workspace paths for portability across workstations

### Documentation
- Updated README.md with migration notes and new setup instructions
- Added troubleshooting section to MIGRATION.md

## [1.0] - 2025-06-18

### Summary
Frozen version corresponding to manuscript submission/publication. Tagged at commit `d6ec164`.

### Included Features
- Complete bulk RNA-seq analysis pipeline for DRP1 mutations (G32A, R403C)
- Differential expression analysis using edgeR/limma voom
- Gene set enrichment analysis (GSEA) with multiple databases
- SynGO synaptic ontology enrichment
- Calcium signaling gene analysis
- Comprehensive visualization suite:
  - Volcano plots (standard and vertical comparison)
  - MA plots
  - Heatmaps (correlation, gene expression)
  - MDS plots
  - GSEA enrichment plots

### Analysis Details
- **Samples:** 26 (3 genotypes × 2 timepoints × biological replicates)
- **Contrasts:** 9 statistical comparisons
  - 4 pairwise (mutation vs control at each timepoint)
  - 3 maturation (time effects within genotype)
  - 2 interaction (mutation-specific maturation changes)
- **Normalization:** TMM (Trimmed Mean of M-values)
- **Statistical method:** limma-voom with empirical Bayes moderation

### Data Files
- Counts matrix: `03_Results/01_Preprocessing/04_FeatureCounts/count_matrices_fc/sorted_counts_matrix.txt`
- Metadata: `03_Results/01_Preprocessing/04_FeatureCounts/count_matrices_fc/metadata.csv`
- SynGO database: `00_Data/SynGO_bulk_20231201/`

---

## Future Sections (Template for Reviewer Comments)

Use the sections below to track changes made in response to reviewer feedback:

## [1.1] - YYYY-MM-DD (In Progress)

### Added - Reviewer Response
<!-- Example:
- New analysis of mitochondrial genes per Reviewer #2 comment
- Additional GSEA database (WikiPathways) per Reviewer #3 suggestion
-->

### Changed - Reviewer Response
<!-- Example:
- Updated FDR threshold from 0.05 to 0.1 per Reviewer #1
- Re-ran GSEA with updated gene sets
- Modified volcano plot color scheme for better visibility
-->

### Figures - Reviewer Response
<!-- Example:
- Regenerated Figure 3 with higher resolution
- Added supplementary figure showing cell type markers
- Updated Figure 5 to include NNAT (Neuronatin) as requested
-->

---

## Notes for Maintainers

### Version Numbering
- **Major version (X.0):** Substantial changes to analysis pipeline or methodology
- **Minor version (1.X):** New analyses, figure updates, or reviewer-requested changes
- **Patch version (1.0.X):** Bug fixes, documentation updates, minor corrections

### Commit Message Convention
When making changes for reviewer comments, use descriptive commit messages:
```
[Reviewer 2] Add mitochondrial gene enrichment analysis

- Perform GSEA on mitochondrial GO terms
- Create supplementary figure S5
- Update methods section accordingly

Addresses Reviewer #2, Comment #3
```

### Tracking Reviewer Comments
For each reviewer comment being addressed:
1. Add entry to CHANGELOG.md under appropriate section
2. Reference the specific reviewer and comment number
3. Describe what was changed and why
4. Note any new files or figures generated
5. Commit with descriptive message linking to reviewer comment

### Reproducibility Notes
- Always document the exact version of code used for each submission
- Tag git commits corresponding to manuscript revisions (e.g., `v1.1-revision1`)
- Include R session info in analysis outputs when reporting new results
- Document any changes to statistical thresholds or parameters
