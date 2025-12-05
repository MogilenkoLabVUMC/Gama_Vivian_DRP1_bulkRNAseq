# Repository Migration Documentation

## Overview

This document describes the migration of the GVDRP1 bulk RNA-seq analysis repository from the original workstation to a new workstation environment, completed on 2025-11-19.

## Frozen Version

**Git Tag:** `v1.0`
**Commit:** `d6ec164b534c79e03c5d5f204522af8ccea02a6c`
**Description:** This tag marks the frozen state of the repository corresponding to the analysis for the manuscript undergoing review.

## Migration Context

### Why Migration Was Needed

1. **Workstation Change:** Project moved to new computational environment
2. **Path Updates:** Original absolute paths no longer valid on new system
3. **Analysis Focus Shift:** Post-migration work focuses on counts matrix analysis for reviewer comments, not raw sequencing data processing
4. **Container Standardization:** Adoption of `scdock-r-dev:v0.5.1` container across projects

### What Changed

#### 1. Development Container Configuration

**Before (Original Workstation):**
- Direct Docker image mount
- Hardcoded external data mounts:
  - `/home/karijolichlab/data/Mogilenko_lab/data/Gama_Vivian_DRP1/raw` → `/data/reads`
  - `/home/karijolichlab/data/Mogilenko_lab/data/GRCh38_hs_genome` → `/data/human_ref`
  - `/home/karijolichlab/data/Mogilenko_lab/pipeline/R_GSEA_visualisations` → GSEA module
  - `/home/karijolichlab/data/Mogilenko_lab/pipeline/bulkRNAseq_scripts` → Bash scripts

**After (New Workstation):**
- Docker Compose-based configuration
- No raw data mounts (preprocessing complete, counts matrices available)
- Updated RNAseq-toolkit mount: `/data1/users/antonz/pipeline/RNAseq-toolkit`
- Environment variable support via `.env` file

#### 2. Code Portability Improvements

**Path Standardization:**
- Replaced hardcoded `/workspaces/GVDRP1/` paths with relative paths using `here::here()`
- Affected file: `02_Analysis/Analysis_pipeline_fin.R` → `1.1.main_pipeline.R`(+ 7 path corrections)
- All analysis scripts now use workspace-relative paths

**Git Submodule Update:**
- RNAseq-toolkit submodule path updated from old location to `/data1/users/antonz/pipeline/RNAseq-toolkit`
- Submodule now properly initialized and accessible

#### 3. Documentation Enhancements

**New Files:**
- `MIGRATION.md` (this file) - Migration documentation
- `CHANGELOG.md` - Tracking changes for reviewer responses
- Updated `README.md` - Added migration notes and new setup instructions

**Backup Files:**
- `.devcontainer/devcontaine.back.json` - Preserved original container config (not tracked in git)

## What Was Removed

### Raw Data Dependencies (No Longer Needed)

The following external data sources are **not required** for continued analysis:

1. **Raw sequencing reads** (`/data/reads`)
   - FASTQ files from RNA-seq experiments
   - Only needed for preprocessing steps (STAR alignment, featureCounts)

2. **Reference genome** (`/data/human_ref`)
   - GRCh38 genome assembly and annotations
   - Only needed for read alignment

3. **Bash preprocessing scripts** (`01_Scripts/bash/`)
   - Scripts for quality control, alignment, quantification
   - Preprocessing is complete; results preserved in counts matrices

### Why This Is Safe

**All downstream analysis depends only on:**
- ✅ Counts matrix: `03_Results/01_Preprocessing/04_FeatureCounts/count_matrices_fc/sorted_counts_matrix.txt`
- ✅ Sample metadata: `03_Results/01_Preprocessing/04_FeatureCounts/count_matrices_fc/metadata.csv`
- ✅ Reference databases: `00_Data/SynGO_bulk_20231201/` (included in repo)
- ✅ Analysis scripts: `02_Analysis/` and `01_Scripts/R_scripts/`

The complete differential expression analysis, GSEA, SynGO enrichment, and visualization pipeline can be executed without raw sequencing data.

## Setup Instructions for New Environment

### Prerequisites

1. Docker installed and running
2. VS Code with Dev Containers extension
3. Access to `/data1/users/antonz/pipeline/RNAseq-toolkit/` (RNAseq-toolkit submodule location)

### Setup Steps

1. **Clone repository:**
```bash
git clone <repository-url>
cd Gama_Vivian_DRP1_bulkRNAseq
```

2. **Initialize git submodule:**
```bash
git submodule update --init --recursive
```

3. **Create environment file:**
```bash
cp test/.env.example .devcontainer/.env
# Edit .devcontainer/.env to set LOCAL_UID and LOCAL_GID
```

4. **Open in VS Code:**
```bash
code .
```

5. **Rebuild container:**
   - Command Palette (Cmd/Ctrl+Shift+P)
   - "Dev Containers: Rebuild Container"

6. **Verify setup:**
```r
# In R console within container
library(here)
here::here()  # Should show /workspaces/Gama_Vivian_DRP1_bulkRNAseq

# Check counts matrix
counts <- read.table(here::here("03_Results/01_Preprocessing/04_FeatureCounts/count_matrices_fc/sorted_counts_matrix.txt"),
                     header = TRUE, row.names = 1, sep = "\t")
dim(counts)  # Should show gene count × 25 samples
```

## Analysis Capabilities Post-Migration

### What You CAN Do

- ✅ Re-run complete differential expression analysis
- ✅ Modify statistical contrasts and thresholds
- ✅ Generate all publication figures
- ✅ Perform GSEA with different databases
- ✅ Run SynGO enrichment analysis
- ✅ Explore new gene sets (e.g., calcium signaling)
- ✅ Create custom visualizations
- ✅ Respond to reviewer comments requiring re-analysis

### What You CANNOT Do

- ❌ Re-align raw sequencing reads
- ❌ Change featureCounts parameters
- ❌ QC raw FASTQ files
- ❌ Re-run preprocessing pipeline

**Note:** If preprocessing needs to be re-run, the original raw data and reference genome would need to be mounted as in the original configuration (see `.devcontainer/devcontaine.back.json` for reference).

## Container Specifications

**Image:** `scdock-r-dev:v0.5.1`
**Source:** [scbio-docker v0.5.1](https://github.com/tony-zhelonkin/scbio-docker/tree/v0.5.1)
**Workspace:** `/workspaces/Gama_Vivian_DRP1_bulkRNAseq`
**User:** `devuser` (UID/GID mapped to host user)

## Troubleshooting

### Issue: "Cannot find RNAseq-toolkit module"

**Solution:** Ensure git submodule is initialized and path is mounted:
```bash
git submodule update --init --recursive
```

### Issue: "Permission denied" errors in container

**Solution:** Check `.env` file has correct `LOCAL_UID` and `LOCAL_GID`:
```bash
id -u  # Get your UID
id -g  # Get your GID
# Update .devcontainer/.env accordingly
```

### Issue: "Cannot find counts matrix"

**Solution:** Verify file path using `here::here()`:
```r
library(here)
list.files(here::here("03_Results/01_Preprocessing/04_FeatureCounts/count_matrices_fc/"))
```

## References

- **Original Analysis:** Git tag `v1.0`
- **Container Documentation:** See test/test/README.md for modern dev container template
- **RNAseq-toolkit:** https://github.com/tony-zhelonkin/RNAseq-toolkit (if public)

## Contact

For questions about the migration or repository setup:
- Anton Zhelonkin <anton.bioinf.md@gmail.com> / <antonz@uchicago.edu>

---
*Migration completed: 2025-11-19*
