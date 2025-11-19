# Quick Setup Guide

This guide provides step-by-step instructions to get started with the migrated GVDRP1 bulk RNA-seq analysis repository.

## Prerequisites

- Docker installed and running
- VS Code with Dev Containers extension
- SSH keys configured for GitHub access (to clone the RNAseq-toolkit submodule)

## Setup Steps

### 1. Initialize Git Submodule

The RNAseq-toolkit is required for GSEA helper functions. Initialize it with:

```bash
git submodule update --init --recursive
```

**Note:** This requires SSH authentication to GitHub. If you encounter permission errors:
- Ensure your SSH keys are set up: `ssh-add -l`
- Test GitHub access: `ssh -T git@github.com`
- If needed, generate SSH keys: `ssh-keygen -t ed25519 -C "your_email@example.com"`
- Add to GitHub: https://github.com/settings/keys

### 2. Verify Environment Configuration

The `.devcontainer/.env` file has been pre-configured with your user ID:

```bash
LOCAL_UID=788715489
LOCAL_GID=788600513
```

**Optional:** Adjust resource limits in `.env` if needed:
- `MAX_CPUS` (default: 20)
- `MAX_MEMORY` (default: 100G)

### 3. Open in VS Code Dev Container

```bash
code /data1/users/antonz/projects/GVDRP1_prj/Gama_Vivian_DRP1_bulkRNAseq
```

Then:
1. Command Palette (`Cmd/Ctrl+Shift+P`)
2. Select: **"Dev Containers: Reopen in Container"**
3. Wait for container to build and start (first time may take a few minutes)

### 4. Verify Setup

Once inside the container, the post-start script will automatically run sanity checks. You should see:

```
==== Devcontainer sanity checks ====
Python: OK
R:      OK
httpgd: OK
...
```

### 5. Test Analysis Pipeline

Verify that the analysis can access all required files:

```bash
# In container terminal
cd /workspaces/Gama_Vivian_DRP1_bulkRNAseq
ls -l 03_Results/01_Preprocessing/04_FeatureCounts/count_matrices_fc/
ls -l 00_Data/SynGO_bulk_20231201/
ls -l 01_Scripts/RNAseq-toolkit/
```

### 6. Run R Analysis (Optional Test)

Start R and test data loading:

```r
# In R console (use 'radian' or 'r-base')
library(here)
here::here()  # Should show /workspaces/Gama_Vivian_DRP1_bulkRNAseq

# Test loading counts matrix
counts_file <- here::here("03_Results/01_Preprocessing/04_FeatureCounts/count_matrices_fc/sorted_counts_matrix.txt")
file.exists(counts_file)  # Should be TRUE

counts <- read.table(counts_file, header = TRUE, row.names = 1, sep = "\t")
dim(counts)  # Should show dimensions of counts matrix
```

## What Changed from Original Setup

### ✅ New in This Setup

- **Portable paths:** All `/workspaces/GVDRP1/` hardcoded paths replaced with `here::here()`
- **Git submodule:** RNAseq-toolkit now version-controlled as submodule
- **Environment variables:** `.env` file for user ID and resource configuration
- **Docker Compose:** Modern container orchestration replacing direct image mounts
- **Documentation:** MIGRATION.md, CHANGELOG.md, and this SETUP.md

### ❌ Removed from Original Setup

- **Raw data mounts:** No longer needed (preprocessing complete)
- **External GSEA module mount:** Now part of RNAseq-toolkit submodule
- **Hardcoded workstation paths:** All replaced with relative paths

## Container Features

### Available Commands

- **R environments:**
  - `radian` - Interactive R console (recommended)
  - `r-base` - Standard R session
  - `R` - Direct R invocation

- **Python environments:**
  - `usepy base|squid|atac|comms` - Switch Python virtual environments
  - `py-base`, `py-squid`, `py-atac`, `py-comms` - Direct environment invocation

### Useful Aliases

- `ll` - List files with details
- Standard git, tmux, vim commands available

## Working with the Analysis

### Main Analysis Script

The primary analysis pipeline is:

```
02_Analysis/Analysis_pipeline_fin.R
```

This script now uses:
- Relative paths via `here::here()`
- `config$helper_root` pointing to `01_Scripts/RNAseq-toolkit`
- All dependencies loaded from the RNAseq-toolkit submodule

### Key Configuration

Edit configuration at the top of `Analysis_pipeline_fin.R`:

```r
config <- list(
  counts_file   = "03_Results/01_Preprocessing/04_FeatureCounts/count_matrices_fc/sorted_counts_matrix.txt",
  metadata_file = "03_Results/01_Preprocessing/04_FeatureCounts/count_matrices_fc/metadata.csv",
  out_root      = "03_Results/02_Analysis",
  helper_root   = "01_Scripts/RNAseq-toolkit",   # git submodule
  p_cutoff      = 0.05,
  fc_cutoff     = 2,
  ...
)
```

## Troubleshooting

### Issue: "Permission denied" when cloning submodule

**Solution:** Set up SSH keys for GitHub (see Step 1 above)

### Issue: Container won't start / permission errors

**Solution:** Check `.env` file has correct UID/GID:
```bash
id -u  # Should match LOCAL_UID
id -g  # Should match LOCAL_GID
```

### Issue: "Cannot find GSEA helpers"

**Solution:** Ensure RNAseq-toolkit submodule is initialized:
```bash
git submodule status
# Should show commit hash, not a minus sign
```

### Issue: httpgd or other R packages not found

Key packages are pre-installed and frozen from modifications/updates with the [scbio-docker v0.5.1](https://github.com/tony-zhelonkin/scbio-docker/tree/v0.5.1), but should any packages necessary for the analysis be missing, they can be installed at runtime.

**Solution:** Install from within container:
```r
install.packages("httpgd")
# Or for Bioconductor packages:
BiocManager::install("package_name")
```

## Next Steps for Analysis

1. **Review reviewer comments** - Document in CHANGELOG.md
2. **Plan analysis updates** - Track in todo list
3. **Run modified analyses** - Execute specific contrasts or gene sets
4. **Generate new figures** - Save to appropriate Results subdirectories
5. **Document changes** - Update CHANGELOG.md with what was modified

## References

- **Container source:** [scbio-docker v0.5.1](https://github.com/tony-zhelonkin/scbio-docker/tree/v0.5.1)
- **RNAseq toolkit:** [RNAseq-toolkit](https://github.com/tony-zhelonkin/RNAseq-toolkit)
- **Frozen analysis:** Git tag `v1.0` (commit `d6ec164`)

## Support

For issues with:
- **Container:** See scbio-docker documentation
- **Analysis pipeline:** Check MIGRATION.md troubleshooting section
- **Git submodules:** `git submodule --help`

---

*Setup guide created: 2025-11-19*
