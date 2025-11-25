# Installation Guide

This guide provides instructions for setting up the environment to reproduce the DRP1 bulk RNA-seq analysis.

## Prerequisites

- **Docker** (20.10+)
- **VS Code** (1.75+) with Dev Containers extension
- **Git** (2.30+) with submodule support
- **Disk space:** ~30 GB (~17 GB container + 5 GB analysis results)
- **RAM:** 8 GB minimum, 16 GB recommended

## Quick Start

```bash
# 1. Build Docker image (scbio-docker v0.5.1)
git clone git@github.com:tony-zhelonkin/scbio-docker.git
cd scbio-docker
git fetch --tags
git checkout tags/v0.5.1
./build-optimized.sh --tag scdock-r-dev:v0.5.1

# 2. Clone this analysis repository
cd ..
git clone <repo-url>
cd Gama_Vivian_DRP1_bulkRNAseq

# 3. Initialize git submodules (RNAseq-toolkit)
git submodule update --init --recursive

# 4. Open in VS Code and launch container
code .
# Command Palette → "Dev Containers: Reopen in Container"

# 5. Inside container, install runtime R packages
Rscript 02_Analysis/0.runtime_installs.R

# 6. Verify installation
Rscript freeze_requirements.R

# 7. Run analysis
Rscript 02_Analysis/1a.Main_pipeline.R
```

**Estimated time:** 30-45 minutes (20-30 min container build + 10-15 min package installation)

## Detailed Setup

### Step 1: Build Docker Image

The analysis was run in the **scbio-docker v0.5.1** container environment. You must build this image first.

**Repository:** https://github.com/tony-zhelonkin/scbio-docker/tree/v0.5.1

```bash
# Clone repository
git clone https://github.com/tony-zhelonkin/scbio-docker.git
cd scbio-docker

# Checkout release tag v0.5.1
git fetch --tags
git checkout tags/v0.5.1

# Build image (~20-30 minutes)
./build-optimized.sh --tag scdock-r-dev:v0.5.1

# Optional: provide GitHub PAT to avoid rate limits
./build-optimized.sh \
  --github-pat ghp_your_token_here \
  --tag scdock-r-dev:v0.5.1

# Verify build
docker images | grep scdock-r-dev
```

**What this provides:**
- R 4.3+ with ~80 pre-installed packages
- Python 3.9+ with scientific libraries
- Bioconductor ecosystem
- System tools (samtools, bcftools, bedtools)

**For details:** See the [scbio-docker README](https://github.com/tony-zhelonkin/scbio-docker/blob/v0.5.1/README.md)

### Step 2: Clone Analysis Repository

```bash
cd ..  # Leave scbio-docker directory
git clone <repo-url>
cd Gama_Vivian_DRP1_bulkRNAseq

# Initialize submodules (RNAseq-toolkit helper functions)
git submodule update --init --recursive

# Verify submodule
ls -la 01_Scripts/RNAseq-toolkit/scripts/
```

### Step 3: Launch Development Container

**Option A: VS Code Dev Containers (Recommended)**

1. Open folder in VS Code:
   ```bash
   code .
   ```

2. Launch container:
   - Command Palette (`Ctrl+Shift+P` / `Cmd+Shift+P`)
   - Select: **"Dev Containers: Reopen in Container"**
   - Wait for container to start (~30 seconds)

3. Verify environment inside container:
   ```bash
   R --version      # Should show R 4.3+
   python3 --version # Should show Python 3.9+
   which radian     # Should show /usr/local/bin/radian
   here::here()     # Should show /workspaces/Gama_Vivian_DRP1_bulkRNAseq
   ```

**Option B: Manual Docker Run**

```bash
docker run -it --rm \
  -v $(pwd):/workspaces/Gama_Vivian_DRP1_bulkRNAseq \
  -w /workspaces/Gama_Vivian_DRP1_bulkRNAseq \
  --memory=16g \
  scdock-r-dev:v0.5.1 \
  bash
```

### Step 4: Install Runtime Dependencies

The container includes most packages, but some must be installed at runtime:

```bash
# Inside container - automated installation
Rscript 02_Analysis/0.runtime_installs.R
```

This installs:
- **WGCNA** (co-expression analysis) # was not used in the paper submission
- **org.Hs.eg.db** (human gene annotations)
- **msigdbr** (MSigDB gene sets)
- **devtools** (GitHub package installation)

**Manual installation (if automated script fails):**

```r
# In R console
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("WGCNA", "devtools", "org.Hs.eg.db"))
install.packages('msigdbr', repos = 'https://cloud.r-project.org')
```

**Python packages** (most are already in container):

```bash
# Optional: ensure all core packages present
pip install -r requirements.txt
```

### Step 5: Verify Installation

```bash
# Generate session info to verify R packages
Rscript freeze_requirements.R

# Check output files created
ls -lh R_session_info.txt R_packages.txt

# Test Python imports
python3 -c "import pandas, numpy, matplotlib, upsetplot; print('✓ Python OK')"

# Quick test of analysis pipeline (~5-10 min with checkpoints)
Rscript 02_Analysis/1a.Main_pipeline.R
```

If the main pipeline runs without errors, your environment is correctly configured.

## Package Requirements

### R Environment

**Package manifest files:**
- `R_session_info.txt` - Complete sessionInfo() output with all versions
- `R_packages.txt` - Simple package list
- `02_Analysis/0.runtime_installs.R` - Runtime installation script

**Core packages** (pre-installed in container):

| Package | Version | Purpose |
|---------|---------|---------|
| edgeR | 3.42+ | Normalization, TMM |
| limma | 3.56+ | Differential expression |
| clusterProfiler | 4.8+ | GSEA wrapper |
| fgsea | 1.26+ | Fast GSEA |
| msigdbr | 7.5+ | MSigDB gene sets |
| org.Hs.eg.db | 3.17+ | Gene annotations |
| WGCNA | 1.72+ | Co-expression |
| ggplot2 | 3.4+ | Visualization |
| dplyr | 1.1+ | Data manipulation |

**See:** `R_session_info.txt` for complete list (generated after setup)

### Python Environment

**Requirement files:**
- `requirements.txt` - Core packages (install with `pip install -r requirements.txt`)
- `python_requirements_freeze.txt` - Exact versions for reproducibility

**Core packages** (pre-installed in container):

| Package | Version | Purpose |
|---------|---------|---------|
| pandas | 1.3+ | Data manipulation |
| numpy | 1.21+ | Numerical computing |
| matplotlib | 3.4+ | Plotting |
| upsetplot | 0.6+ | Intersection plots |

## Troubleshooting

### 1. Docker image not found

**Symptom:** `Error: No such image: scdock-r-dev:v0.5.1`

**Solution:** Build the scbio-docker image first (Step 1)

### 2. Submodule directory empty

**Symptom:** `01_Scripts/RNAseq-toolkit/` is empty

**Solution:**
```bash
git submodule update --init --recursive
```

### 3. Permission denied in container

**Symptom:** Cannot write files

**Solution:** Configure UID/GID in `.devcontainer/.env`:
```bash
echo "HOST_UID=$(id -u)" > .devcontainer/.env
echo "HOST_GID=$(id -g)" >> .devcontainer/.env

# Rebuild container
# Command Palette → "Dev Containers: Rebuild Container"
```

### 4. R package installation fails

**Symptom:** BiocManager packages fail to install

**Solution:**
```r
# Update BiocManager
install.packages("BiocManager")
BiocManager::install(version = "3.17")

# Install packages individually
BiocManager::install("WGCNA")
BiocManager::install("org.Hs.eg.db")
```

### 5. "Installation paths not writeable" warning

**This is NORMAL and EXPECTED.** The warning appears when BiocManager checks for updates to system packages (which are read-only by design for reproducibility). Your packages still install successfully to the user library (`~/R/...`).

To suppress:
```r
BiocManager::install("package", update = FALSE)
```

### 6. Checkpoint cache errors

**Symptom:** Analysis fails loading cached checkpoints

**Solution:** Delete cache and recompute:
```bash
rm -rf 03_Results/02_Analysis/checkpoints/*.rds
Rscript 02_Analysis/1a.Main_pipeline.R
```

Or set `force_recompute = TRUE` in `02_Analysis/1a.Main_pipeline.R`

## Reproducing Exact Environment

For full reproducibility using exact package versions:

```bash
# R packages - generate session info
Rscript freeze_requirements.R

# Compare your output to committed R_session_info.txt
diff R_session_info.txt <(cat <<committed version>>)

# Python packages - install exact frozen versions
pip install -r python_requirements_freeze.txt
```

## Summary

**Minimum setup steps:**
1. Build scbio-docker v0.5.1 image (~20-30 min)
2. Clone this repository with submodules (~2 min)
3. Open in VS Code Dev Container (~1 min)
4. Install runtime R packages (~5-10 min)
5. Verify with test analysis (~5-10 min)

**Total time:** ~30-50 minutes

**Result:** Fully configured environment ready to reproduce analysis

## Additional Resources

- **Container documentation:** [scbio-docker v0.5.1](https://github.com/tony-zhelonkin/scbio-docker/tree/v0.5.1)
- **Analysis documentation:** [CLAUDE.md](CLAUDE.md) - Detailed analysis guide
- **Script reference:** [02_Analysis/SCRIPTS.md](02_Analysis/SCRIPTS.md) - Script inventory
- **Session history:** [docs/SESSION_HISTORY.md](docs/SESSION_HISTORY.md) - Development log