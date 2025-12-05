#!/usr/bin/env Rscript
###############################################################################
## Freeze R Environment - Capture all installed packages and versions       ##
###############################################################################
## This script documents all R packages used in the analysis pipeline       ##
## for reproducibility purposes.                                            ##
###############################################################################

library(here)

message("Freezing R environment...")

# Load all packages used in the analysis
packages_to_load <- c(
  # Core differential expression
  "edgeR", "limma", "DESeq2",

  # Data manipulation
  "dplyr", "tidyr", "reshape2",

  # Visualization
  "ggplot2", "pheatmap", "RColorBrewer", "viridis", "patchwork",
  "VennDiagram", "grid", "UpSetR", "gridExtra",

  # Gene set enrichment
  "clusterProfiler", "fgsea", "msigdbr", "enrichplot",

  # Annotation
  "org.Hs.eg.db", "AnnotationDbi",

  # Co-expression analysis
  #"WGCNA", was not at the moment of publication revision process

  # Utilities
  "here", "readxl", "writexl", "readr",

  # GitHub packages (may not load via library)
  "devtools"
)

# Try to load each package
for (pkg in packages_to_load) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
    message(sprintf("✓ Loaded: %s", pkg))
  } else {
    message(sprintf("✗ Not installed: %s", pkg))
  }
}

# Capture session information
session_file <- here::here("R_session_info.txt")
sink(session_file)
cat("###############################################################################\n")
cat("## R Session Information - Gama_Vivian_DRP1_bulkRNAseq\n")
cat("###############################################################################\n")
cat("## Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("## Purpose: Document R package versions for reproducibility\n")
cat("###############################################################################\n\n")

sessionInfo()

cat("\n\n###############################################################################\n")
cat("## Detailed Package Versions\n")
cat("###############################################################################\n\n")

# Get all installed packages with versions
installed <- as.data.frame(installed.packages()[, c("Package", "Version", "Priority")])
installed <- installed[order(installed$Package), ]

# Separate base and user-installed packages
base_pkgs <- installed[!is.na(installed$Priority), ]
user_pkgs <- installed[is.na(installed$Priority), ]

cat("BASE R PACKAGES (", nrow(base_pkgs), "):\n", sep = "")
cat("----------------------------------------\n")
for (i in 1:nrow(base_pkgs)) {
  cat(sprintf("%-30s %s\n", base_pkgs$Package[i], base_pkgs$Version[i]))
}

cat("\n\nUSER-INSTALLED PACKAGES (", nrow(user_pkgs), "):\n", sep = "")
cat("----------------------------------------\n")
for (i in 1:nrow(user_pkgs)) {
  cat(sprintf("%-30s %s\n", user_pkgs$Package[i], user_pkgs$Version[i]))
}

cat("\n\n###############################################################################\n")
cat("## Runtime Install Commands (from 02_Analysis/0.1.runtime_installs.R)\n")
cat("###############################################################################\n\n")
cat("# Bioconductor packages\n")
cat('BiocManager::install(c("WGCNA", "devtools", "org.Hs.eg.db"))\n\n')
cat("# CRAN packages\n")
cat("install.packages('msigdbr', repos = 'https://cloud.r-project.org')\n\n")
cat("# GitHub packages\n")
cat("# Note: PerseusR may not be actively used in current analysis\n")
cat('devtools::install_github("jdrudolph/PerseusR")\n')

cat("\n\n###############################################################################\n")
cat("## Key Package Sources\n")
cat("###############################################################################\n\n")

key_packages <- c("edgeR", "limma", "clusterProfiler", "org.Hs.eg.db",
                  "msigdbr", "WGCNA", "fgsea", "ggplot2", "dplyr")

for (pkg in key_packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    desc <- packageDescription(pkg)
    cat(sprintf("\n%s:\n", pkg))
    cat(sprintf("  Version: %s\n", desc$Version))
    if (!is.null(desc$Repository)) {
      cat(sprintf("  Repository: %s\n", desc$Repository))
    }
    if (!is.null(desc$GithubRepo)) {
      cat(sprintf("  GitHub: %s/%s\n", desc$GithubUsername, desc$GithubRepo))
    }
    if (!is.null(desc$biocViews)) {
      cat("  Source: Bioconductor\n")
    }
  }
}

sink()

message(sprintf("\n✓ Session information saved to: %s", session_file))

# Also create a simple package list for quick reference
pkg_list_file <- here::here("R_packages.txt")
writeLines(
  c(
    "# R Packages Used in Analysis",
    "# Format: Package (Version)",
    "# Generated: " %+% format(Sys.time(), "%Y-%m-%d"),
    "",
    paste0(user_pkgs$Package, " (", user_pkgs$Version, ")")
  ),
  pkg_list_file
)

message(sprintf("✓ Package list saved to: %s", pkg_list_file))

message("\nDone! Files created:")
message(sprintf("  - %s (detailed session info)", session_file))
message(sprintf("  - %s (simple package list)", pkg_list_file))