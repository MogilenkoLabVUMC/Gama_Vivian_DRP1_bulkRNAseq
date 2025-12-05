#!/usr/bin/env Rscript
###############################################################################
## Calcium Gene Diagnostic Script
## Purpose: Investigate which calcium genes are present in the dataset
##          and their statistical significance across all contrasts
###############################################################################

library(here)
library(limma)
library(dplyr)

# Load checkpoint
fit <- readRDS(here::here("03_Results/02_Analysis/checkpoints/fit_object.rds"))
contrasts <- readRDS(here::here("03_Results/02_Analysis/checkpoints/contrasts_matrix.rds"))

# Calcium genes list (from config)
calcium_genes <- c("NNAT","CACNG3","CACNA1C","CACNA1S","ATP2A1",
                   "RYR1","MYLK3","CASR","VDR","STIM1","STIM2",
                   "ORAI1","CALB1","CALR","PNPO")

cat("###############################################################################\n")
cat("## CALCIUM GENE DIAGNOSTIC REPORT\n")
cat("###############################################################################\n\n")

cat("=== 1. CALCIUM GENE PRESENCE IN COUNT MATRIX ===\n\n")

# Check presence in fit object
present <- calcium_genes %in% rownames(fit)
names(present) <- calcium_genes

cat(sprintf("Total calcium genes checked: %d\n", length(calcium_genes)))
cat(sprintf("Present in dataset: %d\n", sum(present)))
cat(sprintf("Missing from dataset: %d\n\n", sum(!present)))

cat("Present genes:\n")
for (g in calcium_genes[present]) {
  cat(sprintf("  ✓ %s\n", g))
}

if (any(!present)) {
  cat("\nMissing genes:\n")
  for (g in calcium_genes[!present]) {
    cat(sprintf("  ✗ %s\n", g))
  }

  cat("\n=== Searching for Alternative Symbols ===\n\n")
  missing <- calcium_genes[!present]
  for (g in missing) {
    matches <- grep(g, rownames(fit), ignore.case = TRUE, value = TRUE)
    if (length(matches) > 0) {
      cat(sprintf("  %s alternatives: %s\n", g, paste(matches, collapse = ", ")))
    } else {
      cat(sprintf("  %s: NO MATCHES FOUND\n", g))
    }
  }
}

cat("\n\n=== 2. CALCIUM GENE STATISTICS BY CONTRAST ===\n\n")

# Create summary table
summary_list <- list()

for (co in colnames(contrasts)) {
  tbl <- limma::topTable(fit, coef = co, number = Inf, sort.by = "t")
  present_in_contrast <- calcium_genes[calcium_genes %in% rownames(tbl)]

  if (length(present_in_contrast) > 0) {
    # Get stats for present genes
    calcium_stats <- tbl[present_in_contrast, , drop = FALSE]

    # Count by different thresholds
    sig_fdr_2fc <- sum(calcium_stats$adj.P.Val <= 0.1 & abs(calcium_stats$logFC) >= 1)
    sig_fdr <- sum(calcium_stats$adj.P.Val <= 0.1)
    sig_p_2fc <- sum(calcium_stats$P.Value <= 0.05 & abs(calcium_stats$logFC) >= 1)
    sig_p <- sum(calcium_stats$P.Value <= 0.05)

    cat(sprintf("Contrast: %s\n", co))
    cat(sprintf("  Present: %d/%d calcium genes\n", length(present_in_contrast), length(calcium_genes)))
    cat(sprintf("  Significant (FDR ≤ 0.1 & |FC| ≥ 2): %d\n", sig_fdr_2fc))
    cat(sprintf("  Significant (FDR ≤ 0.1 only): %d\n", sig_fdr))
    cat(sprintf("  Significant (p ≤ 0.05 & |FC| ≥ 2): %d\n", sig_p_2fc))
    cat(sprintf("  Significant (p ≤ 0.05 only): %d\n\n", sig_p))

    # Store for detailed report
    summary_list[[co]] <- list(
      present_genes = present_in_contrast,
      sig_fdr_2fc = rownames(calcium_stats)[calcium_stats$adj.P.Val <= 0.1 & abs(calcium_stats$logFC) >= 1],
      sig_p_2fc = rownames(calcium_stats)[calcium_stats$P.Value <= 0.05 & abs(calcium_stats$logFC) >= 1]
    )
  } else {
    cat(sprintf("Contrast: %s\n", co))
    cat("  ⚠️  No calcium genes present in this contrast\n\n")
  }
}

cat("\n=== 3. DETAILED CALCIUM GENE STATISTICS ===\n\n")

# Show top 3 most significant calcium genes per contrast
for (co in colnames(contrasts)) {
  tbl <- limma::topTable(fit, coef = co, number = Inf, sort.by = "t")
  present_in_contrast <- calcium_genes[calcium_genes %in% rownames(tbl)]

  if (length(present_in_contrast) > 0) {
    calcium_stats <- tbl[present_in_contrast, , drop = FALSE]
    calcium_stats <- calcium_stats[order(calcium_stats$P.Value), ]

    cat(sprintf("=== %s ===\n", co))
    cat("Top 3 calcium genes by p-value:\n")

    n_show <- min(3, nrow(calcium_stats))
    for (i in 1:n_show) {
      gene <- rownames(calcium_stats)[i]
      logFC <- calcium_stats$logFC[i]
      pval <- calcium_stats$P.Value[i]
      fdr <- calcium_stats$adj.P.Val[i]

      sig_marker <- ""
      if (fdr <= 0.1 & abs(logFC) >= 1) {
        sig_marker <- " ***"
      } else if (fdr <= 0.1) {
        sig_marker <- " **"
      } else if (pval <= 0.05) {
        sig_marker <- " *"
      }

      cat(sprintf("  %d. %s: logFC=%.3f, p=%.2e, FDR=%.2e%s\n",
                  i, gene, logFC, pval, fdr, sig_marker))
    }
    cat("\n")
  }
}

cat("\n=== 4. RECOMMENDATIONS ===\n\n")

# Check if any calcium genes are missing
if (any(!present)) {
  cat("⚠️  ACTION REQUIRED:\n")
  cat(sprintf("   %d calcium genes are missing from the dataset.\n", sum(!present)))
  cat("   Missing genes:", paste(calcium_genes[!present], collapse = ", "), "\n")
  cat("   → Check for gene symbol updates or aliases\n")
  cat("   → Consider updating calcium_genes list in config\n\n")
} else {
  cat("✓ All calcium genes present in dataset\n\n")
}

# Check for visualization issues
cat("VISUALIZATION CONSIDERATIONS:\n")
cat("• Current implementation mixes calcium genes with regular labels\n")
cat("• Recommendation: Use separate ggrepel layer with max.overlaps=Inf\n")
cat("• This ensures calcium genes are NEVER suppressed by overlap conflicts\n\n")

cat("###############################################################################\n")
cat("## END OF DIAGNOSTIC REPORT\n")
cat("###############################################################################\n")
