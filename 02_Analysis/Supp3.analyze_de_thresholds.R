#!/usr/bin/env Rscript
# Inspecting DE threshold analysis across all contrasts
# Analyzes gene counts at FDR: 0.5, 0.1, 0.05, 0.01

library(here)
library(dplyr)
library(tidyr)

# Configuration
config <- list(
  out_root = here("03_Results/02_Analysis"),
  thresholds = c(0.5, 0.1, 0.05, 0.01)
)

# Define all contrasts in order
all_contrasts <- c(
  "G32A_vs_Ctrl_D35",
  "G32A_vs_Ctrl_D65",
  "R403C_vs_Ctrl_D35",
  "R403C_vs_Ctrl_D65",
  "Time_Ctrl",
  "Time_G32A",
  "Time_R403C",
  "Maturation_G32A_specific",
  "Maturation_R403C_specific"
)

message("📊 Analyzing DE results across FDR thresholds...\n")

# Initialize results storage
threshold_counts <- list()
top_genes_list <- list()

# Process each contrast
for (contrast_name in all_contrasts) {
  message("Processing: ", contrast_name)

  # Load DE results
  de_file <- here(config$out_root, "DE_results",
                  paste0(contrast_name, "_DE_results.csv"))

  if (!file.exists(de_file)) {
    message("  ⚠️  File not found: ", de_file)
    next
  }

  de_results <- read.csv(de_file, row.names = 1)

  # Count genes at each threshold
  counts <- sapply(config$thresholds, function(thresh) {
    sum(de_results$adj.P.Val <= thresh, na.rm = TRUE)
  })

  threshold_counts[[contrast_name]] <- counts

  # Get top 10 genes by adj.P.Val
  top_genes <- de_results %>%
    filter(!is.na(adj.P.Val)) %>%
    arrange(adj.P.Val) %>%
    slice_head(n = 10) %>%
    mutate(
      Contrast = contrast_name,
      Gene = rownames(.),
      logFC = round(logFC, 2),
      AveExpr = round(AveExpr, 2),
      P.Value = signif(P.Value, 3),
      adj.P.Val = signif(adj.P.Val, 3),
      B = round(B, 2)
    ) %>%
    select(Contrast, Gene, logFC, AveExpr, P.Value, adj.P.Val, B)

  top_genes_list[[contrast_name]] <- top_genes

  message("  ✓ Processed")
}

message("\n", paste0(rep("=", 80), collapse = ""), "\n")
message("SUMMARY: Gene Counts Across FDR Thresholds\n")
message(paste0(rep("=", 80), collapse = ""), "\n")

# Create summary table
summary_df <- do.call(rbind, lapply(names(threshold_counts), function(name) {
  data.frame(
    Contrast = name,
    `FDR_0.5` = threshold_counts[[name]][1],
    `FDR_0.1` = threshold_counts[[name]][2],
    `FDR_0.05` = threshold_counts[[name]][3],
    `FDR_0.01` = threshold_counts[[name]][4],
    check.names = FALSE
  )
}))

# Print summary table
print(summary_df, row.names = FALSE)

# Save summary table
summary_file <- here(config$out_root, "DE_threshold_summary.csv")
write.csv(summary_df, summary_file, row.names = FALSE)
message("\n✓ Summary table saved to: ", summary_file)

# Combine all top genes
all_top_genes <- do.call(rbind, top_genes_list)

# Save top genes table
top_genes_file <- here(config$out_root, "Top10_genes_by_contrast.csv")
write.csv(all_top_genes, top_genes_file, row.names = FALSE)
message("✓ Top 10 genes table saved to: ", top_genes_file)

message("\n", paste0(rep("=", 80), collapse = ""), "\n")
message("KEY FINDINGS\n")
message(paste0(rep("=", 80), collapse = ""), "\n")

# Identify contrasts with no significant genes at FDR 0.1
no_sig <- summary_df %>%
  filter(`FDR_0.1` == 0) %>%
  pull(Contrast)

if (length(no_sig) > 0) {
  message("\n Contrasts with NO genes at FDR ≤ 0.1:")
  for (c in no_sig) {
    message("  - ", c)
  }
}

# Identify contrasts with very few genes at FDR 0.1
few_sig <- summary_df %>%
  filter(`FDR_0.1` > 0 & `FDR_0.1` <= 5) %>%
  pull(Contrast)

if (length(few_sig) > 0) {
  message("\n🟡 Contrasts with 1-5 genes at FDR ≤ 0.1:")
  for (c in few_sig) {
    n <- summary_df %>% filter(Contrast == c) %>% pull(`FDR_0.1`)
    message("  - ", c, " (", n, " genes)")
  }
}

# Identify contrasts with many genes
many_sig <- summary_df %>%
  filter(`FDR_0.1` > 100) %>%
  pull(Contrast)

if (length(many_sig) > 0) {
  message("\n Contrasts with >100 genes at FDR ≤ 0.1:")
  for (c in many_sig) {
    n <- summary_df %>% filter(Contrast == c) %>% pull(`FDR_0.1`)
    message("  - ", c, " (", n, " genes)")
  }
}

# Print top genes for key contrasts
message("\n", paste0(rep("=", 80), collapse = ""), "\n")
message("TOP 5 GENES BY CONTRAST (sorted by adj.P.Val)\n")
message(paste0(rep("=", 80), collapse = ""), "\n")

for (contrast_name in all_contrasts) {
  if (contrast_name %in% names(top_genes_list)) {
    message("\n", contrast_name, ":")
    top_5 <- top_genes_list[[contrast_name]] %>% slice_head(n = 5)
    for (i in 1:nrow(top_5)) {
      message(sprintf("  %d. %-12s  logFC=%-6.2f  adj.P.Val=%s",
                      i,
                      top_5$Gene[i],
                      top_5$logFC[i],
                      format(top_5$adj.P.Val[i], scientific = TRUE, digits = 3)))
    }
  }
}

message("\n✅ Analysis complete!")
