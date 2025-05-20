


#------------------------#
# 9. Transcription Factor Analysis
#------------------------#

library(decoupleR)
library(dorothea)

# Create output directory for TF analysis
tf_dir <- "03_Results/02_Analysis/TF_Analysis"
dir.create(tf_dir, recursive = TRUE, showWarnings = FALSE)

# Load dorothea regulons
data(dorothea_hs, package = "dorothea")

# Extract high-confidence regulons (A, B, and C levels)
regulons <- dorothea_hs %>%
  dplyr::filter(confidence %in% c("A", "B", "C"))

# Run TF activity analysis for each contrast
for (i in 1:ncol(contrasts)) {
  contrast_name <- colnames(contrasts)[i]
  
  # Get results for this contrast
  results <- topTable(fit, coef = i, number = Inf)
  
  # Prepare input for decoupleR - vector of t-statistics
  t_stats <- results$t
  names(t_stats) <- rownames(results)
  
  # Run TF activity analysis
  tf_activities <- run_wmean(
    mat = t_stats,
    network = regulons,
    .source = "tf",
    .target = "target",
    .mor = "mor",
    times = 100,  # Number of permutations
    seed = 42
  )
  
  # Extract normalized enrichment scores and p-values
  tf_results <- tf_activities %>%
    dplyr::filter(statistic == "norm_wmean") %>%
    dplyr::select(source, score, pvalue) %>%
    dplyr::rename(tf = source, nes = score, p.value = pvalue) %>%
    dplyr::arrange(p.value)
  
  # Add BH adjusted p-values
  tf_results$p.adj <- p.adjust(tf_results$p.value, method = "BH")
  
  # Save results to file
  write.csv(tf_results, file = file.path(tf_dir, paste0(contrast_name, "_TF_activities.csv")), row.names = FALSE)
  
  # Create barplot of top 20 TFs (by absolute NES)
  top_tfs <- tf_results %>%
    dplyr::arrange(desc(abs(nes))) %>%
    head(20)
  
  pdf(file.path(tf_dir, paste0(contrast_name, "_top_TFs.pdf")), width = 10, height = 8)
  ggplot(top_tfs, aes(x = reorder(tf, nes), y = nes, fill = nes > 0)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_manual(values = c("#1B9E77", "#D95F02"), name = "Direction", labels = c("Repressed", "Activated")) +
    labs(
      title = paste0("Top 20 TF Activities - ", contrast_name),
      x = "Transcription Factor",
      y = "Normalized Enrichment Score"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12)
    )
  dev.off()
  
  # Create volcano plot for TF activities
  pdf(file.path(tf_dir, paste0(contrast_name, "_TF_volcano.pdf")), width = 8, height = 7)
  ggplot(tf_results, aes(x = nes, y = -log10(p.value), label = tf, color = p.adj < 0.05)) +
    geom_point(size = 2, alpha = 0.7) +
    scale_color_manual(values = c("grey", "#D95F02"), name = "Significant", labels = c("No", "Yes")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    geom_text_repel(
      data = subset(tf_results, p.adj < 0.05 | abs(nes) > 2),
      aes(label = tf),
      size = 3,
      max.overlaps = 20
    ) +
    labs(
      title = paste0("TF Activity Volcano Plot - ", contrast_name),
      x = "Normalized Enrichment Score",
      y = "-log10(p-value)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      legend.position = "bottom"
    )
  dev.off()
}

# Find common TFs across contrasts
tf_lists <- list()
for (i in 1:ncol(contrasts)) {
  contrast_name <- colnames(contrasts)[i]
  tf_results_file <- file.path(tf_dir, paste0(contrast_name, "_TF_activities.csv"))
  
  if (file.exists(tf_results_file)) {
    tf_results <- read.csv(tf_results_file)
    tf_lists[[contrast_name]] <- tf_results$tf[tf_results$p.adj < 0.05]
  }
}

# Create upset plot for TFs
if (length(tf_lists) > 0 && sum(sapply(tf_lists, length)) > 0) {
  pdf(file.path(tf_dir, "TF_upset_plot.pdf"), width = 10, height = 8)
  upset(fromList(tf_lists), order.by = "freq", nsets = length(tf_lists))
  dev.off()
}




# Summarize TF activity results
# Collect top TFs from each contrast
top_tfs_summary <- data.frame()

# List all TF results files
tf_files <- list.files(tf_dir, pattern = "_TF_activities\\.csv$", full.names = TRUE)

for (tf_file in tf_files) {
  # Extract contrast name from filename
  contrast_name <- gsub("_TF_activities\\.csv$", "", basename(tf_file))
  
  # Read results
  if (file.exists(tf_file)) {
    tf_results <- read.csv(tf_file)
    
    # Get top 5 TFs by absolute NES
    if (nrow(tf_results) > 0) {
      tf_results$abs_nes <- abs(tf_results$nes)
      top5 <- tf_results[order(tf_results$abs_nes, decreasing = TRUE), ][1:min(5, nrow(tf_results)), ]
      
      # Add contrast info
      top5$Contrast <- contrast_name
      
      # Keep only relevant columns
      top5 <- top5[, c("tf", "nes", "p.adj", "Contrast")]
      
      # Combine with full results
      top_tfs_summary <- rbind(top_tfs_summary, top5)
    }
  }
}

# Save top TFs to file
if (nrow(top_tfs_summary) > 0) {
  write.csv(top_tfs_summary, file = file.path(summary_dir, "top_TFs_summary.csv"), row.names = FALSE)
}

# Summarize WGCNA results
# List significant modules and their trait associations
module_summary <- data.frame()

for (module in significant_modules) {
  module_color <- gsub("ME", "", module)
  
  # Get module-trait correlations
  module_trait_cor_row <- ME_trait_cor[module, ]
  module_trait_pval_row <- ME_trait_pval[module, ]
  
  # Get number of genes in the module
  module_idx <- which(labels2colors(net$colors) == module_color)
  n_genes <- sum(net$colors == module_idx - 1)
  
  # Get top hub genes (top 3)
  top_hub_genes <- paste(hub_genes[[module]][1:min(3, length(hub_genes[[module]]))], collapse = ", ")
  
  # Create row for this module
  module_row <- data.frame(
    Module = module_color,
    N_genes = n_genes,
    Genotype_cor = module_trait_cor_row["Genotype"],
    Genotype_pval = module_trait_pval_row["Genotype"],
    Day_cor = module_trait_cor_row["Day_D65"],
    Day_pval = module_trait_pval_row["Day_D65"],
    Top_hub_genes = top_hub_genes
  )
  
  # Combine with full results
  module_summary <- rbind(module_summary, module_row)
}

# Save module summary to file
if (nrow(module_summary) > 0) {
  write.csv(module_summary, file = file.path(summary_dir, "WGCNA_module_summary.csv"), row.names = FALSE)
}
