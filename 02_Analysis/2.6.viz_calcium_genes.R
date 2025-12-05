###############################################################################
##  Calcium genes visualization script                                       ##
##  - Loads checkpoints from main Analysis_pipeline.R                        ##
##  - Generates calcium-specific heatmaps, boxplots, and volcano plots       ##
###############################################################################

# -------------------------------------------------------------------- #
# 0.  Configuration                                                    #
# -------------------------------------------------------------------- #
config <- list(
  out_root      = "03_Results/02_Analysis",
  calcium_genes = c(
    "NNAT","CACNG3","CACNA1S","ATP2A1",
    "RYR1","MYLK3","VDR","STIM1","STIM2",
    "ORAI1_1","CALB1","CALR","PNPO"
  )
  # Note: CACNA1C, CASR, ORAI1 not present in dataset. ORAI1_1 is the alternative symbol for ORAI1.
)

checkpoint_dir <- here::here(config$out_root, "checkpoints")

# -------------------------------------------------------------------- #
# 1.  Load required packages                                           #
# -------------------------------------------------------------------- #
required_pkgs <- c("ggplot2", "pheatmap", "RColorBrewer",
                   "reshape2", "limma", "here")

for (p in required_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    message(sprintf("• installing %s …", p))
    install.packages(p, repos = "https://cloud.r-project.org")
  }
  library(p, character.only = TRUE)
}

# Load unified color configuration
source(here::here("01_Scripts/R_scripts/color_config.R"))

# -------------------------------------------------------------------- #
# 2.  Load checkpoints from main pipeline                              #
# -------------------------------------------------------------------- #
message("📂 Loading checkpoints from main analysis...")

# Load model objects
fit <- readRDS(file.path(checkpoint_dir, "fit_object.rds"))
DGE <- readRDS(file.path(checkpoint_dir, "DGE_object.rds"))
contrasts <- readRDS(file.path(checkpoint_dir, "contrasts_matrix.rds"))

# Load QC variables
qc_vars <- readRDS(file.path(checkpoint_dir, "qc_variables.rds"))
ordered_samples <- qc_vars$ordered_samples
annot <- qc_vars$annot
ann_colors <- qc_vars$ann_colors
logCPM <- qc_vars$logCPM

message("✓ Checkpoints loaded successfully")

# -------------------------------------------------------------------- #
# 3.  Load helper functions for volcano plots                          #
# -------------------------------------------------------------------- #
helper_root <- "01_Scripts/RNAseq-toolkit"
source_if_present <- function(...) {
  path <- here::here(...)
  if (file.exists(path)) {
    source(path, echo = FALSE)
  } else {
    warning("helper not found → ", path)
  }
}

source_if_present(helper_root, "scripts/DE/plot_standard_volcano.R")

# -------------------------------------------------------------------- #
# 4.  Calcium genes analysis                                           #
# -------------------------------------------------------------------- #

# Create output directory for calcium gene analysis
calcium_dir <- file.path(config$out_root, "Calcium_genes")
dir.create(calcium_dir, recursive = TRUE, showWarnings = FALSE)

# Extract expression data for calcium genes that are present in our dataset
calcium_genes_present <- config$calcium_genes[config$calcium_genes %in% rownames(logCPM)]

message("📊 Analyzing ", length(calcium_genes_present), " calcium genes...")

if (length(calcium_genes_present) > 0) {
  # Extract expression values
  calcium_expr <- logCPM[calcium_genes_present, ordered_samples]

  # -------------------------------------------------------------------- #
  # 4.1  Heatmap of expression across all samples                        #
  # -------------------------------------------------------------------- #
  message("  • Creating expression heatmap...")

  # Use unified color configuration for heatmap annotations
  # This avoids conflicts between diverging gradient and categorical colors
  heatmap_ann_colors <- get_heatmap_ann_colors()

  pdf(file.path(calcium_dir, "calcium_genes_expression_heatmap.pdf"),
      width = 12, height = length(calcium_genes_present) * 0.4 + 3)
  pheatmap(
    calcium_expr,
    main = "Expression of Calcium Signaling Genes",
    annotation_col = annot,
    annotation_colors = heatmap_ann_colors,
    color = get_diverging_palette(100),  # Unified Blue-White-Orange gradient
    border_color = NA,
    fontsize = 10,
    fontsize_row = 10,
    fontsize_col = 9,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    scale = "row"  # Scale by row for better visualization
  )
  dev.off()

  # -------------------------------------------------------------------- #
  # 4.2  Boxplots to compare expression across genotypes and timepoints  #
  # -------------------------------------------------------------------- #
  message("  • Creating boxplots...")
  calcium_expr_melted <- reshape2::melt(
    data.frame(Gene = rownames(calcium_expr), calcium_expr),
    id.vars = "Gene",
    variable.name = "Sample",
    value.name = "Expression"
  )

  # Add sample metadata
  calcium_expr_melted$Genotype <- DGE$samples$genotype[match(as.character(calcium_expr_melted$Sample), rownames(DGE$samples))]
  calcium_expr_melted$Days <- DGE$samples$days[match(as.character(calcium_expr_melted$Sample), rownames(DGE$samples))]

  # Create a multi-panel plot for each gene
  pdf(file.path(calcium_dir, "calcium_genes_boxplots.pdf"), width = 12, height = 10)

  # Set up multi-panel layout
  num_genes <- length(calcium_genes_present)
  ncols <- min(3, num_genes)
  nrows <- ceiling(num_genes / ncols)
  par(mfrow = c(nrows, ncols))

  # Loop through each gene and create boxplot
  for (gene in calcium_genes_present) {
    gene_data <- subset(calcium_expr_melted, Gene == gene)

    # Create boxplot
    boxplot(
      Expression ~ Genotype:Days,
      data = gene_data,
      main = gene,
      xlab = "",
      ylab = "Log2 CPM",
      col = c("#1B9E77", "#D95F02", "#7570B3", "#1B9E77", "#D95F02", "#7570B3"),
      names = c("Ctrl D35", "G32A D35", "R403C D35", "Ctrl D65", "G32A D65", "R403C D65"),
      las = 2  # Rotate x-axis labels
    )

    # Add points for individual samples
    stripchart(
      Expression ~ Genotype:Days,
      data = gene_data,
      vertical = TRUE,
      method = "jitter",
      add = TRUE,
      pch = 19,
      col = "black",
      cex = 0.6
    )
  }

  dev.off()

  # -------------------------------------------------------------------- #
  # 4.3  Test for differential expression of calcium genes               #
  # -------------------------------------------------------------------- #
  message("  • Testing differential expression...")
  calcium_gene_results <- data.frame()

  for (i in 1:ncol(contrasts)) {
    contrast_name <- colnames(contrasts)[i]

    # Get results for this contrast
    results <- limma::topTable(fit, coef = i, number = Inf)

    # Extract results for calcium genes
    calcium_results <- results[rownames(results) %in% calcium_genes_present, ]

    if (nrow(calcium_results) > 0) {
      # Add contrast name
      calcium_results$Contrast <- contrast_name

      # Combine with full results
      calcium_gene_results <- rbind(calcium_gene_results, calcium_results)
    }
  }

  # Save results to file
  write.csv(calcium_gene_results,
            file = file.path(calcium_dir, "calcium_genes_DE_results.csv"),
            row.names = TRUE)

  # -------------------------------------------------------------------- #
  # 4.4  Create volcano plots highlighting calcium genes                 #
  # -------------------------------------------------------------------- #
  message("  • Creating volcano plots...")
  for (i in 1:ncol(contrasts)) {
    contrast_name <- colnames(contrasts)[i]

    # Get results for this contrast
    results <- limma::topTable(fit, coef = i, number = Inf)

    # Create volcano plot highlighting calcium genes
    volcano_plot <- create_standard_volcano(
      results,
      p_cutoff = 0.05,
      fc_cutoff = 1,  # Using lower FC threshold to capture more subtle changes
      label_method = "none",  # Don't label all significant genes
      highlight_gene = calcium_genes_present,  # Highlight calcium genes
      title = paste0(contrast_name, " - Calcium Genes")
    )

    # Save the plot
    ggsave(
      file.path(calcium_dir, paste0(contrast_name, "_calcium_volcano.pdf")),
      volcano_plot,
      width = 8, height = 7
    )
  }

  message("✓ Calcium gene analysis complete!")
  message("  Results saved to: ", calcium_dir)
} else {
  warning("No calcium genes found in dataset!")
}
