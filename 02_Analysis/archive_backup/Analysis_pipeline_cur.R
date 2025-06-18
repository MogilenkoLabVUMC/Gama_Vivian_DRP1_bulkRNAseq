#------------------------#
# 0. Prep env & Source scripts 
#------------------------#

source("01_Scripts/R_scripts/read_count_matrix.R")
library(edgeR)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(limma)
library(ggplot2)
library(msigdbr)
library(fgsea)
library(WGCNA)

library(decoupleR)
library(clusterProfiler)
library(org.Hs.eg.db)


#------------------------#
# 1. Read data 
#------------------------#
library(edgeR)
library(dplyr)

# Source the function
source("process_rnaseq_data.R")

# File paths
counts_file <- "03_Results/01_Preprocessing/04_FeatureCounts/count_matrices_fc/sorted_counts_matrix.txt"
metadata_file <- "03_Results/01_Preprocessing/04_FeatureCounts/count_matrices_fc/metadata.csv"

# Preview count matrix
cat("Preview of count matrix:\n")
head(read.delim(counts_file, nrows = 5, check.names = FALSE))

# Preview metadata
cat("\nPreview of metadata:\n")
head(read.csv(metadata_file, sep = ";", nrows = 5))

# Process the data
DGE <- process_rnaseq_data(counts_file, metadata_file, annotate = FALSE)
head(DGE$counts)

# Examine the DGE object
print(DGE)
print(head(DGE$counts))
print(DGE$samples)


# Calculate log CPM values
logCPM <- cpm(DGE, log = TRUE)

# Calculate sample correlations
sample_cor <- cor(logCPM)

# Create a custom sample order
# First, separate D35 and D65 samples
D35_samples <- rownames(DGE$samples)[DGE$samples$days == "D35"]
D65_samples <- rownames(DGE$samples)[DGE$samples$days == "D65"]

# Within each day, order by genotype: Control, G32A, R403C
D35_Control <- D35_samples[DGE$samples[D35_samples, "genotype"] == "Control"]
D35_G32A <- D35_samples[DGE$samples[D35_samples, "genotype"] == "G32A"]
D35_R403C <- D35_samples[DGE$samples[D35_samples, "genotype"] == "R403C"]

D65_Control <- D65_samples[DGE$samples[D65_samples, "genotype"] == "Control"]
D65_G32A <- D65_samples[DGE$samples[D65_samples, "genotype"] == "G32A"]
D65_R403C <- D65_samples[DGE$samples[D65_samples, "genotype"] == "R403C"]

# Combine in the desired order
ordered_samples <- c(D35_Control, D35_G32A, D35_R403C, D65_Control, D65_G32A, D65_R403C)

# Reorder the correlation matrix
ordered_cor <- sample_cor[ordered_samples, ordered_samples]

# Create annotation dataframe for the heatmap
annotation_df <- data.frame(
  Genotype = DGE$samples$genotype[match(ordered_samples, rownames(DGE$samples))],
  Days = DGE$samples$days[match(ordered_samples, rownames(DGE$samples))],
  row.names = ordered_samples
)

# Define colors for annotations - using a harmonious, color-blind friendly palette
# Using a more harmonious palette from ColorBrewer
ann_colors <- list(
  Genotype = c(Control = "#1B9E77", G32A = "#D95F02", R403C = "#7570B3"),
  Days = c(D35 = "#E7298A", D65 = "#66A61E")
)

# Create a harmonious color palette for the heatmap that works with the annotation colors
# Using YlOrBr which is color-blind friendly and works well with the annotation colors
heatmap_colors <- colorRampPalette(brewer.pal(9, "YlOrBr"))(100)

# Create a prettier heatmap using pheatmap
pdf("03_Results/02_Analysis/Plots/General/sample_correlation_heatmap_ordered.pdf", width = 12, height = 10)
pheatmap(
  ordered_cor,
  main = "Sample Correlation",
  annotation_row = annotation_df,
  annotation_col = annotation_df,
  annotation_colors = ann_colors,
  color = heatmap_colors,
  border_color = NA,
  fontsize = 10,         # Increased from 8
  fontsize_row = 9,      # Increased from 7
  fontsize_col = 9,      # Increased from 7
  cluster_rows = FALSE,
  cluster_cols = FALSE
)
dev.off()

# Create an MDS plot with the same color scheme
pdf("03_Results/02_Analysis/Plots/General/MDS_plot_ordered.pdf", width = 10, height = 8)
# Create a color vector based on genotype and days - using the same colors as in the heatmap
genotype_colors <- c(Control = "#1B9E77", G32A = "#D95F02", R403C = "#7570B3")
days_shapes <- c(D35 = 16, D65 = 17)  # Different shapes for different days

# Get colors and shapes for each sample
sample_colors <- genotype_colors[DGE$samples$genotype]
sample_shapes <- days_shapes[DGE$samples$days]

# Plot MDS
plotMDS(
  DGE,
  col = sample_colors,
  pch = sample_shapes,
  cex = 1.5,        # Increased from 1.2
  main = "MDS Plot of RNA-seq Samples"
)

# Add a legend with larger text
legend(
  "topright",
  legend = c(
    paste("Control", c("D35", "D65")),
    paste("G32A", c("D35", "D65")),
    paste("R403C", c("D35", "D65"))
  ),
  col = rep(genotype_colors, each = 2),
  pch = rep(days_shapes, 3),
  cex = 1.0,        # Increased from 0.8
  pt.cex = 1.5,     # Increased from 1.2
  bty = "n"
)
dev.off()


#------------------------#
# 2. Design experiment matrix & contrasts 
#------------------------#
design <- model.matrix(~ 0 + group, data = DGE$samples) 
colnames(design) <- levels(DGE$samples$group)

# Define contrasts of interest
# 1. Mutation effects at each time point
# Question: What are the effects of G32A and R403C mutations at D35 and D65?
contrasts <- makeContrasts(
  # Mutation vs Control at D35
  G32A_vs_Ctrl_D35 = D35_G32A - D35_Control,
  R403C_vs_Ctrl_D35 = D35_R403C - D35_Control,
  
  # Mutation vs Control at D65
  G32A_vs_Ctrl_D65 = D65_G32A - D65_Control,
  R403C_vs_Ctrl_D65 = D65_R403C - D65_Control,

  # Mutation D65 vs Mutation D35
  G32A_vs_G32A_D65_vs_D35 = D65_G32A - D35_G32A,
  R403C_vs_R403C_D65_vs_D35 = D65_R403C - D35_R403C,
  
  # 2. Maturation effects within each genotype
  # Question: How does maturation affect each genotype?
  Time_Ctrl = D65_Control - D35_Control,
  Time_G32A = D65_G32A - D35_G32A,
  Time_R403C = D65_R403C - D35_R403C,
  Time_Mutations = (D65_G32A + D65_R403C)/2 - (D35_G32A + D35_R403C)/2,
  
  # 3. Interaction effects (difference-in-difference)
  # Question: Do the mutations alter the normal maturation trajectory?
  Int_G32A = (D65_G32A - D35_G32A) - (D65_Control - D35_Control),
  Int_R403C = (D65_R403C - D35_R403C) - (D65_Control - D35_Control),

  # 4. General mutation effect over time (combined mutations vs control)
  # Question: What is the overall effect of mutations on maturation?
  Int_All_Mutations = ((D65_G32A + D65_R403C)/2 - (D35_G32A + D35_R403C)/2) - (D65_Control - D35_Control),
  
  levels = design
)


# List of genes to highlight in volcano plots (from collaborator)
genes_of_interest <- c(
  "NEURONATIN", "CACNG3", "CACNA1C", "CACNA1S", "ATP2A1", 
  "RYR1", "MYLK3", "CASR", "VDR", "STIM1", "STIM2", 
  "ORAI1", "CALB1", "CALR"
)

# Continuing from your previous script sections...

#------------------------#
# 3. Run DE
#------------------------#

# Fitting moderated t-statistics model
fit <- edgeR::voomLmFit(
  DGE, design,
  sample.weights = FALSE)
fit <- contrasts.fit(fit, contrasts)
fit <- eBayes(fit, robust = TRUE)

# Create output directories
de_dir <- "03_Results/02_Analysis/DE_results"
volcano_dir <- "03_Results/02_Analysis/Plots/Volcano"
gsea_dir <- "03_Results/02_Analysis/Plots/GSEA"
heatmap_dir <- "03_Results/02_Analysis/Plots/Heatmaps"

# Create directories if they don't exist
dir.create(de_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(volcano_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(gsea_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(heatmap_dir, recursive = TRUE, showWarnings = FALSE)

#------------------------#
# 4. Visualise DEGs
#------------------------#

# Source necessary scripts for visualization
source("01_Scripts/GSEA_module/scripts/DE/plot_standard_volcano.R")
source("01_Scripts/GSEA_module/scripts/DE/analyzePathVolcanoViz.R")
source("01_Scripts/GSEA_module/scripts/DE/plotPCA.R")
source("01_Scripts/GSEA_module/scripts/DE/volcano_helpers.R")

# Create PCA plot for all samples
pdf("03_Results/02_Analysis/Plots/General/PCA_plot_all.pdf", width = 10, height = 8)
create_pca_plot(DGE, title = "PCA of All Samples")
dev.off()
# Warning messages:
# 1: In create_pca_plot(DGE, title = "PCA of All Samples") :
#   Expected columns 'group' and 'organ' not found in DGE_object$samples. Plotting may fail or look incorrect.
# 2: In create_pca_plot(DGE, title = "PCA of All Samples") :
#   Found unexpected organ types: Unknown. Assigning default shape (circle).

# Count significant DE genes per contrast
# Using decideTests to determine significant genes with FDR < 0.05
de_results <- decideTests(fit, p.value = 0.05)

# Create a data frame of DE gene counts
de_counts <- data.frame(
  Contrast = colnames(contrasts),
  Total = colSums(abs(de_results) > 0),
  Up = colSums(de_results > 0),
  Down = colSums(de_results < 0)
)

# Save the DE counts
write.csv(de_counts, file = file.path(de_dir, "DE_gene_counts.csv"), row.names = FALSE)

# Create a bar plot of DE gene counts
de_counts_long <- reshape2::melt(
  de_counts, 
  id.vars = "Contrast", 
  measure.vars = c("Up", "Down"),
  variable.name = "Direction",
  value.name = "Count"
)

pdf("03_Results/02_Analysis/Plots/General/DE_gene_counts.pdf", width = 10, height = 6)
ggplot(de_counts_long, aes(x = Contrast, y = Count, fill = Direction)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = c("Up" = "#D95F02", "Down" = "#1B9E77")) +
  theme_minimal(base_size = 12) +
  labs(title = "Number of DE genes (FDR < 0.05)",
       y = "Number of genes",
       x = "Contrast") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# Create volcano plots for each contrast
for (i in 1:ncol(contrasts)) {
  contrast_name <- colnames(contrasts)[i]
  dir.create(file.path(volcano_dir, contrast_name), showWarnings = FALSE)
  
  # Get results table for this contrast
  results <- topTable(fit, coef = i, number = Inf)
  
  # Save results to file
  write.csv(results, file = file.path(de_dir, paste0(contrast_name, "_results.csv")))
  
  # Create standard volcano plot
  volcano_plot <- create_standard_volcano(
    results,
    p_cutoff = 0.05,
    fc_cutoff = 2,
    label_method = "sig",
    highlight_gene = genes_of_interest,
    title = contrast_name
  )
  
  # Save the plot
  ggsave(
    file.path(volcano_dir, contrast_name, paste0(contrast_name, "_volcano.pdf")),
    volcano_plot,
    width = 8, height = 7
  )
  
  # Create vertical volcano plot (for multi-panel display)
  vert_volcano <- create_vertical_volcano(
    results,
    fc_cutoff = 2,
    p_cutoff = 0.05,
    title = contrast_name,
    label_method = "sig",
    max.overlaps = 10
  )
  
  # Save to variable for later combination
  if (i == 1) {
    volcano_list <- list()
  }
  volcano_list[[contrast_name]] <- vert_volcano
}

# Combine all vertical volcanos into one figure
if (length(volcano_list) > 0) {
  # For all main contrasts
  main_contrasts <- c("G32A_vs_Ctrl_D35", "R403C_vs_Ctrl_D35", "G32A_vs_Ctrl_D65", "R403C_vs_Ctrl_D65")
  if (all(main_contrasts %in% names(volcano_list))) {
    combined_volcano <- combine_volcano_row(volcano_list[main_contrasts], labels = main_contrasts)
    ggsave(
      file.path(volcano_dir, "Combined_mutation_effects_volcano.pdf"),
      combined_volcano,
      width = 16, height = 6
    )
  }
  
  # For all time contrasts
  time_contrasts <- c("Time_Ctrl", "Time_G32A", "Time_R403C", "Time_Mutations")
  if (all(time_contrasts %in% names(volcano_list))) {
    combined_time_volcano <- combine_volcano_row(volcano_list[time_contrasts], labels = time_contrasts)
    ggsave(
      file.path(volcano_dir, "Combined_time_effects_volcano.pdf"),
      combined_time_volcano,
      width = 12, height = 6
    )
  }
  
  # For interaction contrasts
  int_contrasts <- c("Int_G32A", "Int_R403C", "Int_All_Mutations")
  if (all(int_contrasts %in% names(volcano_list))) {
    combined_int_volcano <- combine_volcano_row(volcano_list[int_contrasts], labels = int_contrasts)
    ggsave(
      file.path(volcano_dir, "Combined_interaction_effects_volcano.pdf"),
      combined_int_volcano,
      width = 8, height = 6
    )
  }
}

# Compare DE results between contrasts (Venn diagrams)
# Install if needed: if (!requireNamespace("VennDiagram", quietly = TRUE)) install.packages("VennDiagram")
library(VennDiagram) #need to install
library(grid)

# Create pairwise comparisons for mutation effects
pdf(file.path(volcano_dir, "Venn_diagrams_mutation_effects.pdf"), width = 10, height = 8)
par(mfrow = c(2, 1))

# D35 comparison (G32A vs R403C)
d35_contrasts <- c("G32A_vs_Ctrl_D35", "R403C_vs_Ctrl_D35")
if (all(d35_contrasts %in% colnames(contrasts))) {
  # Get significant genes for each contrast
  sig_g32a_d35 <- rownames(fit$coefficients)[which(de_results[, "G32A_vs_Ctrl_D35"] != 0)]
  sig_r403c_d35 <- rownames(fit$coefficients)[which(de_results[, "R403C_vs_Ctrl_D35"] != 0)]
  
  # Create Venn diagram
  venn_list <- list(`G32A (D35)` = sig_g32a_d35, `R403C (D35)` = sig_r403c_d35)
  venn_colors <- c("#D95F02", "#7570B3")
  
  venn_plot_d35 <- venn.diagram(
    x = venn_list,
    category.names = names(venn_list),
    filename = NULL,
    output = TRUE,
    fill = venn_colors,
    main = "DE genes at D35"
  )
  grid.draw(venn_plot_d35)
}

# D65 comparison (G32A vs R403C)
d65_contrasts <- c("G32A_vs_Ctrl_D65", "R403C_vs_Ctrl_D65")
if (all(d65_contrasts %in% colnames(contrasts))) {
  # Get significant genes for each contrast
  sig_g32a_d65 <- rownames(fit$coefficients)[which(de_results[, "G32A_vs_Ctrl_D65"] != 0)]
  sig_r403c_d65 <- rownames(fit$coefficients)[which(de_results[, "R403C_vs_Ctrl_D65"] != 0)]
  
  # Create Venn diagram
  venn_list <- list(`G32A (D65)` = sig_g32a_d65, `R403C (D65)` = sig_r403c_d65)
  venn_colors <- c("#D95F02", "#7570B3")
  
  venn_plot_d65 <- venn.diagram(
    x = venn_list,
    category.names = names(venn_list),
    filename = NULL,
    output = TRUE,
    fill = venn_colors,
    main = "DE genes at D65"
  )
  grid.draw(venn_plot_d65)
}
dev.off()

# Create UpSet plot to compare all contrasts
# Install if needed: if (!requireNamespace("UpSetR", quietly = TRUE)) install.packages("UpSetR")
library(UpSetR) # need to install

# Create a list of significant genes for each contrast
sig_genes_list <- list()
for (contrast in colnames(contrasts)) {
  sig_genes_list[[contrast]] <- rownames(fit$coefficients)[which(de_results[, contrast] != 0)]
}

# Create UpSet plot
pdf(file.path(volcano_dir, "UpSet_plot_all_contrasts.pdf"), width = 12, height = 8)
upset(fromList(sig_genes_list), order.by = "freq", nsets = length(sig_genes_list))
dev.off()
# for some reason, saves a pdf with first page empty, second page with the correct plot 


#------------------------#
# 5. Run generic GSEA
#------------------------#

# Source necessary scripts for GSEA
source("01_Scripts/GSEA_module/scripts/GSEA/GSEA_processing/run_gsea.R")
source("01_Scripts/GSEA_module/scripts/GSEA/GSEA_processing/run_gsea_analysis.R")


# Run GSEA for each contrast # NOT WORKING AS EXPECTED
for (i in 1:ncol(contrasts)) {
  contrast_name <- colnames(contrasts)[i]
  
  # Get results table for this contrast
  results <- topTable(fit, coef = i, number = Inf)
  
  # Create sample annotation for heatmaps
  sample_annotation <- data.frame(
    Genotype = DGE$samples$genotype,
    Days = DGE$samples$days,
    row.names = rownames(DGE$samples)
  )
  
  # Run GSEA analysis with multiple databases
  run_gsea_analysis(
    de_table = results,
    analysis_name = contrast_name,
    rank_metric = "t",
    species = "Homo sapiens",  # Assuming human data
    n_pathways = 30,
    padj_cutoff = 0.05,
    output_dir = file.path(gsea_dir, contrast_name),
    save_plots = TRUE,
    sample_annotation = sample_annotation,
    sample_order = ordered_samples,
    helper_root = "01_Scripts/GSEA_module"
  )
}
# NOT WORKING AS EXPECTED. The output in the comment below
# [DEBUG] sourcing 01_Scripts/GSEA_module/R_GSEA_visualisations/scripts/custom_minimal_theme.R
#[run_gsea_analysis] helper not found → 01_Scripts/GSEA_module/R_GSEA_visualisations/scripts/custom_minimal_theme.R
#[DEBUG] sourcing 01_Scripts/GSEA_module/R_GSEA_visualisations/scripts/GSEA/GSEA_plotting/gsea_plotting_utils.R
#[run_gsea_analysis] helper not found → 01_Scripts/GSEA_module/R_GSEA_visualisations/scripts/GSEA/GSEA_plotting/gsea_plotting_utils.R
#[DEBUG] sourcing 01_Scripts/GSEA_module/R_GSEA_visualisations/scripts/GSEA/GSEA_plotting/gsea_dotplot.R
#[run_gsea_analysis] helper not found → 01_Scripts/GSEA_module/R_GSEA_visualisations/scripts/GSEA/GSEA_plotting/gsea_dotplot.R
#[DEBUG] sourcing 01_Scripts/GSEA_module/R_GSEA_visualisations/scripts/GSEA/GSEA_plotting/gsea_dotplot_facet.R
#[run_gsea_analysis] helper not found → 01_Scripts/GSEA_module/R_GSEA_visualisations/scripts/GSEA/GSEA_plotting/gsea_dotplot_facet.R
#[DEBUG] sourcing 01_Scripts/GSEA_module/R_GSEA_visualisations/scripts/GSEA/GSEA_plotting/gsea_barplot.R
#[run_gsea_analysis] helper not found → 01_Scripts/GSEA_module/R_GSEA_visualisations/scripts/GSEA/GSEA_plotting/gsea_barplot.R
#[DEBUG] sourcing 01_Scripts/GSEA_module/R_GSEA_visualisations/scripts/GSEA/GSEA_plotting/gsea_running_sum_plot.R
#[run_gsea_analysis] helper not found → 01_Scripts/GSEA_module/R_GSEA_visualisations/scripts/GSEA/GSEA_plotting/gsea_running_sum_plot.R
#[DEBUG] sourcing 01_Scripts/GSEA_module/R_GSEA_visualisations/scripts/GSEA/GSEA_plotting/gsea_heatmap.R
#[run_gsea_analysis] helper not found → 01_Scripts/GSEA_module/R_GSEA_visualisations/scripts/GSEA/GSEA_plotting/gsea_heatmap.R
#[DEBUG] sourcing 01_Scripts/GSEA_module/R_GSEA_visualisations/scripts/GSEA/GSEA_processing/run_gsea.R
#[run_gsea_analysis] helper not found → 01_Scripts/GSEA_module/R_GSEA_visualisations/scripts/GSEA/GSEA_processing/run_gsea.R
#[DEBUG] sourcing 01_Scripts/GSEA_module/R_GSEA_visualisations/scripts/DE/volcano_helpers.R
#[run_gsea_analysis] helper not found → 01_Scripts/GSEA_module/R_GSEA_visualisations/scripts/DE/volcano_helpers.R
#▶  Hallmark
#Fetching MSigDB sets for species='Homo sapiens', db_species='MM', collection='H', subcollection=''
#Running clusterProfiler::GSEA...
#GSEA completed successfully.
#Error in get_db_plot_params(db_name) : 
#  could not find function "get_db_plot_params"
# NOTA BENE! The data is human, but the script uses db_species='MM' for some reason.
# NOTA BENE! I also need to include the SynGo database for GSEA testing
# (python3.10) root@b24371903d2b:/workspaces/GVDRP1/00_Data/SynGO_bulk_20231201# ll
# total 508
# drwx------ 2 mogilenko_lab mogilenko_lab   4096 May 16 11:03  ./
# drwxrwxr-x 6 mogilenko_lab mogilenko_lab   4096 May 16 11:03  ../
# -rw-r--r-- 1 mogilenko_lab mogilenko_lab    165 May 16 11:03 '~$syngo_annotations.xlsx'
# -rw-r--r-- 1 mogilenko_lab mogilenko_lab    165 May 16 11:03 '~$syngo_ontologies.xlsx'
# -rw-r--r-- 1 mogilenko_lab mogilenko_lab 248715 May 16 11:03  syngo_annotations.xlsx
# -rw-r--r-- 1 mogilenko_lab mogilenko_lab 102741 May 16 11:03  syngo_genes.xlsx
# -rw-r--r-- 1 mogilenko_lab mogilenko_lab   4845 May 16 11:03  SynGO_lookup_codes.json
# -rw-r--r-- 1 mogilenko_lab mogilenko_lab 136873 May 16 11:03  syngo_ontologies.xlsx
# (python3.10) root@b24371903d2b:/workspaces/GVDRP1/00_Data/SynGO_bulk_20231201# 

#------------------------#
# 8. Analyze genes of interest (calcium signaling)
#------------------------#

# Focus on calcium signaling genes from the collaborator
calcium_genes <- c(
  "NEURONATIN", "CACNG3", "CACNA1C", "CACNA1S", "ATP2A1", 
  "RYR1", "MYLK3", "CASR", "VDR", "STIM1", "STIM2", 
  "ORAI1", "CALB1", "CALR"
)

# Create output directory for calcium gene analysis
calcium_dir <- "03_Results/02_Analysis/Calcium_genes"
dir.create(calcium_dir, recursive = TRUE, showWarnings = FALSE)

# Extract expression data for calcium genes that are present in our dataset
calcium_genes_present <- calcium_genes[calcium_genes %in% rownames(logCPM)]

if (length(calcium_genes_present) > 0) {
  # Extract expression values
  calcium_expr <- logCPM[calcium_genes_present, ordered_samples]
  
  # Create heatmap of expression across all samples
  pdf(file.path(calcium_dir, "calcium_genes_expression_heatmap.pdf"), width = 12, height = length(calcium_genes_present) * 0.4 + 3)
  pheatmap(
    calcium_expr,
    main = "Expression of Calcium Signaling Genes",
    annotation_col = annot,
    annotation_colors = ann_colors,
    color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
    border_color = NA,
    fontsize = 10,
    fontsize_row = 10,
    fontsize_col = 9,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    scale = "row"  # Scale by row for better visualization
  )
  dev.off()
  
  # Create boxplots to compare expression across genotypes and timepoints
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
  
  # Test for differential expression of calcium genes in each contrast
  calcium_gene_results <- data.frame()
  
  for (i in 1:ncol(contrasts)) {
    contrast_name <- colnames(contrasts)[i]
    
    # Get results for this contrast
    results <- topTable(fit, coef = i, number = Inf)
    
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
  write.csv(calcium_gene_results, file = file.path(calcium_dir, "calcium_genes_DE_results.csv"), row.names = TRUE)
  
  # Create volcano plots highlighting calcium genes for each contrast
  for (i in 1:ncol(contrasts)) {
    contrast_name <- colnames(contrasts)[i]
    
    # Get results for this contrast
    results <- topTable(fit, coef = i, number = Inf)
    
    # Create volcano plot highlighting calcium genes
    volcano_plot <- create_standard_volcano(
      results,
      p_cutoff = 0.05,
      fc_cutoff = 1,  # Using lower FC threshold to capture more subtle changes
      label_method = "none",  # Don't label all significant genes
      highlight_gene = calcium_genes_present,  # Highlight calcium genes
      # the problem is that the function highlights the gene names with the color of the dots, making text unreadable against the background of same color dots 
      title = paste0(contrast_name, " - Calcium Genes")
    )
    
    # Save the plot
    ggsave(
      file.path(calcium_dir, paste0(contrast_name, "_calcium_volcano.pdf")),
      volcano_plot,
      width = 8, height = 7
    )
  }
}





  

# finished here

#------------------------#
# 10. SynGO GSEA
#------------------------#

# This section requires downloading and processing SynGO gene sets.
# The general approach would be similar to the GSEA section, but using SynGO-specific gene sets.

# Placeholder for SynGO analysis - would require downloading and formatting SynGO database
syngo_dir <- "03_Results/02_Analysis/SynGO"
dir.create(syngo_dir, recursive = TRUE, showWarnings = FALSE)

# Note: Implementation would depend on how SynGO data is structured and available.
# Typically, you would:
# 1. Download SynGO gene sets
# 2. Format them for use with GSEA or overrepresentation analysis
# 3. Run analysis similar to GSEA section
# 4. Create visualizations specific to SynGO (focusing on pre- vs post-synaptic enrichment)

# Once SynGO gene sets are available, code similar to GSEA section can be used

#------------------------#
# 11. Summary of Findings
#------------------------#

# Create a summary of key findings
summary_dir <- "03_Results/02_Analysis/Summary"
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

# Count significant DEGs
de_summary <- data.frame(
  Contrast = colnames(contrasts),
  Total_DEGs = colSums(abs(de_results) > 0),
  Up_regulated = colSums(de_results > 0),
  Down_regulated = colSums(de_results < 0)
)

# Add percentages
de_summary$Up_percent <- round(de_summary$Up_regulated / de_summary$Total_DEGs * 100, 1)
de_summary$Down_percent <- round(de_summary$Down_regulated / de_summary$Total_DEGs * 100, 1)

# Save summary to file
write.csv(de_summary, file = file.path(summary_dir, "DE_summary.csv"), row.names = FALSE)

# Create a bar plot of DEG counts
pdf(file.path(summary_dir, "DE_summary_barplot.pdf"), width = 10, height = 6)
de_summary_long <- reshape2::melt(
  de_summary,
  id.vars = "Contrast",
  measure.vars = c("Up_regulated", "Down_regulated"),
  variable.name = "Direction",
  value.name = "Count"
)

ggplot(de_summary_long, aes(x = Contrast, y = Count, fill = Direction)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = c("Up_regulated" = "#D95F02", "Down_regulated" = "#1B9E77")) +
  geom_text(aes(label = Count), position = position_dodge(width = 0.9), vjust = -0.5) +
  labs(
    title = "Differential Expression Summary",
    x = "Contrast",
    y = "Number of DEGs (FDR < 0.05)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )
dev.off()

# Summarize GSEA results
# Create a table of top enriched pathways for each contrast

# Collect top pathways from each contrast and database
top_pathways <- data.frame()

# List all GSEA result directories
gsea_result_dirs <- list.dirs(gsea_dir, recursive = FALSE)

for (contrast_dir in gsea_result_dirs) {
  contrast_name <- basename(contrast_dir)
  
  # List all database subdirectories
  db_dirs <- list.dirs(contrast_dir, recursive = FALSE)
  
  for (db_dir in db_dirs) {
    db_name <- basename(db_dir)
    
    # Check for hallmark results file
    csv_files <- list.files(db_dir, pattern = ".*_results\\.csv$", full.names = TRUE)
    
    for (csv_file in csv_files) {
      if (file.exists(csv_file)) {
        # Read results
        results <- read.csv(csv_file)
        
        # Get top 5 pathways
        if (nrow(results) > 0) {
          top5 <- results[order(results$pvalue), ][1:min(5, nrow(results)), ]
          
          # Add contrast and database info
          top5$Contrast <- contrast_name
          top5$Database <- db_name
          
          # Combine with full results
          top_pathways <- rbind(top_pathways, top5)
        }
      }
    }
  }
}

# Save top pathways to file
if (nrow(top_pathways) > 0) {
  write.csv(top_pathways, file = file.path(summary_dir, "top_pathways_summary.csv"), row.names = FALSE)
}

# Create an overall summary for results interpretation
cat("
## RNA-seq Analysis Summary

### Overview
- Total samples analyzed: ", nrow(DGE$samples), "
- Genotypes: ", paste(unique(DGE$samples$genotype), collapse = ", "), "
- Time points: ", paste(unique(DGE$samples$days), collapse = ", "), "

### Differential Expression Results
", file = file.path(summary_dir, "analysis_summary.md"))

# Add DE summary
cat("
#### Number of differentially expressed genes (FDR < 0.05):
", file = file.path(summary_dir, "analysis_summary.md"), append = TRUE)

for (i in 1:nrow(de_summary)) {
  cat(sprintf("- %s: %d total DEGs (%d up, %d down)\n", 
              de_summary$Contrast[i], 
              de_summary$Total_DEGs[i], 
              de_summary$Up_regulated[i], 
              de_summary$Down_regulated[i]),
      file = file.path(summary_dir, "analysis_summary.md"), append = TRUE)
}

# Add pathway enrichment summary
if (nrow(top_pathways) > 0) {
  cat("\n### Top Enriched Pathways\n", file = file.path(summary_dir, "analysis_summary.md"), append = TRUE)
  
  for (contrast in unique(top_pathways$Contrast)) {
    cat(sprintf("\n#### %s\n", contrast), file = file.path(summary_dir, "analysis_summary.md"), append = TRUE)
    
    for (db in unique(top_pathways$Database[top_pathways$Contrast == contrast])) {
      cat(sprintf("**%s database:**\n", db), file = file.path(summary_dir, "analysis_summary.md"), append = TRUE)
      
      # Get top pathways for this contrast and database
      top <- top_pathways[top_pathways$Contrast == contrast & top_pathways$Database == db, ]
      top <- top[order(top$pvalue), ][1:min(3, nrow(top)), ]
      
      for (i in 1:nrow(top)) {
        cat(sprintf("- %s (p-value = %.4f, NES = %.2f)\n", 
                    top$Description[i], 
                    top$pvalue[i], 
                    top$NES[i]),
            file = file.path(summary_dir, "analysis_summary.md"), append = TRUE)
      }
    }
  }
}

# Add calcium genes summary
if (length(calcium_genes_present) > 0) {
  cat("\n### Calcium Signaling Genes Analysis\n", file = file.path(summary_dir, "analysis_summary.md"), append = TRUE)
  
  # Get calcium gene results
  calcium_file <- file.path(calcium_dir, "calcium_genes_DE_results.csv")
  if (file.exists(calcium_file)) {
    calcium_results <- read.csv(calcium_file, row.names = 1)
    
    # Summarize significantly DE calcium genes
    sig_calcium <- calcium_results[calcium_results$adj.P.Val < 0.05, ]
    
    if (nrow(sig_calcium) > 0) {
      cat("Differentially expressed calcium signaling genes (FDR < 0.05):\n", file = file.path(summary_dir, "analysis_summary.md"), append = TRUE)
      
      for (contrast in unique(sig_calcium$Contrast)) {
        cat(sprintf("\n#### %s\n", contrast), file = file.path(summary_dir, "analysis_summary.md"), append = TRUE)
        
        # Get significant genes for this contrast
        sig_genes <- sig_calcium[sig_calcium$Contrast == contrast, ]
        
        for (i in 1:nrow(sig_genes)) {
          cat(sprintf("- %s (log2FC = %.2f, adj. p-value = %.4f)\n", 
                      rownames(sig_genes)[i], 
                      sig_genes$logFC[i], 
                      sig_genes$adj.P.Val[i]),
              file = file.path(summary_dir, "analysis_summary.md"), append = TRUE)
        }
      }
    } else {
      cat("No calcium signaling genes were significantly differentially expressed.\n", file = file.path(summary_dir, "analysis_summary.md"), append = TRUE)
    }
  }
}

# Generate HTML report
if (requireNamespace("rmarkdown", quietly = TRUE)) {
  # Create a temporary Rmd file
  rmd_file <- file.path(summary_dir, "analysis_report.Rmd")
  
  # Write YAML header and include the markdown file
  cat("---
title: \"RNA-seq Analysis Report\"
date: \"`r format(Sys.time(), '%d %B, %Y')`\"
output: 
  html_document:
    toc: true
    toc_float: true
    theme: cosmo
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)
