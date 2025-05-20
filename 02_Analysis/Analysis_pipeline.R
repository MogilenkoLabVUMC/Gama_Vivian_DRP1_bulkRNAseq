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
  
  # 3. Interaction effects (difference-in-difference)
  # Question: Do the mutations alter the normal maturation trajectory?
  Int_G32A = (D65_G32A - D35_G32A) - (D65_Control - D35_Control),
  Int_R403C = (D65_R403C - D35_R403C) - (D65_Control - D35_Control),
  
  levels = design
)


# List of genes to highlight in volcano plots (from collaborator)
genes_of_interest <- c(
  "NEURONATIN", "CACNG3", "CACNA1C", "CACNA1S", "ATP2A1", 
  "RYR1", "MYLK3", "CASR", "VDR", "STIM1", "STIM2", 
  "ORAI1", "CALB1", "CALR"
)


#------------------------#
# 3. Run DE
#------------------------#


# Fitting moderated t-statistics model
fit <- edgeR::voomLmFit(
  DGE, design,
  sample.weights = FALSE)
fit <- contrasts.fit(fit, contrasts)
fit <- eBayes(fit, robust = TRUE)

#------------------------#
# 4. Visualise DEGs
#------------------------#




#------------------------#
# 5. Run generic GSEA
#------------------------#


#------------------------#
# 6. Run SynGO GSEA
#------------------------#


#------------------------#
# 7. Run WGCA
#------------------------#


