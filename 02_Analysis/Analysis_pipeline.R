#------------------------#
# 0. Prep env & Source scripts 
#------------------------#

source("01_Scripts/R_scripts/read_count_matrix.R")


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


# Create MDS plot to visualize sample relationships
pdf("03_Results/02_Analysis/Plots/General/PlotsMDS_plot.pdf")
plotMDS(DGE, col = as.numeric(DGE$samples$group), labels = paste(DGE$samples$genotype, DGE$samples$days, sep="_"))
legend("topright", legend = levels(DGE$samples$group), col = 1:length(levels(DGE$samples$group)), pch = 16)
dev.off()

# Create a heatmap of sample correlations
logCPM <- cpm(DGE, log = TRUE)
sample_cor <- cor(logCPM)
pdf("03_Results/01_Preprocessing/04_FeatureCounts/count_matrices_fc/sample_correlation_heatmap.pdf")
heatmap(sample_cor, main = "Sample Correlation")
dev.off()



#------------------------#
# 2. Design experiment matrix & conntrasts 
#------------------------#



#------------------------#
# 3. Run DE
#------------------------#


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


