library(edgeR)
library(dplyr)

process_rnaseq_data <- function(counts_file, metadata_file, annotate = FALSE) {
  # Read count matrix with proper header handling
  counts <- read.delim(counts_file, check.names = FALSE, stringsAsFactors = FALSE, header = TRUE)
  
  # Set gene IDs as rownames and remove the Geneid column
  rownames(counts) <- counts$Geneid
  counts <- counts[, -which(colnames(counts) == "Geneid"), drop = FALSE]
  
  # Convert counts to numeric matrix
  counts <- as.matrix(counts)
  mode(counts) <- "numeric"
  
  # Read metadata
  metadata <- read.csv(metadata_file, sep = ";", stringsAsFactors = FALSE)
  
  # Make sample names as rownames in metadata
  rownames(metadata) <- metadata$sample
  
  # Ensure counts and metadata have matching samples
  # First, clean up sample names in counts to match metadata
  count_sample_names <- gsub("_Aligned$", "", colnames(counts))
  colnames(counts) <- count_sample_names
  
  # Find common samples
  common_samples <- intersect(colnames(counts), rownames(metadata))
  
  if(length(common_samples) == 0) {
    stop("No matching sample names found between count matrix and metadata")
  }
  
  counts <- counts[, common_samples, drop = FALSE]
  metadata <- metadata[common_samples, , drop = FALSE]
  
  # Print some debugging info
  cat("Number of samples in count matrix:", ncol(counts), "\n")
  cat("Number of genes in count matrix:", nrow(counts), "\n")
  cat("First few sample names:", head(colnames(counts)), "\n")
  
  # Create DGEList object
  DGE <- DGEList(counts = counts)
  
  # Add sample information
  DGE$samples$genotype <- factor(metadata$genotype)
  DGE$samples$days <- factor(metadata$days)
  DGE$samples$rep <- factor(metadata$rep)
  DGE$samples$cell_line <- factor(metadata$cell_line)
  DGE$samples$cell_type <- factor(metadata$cell_type)
  
  # Create a group factor combining genotype and days
  DGE$samples$group <- factor(paste(metadata$days, metadata$genotype, sep="_"))
  
  # Filter low expressed genes
  keep <- filterByExpr(DGE)
  DGE <- DGE[keep, , keep.lib.sizes = FALSE]
  
  # Normalize library sizes
  DGE <- calcNormFactors(DGE, method = "TMM")
  
  # Annotate genes if requested
  if (annotate) {
    require(biomaRt)
    ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    gene_ids <- rownames(DGE)
    
    annotations <- getBM(
      attributes = c("ensembl_gene_id", "external_gene_name", "description"),
      filters = "ensembl_gene_id",
      values = gene_ids,
      mart = ensembl
    )
    
    # Keep only unique gene names and handle empty mappings
    mapping_table <- annotations %>%
      filter(external_gene_name != "") %>%
      distinct(external_gene_name, .keep_all = TRUE) %>%
      select(ensembl_gene_id, external_gene_name)
    
    # Match and rename genes, keeping original names if mapping fails
    matched_idx <- match(rownames(DGE), mapping_table$ensembl_gene_id)
    rownames(DGE) <- ifelse(!is.na(matched_idx) & mapping_table$external_gene_name[matched_idx] != "",
                            mapping_table$external_gene_name[matched_idx],
                            rownames(DGE))
  }
  
  return(DGE)
}
