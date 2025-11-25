###############################################################################
##  Mechanistic Cascade: Energy Crisis → Translation Failure                ##
##  Gene sets extracted from REAL significant GSEA pathways                 ##
###############################################################################
##  STORY: Mt Central Dogma ↑ → Mt Ribosomes ↑ → OXPHOS ↑ (compensation)   ##
##         BUT: Synaptic Ribosomes ↓ → Failed local translation → Ca2+     ##
###############################################################################

library(here)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(grid)

message("📂 Loading checkpoints...")
checkpoint_dir <- here("03_Results/02_Analysis/checkpoints")
all_gsea_results <- readRDS(file.path(checkpoint_dir, "all_gsea_results.rds"))
syngo_gsea_results <- readRDS(file.path(checkpoint_dir, "syngo_gsea_results.rds"))
mitocarta_gsea_results <- readRDS(file.path(checkpoint_dir, "mitocarta_gsea_results.rds"))
fit <- readRDS(file.path(checkpoint_dir, "fit_object.rds"))
de_results <- readRDS(file.path(checkpoint_dir, "de_results.rds"))

# Output directory
out_dir <- here("03_Results/02_Analysis/Plots/Mito_translation_cascade")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("✓ Checkpoints loaded\n")

###############################################################################
##  Configuration                                                            ##
###############################################################################

# Trajectory framework contrasts (matching cross-database validation)
all_contrasts <- c(
  "G32A_vs_Ctrl_D35", "Maturation_G32A_specific", "G32A_vs_Ctrl_D65",
  "R403C_vs_Ctrl_D35", "Maturation_R403C_specific", "R403C_vs_Ctrl_D65"
)

contrast_labels <- c(
  "G32A\nEarly", "G32A\nTrajDev", "G32A\nLate",
  "R403C\nEarly", "R403C\nTrajDev", "R403C\nLate"
)

# Max genes per module
MAX_GENES_PER_MODULE <- 20

###############################################################################
##  Helper Functions                                                         ##
###############################################################################

#' Extract core enrichment genes from a GSEA result
#' @param gsea_result A gseaResult object
#' @param pathway_name Pattern or exact name to match
#' @param exact_match If TRUE, use exact string match
#' @return Character vector of gene symbols
extract_pathway_genes <- function(gsea_result, pathway_name, exact_match = FALSE) {
  if (is.null(gsea_result)) {
    message(sprintf("    WARNING: GSEA result is NULL for '%s'", pathway_name))
    return(character(0))
  }

  df <- gsea_result@result

  if (exact_match) {
    idx <- df$Description == pathway_name
  } else {
    idx <- grepl(pathway_name, df$Description, ignore.case = TRUE)
  }

  if (!any(idx)) {
    message(sprintf("    WARNING: No pathway matching '%s'", pathway_name))
    return(character(0))
  }

  # Get the most significant matching pathway
  matched <- df[idx, ]
  best <- matched[which.min(matched$p.adjust), ]

  message(sprintf("    Found: '%s' (NES=%.2f, p.adj=%.2e, n=%d)",
                  best$Description, best$NES, best$p.adjust, best$setSize))

  # Extract core enrichment genes (forward-slash delimited)
  genes <- unlist(strsplit(best$core_enrichment, "/"))
  return(genes)
}

#' Limit gene set to top N by absolute logFC
#' @param genes Character vector of gene symbols
#' @param logfc_matrix Matrix with genes as rows, contrasts as columns
#' @param n Maximum number of genes to return
#' @return Character vector of top N genes
limit_genes <- function(genes, logfc_matrix, n = 20) {
  genes_present <- genes[genes %in% rownames(logfc_matrix)]

  if (length(genes_present) == 0) return(character(0))
  if (length(genes_present) <= n) return(genes_present)

  # Rank by mean absolute logFC across TrajDev contrasts
  trajdev_cols <- grep("TrajDev", colnames(logfc_matrix))
  if (length(trajdev_cols) == 0) {
    # Fallback: use all columns
    mean_abs_fc <- rowMeans(abs(logfc_matrix[genes_present, , drop = FALSE]), na.rm = TRUE)
  } else {
    mean_abs_fc <- rowMeans(abs(logfc_matrix[genes_present, trajdev_cols, drop = FALSE]), na.rm = TRUE)
  }

  top_genes <- names(sort(mean_abs_fc, decreasing = TRUE)[1:n])
  return(top_genes)
}

###############################################################################
##  Extract logFC Matrix for All Contrasts                                  ##
###############################################################################

message("📊 Extracting expression data for trajectory framework...\n")

logfc_matrix <- sapply(all_contrasts, function(contrast) {
  coef_idx <- which(colnames(fit$coefficients) == contrast)
  if (length(coef_idx) == 0) {
    warning(sprintf("Contrast '%s' not found in fit object", contrast))
    return(rep(NA, nrow(fit$coefficients)))
  }
  logfc <- fit$coefficients[, coef_idx]
  names(logfc) <- rownames(fit$coefficients)
  return(logfc)
})

colnames(logfc_matrix) <- contrast_labels
message(sprintf("  Extracted logFC for %d contrasts\n", ncol(logfc_matrix)))

###############################################################################
##  Extract Gene Sets from REAL Significant Pathways                        ##
###############################################################################

message("📊 Extracting gene sets from significant GSEA pathways...\n")

## MODULE 1: Mitochondrial Central Dogma (MitoCarta)
message("  Module 1: Mitochondrial Central Dogma (MitoCarta)...")
central_dogma_genes_raw <- extract_pathway_genes(
  mitocarta_gsea_results[["Maturation_G32A_specific"]],
  "Mitochondrial_central_dogma", exact_match = TRUE
)

## MODULE 2: Mitochondrial Ribosomes (MitoCarta)
message("  Module 2: Mitochondrial Ribosomes (MitoCarta)...")
mito_ribosome_genes_raw <- extract_pathway_genes(
  mitocarta_gsea_results[["Maturation_G32A_specific"]],
  "Mitochondrial_ribosome", exact_match = TRUE
)

## MODULE 3: ATP Synthase (GOCC) - Top hit from Fig 4
# Replaces generic OXPHOS
message("  Module 3: ATP Synthase (GOCC - Top hit from Fig 4)...")
oxphos_genes_raw <- extract_pathway_genes(
  all_gsea_results[["Maturation_G32A_specific"]][["gocc"]],
  "GOCC_ATPASE_COMPLEX", exact_match = TRUE
)

## MODULE 4 & 5: Synaptic Ribosomes (SynGO) - split into Common vs Postsynaptic-only
message("  Module 4 & 5: Synaptic Ribosomes (SynGO)...")
message("    Extracting presynaptic ribosome...")
presynaptic_ribo_genes <- extract_pathway_genes(
  syngo_gsea_results[["Maturation_G32A_specific"]],
  "presynaptic ribosome"
)

message("    Extracting postsynaptic ribosome...")
postsynaptic_ribo_genes <- extract_pathway_genes(
  syngo_gsea_results[["Maturation_G32A_specific"]],
  "postsynaptic ribosome"
)

# Split into Common (intersection) and Postsynaptic-only (setdiff)
synaptic_common_genes_raw <- intersect(presynaptic_ribo_genes, postsynaptic_ribo_genes)
postsynaptic_only_genes_raw <- setdiff(postsynaptic_ribo_genes, presynaptic_ribo_genes)

message(sprintf("    Common (pre ∩ post): %d genes", length(synaptic_common_genes_raw)))
message(sprintf("    Postsynaptic-only: %d genes", length(postsynaptic_only_genes_raw)))

## MODULE 6: Calcium Signaling (Reactome) - Top hit from Fig 4
# Replaces hardcoded config genes
# NOTE: Using G32A Early contrast where this pathway is most significant (NES=1.78)
message("  Module 6: Calcium Signaling (Reactome - Top hit from Fig 4)...")
calcium_genes_raw <- extract_pathway_genes(
  all_gsea_results[["G32A_vs_Ctrl_D35"]][["reactome"]],
  "REACTOME_ELEVATION_OF_CYTOSOLIC_CA2_LEVELS", exact_match = TRUE
)

###############################################################################
##  Limit Genes Per Module (Top 20 by |logFC|)                              ##
###############################################################################

message(sprintf("📊 Limiting to top %d genes per module by |logFC|...\n", MAX_GENES_PER_MODULE))

central_dogma_genes <- limit_genes(central_dogma_genes_raw, logfc_matrix, MAX_GENES_PER_MODULE)
mito_ribosome_genes <- limit_genes(mito_ribosome_genes_raw, logfc_matrix, MAX_GENES_PER_MODULE)
oxphos_genes <- limit_genes(oxphos_genes_raw, logfc_matrix, MAX_GENES_PER_MODULE)
synaptic_common_genes <- limit_genes(synaptic_common_genes_raw, logfc_matrix, MAX_GENES_PER_MODULE)
postsynaptic_only_genes <- limit_genes(postsynaptic_only_genes_raw, logfc_matrix, MAX_GENES_PER_MODULE)
calcium_genes <- limit_genes(calcium_genes_raw, logfc_matrix, MAX_GENES_PER_MODULE)

message(sprintf("  Mt Central Dogma: %d → %d genes", length(central_dogma_genes_raw), length(central_dogma_genes)))
message(sprintf("  Mt Ribosomes: %d → %d genes", length(mito_ribosome_genes_raw), length(mito_ribosome_genes)))
message(sprintf("  ATP Synthase: %d → %d genes", length(oxphos_genes_raw), length(oxphos_genes)))
message(sprintf("  Synaptic Common: %d → %d genes", length(synaptic_common_genes_raw), length(synaptic_common_genes)))
message(sprintf("  Postsynaptic Only: %d → %d genes", length(postsynaptic_only_genes_raw), length(postsynaptic_only_genes)))
message(sprintf("  Calcium Signaling: %d → %d genes\n", length(calcium_genes_raw), length(calcium_genes)))

###############################################################################
##  Combine Gene Sets into Modules                                          ##
###############################################################################

all_module_genes <- list(
  "1. Mt Central Dogma" = central_dogma_genes,
  "2. Mt Ribosomes" = mito_ribosome_genes,
  "3. ATP Synthase (Complex V)" = oxphos_genes,
  "4. Synaptic Ribo (Common)" = synaptic_common_genes,
  "5. Postsynaptic Ribo (Only)" = postsynaptic_only_genes,
  "6. Calcium Signaling" = calcium_genes
)

# Build cascade data matrix
cascade_data <- data.frame()
module_annotation <- c()

for (module_name in names(all_module_genes)) {
  genes <- all_module_genes[[module_name]]
  genes_present <- genes[genes %in% rownames(logfc_matrix)]

  if (length(genes_present) > 0) {
    module_data <- logfc_matrix[genes_present, , drop = FALSE]
    cascade_data <- rbind(cascade_data, module_data)
    module_annotation <- c(module_annotation, rep(module_name, length(genes_present)))
  }
}

message(sprintf("  Total genes for heatmap: %d\n", nrow(cascade_data)))

###############################################################################
##  Create Cascade Heatmap                                                  ##
###############################################################################

message("📊 Creating mechanistic cascade heatmap...\n")

# Color scheme for logFC
col_fun <- colorRamp2(
  c(-1.5, -0.75, 0, 0.75, 1.5),
  c("#0571b0", "#92c5de", "white", "#f4a582", "#ca0020")
)

# Module colors (Updated to match Fig 4 Semantic Colors where possible)
module_colors <- c(
  "1. Mt Central Dogma" = "#CC79A7",
  "2. Mt Ribosomes" = "#E69F00",
  "3. ATP Synthase (Complex V)" = "#DAA520", # Goldenrod (matches Fig 4)
  "4. Synaptic Ribo (Common)" = "#009E73",
  "5. Postsynaptic Ribo (Only)" = "#D55E00",
  "6. Calcium Signaling" = "#DC143C"       # Crimson (matches Fig 4)
)

# Row annotation
ha_module <- rowAnnotation(
  Module = module_annotation,
  col = list(Module = module_colors),
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    Module = list(
      title = "Module",
      title_gp = gpar(fontface = "bold", fontsize = 9),
      labels_gp = gpar(fontsize = 8)
    )
  ),
  width = unit(0.4, "cm")
)

# Column annotation
mutation_groups <- c(rep("G32A", 3), rep("R403C", 3))
trajectory_groups <- rep(c("Early", "TrajDev", "Late"), 2)

ha_col <- HeatmapAnnotation(
  Mutation = mutation_groups,
  Trajectory = trajectory_groups,
  col = list(
    Mutation = c("G32A" = "#E41A1C", "R403C" = "#377EB8"),
    Trajectory = c("Early" = "#FEE0D2", "TrajDev" = "#FC9272", "Late" = "#DE2D26")
  ),
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 9),
  simple_anno_size = unit(0.3, "cm"),
  annotation_legend_param = list(
    Mutation = list(
      title = "Mutation",
      title_gp = gpar(fontface = "bold", fontsize = 9),
      labels_gp = gpar(fontsize = 8)
    ),
    Trajectory = list(
      title = "Stage",
      title_gp = gpar(fontface = "bold", fontsize = 9),
      labels_gp = gpar(fontsize = 8)
    )
  )
)

# Figure dimensions
n_genes <- nrow(cascade_data)
fig_height <- max(8, n_genes * 0.25 + 3)
fig_width <- 7

message(sprintf("  Figure dimensions: %.1f x %.1f inches (%d genes)\n",
                fig_width, fig_height, n_genes))

# Create heatmap
pdf(file.path(out_dir, "Mechanistic_Cascade_Heatmap.pdf"),
    width = fig_width, height = fig_height)

ht <- Heatmap(
  as.matrix(cascade_data),
  name = "logFC",
  col = col_fun,

  # Row settings
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 7),
  cluster_rows = FALSE,
  show_row_dend = FALSE,
  row_split = factor(module_annotation, levels = names(all_module_genes)),
  row_title = gsub("^[0-9]+\\. ", "", names(all_module_genes)),
  row_title_gp = gpar(fontface = "bold", fontsize = 9),
  row_title_rot = 0,
  row_gap = unit(2, "mm"),

  # Column settings
  cluster_columns = FALSE,
  show_column_dend = FALSE,
  column_names_gp = gpar(fontsize = 8),
  column_names_rot = 0,
  column_split = factor(mutation_groups, levels = c("G32A", "R403C")),
  column_title = "Mechanistic Cascade: Pathway Core Genes (Top 20 by |logFC|)",
  column_title_gp = gpar(fontface = "bold", fontsize = 11),

  # Annotations
  left_annotation = ha_module,
  top_annotation = ha_col,

  # Appearance
  border = TRUE,
  rect_gp = gpar(col = "gray90", lwd = 0.2),
  width = unit(3.5, "cm"),

  # Legend
  heatmap_legend_param = list(
    title = "logFC",
    at = c(-1.5, -0.75, 0, 0.75, 1.5),
    labels = c("-1.5", "-0.75", "0", "0.75", "1.5"),
    legend_height = unit(3, "cm"),
    title_gp = gpar(fontface = "bold", fontsize = 9),
    labels_gp = gpar(fontsize = 8)
  )
)

draw(ht,
     heatmap_legend_side = "right",
     annotation_legend_side = "right",
     merge_legend = TRUE)

dev.off()

message("✓ Cascade heatmap complete\n")

###############################################################################
##  Create Module Summary Heatmap                                           ##
###############################################################################

message("📊 Creating module summary heatmap...\n")

# Calculate mean logFC for each module × contrast
module_stats <- data.frame()

for (module_name in names(all_module_genes)) {
  genes <- all_module_genes[[module_name]]
  genes_present <- genes[genes %in% rownames(logfc_matrix)]

  if (length(genes_present) > 0) {
    means <- colMeans(logfc_matrix[genes_present, , drop = FALSE], na.rm = TRUE)

    module_stats <- rbind(module_stats, data.frame(
      Module = module_name,
      Contrast = names(means),
      Mean_logFC = as.numeric(means),
      N_genes = length(genes_present)
    ))
  }
}

# Clean module names
module_stats$Module_Clean <- gsub("^[0-9]+\\. ", "", module_stats$Module)

# Create wide format
module_wide <- module_stats %>%
  select(Module_Clean, Contrast, Mean_logFC) %>%
  pivot_wider(names_from = Contrast, values_from = Mean_logFC)

module_mat <- as.matrix(module_wide[, -1])
rownames(module_mat) <- module_wide$Module_Clean

# Reorder rows
module_order <- gsub("^[0-9]+\\. ", "", names(all_module_genes))
module_mat <- module_mat[module_order[module_order %in% rownames(module_mat)], ]

# Summary heatmap
col_fun_summary <- colorRamp2(
  c(-0.6, -0.3, 0, 0.3, 0.6),
  c("#0571b0", "#92c5de", "white", "#f4a582", "#ca0020")
)

pdf(file.path(out_dir, "Module_Summary_Heatmap.pdf"), width = 6, height = 5)

ht_summary <- Heatmap(
  module_mat,
  name = "Mean\nlogFC",
  col = col_fun_summary,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 10),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_names_gp = gpar(fontsize = 8),
  column_names_rot = 0,
  column_split = factor(c(rep("G32A", 3), rep("R403C", 3)),
                        levels = c("G32A", "R403C")),
  column_title = "Module-Level Summary (Mean logFC)",
  column_title_gp = gpar(fontface = "bold", fontsize = 11),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.2f", module_mat[i, j]),
              x, y, gp = gpar(fontsize = 7))
  },
  border = TRUE,
  heatmap_legend_param = list(
    title = "Mean\nlogFC",
    legend_height = unit(2.5, "cm"),
    title_gp = gpar(fontface = "bold", fontsize = 9),
    labels_gp = gpar(fontsize = 8)
  )
)

draw(ht_summary)
dev.off()

message("✓ Module summary heatmap complete\n")

###############################################################################
##  Summary Statistics                                                       ##
###############################################################################

message("📝 Generating summary...\n")

cat("\n=== Mechanistic Cascade: Real Pathway Gene Sets ===\n\n")

cat("Modules extracted from SIGNIFICANT GSEA pathways:\n")
cat("  1. Mt Central Dogma: MitoCarta (p=2.14e-13)\n")
cat("  2. Mt Ribosomes: MitoCarta (p=7.98e-04)\n")
cat("  3. ATP Synthase: GOCC_ATPASE_COMPLEX (p=0.002)\n")
cat("  4. Synaptic Ribo Common: SynGO pre∩post (p<1e-12)\n")
cat("  5. Postsynaptic Ribo Only: SynGO post-pre\n")
cat("  6. Calcium Signaling: REACTOME_ELEVATION_OF_CYTOSOLIC_CA2_LEVELS (p=0.03 in Early)\n\n")

cat("Gene counts per module:\n")
for (mod in names(all_module_genes)) {
  cat(sprintf("  %s: %d genes\n", mod, length(all_module_genes[[mod]])))
}

cat("\nModule Mean logFC Summary:\n\n")
print(module_stats %>%
        select(Module_Clean, Contrast, Mean_logFC) %>%
        pivot_wider(names_from = Contrast, values_from = Mean_logFC) %>%
        as.data.frame(),
      row.names = FALSE)

cat("\n")
message("✅ Mechanistic cascade visualization complete!")
message(sprintf("📁 Output directory: %s", out_dir))
message("\n🎯 KEY OUTPUTS:")
message("  1. Mechanistic_Cascade_Heatmap.pdf - Gene-level heatmap (6 contrasts)")
message("  2. Module_Summary_Heatmap.pdf - Module-level mean summary")
