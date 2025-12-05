###############################################################################
##  Mitochondrial & Calcium Signaling Cascade Heatmap                       ##
##  Focused visualization: Complex V, Calcium Genes, Mt Central Dogma       ##
##                                                                           ##
##  Style: Matches Panel_C_Expression_Heatmap.pdf (synaptic ribosomes)      ##
##  - Hierarchical clustering within modules                                 ##
##  - Narrower logFC scale (-0.6 to 0.6) to reveal subtle patterns          ##
##  - Dendrogram showing gene co-regulation                                  ##
###############################################################################

library(here)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(grid)

# Load unified color configuration
source(here("01_Scripts/R_scripts/color_config.R"))

message("📂 Loading checkpoints...")
checkpoint_dir <- here("03_Results/02_Analysis/checkpoints")
mitocarta_gsea_results <- readRDS(file.path(checkpoint_dir, "mitocarta_gsea_results.rds"))
fit <- readRDS(file.path(checkpoint_dir, "fit_object.rds"))

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
  "Early", "TrajDev", "Late",
  "Early", "TrajDev", "Late"
)

# Calcium genes from main pipeline config - MUST include NNAT and PNPO
# These are genes of specific interest for calcium signaling in the study
CALCIUM_GENES_CONFIG <- c(
  "NNAT",    # Neuronatin - key imprinted gene, calcium regulation

  "CACNG3",  # Voltage-dependent calcium channel gamma-3
  "CACNA1S", # Voltage-dependent L-type calcium channel alpha-1S
  "ATP2A1",  # Sarcoplasmic/endoplasmic reticulum calcium ATPase 1
  "RYR1",    # Ryanodine receptor 1 (skeletal muscle type)
  "MYLK3",   # Myosin light chain kinase 3
  "VDR",     # Vitamin D receptor (regulates calcium homeostasis)
  "STIM1",   # Stromal interaction molecule 1 (ER calcium sensor)
  "STIM2",   # Stromal interaction molecule 2
  "ORAI1",   # Calcium release-activated calcium modulator 1
  "CALB1",   # Calbindin 1 (calcium-binding protein)
  "CALR",    # Calreticulin (calcium-binding ER protein)
  "PNPO"     # Pyridoxamine 5'-phosphate oxidase (B6 metabolism, affects calcium)
)

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
##  Extract Gene Sets: 3 FOCUSED MODULES                                    ##
##  1. ATP Synthase (Complex V) - ATP5* genes from KEGG OXPHOS pathway      ##
##  2. Calcium Signaling - Config genes (NNAT, PNPO, etc.)                  ##
##  3. Mitochondrial Central Dogma - MitoCarta pathway                      ##
###############################################################################

message("📊 Extracting gene sets for 3 focused modules...\n")

## MODULE 1: ATP Synthase (Complex V) - TRUE ATP5* genes
## The GOCC_ATPASE_COMPLEX contains chromatin remodeling ATPases (wrong!)
## We use curated ATP synthase subunit genes from KEGG OXPHOS pathway
message("  Module 1: ATP Synthase (Complex V) - Curated ATP5* subunits...")

# ATP synthase F1 subunits (catalytic core)
atp_f1 <- c("ATP5F1A", "ATP5F1B", "ATP5F1C", "ATP5F1D", "ATP5F1E")
# ATP synthase FO subunits (proton channel)
atp_fo <- c("ATP5PB", "ATP5MC1", "ATP5MC2", "ATP5MC3", "ATP5PD", "ATP5PF",
            "ATP5MF", "ATP5MG", "ATP5PO", "ATP5ME", "ATP5MD")
# ATP synthase peripheral stalk
atp_stalk <- c("ATP5MJ", "ATP5MK", "ATP5ML", "ATP5MF")
# Coupling factors and assembly factors
atp_assembly <- c("ATPAF1", "ATPAF2", "TMEM70", "ATP5IF1")

atp_synthase_genes_curated <- unique(c(atp_f1, atp_fo, atp_stalk, atp_assembly))
atp_synthase_genes <- atp_synthase_genes_curated[atp_synthase_genes_curated %in% rownames(logfc_matrix)]
message(sprintf("    ATP Synthase: %d/%d curated genes found in expression data",
                length(atp_synthase_genes), length(atp_synthase_genes_curated)))

## MODULE 2: Calcium Signaling - Config genes (prioritized, curated list)
message("  Module 2: Calcium Signaling - Config genes (incl. NNAT, PNPO)...")
calcium_genes <- CALCIUM_GENES_CONFIG[CALCIUM_GENES_CONFIG %in% rownames(logfc_matrix)]
message(sprintf("    Calcium genes: %d/%d config genes found in expression data",
                length(calcium_genes), length(CALCIUM_GENES_CONFIG)))

# Report which config genes are missing
missing_calcium <- setdiff(CALCIUM_GENES_CONFIG, rownames(logfc_matrix))
if (length(missing_calcium) > 0) {
  message(sprintf("    Missing: %s", paste(missing_calcium, collapse = ", ")))
}

## MODULE 3: Mitochondrial Central Dogma (MitoCarta)
## Extract core enrichment genes, then limit to top 25 by mean |logFC| across TrajDev
message("  Module 3: Mitochondrial Central Dogma (MitoCarta)...")
central_dogma_genes_all <- extract_pathway_genes(
  mitocarta_gsea_results[["Maturation_G32A_specific"]],
  "Mitochondrial_central_dogma", exact_match = TRUE
)
central_dogma_genes_present <- central_dogma_genes_all[central_dogma_genes_all %in% rownames(logfc_matrix)]

# Limit to top 25 by mean |logFC| across TrajDev contrasts (like Panel_C uses focused set)
MAX_CENTRAL_DOGMA_GENES <- 25
if (length(central_dogma_genes_present) > MAX_CENTRAL_DOGMA_GENES) {
  trajdev_cols <- grep("TrajDev", colnames(logfc_matrix))
  mean_abs_fc <- rowMeans(abs(logfc_matrix[central_dogma_genes_present, trajdev_cols, drop = FALSE]), na.rm = TRUE)
  central_dogma_genes <- names(sort(mean_abs_fc, decreasing = TRUE)[1:MAX_CENTRAL_DOGMA_GENES])
  message(sprintf("    Limiting to top %d genes by |logFC| (from %d total)",
                  MAX_CENTRAL_DOGMA_GENES, length(central_dogma_genes_present)))
} else {
  central_dogma_genes <- central_dogma_genes_present
}

###############################################################################
##  Combine Gene Sets into 3 Modules                                        ##
###############################################################################

message("\n📊 Building heatmap data matrix...\n")

# Define modules with clear, descriptive names
all_module_genes <- list(
  "ATP Synthase (Complex V)" = atp_synthase_genes,
  "Calcium Signaling" = calcium_genes,
  "Mt Central Dogma" = central_dogma_genes
)

# Report gene counts
for (mod in names(all_module_genes)) {
  message(sprintf("  %s: %d genes", mod, length(all_module_genes[[mod]])))
}

# Build cascade data matrix with module annotation
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

message(sprintf("\n  Total genes for heatmap: %d\n", nrow(cascade_data)))

###############################################################################
##  Create Cascade Heatmap (Panel_C style: clustering + dendrogram)         ##
###############################################################################

message("📊 Creating mechanistic cascade heatmap (Panel_C style)...\n")

# Color scheme for logFC - NARROWER SCALE to reveal subtle patterns
# Using -0.6 to 0.6 like Panel_C (synaptic ribosomes)
col_fun <- colorRamp2(c(-0.6, 0, 0.6), c("#2166AC", "#F7F7F7", "#B35806"))

# Module colors - simplified for 3 modules
module_colors <- c(
 "ATP Synthase (Complex V)" = "#2E8B57",  # Sea green
 "Calcium Signaling"        = "#DC143C",  # Crimson
 "Mt Central Dogma"         = "#9467BD"   # Purple
)

# Convert module annotation to factor with explicit order
module_factor <- factor(module_annotation, levels = names(all_module_genes))

# Row annotation (right side, like Panel_C)
ha_module <- rowAnnotation(
  Module = module_factor,
  col = list(Module = module_colors),
  show_annotation_name = TRUE,
  annotation_name_side = "top",
  annotation_legend_param = list(
    Module = list(
      title = "Functional Module",
      title_gp = gpar(fontface = "bold", fontsize = 10),
      labels_gp = gpar(fontsize = 9)
    )
  )
)

# Column split for mutations
column_split <- factor(rep(c("G32A", "R403C"), each = 3), levels = c("G32A", "R403C"))

# Figure dimensions - taller to accommodate genes with clustering
n_genes <- nrow(cascade_data)
fig_height <- max(10, n_genes * 0.18 + 2)  # Adjusted for better gene label visibility
fig_width <- 6

message(sprintf("  Figure dimensions: %.1f x %.1f inches (%d genes)\n",
                fig_width, fig_height, n_genes))

# Create heatmap with CLUSTERING WITHIN MODULES (like Panel_C)
pdf(file.path(out_dir, "Mechanistic_Cascade_Heatmap.pdf"),
    width = fig_width, height = fig_height)

ht <- Heatmap(
  as.matrix(cascade_data),
  name = "logFC",
  col = col_fun,

  # Row settings - ENABLE CLUSTERING within each module (key difference from before!)
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 8),
  cluster_rows = TRUE,              # <-- ENABLE clustering within modules
  cluster_row_slices = FALSE,       # Don't reorder module slices
  show_row_dend = TRUE,             # <-- SHOW dendrogram
  row_dend_width = unit(10, "mm"),  # Dendrogram width
  row_split = module_factor,
  row_title_gp = gpar(fontface = "bold", fontsize = 10),
  row_title_rot = 0,
  row_gap = unit(2, "mm"),

  # Column settings
  cluster_columns = FALSE,
  show_column_dend = FALSE,
  column_names_gp = gpar(fontsize = 9),
  column_split = column_split,
  column_title_gp = gpar(fontface = "bold", fontsize = 11),
  column_gap = unit(2, "mm"),

  # Annotations
  right_annotation = ha_module,

  # Appearance
  border = TRUE,
  rect_gp = gpar(col = "gray80", lwd = 0.5),

  # Legend
  heatmap_legend_param = list(
    title = "logFC",
    at = c(-0.6, -0.3, 0, 0.3, 0.6),
    labels = c("-0.6", "-0.3", "0", "0.3", "0.6"),
    legend_height = unit(4, "cm"),
    title_gp = gpar(fontface = "bold", fontsize = 10)
  )
)

draw(ht,
     heatmap_legend_side = "right",
     column_title = "Mitochondrial & Calcium Gene Expression Trajectories")

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
      Contrast = colnames(logfc_matrix),
      Mean_logFC = as.numeric(means),
      N_genes = length(genes_present)
    ))
  }
}

# Create wide format for heatmap - use unique column names
# First, add mutation info to make columns unique
module_stats$ColName <- paste0(
  rep(c("G32A", "R403C"), each = 3, times = length(unique(module_stats$Module))),
  "_",
  module_stats$Contrast
)[1:nrow(module_stats)]  # Handle potential length mismatch

# Actually let's do this more carefully
# The contrast names are "Early", "TrajDev", "Late" repeated twice
# We need to make them unique before pivoting
module_stats_unique <- module_stats %>%
  mutate(
    Mutation = rep(rep(c("G32A", "R403C"), each = 3), times = length(unique(Module))),
    ColName = paste0(Mutation, "\n", Contrast)
  ) %>%
  select(Module, ColName, Mean_logFC) %>%
  distinct()

module_wide <- module_stats_unique %>%
  pivot_wider(names_from = ColName, values_from = Mean_logFC)

module_mat <- as.matrix(module_wide[, -1])
rownames(module_mat) <- module_wide$Module

# Ensure correct column order
expected_cols <- c("G32A\nEarly", "G32A\nTrajDev", "G32A\nLate",
                   "R403C\nEarly", "R403C\nTrajDev", "R403C\nLate")
module_mat <- module_mat[, expected_cols[expected_cols %in% colnames(module_mat)]]

# Summary heatmap with same color scale
col_fun_summary <- colorRamp2(c(-0.3, 0, 0.3), c("#2166AC", "#F7F7F7", "#B35806"))

pdf(file.path(out_dir, "Module_Summary_Heatmap.pdf"), width = 7, height = 4)

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
  column_split = factor(rep(c("G32A", "R403C"), each = 3),
                        levels = c("G32A", "R403C")),
  column_title = "Module-Level Summary (Mean logFC)",
  column_title_gp = gpar(fontface = "bold", fontsize = 11),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.2f", module_mat[i, j]),
              x, y, gp = gpar(fontsize = 8))
  },
  border = TRUE,
  heatmap_legend_param = list(
    title = "Mean\nlogFC",
    at = c(-0.3, -0.15, 0, 0.15, 0.3),
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

cat("\n=== Mitochondrial & Calcium Gene Expression Analysis ===\n\n")

cat("THREE FOCUSED MODULES:\n")
cat("  1. ATP Synthase (Complex V): Curated ATP5* subunit genes from KEGG hsa00190\n")
cat("     - F1 subunits (catalytic): ATP5F1A/B/C/D/E\n")
cat("     - FO subunits (proton channel): ATP5PB, ATP5MC1-3, ATP5PD/PF/MF/MG/PO/ME/MD\n")
cat("     - Assembly factors: ATPAF1/2, TMEM70\n\n")
cat("  2. Calcium Signaling: Study-prioritized genes\n")
cat("     - Key targets: NNAT (neuronatin), PNPO\n")
cat("     - Channels: CACNG3, CACNA1S, RYR1\n")
cat("     - ER calcium sensors: STIM1, STIM2, CALR\n")
cat("     - Other: ATP2A1, MYLK3, VDR, CALB1\n\n")
cat("  3. Mt Central Dogma: MitoCarta pathway (GSEA p.adj=2.14e-13)\n")
cat("     - Top 25 genes by |logFC| from core enrichment\n\n")

cat("Gene counts per module:\n")
for (mod in names(all_module_genes)) {
  cat(sprintf("  %s: %d genes\n", mod, length(all_module_genes[[mod]])))
}

cat("\nModule Mean logFC Summary:\n\n")
print(module_stats %>%
        select(Module, Contrast, Mean_logFC) %>%
        pivot_wider(names_from = Contrast, values_from = Mean_logFC) %>%
        as.data.frame(),
      row.names = FALSE)

# List actual genes included
cat("\n\nGenes included in each module:\n")
for (mod in names(all_module_genes)) {
  genes <- all_module_genes[[mod]]
  cat(sprintf("\n%s (%d genes):\n", mod, length(genes)))
  cat(sprintf("  %s\n", paste(sort(genes), collapse = ", ")))
}

cat("\n")
message("✅ Mitochondrial & Calcium visualization complete!")
message(sprintf("📁 Output directory: %s", out_dir))
message("\n🎯 KEY OUTPUTS:")
message("  1. Mechanistic_Cascade_Heatmap.pdf - Gene-level heatmap with clustering")
message("  2. Module_Summary_Heatmap.pdf - Module-level mean summary")
message("\n🎯 STYLE NOTES:")
message("  - Matches Panel_C_Expression_Heatmap.pdf (synaptic ribosomes)")
message("  - Hierarchical clustering within modules shows co-regulated genes")
message("  - Narrower logFC scale (-0.6 to 0.6) reveals subtle patterns")
message("  - Dendrogram indicates gene co-regulation relationships")
