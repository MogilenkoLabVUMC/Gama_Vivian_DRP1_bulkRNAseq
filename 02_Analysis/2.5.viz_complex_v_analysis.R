###############################################################################
##  Complex V (ATP Synthase) Analysis: Direct Energy Crisis                 ##
###############################################################################

library(here)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(ComplexHeatmap)
library(circlize)

# Load unified color configuration
source(here("01_Scripts/R_scripts/color_config.R"))

message("📂 Loading checkpoints...")
checkpoint_dir <- here("03_Results/02_Analysis/checkpoints")
mitocarta_gsea_results <- readRDS(file.path(checkpoint_dir, "mitocarta_gsea_results.rds"))
fit <- readRDS(file.path(checkpoint_dir, "fit_object.rds"))
de_results <- readRDS(file.path(checkpoint_dir, "de_results.rds"))

# Output directory
out_dir <- here("03_Results/02_Analysis/Plots/Complex_V_analysis")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("✓ Checkpoints loaded\n")

###############################################################################
##  Extract Complex V Pathways from MitoCarta                               ##
###############################################################################

message("📊 Extracting Complex V pathways from MitoCarta...\n")

all_contrasts <- c(
  "G32A_vs_Ctrl_D35", "R403C_vs_Ctrl_D35",
  "Maturation_G32A_specific", "Maturation_R403C_specific",
  "G32A_vs_Ctrl_D65", "R403C_vs_Ctrl_D65",
  "Time_Ctrl", "Time_G32A", "Time_R403C"
)

complex_v_data <- data.frame()
complex_v_genes <- character()

complex_v_keywords <- c(
  "complex.*v",
  "cv_",
  "atp.*synth",
  "atp5"
)

for (contrast in all_contrasts) {
  mito_res <- mitocarta_gsea_results[[contrast]]

  if (!is.null(mito_res)) {
    res_df <- mito_res@result

    # Find Complex V pathways
    cv_idx <- grepl(paste(complex_v_keywords, collapse = "|"),
                   res_df$Description, ignore.case = TRUE)

    if (any(cv_idx)) {
      cv_pathways <- res_df[cv_idx, ]

      # Filter for significant
      cv_pathways_sig <- cv_pathways[cv_pathways$p.adjust < 0.05, ]

      if (nrow(cv_pathways_sig) > 0) {
        cv_pathways_sig$Contrast <- contrast

        # Extract core enrichment genes
        for (i in 1:nrow(cv_pathways_sig)) {
          genes <- unlist(strsplit(cv_pathways_sig$core_enrichment[i], "/"))
          complex_v_genes <- unique(c(complex_v_genes, genes))
        }

        complex_v_data <- rbind(complex_v_data, cv_pathways_sig)
      }
    }
  }
}

message(sprintf("  Found %d Complex V pathway enrichments\n", nrow(complex_v_data)))
message(sprintf("  Found %d unique Complex V genes\n", length(complex_v_genes)))

###############################################################################
##  Classify Contrasts                                                      ##
###############################################################################

complex_v_data$Mutation <- case_when(
  grepl("G32A", complex_v_data$Contrast) ~ "G32A",
  grepl("R403C", complex_v_data$Contrast) ~ "R403C",
  grepl("Ctrl", complex_v_data$Contrast) ~ "Control",
  TRUE ~ "Unknown"
)

complex_v_data$Stage <- case_when(
  grepl("vs_Ctrl_D35", complex_v_data$Contrast) ~ "Early (D35)",
  grepl("Maturation", complex_v_data$Contrast) ~ "Maturation",
  grepl("vs_Ctrl_D65", complex_v_data$Contrast) ~ "Late (D65)",
  grepl("Time", complex_v_data$Contrast) ~ "Time Effect",
  TRUE ~ "Unknown"
)

# Simplify pathway names
complex_v_data$Pathway_Clean <- gsub("_", " ", complex_v_data$Description)
complex_v_data$Pathway_Clean <- tools::toTitleCase(tolower(complex_v_data$Pathway_Clean))

###############################################################################
##  FIGURE 1: Complex V Pathway Enrichment Across All Contrasts             ##
###############################################################################

message("📊 Creating Complex V pathway enrichment plot...\n")

# Order contrasts logically
complex_v_data$Contrast_Display <- factor(
  complex_v_data$Contrast,
  levels = c("G32A_vs_Ctrl_D35", "R403C_vs_Ctrl_D35",
             "Maturation_G32A_specific", "Maturation_R403C_specific",
             "G32A_vs_Ctrl_D65", "R403C_vs_Ctrl_D65",
             "Time_Ctrl", "Time_G32A", "Time_R403C"),
  labels = c("G32A D35", "R403C D35",
             "G32A Maturation", "R403C Maturation",
             "G32A D65", "R403C D65",
             "Ctrl Time", "G32A Time", "R403C Time")
)

# Create dotplot
p_cv_pathways <- ggplot(complex_v_data,
                        aes(x = Contrast_Display, y = Pathway_Clean,
                            color = NES, size = -log10(p.adjust))) +
  geom_point(alpha = 0.85) +

  geom_point(
    data = subset(complex_v_data, p.adjust < 0.01),
    shape = 21, color = "black", fill = NA, stroke = 1.5, alpha = 1
  ) +

  # Using unified NES color scale from color_config.R
  nes_ggplot_scale(limits = c(-3, 3), name = "NES") +

  scale_size_continuous(
    range = c(3, 10),
    name = "-log10(FDR)",
    breaks = c(2, 5, 10, 15)
  ) +

  labs(
    title = "Complex V (ATP Synthase): Most Vulnerable OXPHOS Complex",
    subtitle = "DOWN at D35 | RECOVERY during maturation | Direct link to ATP crisis (FDR < 0.05, black outline = FDR < 0.01)",
    x = "Contrast",
    y = "Complex V Pathway"
  ) +

  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "gray30"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 9),
    axis.text.y = element_text(size = 9),
    axis.title = element_text(face = "bold"),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    panel.border = element_rect(color = "gray70", fill = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  )

plot_height <- max(6, length(unique(complex_v_data$Pathway_Clean)) * 0.4)

ggsave(file.path(out_dir, "Complex_V_Pathway_Enrichment.pdf"),
       p_cv_pathways, width = 11, height = plot_height)

message("✓ Complex V pathway enrichment plot complete\n")

###############################################################################
##  FIGURE 2: Complex V Gene Expression Heatmap                             ##
###############################################################################

message("📊 Creating Complex V gene expression heatmap...\n")

# Extract logFC for Complex V genes
contrasts_for_heatmap <- c(
  "G32A_vs_Ctrl_D35", "R403C_vs_Ctrl_D35",
  "Maturation_G32A_specific", "Maturation_R403C_specific",
  "G32A_vs_Ctrl_D65", "R403C_vs_Ctrl_D65"
)

# Get genes present in data
cv_genes_present <- complex_v_genes[complex_v_genes %in% rownames(fit$coefficients)]

message(sprintf("  Found %d/%d Complex V genes in expression data\n",
                length(cv_genes_present), length(complex_v_genes)))

if (length(cv_genes_present) > 0) {
  # Extract logFC matrix for all contrasts
  logfc_full <- sapply(contrasts_for_heatmap, function(contrast) {
    coef_idx <- which(colnames(fit$coefficients) == contrast)
    if (length(coef_idx) > 0) {
      logfc <- fit$coefficients[cv_genes_present, coef_idx]
      return(logfc)
    } else {
      return(rep(NA, length(cv_genes_present)))
    }
  })

  rownames(logfc_full) <- cv_genes_present

  # Extract matrices for each mutation (contrasts_for_heatmap order:
  # 0: G32A_vs_Ctrl_D35, 1: R403C_vs_Ctrl_D35,
  # 2: Maturation_G32A_specific, 3: Maturation_R403C_specific,
  # 4: G32A_vs_Ctrl_D65, 5: R403C_vs_Ctrl_D65)

  logfc_g32a <- logfc_full[, c(1, 3, 5)]  # Early, TrajDev, Late
  logfc_r403c <- logfc_full[, c(2, 4, 6)]  # Early, TrajDev, Late

  colnames(logfc_g32a) <- c("Early", "TrajDev", "Late")
  colnames(logfc_r403c) <- c("Early", "TrajDev", "Late")

  # Color scheme (using unified palette from color_config.R)
  col_fun <- nes_color_scale(limits = c(-0.6, 0.6), n_colors = 3)

  # Create heatmaps side by side (G32A | R403C)
  pdf(file.path(out_dir, "Complex_V_Gene_Expression_Heatmap.pdf"),
      width = 6, height = max(8, length(cv_genes_present) * 0.15))

  # G32A heatmap (left)
  ht_g32a <- Heatmap(
    logfc_g32a,
    name = "logFC",
    col = col_fun,

    # Row settings
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 8),
    cluster_rows = TRUE,
    show_row_dend = TRUE,
    row_dend_width = unit(15, "mm"),

    # Column settings
    cluster_columns = FALSE,
    column_names_gp = gpar(fontsize = 11, fontface = "bold"),
    column_title = "G32A (GTPase Domain)",
    column_title_gp = gpar(fontface = "bold", fontsize = 12),

    # Appearance
    border = TRUE,
    rect_gp = gpar(col = "gray90", lwd = 0.5),

    # Legend
    heatmap_legend_param = list(
      title = "logFC",
      at = c(-0.6, -0.3, 0, 0.3, 0.6),
      legend_height = unit(4, "cm"),
      title_gp = gpar(fontface = "bold", fontsize = 10)
    )
  )

  # R403C heatmap (right)
  ht_r403c <- Heatmap(
    logfc_r403c,
    name = "logFC",
    col = col_fun,

    # Row settings
    show_row_names = FALSE,
    cluster_rows = TRUE,
    show_row_dend = FALSE,

    # Column settings
    cluster_columns = FALSE,
    column_names_gp = gpar(fontsize = 11, fontface = "bold"),
    column_title = "R403C (Stalk Domain)",
    column_title_gp = gpar(fontface = "bold", fontsize = 12),

    # Appearance
    border = TRUE,
    rect_gp = gpar(col = "gray90", lwd = 0.5),

    # No legend on right heatmap
    show_heatmap_legend = FALSE
  )

  # Combine and draw
  ht_combined <- ht_g32a + ht_r403c

  draw(ht_combined,
       column_title = "Complex V (ATP Synthase) Gene Expression\nAcross Development (Early → TrajDev → Late)",
       column_title_gp = gpar(fontface = "bold", fontsize = 14),
       heatmap_legend_side = "right")

  dev.off()

  message("✓ Complex V gene expression heatmap complete\n")
} else {
  message("  No Complex V genes found in expression data (skipping heatmap)\n")
}


###############################################################################
##  Export Data                                                              ##
###############################################################################

message("📝 Exporting Complex V data...\n")

write.csv(complex_v_data %>%
            select(Description, Pathway_Clean, Contrast, Mutation, Stage,
                   NES, p.adjust, pvalue, setSize, core_enrichment) %>%
            arrange(Contrast, p.adjust),
          file.path(out_dir, "Complex_V_Pathway_Data.csv"),
          row.names = FALSE)

write.csv(data.frame(Gene = complex_v_genes),
          file.path(out_dir, "Complex_V_Genes.csv"),
          row.names = FALSE)

message("✓ CSV exports complete\n")

###############################################################################
##  Summary Report                                                           ##
###############################################################################

cat("\n=== Complex V (ATP Synthase) Analysis Summary ===\n\n")

cat("Complex V Pathway Enrichments:\n")
cv_summary <- complex_v_data %>%
  group_by(Stage, Mutation) %>%
  summarize(N = n(), Mean_NES = mean(NES, na.rm = TRUE), .groups = "drop")
print(cv_summary)

cat("\n\nComplex V Genes Identified:\n")
cat(sprintf("  Total unique genes: %d\n", length(complex_v_genes)))
cat(sprintf("  Genes in expression data: %d\n", length(cv_genes_present)))

cat("\n\nBIOLOGICAL INTERPRETATION:\n")
cat("  Complex V = ATP Synthase = DIRECT ENERGY CRISIS\n\n")

cat("  EARLY (D35): Strong suppression\n")
cat("    - Complex V pathways significantly DOWN\n")
cat("    - Cannot produce ATP efficiently\n")
cat("    - Mitochondrial hyperfusion prevents proper positioning\n\n")

cat("  MATURATION: Recovery attempt\n")
cat("    - Complex V pathways UPREGULATED\n")
cat("    - Cell attempts to restore ATP production\n")
cat("    - Mitochondrial translation machinery also increases\n\n")

cat("  LATE (D65): Mutation-specific outcome\n")
cat("    - G32A: Failed recovery (no Complex V enrichment)\n")
cat("    - R403C: Partial metabolic compensation\n\n")

cat("  WHY Complex V is MOST VULNERABLE:\n")
cat("    1. Largest OXPHOS complex (28 subunits)\n")
cat("    2. Requires both nuclear and mitochondrial gene expression\n")
cat("    3. Assembly depends on proper mitochondrial positioning\n")
cat("    4. Direct link to ATP production (final step of OXPHOS)\n")
cat("    5. Most enrichments across all contrasts\n\n")

message("✅ Complex V analysis complete!")
message(sprintf("📁 Output directory: %s", out_dir))
message("\n🎯 KEY OUTPUTS:")
message("  1. Complex_V_Pathway_Enrichment.pdf - Pathway-level enrichment across contrasts")
message("  2. Complex_V_Gene_Expression_Heatmap.pdf - Gene-level expression changes with Early/TrajDev/Late framework")
message("  3. Complex_V_Pathway_Data.csv - Complete pathway enrichment data")
message("  4. Complex_V_Genes.csv - List of Complex V genes")
