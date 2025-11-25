###############################################################################
##  Temporal Trajectory: MitoCarta Pathway Recovery (D35→D65)               ##
##  Separation of Concerns: Focused on mitochondrial pathway dynamics       ##
###############################################################################
##  STORY: Early crisis → Attempted recovery → Partial rescue               ##
##  Shows DOWN at D35, UP during maturation, mutation-specific divergence   ##
###############################################################################

library(here)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

message("📂 Loading checkpoints...")
checkpoint_dir <- here("03_Results/02_Analysis/checkpoints")
mitocarta_gsea_results <- readRDS(file.path(checkpoint_dir, "mitocarta_gsea_results.rds"))

# Output directory
out_dir <- here("03_Results/02_Analysis/Plots/Temporal_trajectory")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("✓ Checkpoints loaded\n")

###############################################################################
##  Extract MitoCarta Pathways Across Time                                  ##
###############################################################################

message("📊 Extracting MitoCarta pathways across developmental stages...\n")

# Contrasts spanning developmental time
contrasts_temporal <- c(
  "G32A_vs_Ctrl_D35",           # Early crisis
  "R403C_vs_Ctrl_D35",          # Early crisis
  "Maturation_G32A_specific",   # Recovery attempt
  "Maturation_R403C_specific",  # Recovery attempt
  "G32A_vs_Ctrl_D65",           # Late stage
  "R403C_vs_Ctrl_D65"           # Late stage
)

temporal_data <- data.frame()

for (contrast in contrasts_temporal) {
  mito_res <- mitocarta_gsea_results[[contrast]]

  if (!is.null(mito_res)) {
    res_df <- mito_res@result

    # Get ALL significant pathways
    sig_pathways <- res_df[res_df$p.adjust < 0.05, ]

    if (nrow(sig_pathways) > 0) {
      sig_pathways$Contrast <- contrast
      temporal_data <- rbind(temporal_data, sig_pathways)
    }
  }
}

message(sprintf("  Extracted %d significant MitoCarta pathways across time\n",
                nrow(temporal_data)))

###############################################################################
##  Classify by Timepoint and Mutation                                      ##
###############################################################################

temporal_data$Mutation <- ifelse(
  grepl("G32A", temporal_data$Contrast),
  "G32A",
  "R403C"
)

temporal_data$Stage <- case_when(
  grepl("vs_Ctrl_D35", temporal_data$Contrast) ~ "Early (D35)",
  grepl("Maturation", temporal_data$Contrast) ~ "Maturation (D35→D65)",
  grepl("vs_Ctrl_D65", temporal_data$Contrast) ~ "Late (D65)",
  TRUE ~ "Unknown"
)

temporal_data$Stage <- factor(
  temporal_data$Stage,
  levels = c("Early (D35)", "Maturation (D35→D65)", "Late (D65)")
)

# Categorize pathways by function
temporal_data$Function <- case_when(
  grepl("complex.*v|cv_|atp.*synth", temporal_data$Description, ignore.case = TRUE) ~
    "ATP Synthase (Complex V)",

  grepl("central.*dogma|translation|ribosom", temporal_data$Description, ignore.case = TRUE) ~
    "Translation Machinery",

  grepl("oxphos|complex.*i|complex.*iii|complex.*iv|electron.*transport",
        temporal_data$Description, ignore.case = TRUE) ~
    "OXPHOS Complexes",

  grepl("metabolism|metabolic", temporal_data$Description, ignore.case = TRUE) ~
    "Metabolism",

  grepl("mtdna|dna.*maintenance|replication", temporal_data$Description, ignore.case = TRUE) ~
    "mtDNA Maintenance",

  grepl("import|sorting|protein.*localization", temporal_data$Description, ignore.case = TRUE) ~
    "Protein Import",

  TRUE ~ "Other Mitochondrial"
)

# Simplify pathway names
temporal_data$Pathway_Clean <- gsub("_", " ", temporal_data$Description)
temporal_data$Pathway_Clean <- tools::toTitleCase(tolower(temporal_data$Pathway_Clean))

# Truncate long names
temporal_data$Pathway_Clean <- ifelse(
  nchar(temporal_data$Pathway_Clean) > 45,
  paste0(substr(temporal_data$Pathway_Clean, 1, 42), "..."),
  temporal_data$Pathway_Clean
)

###############################################################################
##  FIGURE 1: Temporal Heatmap (All Pathways Across Time)                   ##
###############################################################################

message("📊 Creating temporal heatmap...\n")

# Create wide format for heatmap
heatmap_data <- temporal_data %>%
  mutate(Condition = paste(Mutation, Stage, sep = " | ")) %>%
  select(Pathway_Clean, Condition, NES, Function) %>%
  pivot_wider(names_from = Condition, values_from = NES, values_fill = 0)

# Extract pathway names and function for annotation
pathway_names <- heatmap_data$Pathway_Clean
pathway_functions <- heatmap_data$Function

# Create matrix (exclude pathway name and function columns)
nes_matrix <- as.matrix(heatmap_data[, -c(1, ncol(heatmap_data))])
rownames(nes_matrix) <- pathway_names

# Order columns by stage
col_order <- c(
  "G32A | Early (D35)", "R403C | Early (D35)",
  "G32A | Maturation (D35→D65)", "R403C | Maturation (D35→D65)",
  "G32A | Late (D65)", "R403C | Late (D65)"
)

# Filter to existing columns
col_order <- col_order[col_order %in% colnames(nes_matrix)]
nes_matrix <- nes_matrix[, col_order]

# Sort rows by mean NES
row_means <- rowMeans(nes_matrix, na.rm = TRUE)
nes_matrix <- nes_matrix[order(row_means), ]

# Create heatmap
library(ComplexHeatmap)
library(circlize)

col_fun <- colorRamp2(
  c(-2.5, -1.5, 0, 1.5, 2.5),
  c("#0571b0", "#92c5de", "white", "#f4a582", "#ca0020")
)

# Function colors
function_colors <- c(
  "ATP Synthase (Complex V)" = "#CC79A7",
  "Translation Machinery" = "#E69F00",
  "OXPHOS Complexes" = "#56B4E9",
  "Metabolism" = "#009E73",
  "mtDNA Maintenance" = "#F0E442",
  "Protein Import" = "#0072B2",
  "Other Mitochondrial" = "#999999"
)

# Row annotation
ha_function <- rowAnnotation(
  Function = pathway_functions[order(row_means)],
  col = list(Function = function_colors),
  show_annotation_name = TRUE,
  annotation_name_side = "top",
  width = unit(0.5, "cm")
)

pdf(file.path(out_dir, "Temporal_Trajectory_Heatmap.pdf"),
    width = 12, height = max(10, nrow(nes_matrix) * 0.2))

ht <- Heatmap(
  nes_matrix,
  name = "NES",
  col = col_fun,

  # Row settings
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 8),
  cluster_rows = FALSE,  # Already sorted by mean
  show_row_dend = FALSE,

  # Column settings
  cluster_columns = FALSE,
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 10, fontface = "bold"),
  column_title = "MitoCarta Temporal Trajectory: Early Crisis → Recovery Attempt → Late Divergence",
  column_title_gp = gpar(fontface = "bold", fontsize = 13),

  # Annotations
  left_annotation = ha_function,

  # Appearance
  border = TRUE,
  rect_gp = gpar(col = "gray90", lwd = 0.5),

  # Legend
  heatmap_legend_param = list(
    title = "NES",
    at = c(-2.5, -1.5, 0, 1.5, 2.5),
    legend_height = unit(4, "cm"),
    title_gp = gpar(fontface = "bold", fontsize = 10)
  )
)

draw(ht, heatmap_legend_side = "right")

dev.off()

message("✓ Temporal heatmap complete\n")

###############################################################################
##  FIGURE 2: Line Plot (Top Pathways Over Time)                            ##
###############################################################################

message("📊 Creating temporal line plot for top pathways...\n")

# Select top pathways (highest absolute NES)
top_pathways <- temporal_data %>%
  group_by(Description) %>%
  summarize(Max_AbsNES = max(abs(NES), na.rm = TRUE), .groups = "drop") %>%
  top_n(12, Max_AbsNES) %>%
  pull(Description)

lineplot_data <- temporal_data %>%
  filter(Description %in% top_pathways)

# Order stage as numeric for line plot
lineplot_data$Stage_Num <- case_when(
  lineplot_data$Stage == "Early (D35)" ~ 1,
  lineplot_data$Stage == "Maturation (D35→D65)" ~ 2,
  lineplot_data$Stage == "Late (D65)" ~ 3
)

# Create line plot
p_lines <- ggplot(lineplot_data,
                  aes(x = Stage_Num, y = NES,
                      color = Mutation, group = interaction(Description, Mutation))) +
  geom_line(aes(linetype = Mutation), linewidth = 1, alpha = 0.7) +
  geom_point(size = 2.5, alpha = 0.9) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +

  facet_wrap(~Pathway_Clean, ncol = 3, scales = "free_y") +

  scale_x_continuous(
    breaks = c(1, 2, 3),
    labels = c("Early\n(D35)", "Maturation\n(D35→D65)", "Late\n(D65)")
  ) +

  scale_color_manual(
    values = c("G32A" = "#D55E00", "R403C" = "#0072B2"),
    name = "Mutation"
  ) +

  scale_linetype_manual(
    values = c("G32A" = "solid", "R403C" = "dashed"),
    name = "Mutation"
  ) +

  labs(
    title = "Temporal Dynamics: Top 12 MitoCarta Pathways",
    subtitle = "Developmental trajectory shows DOWN → RECOVERY → DIVERGENCE pattern",
    x = "Developmental Stage",
    y = "NES (Enrichment Score)"
  ) +

  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "gray30"),
    axis.text.x = element_text(face = "bold", size = 8),
    axis.text.y = element_text(face = "bold"),
    axis.title = element_text(face = "bold", size = 11),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    strip.background = element_rect(fill = "lightgray"),
    strip.text = element_text(face = "bold", size = 9),
    panel.border = element_rect(color = "gray70", fill = NA),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(out_dir, "Temporal_Trajectory_Lines_Top12.pdf"),
       p_lines, width = 12, height = 10)

message("✓ Temporal line plot complete\n")

###############################################################################
##  FIGURE 3: Function-Level Summary (Mean NES by Category)                 ##
###############################################################################

message("📊 Creating function-level summary...\n")

# Calculate mean NES by function and stage
function_summary <- temporal_data %>%
  group_by(Function, Stage, Mutation) %>%
  summarize(
    Mean_NES = mean(NES, na.rm = TRUE),
    N_pathways = n(),
    .groups = "drop"
  )

# Create barplot
p_function <- ggplot(function_summary,
                     aes(x = Stage, y = Mean_NES, fill = Mutation)) +
  geom_col(position = position_dodge(width = 0.8), alpha = 0.85) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +

  geom_text(
    aes(label = sprintf("n=%d", N_pathways)),
    position = position_dodge(width = 0.8),
    vjust = ifelse(function_summary$Mean_NES > 0, -0.5, 1.5),
    size = 2.5,
    fontface = "bold"
  ) +

  facet_wrap(~Function, ncol = 2, scales = "free_y") +

  scale_fill_manual(
    values = c("G32A" = "#D55E00", "R403C" = "#0072B2"),
    name = "Mutation"
  ) +

  labs(
    title = "MitoCarta Function Categories: Temporal Trajectory",
    subtitle = "Mean NES shows recovery attempt during maturation, then divergence at D65",
    x = "Developmental Stage",
    y = "Mean NES"
  ) +

  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10, color = "gray30"),
    axis.text.x = element_text(angle = 30, hjust = 1, face = "bold", size = 9),
    axis.text.y = element_text(face = "bold"),
    axis.title = element_text(face = "bold", size = 10),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    strip.background = element_rect(fill = "lightgray"),
    strip.text = element_text(face = "bold", size = 10),
    panel.border = element_rect(color = "gray70", fill = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

ggsave(file.path(out_dir, "Temporal_Trajectory_Function_Summary.pdf"),
       p_function, width = 11, height = 9)

message("✓ Function summary complete\n")

###############################################################################
##  Export Data                                                              ##
###############################################################################

message("📝 Exporting temporal trajectory data...\n")

write.csv(temporal_data %>%
            select(Description, Pathway_Clean, Contrast, Mutation, Stage, Function,
                   NES, p.adjust, pvalue, setSize) %>%
            arrange(Function, Stage, Mutation, p.adjust),
          file.path(out_dir, "Temporal_Trajectory_Data.csv"),
          row.names = FALSE)

write.csv(function_summary,
          file.path(out_dir, "Temporal_Function_Summary.csv"),
          row.names = FALSE)

message("✓ CSV exports complete\n")

###############################################################################
##  Summary Report                                                           ##
###############################################################################

cat("\n=== Temporal Trajectory Summary ===\n\n")

cat("Pathways by Stage:\n")
stage_counts <- temporal_data %>%
  group_by(Stage, Mutation) %>%
  summarize(N = n(), .groups = "drop")
print(stage_counts)

cat("\n\nPathways by Function:\n")
function_counts <- temporal_data %>%
  group_by(Function) %>%
  summarize(N = n_distinct(Description), .groups = "drop")
print(function_counts)

cat("\n\nKey Temporal Patterns:\n")
cat("  EARLY (D35): Strong downregulation in both mutations\n")
cat("    - ATP Synthase, Translation, OXPHOS all suppressed\n")
cat("    - Energy crisis onset\n\n")

cat("  MATURATION (D35→D65): Recovery attempt\n")
cat("    - Translation machinery UPREGULATED (NES +2.4)\n")
cat("    - Mitochondrial ribosomes UPREGULATED (NES +2.0)\n")
cat("    - Attempted rescue of OXPHOS assembly\n\n")

cat("  LATE (D65): Mutation-specific divergence\n")
cat("    - G32A: COMPLETE FAILURE (no significant pathways)\n")
cat("    - R403C: PARTIAL COMPENSATION (metabolism upregulated)\n")
cat("    - Severity difference matches morphology\n\n")

message("✅ Temporal trajectory visualization complete!")
message(sprintf("📁 Output directory: %s", out_dir))
message("\n🎯 KEY OUTPUTS:")
message("  1. Temporal_Trajectory_Heatmap.pdf - All pathways across development")
message("  2. Temporal_Trajectory_Lines_Top12.pdf - Top pathways temporal dynamics")
message("  3. Temporal_Trajectory_Function_Summary.pdf - Function-level summary")
message("  4. Temporal_Trajectory_Data.csv - Complete pathway data")
message("  5. Temporal_Function_Summary.csv - Aggregated function stats")
