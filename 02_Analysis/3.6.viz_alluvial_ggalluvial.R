#!/usr/bin/env Rscript
###############################################################################
## Script: 3.6.viz_alluvial_ggalluvial.R
## Purpose: Classical alluvial diagram showing trajectory flow from early
##          defects through mechanisms to late outcomes
## Output: 03_Results/02_Analysis/Plots/Trajectory_Flow/
## Runtime: ~1-2 minutes
## Dependencies: ggalluvial, ggplot2, dplyr, tidyr, patchwork
###############################################################################

# ============================================================================
# Setup
# ============================================================================

library(here)
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(ggalluvial)
library(patchwork)

# Output directory (alluvial subfolder)
out_dir <- here("03_Results/02_Analysis/Plots/Trajectory_Flow/alluvial")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

cat("================================================================================\n")
cat("CLASSICAL ALLUVIAL DIAGRAM - ggalluvial\n")
cat("================================================================================\n")
cat("Output directory:", out_dir, "\n\n")

# ============================================================================
# Load and Prepare Data
# ============================================================================

cat("Loading data from gsea_results_wide.csv...\n")

# Load wide format data with NES and p.adjust per contrast
wide_file <- here("03_Results/02_Analysis/Python_exports/gsea_results_wide.csv")
df_wide <- read_csv(wide_file, show_col_types = FALSE)

# Also load pattern classifications from master table
master_file <- here("03_Results/02_Analysis/master_gsea_table.csv")
df_patterns <- read_csv(master_file, show_col_types = FALSE) %>%
  select(pathway_id, database, Description,
         starts_with("Pattern_"), starts_with("Confidence_")) %>%
  distinct(pathway_id, database, .keep_all = TRUE)

# Merge patterns into wide data
df_wide <- df_wide %>%
  left_join(df_patterns, by = c("pathway_id", "database", "Description"))

cat("  Loaded", nrow(df_wide), "unique pathways\n\n")

# ============================================================================
# Compute Derived Columns
# ============================================================================

# Define thresholds (from pattern_definitions.py)
PADJ_SIGNIFICANT <- 0.05
NES_EFFECT <- 0.5
IMPROVEMENT_RATIO <- 0.7
WORSENING_RATIO <- 1.3

compute_alluvial_data <- function(df, mutation) {
  # Column names from gsea_results_wide.csv
  # NES columns: NES_{contrast}
  # p.adjust columns: p.adjust_{contrast}
  if (mutation == "G32A") {
    early_nes_col <- "NES_G32A_vs_Ctrl_D35"
    trajdev_nes_col <- "NES_Maturation_G32A_specific"
    late_nes_col <- "NES_G32A_vs_Ctrl_D65"
    early_padj_col <- "p.adjust_G32A_vs_Ctrl_D35"
    trajdev_padj_col <- "p.adjust_Maturation_G32A_specific"
    late_padj_col <- "p.adjust_G32A_vs_Ctrl_D65"
  } else {
    early_nes_col <- "NES_R403C_vs_Ctrl_D35"
    trajdev_nes_col <- "NES_Maturation_R403C_specific"
    late_nes_col <- "NES_R403C_vs_Ctrl_D65"
    early_padj_col <- "p.adjust_R403C_vs_Ctrl_D35"
    trajdev_padj_col <- "p.adjust_Maturation_R403C_specific"
    late_padj_col <- "p.adjust_R403C_vs_Ctrl_D65"
  }
  pattern_col <- paste0("Pattern_", mutation)

  # Check which columns exist
  available_cols <- colnames(df)

  # Use available columns
  result <- df %>%
    mutate(
      # Extract NES values
      early_nes = .data[[early_nes_col]],
      trajdev_nes = .data[[trajdev_nes_col]],
      late_nes = .data[[late_nes_col]],
      pattern = .data[[pattern_col]],

      # Extract p.adjust values (with fallback to 1 if not available)
      early_padj = if (early_padj_col %in% available_cols) .data[[early_padj_col]] else 1,
      trajdev_padj = if (trajdev_padj_col %in% available_cols) .data[[trajdev_padj_col]] else 1,
      late_padj = if (late_padj_col %in% available_cols) .data[[late_padj_col]] else 1,

      # Has Early Defect: significant AND substantial effect
      has_early_defect = !is.na(early_padj) & early_padj < PADJ_SIGNIFICANT &
                         !is.na(early_nes) & abs(early_nes) > NES_EFFECT,

      # Mechanism: Active if TrajDev is significant
      is_active = !is.na(trajdev_padj) & trajdev_padj < PADJ_SIGNIFICANT &
                  !is.na(trajdev_nes) & abs(trajdev_nes) > NES_EFFECT,

      # Late Outcome based on |Late|/|Early| ratio
      late_early_ratio = ifelse(abs(early_nes) > 0.1, abs(late_nes) / abs(early_nes), NA),

      # Determine outcome
      outcome = case_when(
        pattern == "Late_onset" ~ "Late_onset",
        !has_early_defect ~ "No_Early_Defect",
        is.na(late_early_ratio) ~ "Improved",
        late_early_ratio < IMPROVEMENT_RATIO ~ "Improved",
        late_early_ratio > WORSENING_RATIO ~ "Worsened",
        abs(late_nes) < 0.3 & abs(early_nes) > NES_EFFECT ~ "Resolved",
        TRUE ~ "Improved"
      ),

      # Mechanism label
      mechanism = case_when(
        pattern == "Late_onset" ~ "Late_onset",
        !has_early_defect ~ "No_Early_Defect",
        is_active ~ "Active",
        TRUE ~ "Passive"
      ),

      # Early status
      early_status = case_when(
        pattern == "Late_onset" ~ "No_Early_Defect",
        has_early_defect ~ "Early_Defect",
        TRUE ~ "No_Early_Defect"
      )
    )

  return(result)
}

# ============================================================================
# Create Alluvial Data for Both Mutations
# ============================================================================

cat("Computing derived columns for G32A...\n")
df_g32a <- compute_alluvial_data(df_wide, "G32A")

cat("Computing derived columns for R403C...\n")
df_r403c <- compute_alluvial_data(df_wide, "R403C")

# Print summaries
cat("\n--- G32A Summary ---\n")
cat("Early Defect:", sum(df_g32a$has_early_defect, na.rm = TRUE), "\n")
cat("Active Mechanism:", sum(df_g32a$mechanism == "Active", na.rm = TRUE), "\n")
cat("Passive Mechanism:", sum(df_g32a$mechanism == "Passive", na.rm = TRUE), "\n")
cat("Late_onset:", sum(df_g32a$pattern == "Late_onset", na.rm = TRUE), "\n")
cat("Outcomes:\n")
print(table(df_g32a$outcome))

cat("\n--- R403C Summary ---\n")
cat("Early Defect:", sum(df_r403c$has_early_defect, na.rm = TRUE), "\n")
cat("Active Mechanism:", sum(df_r403c$mechanism == "Active", na.rm = TRUE), "\n")
cat("Passive Mechanism:", sum(df_r403c$mechanism == "Passive", na.rm = TRUE), "\n")
cat("Late_onset:", sum(df_r403c$pattern == "Late_onset", na.rm = TRUE), "\n")
cat("Outcomes:\n")
print(table(df_r403c$outcome))

# ============================================================================
# Prepare Data for ggalluvial (Focus on Early Defects + Late_onset)
# ============================================================================

prepare_alluvial_df <- function(df, mutation_name) {
  # Filter to pathways of interest: those with early defects OR Late_onset
  df_filtered <- df %>%
    filter(has_early_defect | pattern == "Late_onset") %>%
    select(pathway_id, early_status, mechanism, outcome) %>%
    mutate(
      # Clean up labels for display
      Early = case_when(
        early_status == "Early_Defect" ~ "Early Defect",
        TRUE ~ "No Early Defect"
      ),
      Mechanism = case_when(
        mechanism == "Active" ~ "Active\nCompensation",
        mechanism == "Passive" ~ "Passive\nBuffering",
        mechanism == "Late_onset" ~ "Late-onset\nPathway",
        TRUE ~ mechanism
      ),
      Outcome = case_when(
        outcome == "Improved" ~ "Improved",
        outcome == "Resolved" ~ "Fully\nResolved",
        outcome == "Worsened" ~ "Worsened",
        outcome == "Late_onset" ~ "Late-onset\nDefect",
        TRUE ~ outcome
      ),
      mutation = mutation_name
    ) %>%
    filter(!is.na(Mechanism) & !is.na(Outcome))

  return(df_filtered)
}

alluvial_g32a <- prepare_alluvial_df(df_g32a, "G32A")
alluvial_r403c <- prepare_alluvial_df(df_r403c, "R403C")

cat("\n--- Alluvial Data Prepared ---\n")
cat("G32A pathways in diagram:", nrow(alluvial_g32a), "\n")
cat("R403C pathways in diagram:", nrow(alluvial_r403c), "\n")

# ============================================================================
# Create Alluvial Plot Function
# ============================================================================

create_alluvial_plot <- function(df, mutation_name) {
  # Define color palette
  mechanism_colors <- c(
    "Active\nCompensation" = "#2166AC",    # Blue
    "Passive\nBuffering" = "#B2182B",      # Red
    "Late-onset\nPathway" = "#7570B3"      # Purple
  )

  # Count data for annotations
  n_total <- nrow(df)
  n_early <- sum(df$Early == "Early Defect")
  n_active <- sum(df$Mechanism == "Active\nCompensation")
  n_passive <- sum(df$Mechanism == "Passive\nBuffering")
  n_lateonset <- sum(df$Mechanism == "Late-onset\nPathway")
  n_improved <- sum(df$Outcome == "Improved")
  n_resolved <- sum(df$Outcome == "Fully\nResolved")
  n_worsened <- sum(df$Outcome == "Worsened")

  # Create the plot
  p <- ggplot(df, aes(axis1 = Early, axis2 = Mechanism, axis3 = Outcome)) +
    geom_alluvium(aes(fill = Mechanism), width = 1/8, alpha = 0.8) +
    geom_stratum(width = 1/8, fill = "gray95", color = "gray40", linewidth = 0.5) +
    geom_text(stat = "stratum",
              aes(label = after_stat(stratum)),
              size = 3.5, fontface = "bold") +
    scale_fill_manual(values = mechanism_colors, name = "Mechanism") +
    scale_x_discrete(limits = c("Early\nStatus", "Mechanism", "Late\nOutcome"),
                     expand = c(0.15, 0.05)) +
    labs(
      title = paste0(mutation_name, " Trajectory Flow"),
      subtitle = paste0("n = ", n_total, " pathways (Early defects: ", n_early,
                       ", Late-onset: ", n_lateonset, ")"),
      y = "Number of Pathways"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40"),
      axis.text.x = element_text(face = "bold", size = 11),
      axis.title.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank()
    )

  return(p)
}

# ============================================================================
# Generate Individual Plots
# ============================================================================

cat("\n========================================\n")
cat("Generating alluvial plots...\n")
cat("========================================\n")

# G32A plot
cat("\n1. Creating G32A alluvial plot...\n")
p_g32a <- create_alluvial_plot(alluvial_g32a, "G32A")

ggsave(file.path(out_dir, "alluvial_ggalluvial_G32A.pdf"),
       p_g32a, width = 10, height = 8, dpi = 300)
ggsave(file.path(out_dir, "alluvial_ggalluvial_G32A.png"),
       p_g32a, width = 10, height = 8, dpi = 300)
cat("   Saved: alluvial_ggalluvial_G32A.pdf/png\n")

# R403C plot
cat("\n2. Creating R403C alluvial plot...\n")
p_r403c <- create_alluvial_plot(alluvial_r403c, "R403C")

ggsave(file.path(out_dir, "alluvial_ggalluvial_R403C.pdf"),
       p_r403c, width = 10, height = 8, dpi = 300)
ggsave(file.path(out_dir, "alluvial_ggalluvial_R403C.png"),
       p_r403c, width = 10, height = 8, dpi = 300)
cat("   Saved: alluvial_ggalluvial_R403C.pdf/png\n")

# ============================================================================
# Generate Combined Figure
# ============================================================================

cat("\n3. Creating combined figure...\n")

p_combined <- p_g32a + p_r403c +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Pathway Trajectory Flow: Early Defects to Late Outcomes",
    subtitle = "DRP1 mutations show predominantly compensatory responses during neuronal maturation",
    theme = theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40")
    )
  ) &
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "alluvial_combined.pdf"),
       p_combined, width = 16, height = 8, dpi = 300)
ggsave(file.path(out_dir, "alluvial_combined.png"),
       p_combined, width = 16, height = 8, dpi = 300)
cat("   Saved: alluvial_combined.pdf/png\n")

# ============================================================================
# Summary Statistics
# ============================================================================

cat("\n================================================================================\n")
cat("SUMMARY STATISTICS\n")
cat("================================================================================\n")

for (mutation in c("G32A", "R403C")) {
  df_mut <- if (mutation == "G32A") alluvial_g32a else alluvial_r403c

  cat("\n", mutation, ":\n", sep = "")
  cat("  Total pathways in diagram: ", nrow(df_mut), "\n", sep = "")

  # Early status
  cat("  Early Status:\n")
  for (status in unique(df_mut$Early)) {
    n <- sum(df_mut$Early == status)
    pct <- round(100 * n / nrow(df_mut), 1)
    cat("    ", status, ": ", n, " (", pct, "%)\n", sep = "")
  }

  # Mechanism breakdown
  cat("  Mechanism:\n")
  for (mech in unique(df_mut$Mechanism)) {
    n <- sum(df_mut$Mechanism == mech)
    pct <- round(100 * n / nrow(df_mut), 1)
    cat("    ", gsub("\n", " ", mech), ": ", n, " (", pct, "%)\n", sep = "")
  }

  # Outcome breakdown
  cat("  Outcome:\n")
  for (outcome in unique(df_mut$Outcome)) {
    n <- sum(df_mut$Outcome == outcome)
    pct <- round(100 * n / nrow(df_mut), 1)
    cat("    ", gsub("\n", " ", outcome), ": ", n, " (", pct, "%)\n", sep = "")
  }
}

cat("\n================================================================================\n")
cat("COMPLETE!\n")
cat("================================================================================\n")
cat("\nOutput files in:", out_dir, "\n")
cat("  - alluvial_ggalluvial_G32A.pdf/png\n")
cat("  - alluvial_ggalluvial_R403C.pdf/png\n")
cat("  - alluvial_combined.pdf/png\n")
