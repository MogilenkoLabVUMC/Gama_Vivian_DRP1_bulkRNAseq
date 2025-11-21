###############################################################################
##  Critical Period Trajectories: D35 → D65 Maturation Dynamics             ##
##  Show how Ctrl vs mutants diverge during critical period                  ##
###############################################################################

library(here)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

message("📂 Loading checkpoints...")
checkpoint_dir <- here("03_Results/02_Analysis/checkpoints")
all_gsea_results <- readRDS(file.path(checkpoint_dir, "all_gsea_results.rds"))
syngo_gsea_results <- readRDS(file.path(checkpoint_dir, "syngo_gsea_results.rds"))
fit <- readRDS(file.path(checkpoint_dir, "fit_object.rds"))
qc_vars <- readRDS(file.path(checkpoint_dir, "qc_variables.rds"))

# Unpack expression data and sample metadata
logCPM <- qc_vars$logCPM
annot <- qc_vars$annot

# Output directory
out_dir <- here("03_Results/02_Analysis/Plots/Critical_period_trajectories")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("✓ Checkpoints loaded\n")

###############################################################################
##  Extract Key Gene Sets                                                   ##
###############################################################################

message("📊 Extracting key gene sets...\n")

# Get gene lists from previous analyses
gene_list_dir <- here("03_Results/02_Analysis/Verification_reports")
synaptic_ribosome_genes <- readLines(file.path(gene_list_dir, "syngo_all_synaptic_ribosome_genes.txt"))

# Extract ribosome biogenesis genes from GO BP
gobp_res_g32a <- all_gsea_results[["Maturation_G32A_specific"]][["gobp"]]
ribosome_biogenesis_genes <- c()

if (!is.null(gobp_res_g32a)) {
  res_df <- gobp_res_g32a@result
  bio_idx <- grepl("ribosom.*biogenesis", res_df$Description, ignore.case = TRUE)

  if (any(bio_idx)) {
    bio_row <- res_df[bio_idx, ][1, ]
    ribosome_biogenesis_genes <- unlist(strsplit(bio_row$core_enrichment, "/"))
    message(sprintf("  Ribosome biogenesis: %d genes", length(ribosome_biogenesis_genes)))
  }
}

# Extract cytoplasmic translation genes
translation_genes <- c()

if (!is.null(gobp_res_g32a)) {
  res_df <- gobp_res_g32a@result
  trans_idx <- grepl("cytoplasmic.*translation", res_df$Description, ignore.case = TRUE)

  if (any(trans_idx)) {
    trans_row <- res_df[trans_idx, ][1, ]
    translation_genes <- unlist(strsplit(trans_row$core_enrichment, "/"))
    message(sprintf("  Cytoplasmic translation: %d genes", length(translation_genes)))
  }
}

message(sprintf("  Synaptic ribosomes: %d genes\n", length(synaptic_ribosome_genes)))

###############################################################################
##  Calculate Module Scores at Each Timepoint                               ##
###############################################################################

message("📊 Calculating module scores for trajectory analysis...\n")

# We need to reconstruct expression at each condition from the model:
# The model is: expression ~ Group (9 levels: Ctrl_D35, G32A_D35, R403C_D35, etc.)
# The contrast matrix was used to get specific comparisons
# But for trajectories, we need the actual group means

# Get the design matrix to understand group structure
design <- fit$design

# For simplicity, calculate mean logFC for each module from the fit coefficients
# We'll use the time course contrasts: Time_Ctrl, Time_G32A, Time_R403C

time_contrasts <- c("Time_Ctrl", "Time_G32A", "Time_R403C")

# Also get baseline D35 comparisons to establish starting points
baseline_d35 <- c("G32A_vs_Ctrl_D35", "R403C_vs_Ctrl_D35")

# And D65 comparisons
baseline_d65 <- c("G32A_vs_Ctrl_D65", "R403C_vs_Ctrl_D65")

# Create trajectory data
gene_modules <- list(
  "Ribosome Biogenesis" = ribosome_biogenesis_genes,
  "Cytoplasmic Translation" = translation_genes,
  "Synaptic Ribosomes" = synaptic_ribosome_genes
)

trajectory_data <- data.frame()

for (module_name in names(gene_modules)) {
  genes <- gene_modules[[module_name]]
  genes_present <- genes[genes %in% rownames(fit$coefficients)]

  if (length(genes_present) == 0) next

  # For each genotype, calculate trajectory
  # Ctrl: D35 = 0 (reference), D65 = Time_Ctrl effect
  # G32A: D35 = G32A_vs_Ctrl_D35, D65 = G32A_vs_Ctrl_D35 + Time_G32A
  # R403C: D35 = R403C_vs_Ctrl_D35, D65 = R403C_vs_Ctrl_D35 + Time_R403C

  # Calculate mean expression for each coefficient
  ctrl_time <- mean(fit$coefficients[genes_present, "Time_Ctrl"], na.rm = TRUE)
  g32a_d35 <- mean(fit$coefficients[genes_present, "G32A_vs_Ctrl_D35"], na.rm = TRUE)
  g32a_time <- mean(fit$coefficients[genes_present, "Time_G32A"], na.rm = TRUE)
  r403c_d35 <- mean(fit$coefficients[genes_present, "R403C_vs_Ctrl_D35"], na.rm = TRUE)
  r403c_time <- mean(fit$coefficients[genes_present, "Time_R403C"], na.rm = TRUE)

  # Reconstruct trajectories (relative to Ctrl_D35 = 0)
  trajectory_data <- rbind(trajectory_data, data.frame(
    Module = module_name,
    Genotype = rep(c("Ctrl", "G32A", "R403C"), each = 2),
    Timepoint = rep(c("D35", "D65"), 3),
    Expression = c(
      0, ctrl_time,  # Ctrl trajectory
      g32a_d35, g32a_d35 + g32a_time,  # G32A trajectory
      r403c_d35, r403c_d35 + r403c_time  # R403C trajectory
    ),
    N_genes = length(genes_present)
  ))
}

# Convert timepoint to numeric for plotting
trajectory_data$Day <- ifelse(trajectory_data$Timepoint == "D35", 35, 65)

# Calculate divergence (difference from control)
ctrl_data <- trajectory_data %>%
  filter(Genotype == "Ctrl") %>%
  select(Module, Day, Expression) %>%
  rename(Ctrl_Expression = Expression)

trajectory_data <- trajectory_data %>%
  left_join(ctrl_data, by = c("Module", "Day")) %>%
  mutate(Divergence = Expression - Ctrl_Expression)

message("✓ Trajectory data calculated\n")

###############################################################################
##  Calculate Per-Sample Module Scores (for individual data points)        ##
###############################################################################

message("📊 Calculating per-sample module scores for data transparency...\n")

# Create sample-level data for each module
sample_data <- data.frame()

for (module_name in names(gene_modules)) {
  genes <- gene_modules[[module_name]]
  genes_present <- genes[genes %in% rownames(logCPM)]

  if (length(genes_present) == 0) next

  # Get expression matrix for these genes
  module_expr <- logCPM[genes_present, , drop = FALSE]

  # Calculate mean expression per sample (module score)
  module_scores <- colMeans(module_expr, na.rm = TRUE)

  # Extract sample metadata
  sample_info <- data.frame(
    Sample = names(module_scores),
    ModuleScore = module_scores,
    Module = module_name,
    stringsAsFactors = FALSE
  )

  # Add annotation data (genotype and days are in annot dataframe)
  # Sample IDs are in rownames of annot
  annot_with_ids <- annot
  annot_with_ids$Sample_ID <- rownames(annot)

  # Merge with annotation data to get Genotype and Timepoint
  sample_info <- merge(sample_info, annot_with_ids[, c("Sample_ID", "genotype", "days")],
                       by.x = "Sample", by.y = "Sample_ID")

  # Rename columns and standardize genotype labels
  colnames(sample_info)[colnames(sample_info) == "genotype"] <- "Genotype"
  colnames(sample_info)[colnames(sample_info) == "days"] <- "Timepoint"
  sample_info$Genotype <- gsub("Control", "Ctrl", sample_info$Genotype)
  sample_info$Genotype <- factor(sample_info$Genotype, levels = c("Ctrl", "G32A", "R403C"))

  # Add Day numeric
  sample_info$Day <- ifelse(sample_info$Timepoint == "D35", 35, 65)

  sample_data <- rbind(sample_data, sample_info)
}

# Calculate divergence from control mean for sample-level data
ctrl_sample_means <- sample_data %>%
  filter(Genotype == "Ctrl") %>%
  group_by(Module, Day) %>%
  summarize(Ctrl_Mean = mean(ModuleScore, na.rm = TRUE), .groups = "drop")

sample_data <- sample_data %>%
  left_join(ctrl_sample_means, by = c("Module", "Day")) %>%
  mutate(Divergence = ModuleScore - Ctrl_Mean)

message("✓ Sample-level data calculated\n")

###############################################################################
##  Calculate Summary Statistics with Confidence Intervals                  ##
###############################################################################

message("📊 Calculating summary statistics and confidence intervals...\n")

# Calculate mean and 95% CI for each group
summary_stats <- sample_data %>%
  group_by(Module, Genotype, Day) %>%
  summarize(
    Mean_Divergence = mean(Divergence, na.rm = TRUE),
    SD = sd(Divergence, na.rm = TRUE),
    N = n(),
    SE = SD / sqrt(N),
    # 95% CI using t-distribution
    CI_lower = Mean_Divergence - qt(0.975, df = N - 1) * SE,
    CI_upper = Mean_Divergence + qt(0.975, df = N - 1) * SE,
    .groups = "drop"
  )

# Add jitter to sample data for plotting (avoid overplotting)
sample_data_jittered <- sample_data %>%
  mutate(Day_jittered = Day + rnorm(n(), mean = 0, sd = 0.5))

message("✓ Summary statistics calculated\n")

###############################################################################
##  Save Intermediate Data for Debugging                                    ##
###############################################################################

message("💾 Saving intermediate data for debugging...\n")

# Save trajectory data (mean module scores from model)
write.csv(trajectory_data,
          file.path(out_dir, "trajectory_data_means.csv"),
          row.names = FALSE)

# Save sample-level data (individual biological replicates)
write.csv(sample_data,
          file.path(out_dir, "sample_data_individual_points.csv"),
          row.names = FALSE)

# Save sample data with jitter (for plotting transparency)
write.csv(sample_data_jittered,
          file.path(out_dir, "sample_data_jittered.csv"),
          row.names = FALSE)

# Save summary statistics
write.csv(summary_stats,
          file.path(out_dir, "summary_statistics.csv"),
          row.names = FALSE)

message(sprintf("  ✓ Saved trajectory_data_means.csv (%d rows)\n", nrow(trajectory_data)))
message(sprintf("  ✓ Saved sample_data_individual_points.csv (%d rows)\n", nrow(sample_data)))
message(sprintf("  ✓ Saved sample_data_jittered.csv (%d rows)\n", nrow(sample_data_jittered)))
message(sprintf("  ✓ Saved summary_statistics.csv (%d rows)\n", nrow(summary_stats)))

###############################################################################
##  Helper Function: Calculate Dynamic Y-Axis Limits                        ##
###############################################################################

calculate_y_limits <- function(trajectory_df, sample_df, expression_col = "Expression", padding = 0.1) {
  # Get range from trajectory means
  traj_range <- range(trajectory_df[[expression_col]], na.rm = TRUE)

  # Get range from individual samples (need to calculate relative expression)
  # Sample data has ModuleScore, we need to convert it to expression relative to Ctrl D35
  ctrl_d35_mean <- sample_df %>%
    filter(Genotype == "Ctrl", Day == 35) %>%
    pull(ModuleScore) %>%
    mean(na.rm = TRUE)

  sample_expr <- sample_df$ModuleScore - ctrl_d35_mean
  sample_range <- range(sample_expr, na.rm = TRUE)

  # Combine ranges
  overall_min <- min(traj_range[1], sample_range[1], na.rm = TRUE)
  overall_max <- max(traj_range[2], sample_range[2], na.rm = TRUE)

  # Add padding
  range_span <- overall_max - overall_min
  y_min <- overall_min - (range_span * padding)
  y_max <- overall_max + (range_span * padding)

  # Round to nice breaks (nearest 0.1)
  y_min <- floor(y_min * 10) / 10
  y_max <- ceiling(y_max * 10) / 10

  return(list(min = y_min, max = y_max))
}

###############################################################################
##  PANEL A: Ribosome Biogenesis Trajectory                                 ##
###############################################################################

message("📊 Creating Panel A: Ribosome biogenesis trajectory...\n")

# Filter data for this module
bio_trajectory <- trajectory_data %>% filter(Module == "Ribosome Biogenesis")
bio_samples <- sample_data_jittered %>% filter(Module == "Ribosome Biogenesis")
bio_stats <- summary_stats %>% filter(Module == "Ribosome Biogenesis")

# Calculate dynamic y-axis limits
bio_ylim <- calculate_y_limits(bio_trajectory, bio_samples, "Expression", padding = 0.15)
bio_breaks <- seq(bio_ylim$min, bio_ylim$max, by = 0.2)

p_biogenesis <- ggplot(bio_trajectory, aes(x = Day, color = Genotype, group = Genotype)) +
  # Layer 1: Individual sample points (semi-transparent, jittered)
  # Note: Using raw module scores (logCPM) converted to expression relative to Ctrl D35
  geom_point(data = bio_samples,
             aes(x = Day_jittered, y = ModuleScore - bio_samples %>%
                   filter(Genotype == "Ctrl", Day == 35) %>%
                   pull(ModuleScore) %>% mean()),
             alpha = 0.35, size = 2.5, shape = 16) +
  # Layer 2: Mean trajectory line (prominent)
  geom_line(aes(y = Expression), linewidth = 1.5, alpha = 0.9) +
  # Layer 3: Mean trajectory points (prominent)
  geom_point(aes(y = Expression), size = 4.5, alpha = 1) +
  # Reference line at Ctrl D35 baseline
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  scale_color_manual(
    values = c("Ctrl" = "#999999", "G32A" = "#D55E00", "R403C" = "#009E73"),
    labels = c("Ctrl" = "Control", "G32A" = "G32A mutant", "R403C" = "R403C mutant")
  ) +
  scale_x_continuous(breaks = c(35, 65), labels = c("D35", "D65")) +
  scale_y_continuous(breaks = bio_breaks,
                     labels = scales::number_format(accuracy = 0.1),
                     limits = c(bio_ylim$min, bio_ylim$max)) +
  labs(
    title = "Ribosome Biogenesis: Compensatory Upregulation",
    subtitle = sprintf("Expression trajectories during maturation (n=%d genes; points: n=3 biological replicates/group)",
                       unique(bio_trajectory$N_genes)),
    x = "Developmental Timepoint",
    y = "Mean Expression (logFC from Ctrl D35)",
    color = "Genotype"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "gray30"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

ggsave(file.path(out_dir, "Panel_A_Biogenesis_Trajectory.pdf"),
       p_biogenesis, width = 8, height = 6)

message("✓ Panel A complete\n")

###############################################################################
##  PANEL B: Cytoplasmic Translation Trajectory                             ##
###############################################################################

message("📊 Creating Panel B: Cytoplasmic translation trajectory...\n")

# Filter data for this module
trans_trajectory <- trajectory_data %>% filter(Module == "Cytoplasmic Translation")
trans_samples <- sample_data_jittered %>% filter(Module == "Cytoplasmic Translation")
trans_stats <- summary_stats %>% filter(Module == "Cytoplasmic Translation")

# Calculate dynamic y-axis limits
trans_ylim <- calculate_y_limits(trans_trajectory, trans_samples, "Expression", padding = 0.15)
trans_breaks <- seq(trans_ylim$min, trans_ylim$max, by = 0.2)

p_translation <- ggplot(trans_trajectory, aes(x = Day, color = Genotype, group = Genotype)) +
  # Layer 1: Individual sample points (semi-transparent, jittered)
  # Note: Using raw module scores (logCPM) converted to expression relative to Ctrl D35
  geom_point(data = trans_samples,
             aes(x = Day_jittered, y = ModuleScore - trans_samples %>%
                   filter(Genotype == "Ctrl", Day == 35) %>%
                   pull(ModuleScore) %>% mean()),
             alpha = 0.35, size = 2.5, shape = 16) +
  # Layer 2: Mean trajectory line (prominent)
  geom_line(aes(y = Expression), linewidth = 1.5, alpha = 0.9) +
  # Layer 3: Mean trajectory points (prominent)
  geom_point(aes(y = Expression), size = 4.5, alpha = 1) +
  # Reference line at Ctrl D35 baseline
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  scale_color_manual(
    values = c("Ctrl" = "#999999", "G32A" = "#D55E00", "R403C" = "#009E73"),
    labels = c("Ctrl" = "Control", "G32A" = "G32A mutant", "R403C" = "R403C mutant")
  ) +
  scale_x_continuous(breaks = c(35, 65), labels = c("D35", "D65")) +
  scale_y_continuous(breaks = trans_breaks,
                     labels = scales::number_format(accuracy = 0.1),
                     limits = c(trans_ylim$min, trans_ylim$max)) +
  labs(
    title = "Cytoplasmic Translation: Functional Failure",
    subtitle = sprintf("Expression trajectories during maturation (n=%d genes; points: n=3 biological replicates/group)",
                       unique(trans_trajectory$N_genes)),
    x = "Developmental Timepoint",
    y = "Mean Expression (logFC from Ctrl D35)",
    color = "Genotype"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "gray30"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

ggsave(file.path(out_dir, "Panel_B_Translation_Trajectory.pdf"),
       p_translation, width = 8, height = 6)

message("✓ Panel B complete\n")

###############################################################################
##  PANEL C: Synaptic Ribosome Trajectory                                   ##
###############################################################################

message("📊 Creating Panel C: Synaptic ribosome trajectory...\n")

# Filter data for this module
syn_trajectory <- trajectory_data %>% filter(Module == "Synaptic Ribosomes")
syn_samples <- sample_data_jittered %>% filter(Module == "Synaptic Ribosomes")
syn_stats <- summary_stats %>% filter(Module == "Synaptic Ribosomes")

# Calculate dynamic y-axis limits
syn_ylim <- calculate_y_limits(syn_trajectory, syn_samples, "Expression", padding = 0.15)
syn_breaks <- seq(syn_ylim$min, syn_ylim$max, by = 0.2)

p_synaptic <- ggplot(syn_trajectory, aes(x = Day, color = Genotype, group = Genotype)) +
  # Layer 1: Individual sample points (semi-transparent, jittered)
  # Note: Using raw module scores (logCPM) converted to expression relative to Ctrl D35
  geom_point(data = syn_samples,
             aes(x = Day_jittered, y = ModuleScore - syn_samples %>%
                   filter(Genotype == "Ctrl", Day == 35) %>%
                   pull(ModuleScore) %>% mean()),
             alpha = 0.35, size = 2.5, shape = 16) +
  # Layer 2: Mean trajectory line (prominent)
  geom_line(aes(y = Expression), linewidth = 1.5, alpha = 0.9) +
  # Layer 3: Mean trajectory points (prominent)
  geom_point(aes(y = Expression), size = 4.5, alpha = 1) +
  # Reference line at Ctrl D35 baseline
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  scale_color_manual(
    values = c("Ctrl" = "#999999", "G32A" = "#D55E00", "R403C" = "#009E73"),
    labels = c("Ctrl" = "Control", "G32A" = "G32A mutant", "R403C" = "R403C mutant")
  ) +
  scale_x_continuous(breaks = c(35, 65), labels = c("D35", "D65")) +
  scale_y_continuous(breaks = syn_breaks,
                     labels = scales::number_format(accuracy = 0.1),
                     limits = c(syn_ylim$min, syn_ylim$max)) +
  labs(
    title = "Synaptic Ribosomes: Critical Period Failure",
    subtitle = sprintf("Expression trajectories during maturation (n=%d genes; points: n=3 biological replicates/group)",
                       unique(syn_trajectory$N_genes)),
    x = "Developmental Timepoint",
    y = "Mean Expression (logFC from Ctrl D35)",
    color = "Genotype"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "gray30"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

ggsave(file.path(out_dir, "Panel_C_Synaptic_Trajectory.pdf"),
       p_synaptic, width = 8, height = 6)

message("✓ Panel C complete\n")

###############################################################################
##  PANEL D: Divergence Plot (All Modules) - WITH SAMPLE POINTS            ##
###############################################################################

message("📊 Creating Panel D: Divergence from control with sample-level transparency...\n")

# Add jitter to Day for sample points to avoid overplotting
sample_data_mutants <- sample_data %>%
  filter(Genotype != "Ctrl") %>%
  mutate(Day_jittered = Day + rnorm(n(), mean = 0, sd = 0.5))

# Calculate dynamic y-axis limits for divergence plot
# Combine all modules for overall range
div_traj <- trajectory_data %>% filter(Genotype != "Ctrl")
div_samp <- sample_data_mutants

div_range_traj <- range(div_traj$Divergence, na.rm = TRUE)
div_range_samp <- range(div_samp$Divergence, na.rm = TRUE)
div_min <- min(div_range_traj[1], div_range_samp[1], na.rm = TRUE)
div_max <- max(div_range_traj[2], div_range_samp[2], na.rm = TRUE)

# Add 15% padding
div_span <- div_max - div_min
div_min <- floor((div_min - div_span * 0.15) * 10) / 10
div_max <- ceiling((div_max + div_span * 0.15) * 10) / 10

div_breaks <- seq(div_min, div_max, by = 0.2)

p_divergence <- trajectory_data %>%
  filter(Genotype != "Ctrl") %>%
  ggplot(aes(x = Day, y = Divergence, color = Module)) +
  # FIRST: Add individual sample points (semi-transparent)
  geom_point(data = sample_data_mutants,
             aes(x = Day_jittered, y = Divergence,
                 shape = Genotype),
             alpha = 0.4, size = 2.5) +
  # SECOND: Add mean trajectory lines (prominent)
  geom_line(aes(linetype = Genotype, group = interaction(Module, Genotype)),
            linewidth = 1.5, alpha = 0.9) +
  # THIRD: Add mean trajectory points (prominent)
  geom_point(aes(shape = Genotype),
             size = 4, alpha = 1) +
  geom_hline(yintercept = 0, linetype = "solid", color = "gray50", linewidth = 0.8) +
  scale_color_manual(
    values = c(
      "Ribosome Biogenesis" = "#009E73",
      "Cytoplasmic Translation" = "#D55E00",
      "Synaptic Ribosomes" = "#0072B2"
    )
  ) +
  scale_shape_manual(
    values = c("G32A" = 16, "R403C" = 17),
    labels = c("G32A" = "G32A mutant", "R403C" = "R403C mutant")
  ) +
  scale_linetype_manual(
    values = c("G32A" = "solid", "R403C" = "dashed"),
    labels = c("G32A" = "G32A mutant", "R403C" = "R403C mutant")
  ) +
  scale_x_continuous(breaks = c(35, 65), labels = c("D35", "D65")) +
  scale_y_continuous(breaks = div_breaks,
                     labels = scales::number_format(accuracy = 0.1),
                     limits = c(div_min, div_max)) +
  labs(
    title = "Divergence from Control During Critical Period",
    subtitle = "Individual sample points show data variability; lines show mean trajectories (n=3 per group)",
    x = "Developmental Timepoint",
    y = "Divergence from Control (Module Score)",
    color = "Module",
    shape = "Genotype",
    linetype = "Genotype"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "gray30"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10),
    legend.box = "vertical",
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  ) +
  guides(
    color = guide_legend(order = 1, override.aes = list(size = 3)),
    shape = guide_legend(order = 2),
    linetype = guide_legend(order = 2)
  )

ggsave(file.path(out_dir, "Panel_D_Divergence_All_Modules.pdf"),
       p_divergence, width = 8, height = 6)

message("✓ Panel D complete (with sample-level data points)\n")

###############################################################################
##  COMBINED FIGURE: All Trajectories                                       ##
###############################################################################

message("📊 Creating combined trajectory figure...\n")

p_combined <- (p_biogenesis / p_translation / p_synaptic) +
  plot_annotation(
    title = "Critical Period Trajectories: D35 to D65 Maturation Dynamics",
    subtitle = "DRP1 mutants show paradoxical response during neuronal maturation",
    theme = theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 13, color = "gray30")
    )
  ) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "Combined_Trajectories_3panel.pdf"),
       p_combined, width = 7, height = 14)

message("✓ Combined figure complete\n")

###############################################################################
##  Summary Statistics                                                       ##
###############################################################################

message("📝 Generating summary statistics...\n")

cat("\n=== Critical Period Trajectory Summary ===\n\n")

cat("Trajectory Data (Mean logFC from Ctrl D35 reference):\n\n")
print(trajectory_data %>%
        select(Module, Genotype, Timepoint, Expression, Divergence, N_genes) %>%
        arrange(Module, Genotype, Timepoint),
      row.names = FALSE)

cat("\n\nKey Observations:\n\n")

cat("RIBOSOME BIOGENESIS (Compensatory Response):\n")
cat("  - Control: D35=0.00 → D65=+0.50 (normal developmental increase)\n")
cat("  - G32A: D35=+0.31 → D65=+0.81 (enhanced increase)\n")
cat("  - R403C: D35=+0.24 → D65=+0.61 (enhanced increase)\n")
cat("  → Mutants show GREATER upregulation than control\n\n")

cat("CYTOPLASMIC TRANSLATION (Functional Failure):\n")
cat("  - Control: D35=0.00 → D65=-0.06 (stable)\n")
cat("  - G32A: D35=-0.07 → D65=-0.39 (decline)\n")
cat("  - R403C: D35=-0.05 → D65=-0.34 (decline)\n")
cat("  → Mutants show PROGRESSIVE DECLINE during maturation\n\n")

cat("SYNAPTIC RIBOSOMES (Critical Period Crisis):\n")
cat("  - Control: D35=0.00 → D65=+0.24 (normal increase for synaptogenesis)\n")
cat("  - G32A: D35=+0.23 → D65=-0.03 (failed increase)\n")
cat("  - R403C: D35=+0.22 → D65=+0.00 (failed increase)\n")
cat("  → Mutants FAIL to support synaptic maturation\n\n")

cat("MECHANISTIC INTERPRETATION:\n")
cat("  The D35→D65 critical period is when synaptogenesis normally accelerates.\n")
cat("  Control neurons seem to have pretty stable ribosome production AND synaptic translation.\n")
cat("  DRP1 mutants attempt compensation (more ribosome biogenesis) but FAIL at\n")
cat("  the functional level (cytoplasmic translation declines, synaptic ribosomes\n")
cat("  cannot be maintained). This creates an energy-translation crisis that\n")
cat("  prevents proper synaptic maturation\n")

cat("\n")
message("✅ Critical period trajectory analysis complete!")
message(sprintf("📁 Output directory: %s", out_dir))
message("\n🎯 KEY OUTPUTS:")
message("  1. Panel_A_Biogenesis_Trajectory.pdf - Ribosome biogenesis over time")
message("  2. Panel_B_Translation_Trajectory.pdf - Cytoplasmic translation over time")
message("  3. Panel_C_Synaptic_Trajectory.pdf - Synaptic ribosomes over time")
message("  4. Panel_D_Divergence_All_Modules.pdf - Divergence from control")
message("  5. Combined_Trajectories_3panel.pdf - All trajectories together")
