###############################################################################
##  Critical Period Trajectories: D35 → D65 Maturation Dynamics             ##
##  Extended with MitoCarta Compensation Pathways                            ##
###############################################################################
##  PURPOSE: Visualize how cytoplasmic, synaptic, and mitochondrial gene    ##
##           expression changes during the critical developmental period     ##
##           (D35 to D65) in control vs. DRP1 mutant neurons.               ##
##                                                                           ##
##  STORY: The paradoxical response - mutants increase ribosome production  ##
##         but fail at functional translation, while mitochondrial pathways ##
##         show compensatory upregulation.                                   ##
###############################################################################

library(here)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

message("=" |> rep(78) |> paste(collapse = ""))
message("Critical Period Trajectories with MitoCarta Pathways")
message("=" |> rep(78) |> paste(collapse = ""))

# ============================================================================ #
# 1. Load Checkpoint Data                                                      #
# ============================================================================ #

message("\n\U0001F4C2 Loading checkpoints...")
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

message("\U2713 Checkpoints loaded\n")

# ============================================================================ #
# 2. Load MitoCarta Gene Sets                                                  #
# ============================================================================ #

message("\U0001F9EC Loading MitoCarta pathways...")

# Source MitoCarta parser
source(here("01_Scripts/R_scripts/parse_mitocarta_gmx.R"))

# Parse MitoCarta GMX file
mitocarta_path <- here("00_Data/MitoCarta_3.0/MitoPathways3.0.gmx")
mitocarta_lists <- mitocarta_gmt(mitocarta_path)
T2G <- mitocarta_lists$T2G

message("\U2713 MitoCarta loaded\n")

# ============================================================================ #
# 3. Extract Gene Lists for All 7 Modules                                      #
# ============================================================================ #

message("\U0001F4CA Extracting gene sets for all modules...\n")

# --- Existing modules (from GSEA enrichment) ---

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

message(sprintf("  Synaptic ribosomes: %d genes", length(synaptic_ribosome_genes)))

# --- MitoCarta modules (from GMX file) ---

# Extract gene lists for the 4 core compensation pathways
# Note: MitoCarta uses hierarchical names like "MitoCarta3.0_MitoPathways.Mitochondrial_ribosome"

# Helper function to extract genes for a pathway
extract_mitocarta_genes <- function(T2G, pathway_pattern) {
  matching_rows <- grepl(pathway_pattern, T2G$gs_name, ignore.case = TRUE)
  genes <- unique(T2G$gene_symbol[matching_rows])
  return(genes)
}

# Extract the 4 core MitoCarta pathways
# Note: MitoCarta uses hierarchical names, but some top-level pathways are simple names
mito_ribosome_genes <- extract_mitocarta_genes(T2G, "Mitochondrial_ribosome$")
mito_ribo_assembly_genes <- extract_mitocarta_genes(T2G, "Mitochondrial_ribosome_assembly")
mtdna_maintenance_genes <- extract_mitocarta_genes(T2G, "mtDNA_maintenance")
oxphos_genes <- extract_mitocarta_genes(T2G, "^OXPHOS$")  # Top-level OXPHOS pathway

message(sprintf("  Mitochondrial ribosome: %d genes", length(mito_ribosome_genes)))
message(sprintf("  Mito ribosome assembly: %d genes", length(mito_ribo_assembly_genes)))
message(sprintf("  mtDNA maintenance: %d genes", length(mtdna_maintenance_genes)))
message(sprintf("  OXPHOS: %d genes", length(oxphos_genes)))

# ============================================================================ #
# 4. Define All Gene Modules                                                   #
# ============================================================================ #

gene_modules <- list(
  # Existing modules
  "Ribosome Biogenesis" = ribosome_biogenesis_genes,
  "Cytoplasmic Translation" = translation_genes,
  "Synaptic Ribosomes" = synaptic_ribosome_genes,
  # MitoCarta modules
  "Mitochondrial Ribosome" = mito_ribosome_genes,
  "Mito Ribosome Assembly" = mito_ribo_assembly_genes,
  "mtDNA Maintenance" = mtdna_maintenance_genes,
  "OXPHOS" = oxphos_genes
)

message("\n\U2713 All 7 gene modules defined\n")

# ============================================================================ #
# 5. Calculate Module Scores at Each Timepoint                                 #
# ============================================================================ #

message("\U0001F4CA Calculating module scores for trajectory analysis...\n")

# Reconstruction logic:
# Ctrl: D35 = 0 (reference), D65 = Time_Ctrl effect
# G32A: D35 = G32A_vs_Ctrl_D35, D65 = G32A_vs_Ctrl_D35 + Time_G32A
# R403C: D35 = R403C_vs_Ctrl_D35, D65 = R403C_vs_Ctrl_D35 + Time_R403C

trajectory_data <- data.frame()

for (module_name in names(gene_modules)) {
  genes <- gene_modules[[module_name]]
  genes_present <- genes[genes %in% rownames(fit$coefficients)]

  if (length(genes_present) == 0) {
    message(sprintf("  WARNING: No genes found for %s", module_name))
    next
  }

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

# Calculate divergence (difference from control at same timepoint)
ctrl_data <- trajectory_data %>%
  filter(Genotype == "Ctrl") %>%
  select(Module, Day, Expression) %>%
  rename(Ctrl_Expression = Expression)

trajectory_data <- trajectory_data %>%
  left_join(ctrl_data, by = c("Module", "Day")) %>%
  mutate(Divergence = Expression - Ctrl_Expression)

message("\U2713 Trajectory data calculated\n")

# ============================================================================ #
# 6. Calculate Per-Sample Module Scores (for individual data points)           #
# ============================================================================ #

message("\U0001F4CA Calculating per-sample module scores for data transparency...\n")

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

  # Add annotation data
  annot_with_ids <- annot
  annot_with_ids$Sample_ID <- rownames(annot)

  sample_info <- merge(sample_info, annot_with_ids[, c("Sample_ID", "genotype", "days")],
                       by.x = "Sample", by.y = "Sample_ID")

  # Rename and standardize
  colnames(sample_info)[colnames(sample_info) == "genotype"] <- "Genotype"
  colnames(sample_info)[colnames(sample_info) == "days"] <- "Timepoint"
  sample_info$Genotype <- gsub("Control", "Ctrl", sample_info$Genotype)
  sample_info$Genotype <- factor(sample_info$Genotype, levels = c("Ctrl", "G32A", "R403C"))
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

# Add jitter for plotting
sample_data_jittered <- sample_data %>%
  mutate(Day_jittered = Day + rnorm(n(), mean = 0, sd = 0.5))

message("\U2713 Sample-level data calculated\n")

# ============================================================================ #
# 7. Helper Functions for Plotting                                             #
# ============================================================================ #

# Calculate dynamic y-axis limits
calculate_y_limits <- function(trajectory_df, sample_df, expression_col = "Expression", padding = 0.1) {
  traj_range <- range(trajectory_df[[expression_col]], na.rm = TRUE)

  ctrl_d35_mean <- sample_df %>%
    filter(Genotype == "Ctrl", Day == 35) %>%
    pull(ModuleScore) %>%
    mean(na.rm = TRUE)

  sample_expr <- sample_df$ModuleScore - ctrl_d35_mean
  sample_range <- range(sample_expr, na.rm = TRUE)

  overall_min <- min(traj_range[1], sample_range[1], na.rm = TRUE)
  overall_max <- max(traj_range[2], sample_range[2], na.rm = TRUE)

  range_span <- overall_max - overall_min
  y_min <- overall_min - (range_span * padding)
  y_max <- overall_max + (range_span * padding)

  y_min <- floor(y_min * 10) / 10
  y_max <- ceiling(y_max * 10) / 10

  return(list(min = y_min, max = y_max))
}

# Create a single trajectory panel
create_trajectory_panel <- function(module_name, title, subtitle_prefix = "Expression trajectories") {
  mod_trajectory <- trajectory_data %>% filter(Module == module_name)
  mod_samples <- sample_data_jittered %>% filter(Module == module_name)

  if (nrow(mod_trajectory) == 0 || nrow(mod_samples) == 0) {
    message(sprintf("  WARNING: No data for module %s", module_name))
    return(NULL)
  }

  # Calculate y-axis limits
  ylim <- calculate_y_limits(mod_trajectory, mod_samples, "Expression", padding = 0.15)
  breaks <- seq(ylim$min, ylim$max, by = 0.2)

  # Calculate baseline for sample points
  ctrl_d35_mean <- mod_samples %>%
    filter(Genotype == "Ctrl", Day == 35) %>%
    pull(ModuleScore) %>%
    mean(na.rm = TRUE)

  p <- ggplot(mod_trajectory, aes(x = Day, color = Genotype, group = Genotype)) +
    # Layer 1: Individual sample points
    geom_point(data = mod_samples,
               aes(x = Day_jittered, y = ModuleScore - ctrl_d35_mean),
               alpha = 0.35, size = 2.5, shape = 16) +
    # Layer 2: Mean trajectory line
    geom_line(aes(y = Expression), linewidth = 1.5, alpha = 0.9) +
    # Layer 3: Mean trajectory points
    geom_point(aes(y = Expression), size = 4.5, alpha = 1) +
    # Reference line
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
    scale_color_manual(
      values = c("Ctrl" = "#999999", "G32A" = "#D55E00", "R403C" = "#009E73"),
      labels = c("Ctrl" = "Control", "G32A" = "G32A mutant", "R403C" = "R403C mutant")
    ) +
    scale_x_continuous(breaks = c(35, 65), labels = c("D35", "D65")) +
    scale_y_continuous(breaks = breaks,
                       labels = scales::number_format(accuracy = 0.1),
                       limits = c(ylim$min, ylim$max)) +
    labs(
      title = title,
      subtitle = sprintf("%s (n=%d genes; points: n=3 biological replicates/group)",
                         subtitle_prefix, unique(mod_trajectory$N_genes)),
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

  return(p)
}

# ============================================================================ #
# 8. Create Individual Panel Plots                                             #
# ============================================================================ #

message("\U0001F4CA Creating individual trajectory panels...\n")

# Panel A: Ribosome Biogenesis
p_biogenesis <- create_trajectory_panel(
  "Ribosome Biogenesis",
  "Ribosome Biogenesis: Compensatory Upregulation"
)
ggsave(file.path(out_dir, "Panel_A_Biogenesis_Trajectory.pdf"),
       p_biogenesis, width = 8, height = 6)
message("  \U2713 Panel A: Ribosome Biogenesis")

# Panel B: Cytoplasmic Translation
p_translation <- create_trajectory_panel(
  "Cytoplasmic Translation",
  "Cytoplasmic Translation: Functional Failure"
)
ggsave(file.path(out_dir, "Panel_B_Translation_Trajectory.pdf"),
       p_translation, width = 8, height = 6)
message("  \U2713 Panel B: Cytoplasmic Translation")

# Panel C: Synaptic Ribosomes
p_synaptic <- create_trajectory_panel(
  "Synaptic Ribosomes",
  "Synaptic Ribosomes: Critical Period Failure"
)
ggsave(file.path(out_dir, "Panel_C_Synaptic_Trajectory.pdf"),
       p_synaptic, width = 8, height = 6)
message("  \U2713 Panel C: Synaptic Ribosomes")

# Panel D: Mitochondrial Ribosome
p_mito_ribo <- create_trajectory_panel(
  "Mitochondrial Ribosome",
  "Mitochondrial Ribosome: Compensation Response"
)
ggsave(file.path(out_dir, "Panel_D_MitoRibosome_Trajectory.pdf"),
       p_mito_ribo, width = 8, height = 6)
message("  \U2713 Panel D: Mitochondrial Ribosome")

# Panel E: Mito Ribosome Assembly
p_mito_assembly <- create_trajectory_panel(
  "Mito Ribosome Assembly",
  "Mito Ribosome Assembly: Early Stress Recovery"
)
ggsave(file.path(out_dir, "Panel_E_MitoRibosomeAssembly_Trajectory.pdf"),
       p_mito_assembly, width = 8, height = 6)
message("  \U2713 Panel E: Mito Ribosome Assembly")

# Panel F: mtDNA Maintenance
p_mtdna <- create_trajectory_panel(
  "mtDNA Maintenance",
  "mtDNA Maintenance: Shared Compensation"
)
ggsave(file.path(out_dir, "Panel_F_mtDNA_Maintenance_Trajectory.pdf"),
       p_mtdna, width = 8, height = 6)
message("  \U2713 Panel F: mtDNA Maintenance")

# Panel G: OXPHOS
p_oxphos <- create_trajectory_panel(
  "OXPHOS",
  "OXPHOS: Energy Pathway Recovery"
)
ggsave(file.path(out_dir, "Panel_G_OXPHOS_Trajectory.pdf"),
       p_oxphos, width = 8, height = 6)
message("  \U2713 Panel G: OXPHOS")

message("\n\U2713 All individual panels complete\n")

# ============================================================================ #
# 9. Create Combined Figures                                                   #
# ============================================================================ #

message("\U0001F4CA Creating combined figures...\n")

# --- Original 3-panel (cytoplasmic/synaptic only) ---
p_combined_3 <- (p_biogenesis / p_translation / p_synaptic) +
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
       p_combined_3, width = 7, height = 14)
message("  \U2713 Combined 3-panel (original)")

# --- 7-panel vertical stack ---
p_combined_7_vertical <- (p_biogenesis / p_translation / p_synaptic /
                          p_mito_ribo / p_mito_assembly / p_mtdna / p_oxphos) +
  plot_annotation(
    title = "Critical Period Trajectories: Cytoplasmic & Mitochondrial Dynamics",
    subtitle = "DRP1 mutants show paradoxical response with mitochondrial compensation",
    theme = theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 13, color = "gray30")
    )
  ) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "Combined_Trajectories_7panel_vertical.pdf"),
       p_combined_7_vertical, width = 8, height = 32)
message("  \U2713 Combined 7-panel vertical")

# --- 2-column grid layout ---
# Left column: Cytoplasmic/Synaptic (3 panels)
col_left <- (p_biogenesis / p_translation / p_synaptic) +
  plot_layout(guides = "collect")

# Right column: MitoCarta (4 panels)
col_right <- (p_mito_ribo / p_mito_assembly / p_mtdna / p_oxphos) +
  plot_layout(guides = "collect")

p_combined_2col <- (col_left | col_right) +
  plot_annotation(
    title = "Critical Period Trajectories: Cytoplasmic & Mitochondrial Dynamics",
    subtitle = "Left: Cytoplasmic/Synaptic modules | Right: MitoCarta compensation pathways",
    theme = theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 12, color = "gray30")
    )
  ) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "Combined_Trajectories_7panel_2col.pdf"),
       p_combined_2col, width = 16, height = 18)
message("  \U2713 Combined 7-panel 2-column")

# --- MitoCarta only 4-panel ---
p_mitocarta_4 <- (p_mito_ribo / p_mito_assembly / p_mtdna / p_oxphos) +
  plot_annotation(
    title = "MitoCarta Compensation Pathways: D35 to D65 Dynamics",
    subtitle = "Mitochondrial pathways show compensatory upregulation during maturation",
    theme = theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 13, color = "gray30")
    )
  ) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "Combined_MitoCarta_4panel.pdf"),
       p_mitocarta_4, width = 8, height = 20)
message("  \U2713 Combined MitoCarta 4-panel")

message("\n\U2713 All combined figures complete\n")

# ============================================================================ #
# 10. Create Divergence Overlay Plots                                          #
# ============================================================================ #

message("\U0001F4CA Creating divergence overlay plots...\n")

# Color scheme for 7 modules
module_colors <- c(
  "Ribosome Biogenesis" = "#009E73",
  "Cytoplasmic Translation" = "#D55E00",
  "Synaptic Ribosomes" = "#0072B2",
  "Mitochondrial Ribosome" = "#CC79A7",
  "Mito Ribosome Assembly" = "#E69F00",
  "mtDNA Maintenance" = "#56B4E9",
  "OXPHOS" = "#F0E442"
)

# --- 7-module divergence (trajectories only) ---
div_data_7 <- trajectory_data %>%
  filter(Genotype != "Ctrl")

p_divergence_7 <- ggplot(div_data_7, aes(x = Day, y = Divergence, color = Module)) +
  geom_line(aes(linetype = Genotype, group = interaction(Module, Genotype)),
            linewidth = 1.5, alpha = 0.9) +
  geom_point(aes(shape = Genotype), size = 4, alpha = 1) +
  geom_hline(yintercept = 0, linetype = "solid", color = "gray50", linewidth = 0.8) +
  scale_color_manual(values = module_colors) +
  scale_shape_manual(
    values = c("G32A" = 16, "R403C" = 17),
    labels = c("G32A" = "G32A mutant", "R403C" = "R403C mutant")
  ) +
  scale_linetype_manual(
    values = c("G32A" = "solid", "R403C" = "dashed"),
    labels = c("G32A" = "G32A mutant", "R403C" = "R403C mutant")
  ) +
  scale_x_continuous(breaks = c(35, 65), labels = c("D35", "D65")) +
  labs(
    title = "Divergence from Control: All 7 Modules",
    subtitle = "Cytoplasmic/Synaptic vs MitoCarta pathways during critical period",
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
    color = guide_legend(order = 1, nrow = 2, override.aes = list(size = 3)),
    shape = guide_legend(order = 2),
    linetype = guide_legend(order = 2)
  )

ggsave(file.path(out_dir, "Panel_Divergence_All_7_Modules.pdf"),
       p_divergence_7, width = 10, height = 7)
message("  \U2713 Divergence: All 7 modules")

# --- MitoCarta only divergence ---
mito_modules_list <- c("Mitochondrial Ribosome", "Mito Ribosome Assembly",
                       "mtDNA Maintenance", "OXPHOS")

div_data_mito <- trajectory_data %>%
  filter(Genotype != "Ctrl", Module %in% mito_modules_list)

p_divergence_mito <- ggplot(div_data_mito, aes(x = Day, y = Divergence, color = Module)) +
  geom_line(aes(linetype = Genotype, group = interaction(Module, Genotype)),
            linewidth = 1.5, alpha = 0.9) +
  geom_point(aes(shape = Genotype), size = 4, alpha = 1) +
  geom_hline(yintercept = 0, linetype = "solid", color = "gray50", linewidth = 0.8) +
  scale_color_manual(values = module_colors[mito_modules_list]) +
  scale_shape_manual(
    values = c("G32A" = 16, "R403C" = 17),
    labels = c("G32A" = "G32A mutant", "R403C" = "R403C mutant")
  ) +
  scale_linetype_manual(
    values = c("G32A" = "solid", "R403C" = "dashed"),
    labels = c("G32A" = "G32A mutant", "R403C" = "R403C mutant")
  ) +
  scale_x_continuous(breaks = c(35, 65), labels = c("D35", "D65")) +
  labs(
    title = "Divergence from Control: MitoCarta Pathways",
    subtitle = "Mitochondrial compensation during critical period",
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
    color = guide_legend(order = 1, nrow = 1, override.aes = list(size = 3)),
    shape = guide_legend(order = 2),
    linetype = guide_legend(order = 2)
  )

ggsave(file.path(out_dir, "Panel_Divergence_MitoCarta_Only.pdf"),
       p_divergence_mito, width = 9, height = 6)
message("  \U2713 Divergence: MitoCarta only")

message("\n\U2713 All divergence plots complete\n")

# ============================================================================ #
# 11. Export Data                                                              #
# ============================================================================ #

message("\U0001F4DD Exporting data files...\n")

# Trajectory data (mean module scores from model)
write.csv(trajectory_data,
          file.path(out_dir, "trajectory_data_with_mitocarta.csv"),
          row.names = FALSE)
message(sprintf("  \U2713 trajectory_data_with_mitocarta.csv (%d rows)", nrow(trajectory_data)))

# Sample-level data (individual biological replicates)
write.csv(sample_data,
          file.path(out_dir, "sample_data_with_mitocarta.csv"),
          row.names = FALSE)
message(sprintf("  \U2713 sample_data_with_mitocarta.csv (%d rows)", nrow(sample_data)))

# Gene module summary
module_summary <- data.frame(
  Module = names(gene_modules),
  N_genes_total = sapply(gene_modules, length),
  N_genes_in_expression = sapply(gene_modules, function(g) sum(g %in% rownames(logCPM)))
)
write.csv(module_summary,
          file.path(out_dir, "module_gene_counts.csv"),
          row.names = FALSE)
message(sprintf("  \U2713 module_gene_counts.csv (%d modules)", nrow(module_summary)))

message("\n\U2713 Data export complete\n")

# ============================================================================ #
# 12. Summary                                                                  #
# ============================================================================ #

message("=" |> rep(78) |> paste(collapse = ""))
message("\U2705 Critical Period Trajectories with MitoCarta - COMPLETE")
message("=" |> rep(78) |> paste(collapse = ""))

message(sprintf("\n\U0001F4C1 Output directory: %s", out_dir))

message("\n\U0001F3AF KEY OUTPUTS:")
message("  Individual panels:")
message("    - Panel_A_Biogenesis_Trajectory.pdf")
message("    - Panel_B_Translation_Trajectory.pdf")
message("    - Panel_C_Synaptic_Trajectory.pdf")
message("    - Panel_D_MitoRibosome_Trajectory.pdf")
message("    - Panel_E_MitoRibosomeAssembly_Trajectory.pdf")
message("    - Panel_F_mtDNA_Maintenance_Trajectory.pdf")
message("    - Panel_G_OXPHOS_Trajectory.pdf")
message("  Combined figures:")
message("    - Combined_Trajectories_3panel.pdf (original)")
message("    - Combined_Trajectories_7panel_vertical.pdf")
message("    - Combined_Trajectories_7panel_2col.pdf")
message("    - Combined_MitoCarta_4panel.pdf")
message("  Divergence overlays:")
message("    - Panel_Divergence_All_7_Modules.pdf")
message("    - Panel_Divergence_MitoCarta_Only.pdf")
message("  Data:")
message("    - trajectory_data_with_mitocarta.csv")
message("    - sample_data_with_mitocarta.csv")
message("    - module_gene_counts.csv")

message("\n\U0001F4CA Module Summary:")
print(module_summary)

message("\n\U2728 Done!")
