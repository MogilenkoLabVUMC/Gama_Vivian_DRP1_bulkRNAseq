###############################################################################
##  Critical Period Trajectories: GSVA-Based Analysis                         ##
##  D35 -> D65 Maturation Dynamics with Proper Pathway Enrichment Scores      ##
###############################################################################
##  PURPOSE: Visualize how gene module enrichment changes during the critical ##
##           developmental period using GSVA (Gene Set Variation Analysis)    ##
##                                                                            ##
##  OUTPUTS: Two complementary views per module:                              ##
##           1. Trajectory: Control visible, shows normal development         ##
##           2. Divergence: Mutants only, shows deviation from control        ##
##                                                                            ##
##  IMPROVEMENTS over previous version:                                       ##
##           - GSVA scores (not simple mean) for proper pathway enrichment    ##
##           - Shared y-axis across panels for direct comparison              ##
##           - Square 6x6 inch panels for GUI grid flexibility                ##
##           - Clean output organization with archiving of old files          ##
###############################################################################

library(here)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(GSVA)

message("=" |> rep(78) |> paste(collapse = ""))
message("Critical Period Trajectories: GSVA-Based Analysis")
message("=" |> rep(78) |> paste(collapse = ""))

# ============================================================================ #
# 0. Configuration                                                             #
# ============================================================================ #

config <- list(
  # Paths
  checkpoint_dir = here("03_Results/02_Analysis/checkpoints"),
  out_root = here("03_Results/02_Analysis/Plots/Critical_period_trajectories"),

  # GSVA parameters
  kcdf = "Gaussian",      # For continuous log-CPM data
  minSize = 10,           # Minimum gene set size
  maxSize = Inf,          # No upper limit

  # Plot dimensions (square for grid flexibility)
  panel_width = 6,

  panel_height = 6,

  # Force recompute GSVA scores
  force_recompute_gsva = FALSE
)

# Create output directories
gsva_dir <- file.path(config$out_root, "gsva")
traj_dir <- file.path(gsva_dir, "trajectory")
div_dir <- file.path(gsva_dir, "divergence")
combined_dir <- file.path(gsva_dir, "combined")
data_dir <- file.path(gsva_dir, "data")

for (d in c(gsva_dir, traj_dir, div_dir, combined_dir, data_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ============================================================================ #
# 1. Archive Existing Outputs                                                  #
# ============================================================================ #

archive_old_outputs <- function(out_dir) {
  archive_date <- format(Sys.Date(), "%Y-%m-%d")
  archive_path <- file.path(out_dir, "archive", archive_date)

  # Find items to archive (everything except archive folder and README)
  existing_items <- list.files(out_dir, full.names = TRUE, include.dirs = TRUE)
  items_to_archive <- existing_items[!grepl("archive$|README\\.md$|gsva$", existing_items)]

  if (length(items_to_archive) == 0) {
    message("  No files to archive")
    return(invisible(NULL))
  }

  dir.create(archive_path, recursive = TRUE, showWarnings = FALSE)

  # Move items
  for (item in items_to_archive) {
    dest <- file.path(archive_path, basename(item))
    file.rename(item, dest)
  }

  # Write manifest
  manifest <- data.frame(
    original_name = basename(items_to_archive),
    archive_date = archive_date,
    reason = "Replaced by GSVA-based implementation"
  )
  write.csv(manifest, file.path(archive_path, "ARCHIVE_MANIFEST.csv"), row.names = FALSE)

  message(sprintf("  Archived %d items to: %s", length(items_to_archive), archive_path))
}

message("\n\U0001F4E6 Archiving existing outputs...")
archive_old_outputs(config$out_root)

# ============================================================================ #
# 2. Load Checkpoint Data                                                      #
# ============================================================================ #

message("\n\U0001F4C2 Loading checkpoints...")

all_gsea_results <- readRDS(file.path(config$checkpoint_dir, "all_gsea_results.rds"))
syngo_gsea_results <- readRDS(file.path(config$checkpoint_dir, "syngo_gsea_results.rds"))
fit <- readRDS(file.path(config$checkpoint_dir, "fit_object.rds"))
qc_vars <- readRDS(file.path(config$checkpoint_dir, "qc_variables.rds"))

# Unpack expression data and sample metadata
logCPM <- qc_vars$logCPM
annot <- qc_vars$annot

message("\U2713 Checkpoints loaded")

# ============================================================================ #
# 3. Load MitoCarta Gene Sets                                                  #
# ============================================================================ #

message("\n\U0001F9EC Loading MitoCarta pathways...")

# Source MitoCarta parser
source(here("01_Scripts/R_scripts/parse_mitocarta_gmx.R"))

# Parse MitoCarta GMX file
mitocarta_path <- here("00_Data/MitoCarta_3.0/MitoPathways3.0.gmx")
mitocarta_lists <- mitocarta_gmt(mitocarta_path)
T2G <- mitocarta_lists$T2G

message("\U2713 MitoCarta loaded")

# ============================================================================ #
# 4. Extract Gene Lists for All 7 Modules                                      #
# ============================================================================ #

message("\n\U0001F4CA Extracting gene sets for all modules...")

# --- Existing modules (from GSEA enrichment) ---

# Get synaptic ribosome genes from verification reports
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

# Helper function to extract genes for a pathway
extract_mitocarta_genes <- function(T2G, pathway_pattern) {
  matching_rows <- grepl(pathway_pattern, T2G$gs_name, ignore.case = TRUE)
  genes <- unique(T2G$gene_symbol[matching_rows])
  return(genes)
}

# Extract the 4 core MitoCarta pathways
mito_ribosome_genes <- extract_mitocarta_genes(T2G, "Mitochondrial_ribosome$")
mito_ribo_assembly_genes <- extract_mitocarta_genes(T2G, "Mitochondrial_ribosome_assembly")
mtdna_maintenance_genes <- extract_mitocarta_genes(T2G, "mtDNA_maintenance")
oxphos_genes <- extract_mitocarta_genes(T2G, "^OXPHOS$")

message(sprintf("  Mitochondrial ribosome: %d genes", length(mito_ribosome_genes)))
message(sprintf("  Mito ribosome assembly: %d genes", length(mito_ribo_assembly_genes)))
message(sprintf("  mtDNA maintenance: %d genes", length(mtdna_maintenance_genes)))
message(sprintf("  OXPHOS: %d genes", length(oxphos_genes)))

# ============================================================================ #
# 5. Define All Gene Modules                                                   #
# ============================================================================ #

# Module definitions with display names and internal keys
module_info <- list(
  list(key = "Ribosome_Biogenesis", display = "Ribosome Biogenesis",
       genes = ribosome_biogenesis_genes, panel = "A"),
  list(key = "Cytoplasmic_Translation", display = "Cytoplasmic Translation",
       genes = translation_genes, panel = "B"),
  list(key = "Synaptic_Ribosomes", display = "Synaptic Ribosomes",
       genes = synaptic_ribosome_genes, panel = "C"),
  list(key = "Mitochondrial_Ribosome", display = "Mitochondrial Ribosome",
       genes = mito_ribosome_genes, panel = "D"),
  list(key = "Mito_Ribosome_Assembly", display = "Mito Ribosome Assembly",
       genes = mito_ribo_assembly_genes, panel = "E"),
  list(key = "mtDNA_Maintenance", display = "mtDNA Maintenance",
       genes = mtdna_maintenance_genes, panel = "F"),
  list(key = "OXPHOS", display = "OXPHOS",
       genes = oxphos_genes, panel = "G")
)

# Create named list for GSVA
gene_modules <- setNames(
  lapply(module_info, function(x) x$genes),
  sapply(module_info, function(x) x$key)
)

message("\n\U2713 All 7 gene modules defined")

# ============================================================================ #
# 6. Compute GSVA Scores                                                       #
# ============================================================================ #

message("\n\U0001F52C Computing GSVA scores...")

gsva_checkpoint_file <- file.path(config$checkpoint_dir, "gsva_module_scores.rds")

if (!config$force_recompute_gsva && file.exists(gsva_checkpoint_file)) {
  message("  Loading cached GSVA scores...")
  gsva_checkpoint <- readRDS(gsva_checkpoint_file)
  gsva_scores <- gsva_checkpoint$scores
  gene_modules_filtered <- gsva_checkpoint$gene_modules_filtered
} else {
  message("  Computing GSVA scores (this may take a moment)...")

  # Filter gene modules to genes present in expression matrix
  gene_modules_filtered <- lapply(gene_modules, function(genes) {
    genes[genes %in% rownames(logCPM)]
  })

  # Remove modules with too few genes
  module_sizes <- sapply(gene_modules_filtered, length)
  valid_modules <- module_sizes >= config$minSize

  if (!all(valid_modules)) {
    warning(sprintf("Removing modules with fewer than %d genes: %s",
                    config$minSize,
                    paste(names(gene_modules_filtered)[!valid_modules], collapse = ", ")))
    gene_modules_filtered <- gene_modules_filtered[valid_modules]
  }

  # Create GSVA parameter object
  gsva_param <- gsvaParam(
    exprData = as.matrix(logCPM),
    geneSets = gene_modules_filtered,
    kcdf = config$kcdf,
    minSize = config$minSize,
    maxSize = config$maxSize
  )

  # Run GSVA
  gsva_scores <- gsva(gsva_param, verbose = TRUE)

  # Save checkpoint
  gsva_checkpoint <- list(
    scores = gsva_scores,
    gene_modules = gene_modules,
    gene_modules_filtered = gene_modules_filtered,
    parameters = list(
      kcdf = config$kcdf,
      minSize = config$minSize,
      timestamp = Sys.time()
    )
  )
  saveRDS(gsva_checkpoint, gsva_checkpoint_file)
  message(sprintf("  \U2713 GSVA scores cached to: %s", basename(gsva_checkpoint_file)))
}

message(sprintf("\U2713 GSVA scores computed: %d modules x %d samples",
                nrow(gsva_scores), ncol(gsva_scores)))

# ============================================================================ #
# 7. Transform GSVA Scores to Long Format                                      #
# ============================================================================ #

message("\n\U0001F4CA Transforming GSVA scores...")

# Convert to long format
gsva_long <- gsva_scores %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Module") %>%
  tidyr::pivot_longer(
    cols = -Module,
    names_to = "Sample",
    values_to = "GSVA_Score"
  )

# Add sample metadata
annot_df <- annot %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Sample")

gsva_long <- gsva_long %>%
  left_join(annot_df, by = "Sample") %>%
  mutate(
    Genotype = gsub("Control", "Ctrl", genotype),
    Genotype = factor(Genotype, levels = c("Ctrl", "G32A", "R403C")),
    Day = ifelse(days == "D35", 35, 65)
  )

# Calculate group means
gsva_means <- gsva_long %>%
  group_by(Module, Genotype, Day) %>%
  summarize(
    Mean_GSVA = mean(GSVA_Score, na.rm = TRUE),
    SD_GSVA = sd(GSVA_Score, na.rm = TRUE),
    N = n(),
    .groups = "drop"
  )

message("\U2713 Data transformation complete")

# ============================================================================ #
# 8. Create Trajectory View (Relative to Ctrl D35)                             #
# ============================================================================ #

message("\n\U0001F4CA Creating trajectory view (Control visible)...")

# Calculate Ctrl D35 baseline for each module
ctrl_d35_baseline <- gsva_means %>%
  filter(Genotype == "Ctrl", Day == 35) %>%
  select(Module, Baseline = Mean_GSVA)

# Calculate expression relative to Ctrl D35
trajectory_means <- gsva_means %>%
  left_join(ctrl_d35_baseline, by = "Module") %>%
  mutate(Expression = Mean_GSVA - Baseline)

# For sample-level points
sample_trajectory <- gsva_long %>%
  left_join(ctrl_d35_baseline, by = "Module") %>%
  mutate(Expression = GSVA_Score - Baseline) %>%
  mutate(Day_jittered = Day + rnorm(n(), mean = 0, sd = 0.5))

message("\U2713 Trajectory view created")

# ============================================================================ #
# 9. Create Divergence View (Mutants vs Control at Each Timepoint)             #
# ============================================================================ #

message("\n\U0001F4CA Creating divergence view (Mutants only)...")

# Calculate control mean at each timepoint
ctrl_by_day <- gsva_means %>%
  filter(Genotype == "Ctrl") %>%
  select(Module, Day, Ctrl_Mean = Mean_GSVA)

# Calculate divergence for means
divergence_means <- gsva_means %>%
  filter(Genotype != "Ctrl") %>%
  left_join(ctrl_by_day, by = c("Module", "Day")) %>%
  mutate(Divergence = Mean_GSVA - Ctrl_Mean)

# For sample-level points
sample_divergence <- gsva_long %>%
  filter(Genotype != "Ctrl") %>%
  left_join(ctrl_by_day, by = c("Module", "Day")) %>%
  mutate(Divergence = GSVA_Score - Ctrl_Mean) %>%
  mutate(Day_jittered = Day + rnorm(n(), mean = 0, sd = 0.5))

message("\U2713 Divergence view created")

# ============================================================================ #
# 10. Calculate Global Y-Axis Ranges (Shared Across Panels)                    #
# ============================================================================ #

message("\n\U0001F4CF Calculating shared y-axis ranges...")

# For trajectory view
traj_global_range <- range(
  c(trajectory_means$Expression, sample_trajectory$Expression),
  na.rm = TRUE
)
traj_y_limits <- c(
  floor(traj_global_range[1] * 10) / 10 - 0.1,
  ceiling(traj_global_range[2] * 10) / 10 + 0.1
)
traj_y_breaks <- seq(traj_y_limits[1], traj_y_limits[2], by = 0.2)

# For divergence view
div_global_range <- range(
  c(divergence_means$Divergence, sample_divergence$Divergence),
  na.rm = TRUE
)
div_y_limits <- c(
  floor(div_global_range[1] * 10) / 10 - 0.1,
  ceiling(div_global_range[2] * 10) / 10 + 0.1
)
div_y_breaks <- seq(div_y_limits[1], div_y_limits[2], by = 0.2)

message(sprintf("  Trajectory y-range: [%.2f, %.2f]", traj_y_limits[1], traj_y_limits[2]))
message(sprintf("  Divergence y-range: [%.2f, %.2f]", div_y_limits[1], div_y_limits[2]))

# ============================================================================ #
# 11. Define Plotting Functions                                                #
# ============================================================================ #

# Color schemes
genotype_colors <- c("Ctrl" = "#999999", "G32A" = "#D55E00", "R403C" = "#009E73")
genotype_labels <- c("Ctrl" = "Control", "G32A" = "G32A mutant", "R403C" = "R403C mutant")

# Panel plotting function
create_gsva_panel <- function(module_key, display_name,
                               means_df, sample_df,
                               view_type = "trajectory",
                               y_limits, y_breaks,
                               n_genes) {

  # Filter to this module
  mod_means <- means_df %>% filter(Module == module_key)
  mod_samples <- sample_df %>% filter(Module == module_key)

  if (nrow(mod_means) == 0) {
    warning(sprintf("No data for module: %s", module_key))
    return(NULL)
  }

  # Set y variable and label based on view type
  if (view_type == "trajectory") {
    y_var <- "Expression"
    y_label <- "GSVA Score (relative to Ctrl D35)"
  } else {
    y_var <- "Divergence"
    y_label <- "Divergence from Control (GSVA Score)"
  }

  # Build plot
  p <- ggplot(mod_means, aes(x = Day, color = Genotype, group = Genotype)) +
    # Layer 1: Individual sample points (semi-transparent)
    geom_point(data = mod_samples,
               aes(x = Day_jittered, y = .data[[y_var]]),
               alpha = 0.35, size = 2.5, shape = 16) +
    # Layer 2: Mean trajectory line
    geom_line(aes(y = .data[[y_var]]), linewidth = 1.5, alpha = 0.9) +
    # Layer 3: Mean trajectory endpoints
    geom_point(aes(y = .data[[y_var]]), size = 4.5, alpha = 1) +
    # Reference line at zero
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
    # Scales
    scale_color_manual(values = genotype_colors, labels = genotype_labels) +
    scale_x_continuous(breaks = c(35, 65), labels = c("D35", "D65")) +
    scale_y_continuous(limits = y_limits, breaks = y_breaks,
                       labels = scales::number_format(accuracy = 0.1)) +
    # Labels
    labs(
      title = display_name,
      subtitle = sprintf("GSVA enrichment (n=%d genes; points: n=3 replicates/group)", n_genes),
      x = "Developmental Timepoint",
      y = y_label,
      color = "Genotype"
    ) +
    # Theme
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 9, color = "gray40"),
      axis.title = element_text(face = "bold", size = 11),
      axis.text = element_text(face = "bold", size = 10),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      aspect.ratio = 1  # Square panels for grid flexibility
    )

  return(p)
}

# ============================================================================ #
# 12. Generate Individual Trajectory Panels                                    #
# ============================================================================ #

message("\n\U0001F4CA Generating individual trajectory panels...")

trajectory_panels <- list()

for (mod in module_info) {
  n_genes <- length(gene_modules_filtered[[mod$key]])

  p <- create_gsva_panel(
    module_key = mod$key,
    display_name = mod$display,
    means_df = trajectory_means,
    sample_df = sample_trajectory,
    view_type = "trajectory",
    y_limits = traj_y_limits,
    y_breaks = traj_y_breaks,
    n_genes = n_genes
  )

  if (!is.null(p)) {
    trajectory_panels[[mod$key]] <- p

    # Save individual panel
    filename <- sprintf("Panel_%s_%s.pdf", mod$panel, mod$key)
    ggsave(file.path(traj_dir, filename), p,
           width = config$panel_width, height = config$panel_height)
    message(sprintf("  \U2713 %s", filename))
  }
}

message("\U2713 All trajectory panels complete")

# ============================================================================ #
# 13. Generate Individual Divergence Panels                                    #
# ============================================================================ #

message("\n\U0001F4CA Generating individual divergence panels...")

divergence_panels <- list()

for (mod in module_info) {
  n_genes <- length(gene_modules_filtered[[mod$key]])

  p <- create_gsva_panel(
    module_key = mod$key,
    display_name = mod$display,
    means_df = divergence_means,
    sample_df = sample_divergence,
    view_type = "divergence",
    y_limits = div_y_limits,
    y_breaks = div_y_breaks,
    n_genes = n_genes
  )

  if (!is.null(p)) {
    divergence_panels[[mod$key]] <- p

    # Save individual panel
    filename <- sprintf("Panel_%s_%s.pdf", mod$panel, mod$key)
    ggsave(file.path(div_dir, filename), p,
           width = config$panel_width, height = config$panel_height)
    message(sprintf("  \U2713 %s", filename))
  }
}

message("\U2713 All divergence panels complete")

# ============================================================================ #
# 14. Generate Combined Figures                                                #
# ============================================================================ #

message("\n\U0001F4CA Generating combined figures...")

# --- 7-panel grid (2 columns) ---
if (length(trajectory_panels) >= 7) {
  p_grid <- wrap_plots(trajectory_panels, ncol = 2) +
    plot_annotation(
      title = "Critical Period Trajectories: GSVA Enrichment Scores",
      subtitle = "Shared y-axis enables direct comparison across all 7 modules",
      theme = theme(
        plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(size = 12, color = "gray40")
      )
    ) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

  ggsave(file.path(combined_dir, "Trajectory_7panel_grid.pdf"),
         p_grid, width = 13, height = 21)
  message("  \U2713 Trajectory_7panel_grid.pdf")
}

# --- Divergence overlay (all modules) ---
module_colors <- c(
  "Ribosome_Biogenesis" = "#009E73",
  "Cytoplasmic_Translation" = "#D55E00",
  "Synaptic_Ribosomes" = "#0072B2",
  "Mitochondrial_Ribosome" = "#CC79A7",
  "Mito_Ribosome_Assembly" = "#E69F00",
  "mtDNA_Maintenance" = "#56B4E9",
  "OXPHOS" = "#F0E442"
)

p_div_overlay <- ggplot(divergence_means,
                         aes(x = Day, y = Divergence, color = Module)) +
  geom_line(aes(linetype = Genotype, group = interaction(Module, Genotype)),
            linewidth = 1.5, alpha = 0.9) +
  geom_point(aes(shape = Genotype), size = 4, alpha = 1) +
  geom_hline(yintercept = 0, linetype = "solid", color = "gray50", linewidth = 0.8) +
  scale_color_manual(values = module_colors,
                     labels = function(x) gsub("_", " ", x)) +
  scale_shape_manual(values = c("G32A" = 16, "R403C" = 17),
                     labels = c("G32A" = "G32A mutant", "R403C" = "R403C mutant")) +
  scale_linetype_manual(values = c("G32A" = "solid", "R403C" = "dashed"),
                        labels = c("G32A" = "G32A mutant", "R403C" = "R403C mutant")) +
  scale_x_continuous(breaks = c(35, 65), labels = c("D35", "D65")) +
  labs(
    title = "Divergence from Control: All 7 Modules (GSVA)",
    subtitle = "Cytoplasmic/Synaptic vs MitoCarta pathways during critical period",
    x = "Developmental Timepoint",
    y = "Divergence from Control (GSVA Score)",
    color = "Module", shape = "Genotype", linetype = "Genotype"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "gray40"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold"),
    legend.position = "bottom",
    legend.box = "vertical",
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  ) +
  guides(
    color = guide_legend(order = 1, nrow = 2, override.aes = list(size = 3)),
    shape = guide_legend(order = 2),
    linetype = guide_legend(order = 2)
  )

ggsave(file.path(combined_dir, "Divergence_overlay_all.pdf"),
       p_div_overlay, width = 10, height = 7)
message("  \U2713 Divergence_overlay_all.pdf")

# --- Divergence overlay (MitoCarta only) ---
mito_modules <- c("Mitochondrial_Ribosome", "Mito_Ribosome_Assembly",
                  "mtDNA_Maintenance", "OXPHOS")

div_mito <- divergence_means %>% filter(Module %in% mito_modules)

p_div_mito <- ggplot(div_mito, aes(x = Day, y = Divergence, color = Module)) +
  geom_line(aes(linetype = Genotype, group = interaction(Module, Genotype)),
            linewidth = 1.5, alpha = 0.9) +
  geom_point(aes(shape = Genotype), size = 4, alpha = 1) +
  geom_hline(yintercept = 0, linetype = "solid", color = "gray50", linewidth = 0.8) +
  scale_color_manual(values = module_colors[mito_modules],
                     labels = function(x) gsub("_", " ", x)) +
  scale_shape_manual(values = c("G32A" = 16, "R403C" = 17),
                     labels = c("G32A" = "G32A mutant", "R403C" = "R403C mutant")) +
  scale_linetype_manual(values = c("G32A" = "solid", "R403C" = "dashed"),
                        labels = c("G32A" = "G32A mutant", "R403C" = "R403C mutant")) +
  scale_x_continuous(breaks = c(35, 65), labels = c("D35", "D65")) +
  labs(
    title = "Divergence from Control: MitoCarta Pathways (GSVA)",
    subtitle = "Mitochondrial compensation during critical period",
    x = "Developmental Timepoint",
    y = "Divergence from Control (GSVA Score)",
    color = "Module", shape = "Genotype", linetype = "Genotype"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "gray40"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold"),
    legend.position = "bottom",
    legend.box = "vertical",
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  ) +
  guides(
    color = guide_legend(order = 1, nrow = 1, override.aes = list(size = 3)),
    shape = guide_legend(order = 2),
    linetype = guide_legend(order = 2)
  )

ggsave(file.path(combined_dir, "Divergence_overlay_mitocarta.pdf"),
       p_div_mito, width = 9, height = 6)
message("  \U2713 Divergence_overlay_mitocarta.pdf")

message("\n\U2713 All combined figures complete")

# ============================================================================ #
# 15. Export Data Files                                                        #
# ============================================================================ #

message("\n\U0001F4DD Exporting data files...")

# GSVA scores matrix
gsva_scores_df <- gsva_scores %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Module")
write.csv(gsva_scores_df, file.path(data_dir, "gsva_scores_matrix.csv"), row.names = FALSE)
message(sprintf("  \U2713 gsva_scores_matrix.csv (%d modules x %d samples)",
                nrow(gsva_scores), ncol(gsva_scores)))

# Trajectory data (means)
write.csv(trajectory_means, file.path(data_dir, "trajectory_data.csv"), row.names = FALSE)
message(sprintf("  \U2713 trajectory_data.csv (%d rows)", nrow(trajectory_means)))

# Divergence data (means)
write.csv(divergence_means, file.path(data_dir, "divergence_data.csv"), row.names = FALSE)
message(sprintf("  \U2713 divergence_data.csv (%d rows)", nrow(divergence_means)))

# Module gene counts
module_summary <- data.frame(
  Module = names(gene_modules),
  N_genes_defined = sapply(gene_modules, length),
  N_genes_in_expression = sapply(gene_modules, function(g) sum(g %in% rownames(logCPM))),
  N_genes_in_gsva = sapply(names(gene_modules), function(m) {
    if (m %in% names(gene_modules_filtered)) length(gene_modules_filtered[[m]]) else 0
  })
)
write.csv(module_summary, file.path(data_dir, "module_gene_counts.csv"), row.names = FALSE)
message(sprintf("  \U2713 module_gene_counts.csv (%d modules)", nrow(module_summary)))

message("\n\U2713 Data export complete")

# ============================================================================ #
# 16. Summary                                                                  #
# ============================================================================ #

message("\n" |> paste0(rep("=", 78) |> paste(collapse = "")))
message("\U2705 Critical Period Trajectories (GSVA) - COMPLETE")
message(rep("=", 78) |> paste(collapse = ""))

message(sprintf("\n\U0001F4C1 Output directory: %s", gsva_dir))

message("\n\U0001F3AF KEY OUTPUTS:")
message("  Trajectory panels (gsva/trajectory/):")
for (mod in module_info) {
  message(sprintf("    - Panel_%s_%s.pdf", mod$panel, mod$key))
}
message("  Divergence panels (gsva/divergence/):")
for (mod in module_info) {
  message(sprintf("    - Panel_%s_%s.pdf", mod$panel, mod$key))
}
message("  Combined figures (gsva/combined/):")
message("    - Trajectory_7panel_grid.pdf")
message("    - Divergence_overlay_all.pdf")
message("    - Divergence_overlay_mitocarta.pdf")
message("  Data files (gsva/data/):")
message("    - gsva_scores_matrix.csv")
message("    - trajectory_data.csv")
message("    - divergence_data.csv")
message("    - module_gene_counts.csv")

message("\n\U0001F4CA Module Summary:")
print(module_summary)

message("\n\U0001F4DD GSVA Parameters:")
message(sprintf("  - kcdf: %s (for continuous log-CPM data)", config$kcdf))
message(sprintf("  - minSize: %d genes", config$minSize))
message(sprintf("  - Panel dimensions: %dx%d inches (square)", config$panel_width, config$panel_height))
message(sprintf("  - Y-axis: Shared across all panels for direct comparison"))

message("\n\U2728 Done!")
