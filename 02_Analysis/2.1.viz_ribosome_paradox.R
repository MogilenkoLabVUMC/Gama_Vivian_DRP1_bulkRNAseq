###############################################################################
##  Ribosome Paradox: Three Pools, Three Fates                              ##
###############################################################################

library(here)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(tibble)

message("📂 Loading checkpoints...")
checkpoint_dir <- here("03_Results/02_Analysis/checkpoints")
all_gsea_results <- readRDS(file.path(checkpoint_dir, "all_gsea_results.rds"))
syngo_gsea_results <- readRDS(file.path(checkpoint_dir, "syngo_gsea_results.rds"))
mitocarta_gsea_results <- readRDS(file.path(checkpoint_dir, "mitocarta_gsea_results.rds"))

# Output directory
out_dir <- here("03_Results/02_Analysis/Plots/Ribosome_paradox")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("✓ Checkpoints loaded\n")

###############################################################################
##  Extract Three Ribosome Pool Pathways                                    ##
###############################################################################

message("📊 Extracting ribosome pathways from three pools...\n")

# Contrasts: Early (D35), TrajDev (Maturation), and Late (D65)
contrasts_early <- c("G32A_vs_Ctrl_D35", "R403C_vs_Ctrl_D35")
contrasts_trajdev <- c("Maturation_G32A_specific", "Maturation_R403C_specific")
contrasts_late <- c("G32A_vs_Ctrl_D65", "R403C_vs_Ctrl_D65")
all_contrasts <- c(contrasts_early, contrasts_trajdev, contrasts_late)

ribosome_data <- data.frame()

## POOL 1: Cytoplasmic Ribosomal Biogenesis (GO pathways)
message("  Pool 1: Cytoplasmic Ribosomal Biogenesis (GOBP)...")

biogenesis_keywords <- c(
  "ribosomal.*biogenesis",
  "ribosomal.*large.*subunit.*biogenesis",
  "ribosomal.*small.*subunit.*biogenesis",
  "ribosome.*biogenesis"
)

for (contrast in all_contrasts) {
  gobp_res <- all_gsea_results[[contrast]][["gobp"]]

  if (!is.null(gobp_res)) {
    res_df <- if (class(gobp_res) == "gseaResult") {
      gobp_res@result
    } else if ("result" %in% names(gobp_res)) {
      gobp_res$result
    } else {
      NULL
    }

    if (!is.null(res_df) && nrow(res_df) > 0) {
      # Find biogenesis pathways
      biogen_idx <- grepl(paste(biogenesis_keywords, collapse = "|"),
                         res_df$Description, ignore.case = TRUE)

      if (any(biogen_idx)) {
        biogen_pathways <- res_df[biogen_idx, ]

        # Include ALL matching pathways (not just significant)
        # This ensures complete trajectories across development
        # Significance will be indicated visually in plots

        if (nrow(biogen_pathways) > 0) {
          biogen_pathways$Pool <- "1. Cytoplasmic Ribosome Biogenesis"
          biogen_pathways$Contrast <- contrast
          biogen_pathways$Significant <- biogen_pathways$p.adjust < 0.05
          ribosome_data <- rbind(ribosome_data, biogen_pathways)
        }
      }
    }
  }
}

message(sprintf("    Found %d biogenesis pathways\n",
                nrow(ribosome_data[ribosome_data$Pool == "1. Cytoplasmic Ribosome Biogenesis", ])))

## POOL 2: Synaptic Ribosomes (SynGO)
message("  Pool 2: Synaptic Ribosomes (SynGO)...")

syngo_ribosome_keywords <- c(
  "ribosom",
  "translation"
)

for (contrast in all_contrasts) {
  syngo_res <- syngo_gsea_results[[contrast]]

  if (!is.null(syngo_res)) {
    syngo_df <- syngo_res@result

    # Find ribosome/translation pathways
    syngo_ribo_idx <- grepl(paste(syngo_ribosome_keywords, collapse = "|"),
                           syngo_df$Description, ignore.case = TRUE)

    if (any(syngo_ribo_idx)) {
      syngo_ribo <- syngo_df[syngo_ribo_idx, ]

      # Include ALL matching pathways (not just significant)
      # This ensures complete trajectories across development
      # Significance will be indicated visually in plots

      if (nrow(syngo_ribo) > 0) {
        syngo_ribo$Pool <- "2. Synaptic Ribosomes (SynGO)"
        syngo_ribo$Contrast <- contrast
        syngo_ribo$Significant <- syngo_ribo$p.adjust < 0.05
        ribosome_data <- rbind(ribosome_data, syngo_ribo)
      }
    }
  }
}

message(sprintf("    Found %d synaptic ribosome pathways\n",
                nrow(ribosome_data[ribosome_data$Pool == "2. Synaptic Ribosomes (SynGO)", ])))

## POOL 3: Mitochondrial Ribosomes (MitoCarta)
message("  Pool 3: Mitochondrial Ribosomes (MitoCarta)...")

mitocarta_ribosome_keywords <- c(
  "ribosom",
  "translation",
  "central.*dogma"
)

for (contrast in all_contrasts) {
  mito_res <- mitocarta_gsea_results[[contrast]]

  if (!is.null(mito_res)) {
    mito_df <- mito_res@result

    # Find ribosome/translation pathways
    mito_ribo_idx <- grepl(paste(mitocarta_ribosome_keywords, collapse = "|"),
                          mito_df$Description, ignore.case = TRUE)

    if (any(mito_ribo_idx)) {
      mito_ribo <- mito_df[mito_ribo_idx, ]

      # Include ALL matching pathways (not just significant)
      # This ensures complete trajectories across development
      # Significance will be indicated visually in plots

      if (nrow(mito_ribo) > 0) {
        mito_ribo$Pool <- "3. Mitochondrial Ribosomes (MitoCarta)"
        mito_ribo$Contrast <- contrast
        mito_ribo$Significant <- mito_ribo$p.adjust < 0.05
        ribosome_data <- rbind(ribosome_data, mito_ribo)
      }
    }
  }
}

message(sprintf("    Found %d mitochondrial ribosome pathways\n",
                nrow(ribosome_data[ribosome_data$Pool == "3. Mitochondrial Ribosomes (MitoCarta)", ])))

message(sprintf("\n  Total ribosome pathways extracted: %d\n", nrow(ribosome_data)))

###############################################################################
##  Classify Contrasts by Timepoint                                         ##
###############################################################################

ribosome_data$Timepoint <- case_when(
  ribosome_data$Contrast %in% contrasts_early ~ "Early (D35)",
  ribosome_data$Contrast %in% contrasts_trajdev ~ "TrajDev (Maturation)",
  ribosome_data$Contrast %in% contrasts_late ~ "Late (D65)",
  TRUE ~ NA_character_
)

# Factor with correct order for trajectory
ribosome_data$Timepoint <- factor(
  ribosome_data$Timepoint,
  levels = c("Early (D35)", "TrajDev (Maturation)", "Late (D65)")
)

ribosome_data$Mutation <- ifelse(
  grepl("G32A", ribosome_data$Contrast),
  "G32A",
  "R403C"
)

# Create combined factor for plotting
ribosome_data$Condition <- paste(ribosome_data$Mutation, ribosome_data$Timepoint, sep = " - ")

ribosome_data$Condition <- factor(
  ribosome_data$Condition,
  levels = c(
    "G32A - Early (D35)", "G32A - TrajDev (Maturation)", "G32A - Late (D65)",
    "R403C - Early (D35)", "R403C - TrajDev (Maturation)", "R403C - Late (D65)"
  )
)

# Factor pool levels
ribosome_data$Pool <- factor(
  ribosome_data$Pool,
  levels = c("1. Cytoplasmic Ribosome Biogenesis",
             "2. Synaptic Ribosomes (SynGO)",
             "3. Mitochondrial Ribosomes (MitoCarta)")
)

# Simplify pathway names for display
ribosome_data$Pathway_Short <- gsub("^GOBP_|^GOCC_|^GOMF_", "",
                                    ribosome_data$Description)
ribosome_data$Pathway_Short <- gsub("_", " ", ribosome_data$Pathway_Short)
ribosome_data$Pathway_Short <- tools::toTitleCase(tolower(ribosome_data$Pathway_Short))

# Truncate long names
ribosome_data$Pathway_Short <- ifelse(
  nchar(ribosome_data$Pathway_Short) > 50,
  paste0(substr(ribosome_data$Pathway_Short, 1, 47), "..."),
  ribosome_data$Pathway_Short
)

###############################################################################
##  FIGURE 1: Three-Pool Dotplot (Main Story)                               ##
###############################################################################

message("📊 Creating three-pool dotplot (main story)...\n")

# Order pathways by NES within each pool
pathway_order <- ribosome_data %>%
  group_by(Pool, Pathway_Short) %>%
  summarize(Mean_NES = mean(NES, na.rm = TRUE), .groups = "drop") %>%
  arrange(Pool, Mean_NES) %>%
  pull(Pathway_Short)

ribosome_data$Pathway_Short <- factor(ribosome_data$Pathway_Short,
                                      levels = pathway_order)

# Create dotplot with facet by Mutation and Pool
# Separate significant and non-significant for distinct visualization
ribosome_sig <- ribosome_data[ribosome_data$Significant, ]
ribosome_nonsig <- ribosome_data[!ribosome_data$Significant, ]

p_paradox <- ggplot(ribosome_data,
                    aes(x = Timepoint, y = Pathway_Short,
                        color = NES, size = -log10(p.adjust))) +
  # Non-significant points: faded, smaller
  geom_point(
    data = ribosome_nonsig,
    alpha = 0.3, size = 1.5, color = "gray70"
  ) +
  # Significant points: full color and size
  geom_point(
    data = ribosome_sig,
    alpha = 0.85
  ) +

  # Highlight highly significant (p.adjust < 0.01)
  geom_point(
    data = subset(ribosome_data, p.adjust < 0.01 & Significant),
    shape = 21, color = "black", fill = NA, stroke = 1.5, alpha = 1
  ) +

  facet_grid(Pool ~ Mutation, scales = "free_y", space = "free_y") +

  scale_color_gradient2(
    low = "#2166AC", mid = "#F7F7F7", high = "#B35806",
    midpoint = 0,
    name = "NES",
    limits = c(-3.5, 3.5),
    oob = scales::squish
  ) +

  scale_size_continuous(
    range = c(2, 8),
    name = "-log10(FDR)",
    breaks = c(2, 5, 10, 15)
  ) +

  labs(
    title = "Ribosome Paradox: Temporal Trajectory Across Three Pools",
    subtitle = "Complete pathways (faded = non-sig FDR ≥ 0.05) | Colored points = FDR < 0.05 | Black outline = FDR < 0.01",
    x = "Developmental Stage",
    y = "Pathway"
  ) +

  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 9.5, color = "gray30"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 10),
    axis.text.y = element_text(size = 9),
    axis.title = element_text(face = "bold"),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    strip.background = element_rect(fill = "lightgray", color = "gray50"),
    strip.text = element_text(face = "bold", size = 11),
    panel.border = element_rect(color = "gray70", fill = NA),
    panel.spacing = unit(0.5, "lines"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

# Dynamic height
plot_height <- max(12, length(unique(ribosome_data$Pathway_Short)) * 0.25)

ggsave(file.path(out_dir, "Ribosome_Paradox_Three_Pools.pdf"),
       p_paradox, width = 12, height = plot_height, limitsize = FALSE)

message("✓ Three-pool dotplot complete\n")

###############################################################################
##  FIGURE 2: Temporal Trajectory (Early vs Maturation)                     ##
###############################################################################

message("📊 Creating temporal trajectory plot...\n")

# Filter to only Early and Late for comparable trajectory
trajectory_data <- ribosome_data %>%
  filter(Timepoint %in% c("Early (D35)", "Late (D65)")) %>%
  group_by(Pool, Timepoint, Mutation) %>%
  summarize(
    Mean_NES = mean(NES, na.rm = TRUE),
    N_pathways = n(),
    .groups = "drop"
  )

# Create line plot
p_trajectory <- ggplot(trajectory_data,
                       aes(x = Timepoint, y = Mean_NES,
                           color = Pool, group = Pool)) +
  geom_line(linewidth = 1.5, alpha = 0.8) +
  geom_point(size = 4, alpha = 0.9) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +

  facet_wrap(~Mutation, ncol = 2) +

  scale_color_manual(
    values = c(
      "1. Cytoplasmic Ribosome Biogenesis" = "#B35806",
      "2. Synaptic Ribosomes (SynGO)" = "#2166AC",
      "3. Mitochondrial Ribosomes (MitoCarta)" = "#D4A574"
    ),
    name = "Ribosome Pool"
  ) +

  labs(
    title = "Temporal Trajectory: Early vs Late Mutation Effects",
    subtitle = "Mean NES across ribosome pools (comparable mutation effects at D35 and D65)",
    x = "Developmental Stage",
    y = "Mean NES (Enrichment Score)"
  ) +

  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, color = "gray30"),
    axis.text = element_text(face = "bold"),
    axis.title = element_text(face = "bold", size = 11),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    strip.background = element_rect(fill = "lightgray"),
    strip.text = element_text(face = "bold", size = 12),
    panel.border = element_rect(color = "gray70", fill = NA),
    panel.grid.minor = element_blank()
  ) +
  guides(color = guide_legend(nrow = 3, byrow = TRUE))

ggsave(file.path(out_dir, "Ribosome_Temporal_Trajectory.pdf"),
       p_trajectory, width = 10, height = 7)

message("✓ Temporal trajectory complete\n")

###############################################################################
##  Export Data                                                              ##
###############################################################################

message("📝 Exporting ribosome pathway data...\n")

write.csv(ribosome_data %>%
            select(Pool, Pathway_Short, Description, Contrast, Mutation, Timepoint,
                   NES, p.adjust, pvalue, setSize) %>%
            arrange(Pool, Mutation, Timepoint, p.adjust),
          file.path(out_dir, "Ribosome_Paradox_Data.csv"),
          row.names = FALSE)

message("✓ CSV exports complete\n")

###############################################################################
##  Summary Report                                                           ##
###############################################################################

cat("\n=== Ribosome Paradox Summary ===\n\n")

cat("Pathways by Pool:\n")
pool_counts <- ribosome_data %>%
  group_by(Pool) %>%
  summarize(N = n(), .groups = "drop")
print(pool_counts)

cat("\n\nBIOLOGICAL INTERPRETATION (TRAJECTORY FRAMEWORK):\n")
cat("  Early (D35): Initial mutation effects in young neurons\n")
cat("  TrajDev: Mutation-specific maturation divergence from normal trajectory\n")
cat("  Late (D65): Final outcome in mature neurons\n")
cat("  Mathematical relationship: Late = Early + TrajDev\n\n")

cat("  Pool 1 (Cytoplasmic Biogenesis):\n")
cat("    → Early/Late: UPREGULATED (compensatory response)\n")
cat("    → Cells attempt to produce more ribosomes\n\n")

cat("  Pool 2 (Synaptic Ribosomes):\n")
cat("    → Early: Baseline defect\n")
cat("    → TrajDev: Progressive worsening during maturation\n")
cat("    → Late: Severe downregulation (failed compensation)\n\n")

cat("  Pool 3 (Mitochondrial Ribosomes):\n")
cat("    → TrajDev: Compensatory increase during maturation\n")
cat("    → Late: Partial recovery (mitochondria attempt OXPHOS rescue)\n\n")

message("✅ Ribosome paradox visualization complete!")
message(sprintf("📁 Output directory: %s", out_dir))
message("\n🎯 KEY OUTPUTS:")
message("  1. Ribosome_Paradox_Three_Pools.pdf - Main dotplot showing all pathways")
message("  2. Ribosome_Temporal_Trajectory.pdf - Temporal dynamics across development")
message("  3. Ribosome_Paradox_Data.csv - Complete pathway data")
