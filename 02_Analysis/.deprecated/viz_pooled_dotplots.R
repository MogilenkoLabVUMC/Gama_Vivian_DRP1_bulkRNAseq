###############################################################################
##  Cross-Database Pooled Dotplots – Integrated Pathway View                ##
##  Loads checkpoints, NO re-running of GSEA                                 ##
###############################################################################
## WORKFLOW:
##  1. Generate raw CSV files with all pathways (for inspection)
##  2. Apply neuronal maturation-specific filters
##  3. Create visualizations from filtered data
##
## Versions:
##  - Comprehensive: Multiple databases (MSigDB + SynGO)
##  - Focused: Key databases (Hallmark, KEGG, Reactome, SynGO)
##
## Priority contrasts:
##  - Maturation_G32A_specific (highlight ribosome findings)
##  - Maturation_R403C_specific (highlight ribosome findings)
##  - Time_Ctrl (normal maturation reference)
##  - G32A_vs_Ctrl_D65 (mature timepoint baseline)
##  - R403C_vs_Ctrl_D65 (mature timepoint baseline)
###############################################################################

library(here)
library(ggplot2)
library(dplyr)
library(stringr)

# -------------------------------------------------------------------- #
# 0.  Configuration & Checkpoint Loading                               #
# -------------------------------------------------------------------- #
config <- list(
  out_root    = "03_Results/02_Analysis",
  helper_root = "01_Scripts/RNAseq-toolkit"
)

checkpoint_dir <- here::here(config$out_root, "checkpoints")
pooled_dir_comp <- here::here(config$out_root, "Plots/GSEA/Cross_database_pooled/comprehensive")
pooled_dir_focus <- here::here(config$out_root, "Plots/GSEA/Cross_database_pooled/focused")
pooled_csv_dir <- here::here(config$out_root, "Plots/GSEA/Cross_database_pooled/csv_data")
dir.create(pooled_dir_comp, recursive = TRUE, showWarnings = FALSE)
dir.create(pooled_dir_focus, recursive = TRUE, showWarnings = FALSE)
dir.create(pooled_csv_dir, recursive = TRUE, showWarnings = FALSE)

## Load helper functions -------------------------------------------------
source_if_present <- function(...) {
  path <- here::here(...)
  if (file.exists(path)) {
    source(path, echo = FALSE)
  } else {
    warning("helper not found → ", path)
  }
}

source_if_present(config$helper_root, "scripts/GSEA/GSEA_plotting/plot_pooled_contrast_dotplot.R")

## Load checkpoints ------------------------------------------------------
message("📂 Loading checkpoints...")
all_gsea_results   <- readRDS(file.path(checkpoint_dir, "all_gsea_results.rds"))
syngo_gsea_results <- readRDS(file.path(checkpoint_dir, "syngo_gsea_results.rds"))
message("✓ Checkpoints loaded successfully")

# -------------------------------------------------------------------- #
# 1.  Prepare Database Collections                                     #
# -------------------------------------------------------------------- #

## All available databases (comprehensive version)
all_databases <- c("hallmark", "gobp", "gocc", "gomf", "kegg",
                   "reactome", "wiki", "canon", "tf")

## Focused databases (manuscript figures)
focused_databases <- c("hallmark", "kegg", "reactome")

## Priority contrasts for visualization
priority_contrasts <- c(
  "Maturation_G32A_specific",
  "Maturation_R403C_specific",
  "Time_Ctrl",
  "G32A_vs_Ctrl_D65",
  "R403C_vs_Ctrl_D65"
)

message("\n📊 Generating pooled dotplots for ", length(priority_contrasts), " contrasts...")
message("  Databases (comprehensive): ", paste(all_databases, collapse = ", "))
message("  Databases (focused): ", paste(focused_databases, collapse = ", "))

# -------------------------------------------------------------------- #
# 1.5  Generate Raw CSV Files (for inspection before filtering)       #
# -------------------------------------------------------------------- #
message("\n💾 STEP 1: Generating raw CSV files for inspection...")

for (co in priority_contrasts) {
  message(sprintf("\n  Processing %s...", co))

  if (!co %in% names(all_gsea_results)) {
    message(sprintf("    ⚠️  Contrast not found in GSEA results: %s", co))
    next
  }

  # Combine all databases (comprehensive)
  all_pathways_comp <- data.frame()

  for (db in all_databases) {
    if (db %in% names(all_gsea_results[[co]])) {
      db_result <- all_gsea_results[[co]][[db]]

      if (!is.null(db_result)) {
        # Extract result dataframe from S4 object
        if (is(db_result, "gseaResult")) {
          res <- as.data.frame(db_result@result)
        } else if ("result" %in% names(db_result)) {
          res <- db_result$result
        } else {
          res <- as.data.frame(db_result)
        }

        res <- res %>%
          filter(p.adjust < 0.05) %>%
          arrange(p.adjust)

        if (nrow(res) > 0) {
          res$Database <- db
          res$Contrast <- co
          all_pathways_comp <- rbind(all_pathways_comp, res)
        }
      }
    }
  }

  # Add SynGO
  if (co %in% names(syngo_gsea_results)) {
    db_result <- syngo_gsea_results[[co]]
    if (!is.null(db_result)) {
      # Extract result dataframe from S4 object
      if (is(db_result, "gseaResult")) {
        res <- as.data.frame(db_result@result)
      } else if ("result" %in% names(db_result)) {
        res <- db_result$result
      } else {
        res <- as.data.frame(db_result)
      }

      res <- res %>%
        filter(p.adjust < 0.05) %>%
        arrange(p.adjust)

      if (nrow(res) > 0) {
        res$Database <- "syngo"
        res$Contrast <- co
        all_pathways_comp <- rbind(all_pathways_comp, res)
      }
    }
  }

  # Save comprehensive raw CSV
  if (nrow(all_pathways_comp) > 0) {
    raw_csv_file <- file.path(pooled_csv_dir, sprintf("%s_raw_comprehensive.csv", co))
    write.csv(all_pathways_comp, raw_csv_file, row.names = FALSE)
    message(sprintf("    ✓ Saved raw comprehensive CSV: %s (%d pathways)",
                    basename(raw_csv_file), nrow(all_pathways_comp)))

    # Generate focused version (only if we have data)
    all_pathways_focus <- all_pathways_comp %>%
      filter(Database %in% c(focused_databases, "syngo"))

    if (nrow(all_pathways_focus) > 0) {
      raw_csv_file <- file.path(pooled_csv_dir, sprintf("%s_raw_focused.csv", co))
      write.csv(all_pathways_focus, raw_csv_file, row.names = FALSE)
      message(sprintf("    ✓ Saved raw focused CSV: %s (%d pathways)",
                      basename(raw_csv_file), nrow(all_pathways_focus)))
    }
  } else {
    message("    ⚠️  No significant pathways found for this contrast")
  }
}

message("\n✅ Raw CSV files generated.")
message(sprintf("📁 CSV directory: %s", pooled_csv_dir))

# -------------------------------------------------------------------- #
# 1.6  Define Neuronal Maturation Filtering Function                  #
# -------------------------------------------------------------------- #

## Filter pathways to keep only neuronal maturation-relevant terms
filter_neuronal_pathways <- function(pathway_df) {

  if (nrow(pathway_df) == 0) return(pathway_df)

  ## EXCLUSION patterns (non-neuronal tissue-specific)
  exclude_patterns <- c(
    # Tissue-specific (non-neural)
    "PANCREA", "BETA_CELL", "ISLET",
    "CARDIAC", "HEART", "MYOCARDI", "CARDIOMYOCYTE",
    "KIDNEY", "RENAL", "NEPHRON",
    "INTESTIN", "COLON", "GUT", "DIGESTIVE",
    "LUNG", "ALVEOL", "BRONCH",
    "LIVER", "HEPAT",
    "BREAST", "MAMMARY",
    "PROSTAT",
    "MYOGENES",  # Myogenesis (muscle development)
    "ADIPOGEN", "ADIPOCYTE",
    # Metabolic (non-neuronal specific)
    "BILE_ACID",
    "XENOBIOTIC",
    "HEME_METABOL",
    "ANDROGEN_RESPONSE",
    # Blood/vascular (unless brain-specific)
    "COAGULATION",
    "BLOOD_CLOT",
    # Other
    "UV_RESPONSE"
  )

  ## INCLUSION patterns (neuronal/brain-relevant - override exclusions)
  include_patterns <- c(
    # Neuronal-specific
    "NEURO", "NEURAL", "NEURON", "NERVE",
    "SYNAP", "SYNAPT",
    "AXON", "DENDRIT",
    "BRAIN", "CORTEX", "HIPPOCAM", "CEREBR",
    "GLIAL", "ASTROCYT", "OLIGODENDRO",
    "MUSCLE.*NEURO", "NEUROMUSC",  # Keep neuromuscular
    # Neuronal structures/functions
    "RIBOSOM", "TRANSLATION", "PROTEIN_SYNTHES",
    "MITOCHOND", "OXIDATIVE_PHOSPHORYL",
    # Relevant cell biology
    "E2F_TARGET", "G2M_CHECKPOINT",
    "MYC_TARGET",
    "APOPTOSIS", "P53_PATHWAY",
    "HYPOXIA",
    "DNA_REPAIR",
    "TNFA_SIGNALING",
    "INTERFERON", "INFLAMMAT",
    "EPITHELIAL_MESENCHYMAL",
    "MITOTIC_SPINDLE",
    "PROTEIN_SECRETION",
    "IL2_STAT5",
    "GLYCOLYSIS"  # Energy metabolism relevant to neurons
  )

  ## SynGO database: always keep (100% neuronal-relevant)
  syngo_mask <- grepl("syngo", pathway_df$Database, ignore.case = TRUE)

  ## Apply exclusions
  exclude_regex <- paste(exclude_patterns, collapse = "|")
  exclude_mask <- grepl(exclude_regex, pathway_df$Description, ignore.case = TRUE) |
                  grepl(exclude_regex, pathway_df$ID, ignore.case = TRUE)

  ## Apply inclusions (override exclusions)
  include_regex <- paste(include_patterns, collapse = "|")
  include_mask <- grepl(include_regex, pathway_df$Description, ignore.case = TRUE) |
                  grepl(include_regex, pathway_df$ID, ignore.case = TRUE)

  ## Final filter: (not excluded) OR (included) OR (syngo)
  keep_mask <- (!exclude_mask) | include_mask | syngo_mask

  filtered_df <- pathway_df[keep_mask, ]

  ## Report filtering
  n_before <- nrow(pathway_df)
  n_after <- nrow(filtered_df)
  n_removed <- n_before - n_after

  if (n_removed > 0) {
    message(sprintf("    Filtered: %d pathways removed (%d → %d)",
                    n_removed, n_before, n_after))
  }

  return(filtered_df)
}

# -------------------------------------------------------------------- #
# 1.7  Apply Filtering and Save Filtered CSVs                         #
# -------------------------------------------------------------------- #
message("\n🔬 STEP 2: Applying neuronal maturation filters...")

for (co in priority_contrasts) {
  message(sprintf("\n  Processing %s...", co))

  # Load and filter comprehensive version
  comp_csv <- file.path(pooled_csv_dir, sprintf("%s_raw_comprehensive.csv", co))
  if (file.exists(comp_csv)) {
    comp_data <- read.csv(comp_csv, stringsAsFactors = FALSE)
    comp_filtered <- filter_neuronal_pathways(comp_data)

    # Save filtered comprehensive CSV
    filtered_csv <- file.path(pooled_csv_dir, sprintf("%s_filtered_comprehensive.csv", co))
    write.csv(comp_filtered, filtered_csv, row.names = FALSE)
    message(sprintf("    ✓ Saved filtered comprehensive CSV: %s", basename(filtered_csv)))
  }

  # Load and filter focused version
  focus_csv <- file.path(pooled_csv_dir, sprintf("%s_raw_focused.csv", co))
  if (file.exists(focus_csv)) {
    focus_data <- read.csv(focus_csv, stringsAsFactors = FALSE)
    focus_filtered <- filter_neuronal_pathways(focus_data)

    # Save filtered focused CSV
    filtered_csv <- file.path(pooled_csv_dir, sprintf("%s_filtered_focused.csv", co))
    write.csv(focus_filtered, filtered_csv, row.names = FALSE)
    message(sprintf("    ✓ Saved filtered focused CSV: %s", basename(filtered_csv)))
  }
}

message("\n✅ Filtering complete. Filtered CSVs ready for visualization.")

# -------------------------------------------------------------------- #
# 2.  Generate Comprehensive Version (Using Filtered Data)            #
# -------------------------------------------------------------------- #
message("\n📊 STEP 3: Generating comprehensive visualizations...")

for (co in priority_contrasts) {
  message(sprintf("\n  Processing %s...", co))

  # Load filtered CSV data
  filtered_csv <- file.path(pooled_csv_dir, sprintf("%s_filtered_comprehensive.csv", co))

  if (!file.exists(filtered_csv)) {
    message(sprintf("    ⚠️  Filtered CSV not found: %s", basename(filtered_csv)))
    next
  }

  filtered_data <- read.csv(filtered_csv, stringsAsFactors = FALSE)

  if (nrow(filtered_data) == 0) {
    message("    ⚠️  No pathways after filtering")
    next
  }

  message(sprintf("    Loaded %d filtered pathways", nrow(filtered_data)))

  # Create manual dotplot from filtered CSV data
  tryCatch({
    # Select top N pathways per database
    top_pathways <- filtered_data %>%
      group_by(Database) %>%
      arrange(p.adjust) %>%
      slice_head(n = 5) %>%
      ungroup()

    # Determine x-axis limits based on data range
    nes_range <- range(top_pathways$NES)
    x_limits <- c(floor(min(nes_range) * 1.1), ceiling(max(nes_range) * 1.1))

    # Create dotplot
    p <- ggplot(top_pathways,
                aes(x = NES, y = reorder(Description, NES))) +
      geom_point(aes(size = -log10(p.adjust), color = NES), alpha = 0.8) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
      scale_color_gradient2(low = "blue", mid = "gray80", high = "red",
                            midpoint = 0, name = "NES") +
      scale_size_continuous(name = "-log10(FDR)", range = c(3, 10)) +
      scale_x_continuous(limits = x_limits) +
      labs(title = sprintf("Neuronal Maturation Pathways: %s", co),
           subtitle = sprintf("Top 5 pathways per database (n=%d databases, FDR < 0.05, filtered for neuronal relevance)",
                              length(unique(top_pathways$Database))),
           x = "Normalized Enrichment Score (NES)",
           y = NULL) +
      theme_minimal(base_size = 11) +
      theme(
        axis.text.y = element_text(size = 9),
        legend.position = "right",
        panel.grid.major.y = element_line(color = "gray90"),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold")
      )

    # Save plot
    out_file <- file.path(pooled_dir_comp, sprintf("%s_pooled_comprehensive.pdf", co))
    ggsave(out_file, p, width = 10, height = 12)
    message(sprintf("    ✓ Saved: %s", basename(out_file)))

  }, error = function(e) {
    message(sprintf("    ⚠️  Error creating plot: %s", e$message))
  })
}

# -------------------------------------------------------------------- #
# 3.  Generate Focused Version (Using Filtered Data)                  #
# -------------------------------------------------------------------- #
message("\n📌 STEP 4: Generating focused visualizations...")

for (co in priority_contrasts) {
  message(sprintf("\n  Processing %s...", co))

  # Load filtered CSV data
  filtered_csv <- file.path(pooled_csv_dir, sprintf("%s_filtered_focused.csv", co))

  if (!file.exists(filtered_csv)) {
    message(sprintf("    ⚠️  Filtered CSV not found: %s", basename(filtered_csv)))
    next
  }

  filtered_data <- read.csv(filtered_csv, stringsAsFactors = FALSE)

  if (nrow(filtered_data) == 0) {
    message("    ⚠️  No pathways after filtering")
    next
  }

  message(sprintf("    Loaded %d filtered pathways", nrow(filtered_data)))

  # Create manual dotplot from filtered CSV data
  tryCatch({
    # Select top N pathways per database
    top_pathways <- filtered_data %>%
      group_by(Database) %>%
      arrange(p.adjust) %>%
      slice_head(n = 10) %>%
      ungroup()

    # Determine x-axis limits based on data range
    nes_range <- range(top_pathways$NES)
    x_limits <- c(floor(min(nes_range) * 1.1), ceiling(max(nes_range) * 1.1))

    # Create dotplot
    p <- ggplot(top_pathways,
                aes(x = NES, y = reorder(Description, NES))) +
      geom_point(aes(size = -log10(p.adjust), color = NES), alpha = 0.8) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
      scale_color_gradient2(low = "blue", mid = "gray80", high = "red",
                            midpoint = 0, name = "NES") +
      scale_size_continuous(name = "-log10(FDR)", range = c(3, 10)) +
      scale_x_continuous(limits = x_limits) +
      labs(title = sprintf("Neuronal Maturation Pathways: %s", co),
           subtitle = sprintf("Top 10 pathways per database (n=%d databases, FDR < 0.05, filtered for neuronal relevance)",
                              length(unique(top_pathways$Database))),
           x = "Normalized Enrichment Score (NES)",
           y = NULL) +
      theme_minimal(base_size = 11) +
      theme(
        axis.text.y = element_text(size = 9),
        legend.position = "right",
        panel.grid.major.y = element_line(color = "gray90"),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold")
      )

    # Save plot
    out_file <- file.path(pooled_dir_focus, sprintf("%s_pooled_focused.pdf", co))
    ggsave(out_file, p, width = 9, height = 10)
    message(sprintf("    ✓ Saved: %s", basename(out_file)))

  }, error = function(e) {
    message(sprintf("    ⚠️  Error creating plot: %s", e$message))
  })
}

# -------------------------------------------------------------------- #
# 4.  Final Summary Report                                             #
# -------------------------------------------------------------------- #
message("\n📝 Generating summary report...")

sink(file.path(pooled_csv_dir, "../pooled_dotplots_summary.txt"))
cat("===============================================================================\n")
cat("CROSS-DATABASE POOLED DOTPLOT ANALYSIS (NEURONAL MATURATION FILTERED)\n")
cat("===============================================================================\n\n")

cat("ANALYSIS DATE:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("WORKFLOW:\n")
cat("  1. Extract significant pathways from GSEA results (FDR < 0.05)\n")
cat("  2. Apply neuronal maturation filtering (remove non-neuronal tissues)\n")
cat("  3. Generate visualizations from filtered data\n\n")

cat("PRIORITY CONTRASTS:\n")
for (co in priority_contrasts) {
  cat(sprintf("  - %s\n", co))
}
cat("\n")

cat("DATABASE COLLECTIONS:\n")
cat("\nComprehensive version:\n")
for (db in c(all_databases, "syngo")) {
  cat(sprintf("  - %s\n", db))
}
cat("\nFocused version (manuscript):\n")
for (db in c(focused_databases, "syngo")) {
  cat(sprintf("  - %s\n", db))
}

cat("\n\nFILTERING CRITERIA:\n")
cat("EXCLUDED (non-neuronal tissues):\n")
cat("  - Pancreas, beta cells\n")
cat("  - Cardiac, heart, myocytes\n")
cat("  - Kidney, renal\n")
cat("  - Intestinal, colon\n")
cat("  - Lung, respiratory\n")
cat("  - Liver, bile acid\n")
cat("  - Adipose, muscle (non-neural)\n")
cat("  - Xenobiotic, coagulation, heme metabolism\n\n")
cat("INCLUDED (neuronal-relevant):\n")
cat("  - Neuron, synapse, axon, dendrite\n")
cat("  - Brain, cortex, hippocampus, glial\n")
cat("  - Ribosome, mitochondria, oxidative phosphorylation\n")
cat("  - Cell cycle (E2F, G2M), apoptosis, DNA repair\n")
cat("  - All SynGO pathways (100% neuronal)\n")

cat("\n\nOUTPUT DIRECTORIES:\n")
cat(sprintf("  CSV data (raw + filtered): %s\n", pooled_csv_dir))
cat(sprintf("  Comprehensive plots: %s\n", pooled_dir_comp))
cat(sprintf("  Focused plots: %s\n", pooled_dir_focus))

cat("\n\nFILE NAMING CONVENTION:\n")
cat("  {contrast}_raw_comprehensive.csv - All pathways before filtering\n")
cat("  {contrast}_filtered_comprehensive.csv - Neuronal-filtered pathways\n")
cat("  {contrast}_pooled_comprehensive.pdf - Visualization (comprehensive)\n")
cat("  {contrast}_pooled_focused.pdf - Visualization (focused)\n")

cat("\n\nKEY PATHWAYS TO HIGHLIGHT:\n")
cat("  - Ribosome pathways (SynGO: Presynaptic/Postsynaptic Ribosome)\n")
cat("  - Mitochondrial pathways (Hallmark: OxPhos, KEGG: ETC)\n")
cat("  - Cell cycle pathways (Hallmark: E2F, G2M checkpoint)\n")
cat("  - Synaptic pathways (SynGO: Synapse, PSD)\n")

cat("\n===============================================================================\n")
sink()

message("\n✅ Cross-database pooled dotplot analysis complete!")
message(sprintf("📁 CSV data: %s", pooled_csv_dir))
message(sprintf("📁 Comprehensive plots: %s", pooled_dir_comp))
message(sprintf("📁 Focused plots: %s", pooled_dir_focus))
