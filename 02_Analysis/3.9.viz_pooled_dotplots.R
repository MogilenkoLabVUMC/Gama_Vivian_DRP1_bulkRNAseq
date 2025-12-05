###############################################################################
##  3.9 Cross-Database Pooled Dotplots                                       ##
##  Creates pooled GSEA dotplots with:                                        ##
##    - Gene Ratio on X-axis (correctly calculated)                           ##
##    - NES gradient coloring (Blue-White-Orange from color_config.R)         ##
##    - Size proportional to significance (-log10 p.adjust)                   ##
##    - Black outline for significant pathways (padj < 0.05)                  ##
##                                                                             ##
##  Output: 03_Results/02_Analysis/Plots/GSEA/Cross_database_pooled/          ##
###############################################################################

library(here)
library(ggplot2)
library(dplyr)
library(stringr)

# -------------------------------------------------------------------- #
# 0.  Configuration & Helper Loading                                    #
# -------------------------------------------------------------------- #
config <- list(
  out_root    = "03_Results/02_Analysis",
  helper_root = "01_Scripts/RNAseq-toolkit"
)

# Source the shared GSEA dotplot helpers
# (includes color_config.R, calculate_gene_ratio, filter_neuronal_pathways,
#  create_gsea_pooled_dotplot)
source(here::here("01_Scripts/R_scripts/gsea_dotplot_helpers.R"))

# Directory setup
checkpoint_dir <- here::here(config$out_root, "checkpoints")
pooled_dir_comp <- here::here(config$out_root, "Plots/GSEA/Cross_database_pooled/comprehensive")
pooled_dir_focus <- here::here(config$out_root, "Plots/GSEA/Cross_database_pooled/focused")
pooled_csv_dir <- here::here(config$out_root, "Plots/GSEA/Cross_database_pooled/csv_data")
dir.create(pooled_dir_comp, recursive = TRUE, showWarnings = FALSE)
dir.create(pooled_dir_focus, recursive = TRUE, showWarnings = FALSE)
dir.create(pooled_csv_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------- #
# Load GSEA Checkpoints                                                 #
# -------------------------------------------------------------------- #
message("\n", paste(rep("=", 70), collapse = ""))
message("CROSS-DATABASE POOLED DOTPLOTS")
message("Gene Ratio on X-axis | NES gradient | Significance outlines")
message(paste(rep("=", 70), collapse = ""))

message("\n[1/5] Loading checkpoints...")
all_gsea_results   <- readRDS(file.path(checkpoint_dir, "all_gsea_results.rds"))
syngo_gsea_results <- readRDS(file.path(checkpoint_dir, "syngo_gsea_results.rds"))
message("Checkpoints loaded successfully")

# Database collections
all_databases <- c("hallmark", "gobp", "gocc", "gomf", "kegg",
                   "reactome", "wiki", "canon", "tf")
focused_databases <- c("hallmark", "kegg", "reactome")

# Priority contrasts
priority_contrasts <- c(
  "Maturation_G32A_specific",
  "Maturation_R403C_specific",
  "Time_Ctrl",
  "G32A_vs_Ctrl_D65",
  "R403C_vs_Ctrl_D65"
)

# -------------------------------------------------------------------- #
# Step 1: Extract and Save Raw CSV Data                                 #
# -------------------------------------------------------------------- #
message("\n[2/5] Extracting raw pathway data...")

for (co in priority_contrasts) {
  message(sprintf("  Processing %s...", co))

  if (!co %in% names(all_gsea_results)) {
    message(sprintf("    Warning: Contrast not found: %s", co))
    next
  }

  all_pathways <- data.frame()

  # Extract from all databases
  for (db in all_databases) {
    if (db %in% names(all_gsea_results[[co]])) {
      db_result <- all_gsea_results[[co]][[db]]

      if (!is.null(db_result)) {
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
          all_pathways <- rbind(all_pathways, res)
        }
      }
    }
  }

  # Add SynGO
  if (co %in% names(syngo_gsea_results)) {
    db_result <- syngo_gsea_results[[co]]
    if (!is.null(db_result)) {
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
        all_pathways <- rbind(all_pathways, res)
      }
    }
  }

  # Save comprehensive raw CSV
  if (nrow(all_pathways) > 0) {
    write.csv(all_pathways,
              file.path(pooled_csv_dir, sprintf("%s_raw_comprehensive.csv", co)),
              row.names = FALSE)

    # Focused version
    focused_pathways <- all_pathways %>%
      filter(Database %in% c(focused_databases, "syngo"))

    if (nrow(focused_pathways) > 0) {
      write.csv(focused_pathways,
                file.path(pooled_csv_dir, sprintf("%s_raw_focused.csv", co)),
                row.names = FALSE)
    }
  }
}

# -------------------------------------------------------------------- #
# Step 2: Apply Neuronal Filtering (using helper function)              #
# -------------------------------------------------------------------- #
message("\n[3/5] Applying neuronal maturation filters...")

for (co in priority_contrasts) {
  message(sprintf("  Processing %s...", co))

  # Filter comprehensive version
  comp_csv <- file.path(pooled_csv_dir, sprintf("%s_raw_comprehensive.csv", co))
  if (file.exists(comp_csv)) {
    comp_data <- read.csv(comp_csv, stringsAsFactors = FALSE)
    comp_filtered <- filter_neuronal_pathways(comp_data)
    write.csv(comp_filtered,
              file.path(pooled_csv_dir, sprintf("%s_filtered_comprehensive.csv", co)),
              row.names = FALSE)
  }

  # Filter focused version
  focus_csv <- file.path(pooled_csv_dir, sprintf("%s_raw_focused.csv", co))
  if (file.exists(focus_csv)) {
    focus_data <- read.csv(focus_csv, stringsAsFactors = FALSE)
    focus_filtered <- filter_neuronal_pathways(focus_data)
    write.csv(focus_filtered,
              file.path(pooled_csv_dir, sprintf("%s_filtered_focused.csv", co)),
              row.names = FALSE)
  }
}

# -------------------------------------------------------------------- #
# Step 3: Generate Comprehensive Visualizations                        #
# -------------------------------------------------------------------- #
message("\n[4/5] Generating comprehensive visualizations (with significance outlines)...")

for (co in priority_contrasts) {
  message(sprintf("  Processing %s...", co))

  filtered_csv <- file.path(pooled_csv_dir, sprintf("%s_filtered_comprehensive.csv", co))

  if (!file.exists(filtered_csv)) {
    message("    Warning: Filtered CSV not found")
    next
  }

  filtered_data <- read.csv(filtered_csv, stringsAsFactors = FALSE)

  if (nrow(filtered_data) == 0) {
    message("    Warning: No pathways after filtering")
    next
  }

  # Use helper function with significance outlines enabled
  p <- create_gsea_pooled_dotplot(
    filtered_data,
    title = sprintf("Neuronal Maturation Pathways: %s", co),
    subtitle = sprintf("Top 5 pathways per database (n=%d databases, FDR < 0.05, black outline = significant)",
                       length(unique(filtered_data$Database))),
    n_per_database = 5,
    highlight_sig = TRUE,
    sig_cutoff = 0.05
  )

  out_file <- file.path(pooled_dir_comp, sprintf("%s_pooled_comprehensive.pdf", co))
  ggsave(out_file, p, width = 10, height = 12)
  message(sprintf("    Saved: %s", basename(out_file)))
}

# -------------------------------------------------------------------- #
# Step 4: Generate Focused Visualizations                               #
# -------------------------------------------------------------------- #
message("\n[5/5] Generating focused visualizations (with significance outlines)...")

for (co in priority_contrasts) {
  message(sprintf("  Processing %s...", co))

  filtered_csv <- file.path(pooled_csv_dir, sprintf("%s_filtered_focused.csv", co))

  if (!file.exists(filtered_csv)) {
    message("    Warning: Filtered CSV not found")
    next
  }

  filtered_data <- read.csv(filtered_csv, stringsAsFactors = FALSE)

  if (nrow(filtered_data) == 0) {
    message("    Warning: No pathways after filtering")
    next
  }

  # Use helper function with significance outlines enabled
  p <- create_gsea_pooled_dotplot(
    filtered_data,
    title = sprintf("Neuronal Maturation Pathways: %s", co),
    subtitle = sprintf("Top 10 pathways per database (n=%d databases, FDR < 0.05, black outline = significant)",
                       length(unique(filtered_data$Database))),
    n_per_database = 10,
    highlight_sig = TRUE,
    sig_cutoff = 0.05
  )

  out_file <- file.path(pooled_dir_focus, sprintf("%s_pooled_focused.pdf", co))
  ggsave(out_file, p, width = 9, height = 10)
  message(sprintf("    Saved: %s", basename(out_file)))
}

# -------------------------------------------------------------------- #
# Final Summary                                                         #
# -------------------------------------------------------------------- #
message("\n", paste(rep("=", 70), collapse = ""))
message("COMPLETE!")
message(paste(rep("=", 70), collapse = ""))
message(sprintf("\nComprehensive plots: %s", pooled_dir_comp))
message(sprintf("Focused plots: %s", pooled_dir_focus))
message(sprintf("CSV data: %s", pooled_csv_dir))
message("\nKey features:")
message("  - X-axis: Gene Ratio (Leading Edge / Set Size)")
message("  - Color: NES (Blue-White-Orange colorblind-safe gradient)")
message("  - Size: -log10(FDR)")
message("  - Black outline: Significant pathways (padj < 0.05)")
