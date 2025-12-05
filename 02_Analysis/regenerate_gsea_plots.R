###############################################################################
##  Surgical GSEA Plot Regeneration                                          ##
##  Uses checkpointed GSEA results - NO recalculation needed                 ##
###############################################################################
##  Purpose: Regenerate all GSEA plots with updated color scheme             
##  Note: This script uses the unified color palette from color_config.R     ##
###############################################################################

# -------------------------------------------------------------------- #
# 0.  Configuration                                                     #
# -------------------------------------------------------------------- #
config <- list(
  out_root     = "03_Results/02_Analysis",
  helper_root  = "01_Scripts/RNAseq-toolkit",
  n_pathways   = 30,
  padj_cutoff  = 0.05
)

checkpoint_dir <- here::here(config$out_root, "checkpoints")
gsea_root <- here::here(config$out_root, "Plots/GSEA")

# -------------------------------------------------------------------- #
# 1.  Load required packages                                            #
# -------------------------------------------------------------------- #
message("Loading required packages...")
required_pkgs <- c("ggplot2", "dplyr", "here", "methods", "stringr")

for (p in required_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    message(sprintf("Installing %s...", p))
    install.packages(p, repos = "https://cloud.r-project.org")
  }
  library(p, character.only = TRUE)
}

# -------------------------------------------------------------------- #
# 2.  Source helper functions                                           #
# -------------------------------------------------------------------- #
message("\nSourcing helper functions...")

source_if_present <- function(...) {
  path <- here::here(...)
  if (file.exists(path)) {
    source(path, echo = FALSE)
    message("  Loaded: ", basename(path))
  } else {
    warning("Helper not found: ", path)
  }
}

# Load unified color configuration
source(here::here("01_Scripts/R_scripts/color_config.R"))

# GSEA plotting helpers (now use updated colors from color_config.R)
gsea_helpers <- c(
  "scripts/custom_minimal_theme.R",
  "scripts/GSEA/GSEA_plotting/gsea_plotting_utils.R",
  "scripts/GSEA/GSEA_plotting/format_pathway_names.R",
  "scripts/GSEA/GSEA_plotting/gsea_dotplot.R",
  "scripts/GSEA/GSEA_plotting/gsea_dotplot_facet.R",
  "scripts/GSEA/GSEA_plotting/gsea_barplot.R",
  "scripts/GSEA/GSEA_plotting/gsea_running_sum_plot.R",
  "scripts/GSEA/GSEA_plotting/gsea_heatmap.R",
  "scripts/GSEA/GSEA_plotting/plot_all_gsea_results.R"
)

for (h in gsea_helpers) {
  source_if_present(config$helper_root, h)
}

# -------------------------------------------------------------------- #
# 3.  Load checkpointed GSEA results                                    #
# -------------------------------------------------------------------- #
message("\nLoading checkpointed GSEA results...")

# Main GSEA results (MSigDB databases)
all_gsea_file <- file.path(checkpoint_dir, "all_gsea_results.rds")
if (!file.exists(all_gsea_file)) {
  stop("Checkpoint not found: ", all_gsea_file, "\n",
       "Please run the main pipeline first: Rscript 02_Analysis/1.1.main_pipeline.R")
}
all_gsea_results <- readRDS(all_gsea_file)
message("  Loaded all_gsea_results.rds: ", length(all_gsea_results), " contrasts")

# SynGO GSEA results
syngo_gsea_file <- file.path(checkpoint_dir, "syngo_gsea_results.rds")
syngo_gsea_results <- NULL
if (file.exists(syngo_gsea_file)) {
  syngo_gsea_results <- readRDS(syngo_gsea_file)
  message("  Loaded syngo_gsea_results.rds: ", length(syngo_gsea_results), " contrasts")
} else {
  message("  SynGO checkpoint not found (skipping)")
}

# MitoCarta GSEA results
mitocarta_gsea_file <- file.path(checkpoint_dir, "mitocarta_gsea_results.rds")
mitocarta_gsea_results <- NULL
if (file.exists(mitocarta_gsea_file)) {
  mitocarta_gsea_results <- readRDS(mitocarta_gsea_file)
  message("  Loaded mitocarta_gsea_results.rds: ", length(mitocarta_gsea_results), " contrasts")
} else {
  message("  MitoCarta checkpoint not found (skipping)")
}

# -------------------------------------------------------------------- #
# 4.  Regenerate plots for MSigDB GSEA                                  #
# -------------------------------------------------------------------- #
message("\n", paste(rep("=", 70), collapse = ""))
message("REGENERATING GSEA PLOTS (MSigDB databases)")
message(paste(rep("=", 70), collapse = ""))

for (co in names(all_gsea_results)) {
  message("\n>>> Processing contrast: ", co)
  this_res_list <- all_gsea_results[[co]]

  if (is.null(this_res_list) || !length(this_res_list)) {
    message("    Skipping ", co, " (no results)")
    next
  }

  plot_all_gsea_results(
    gsea_list = this_res_list,
    analysis_name = co,
    out_root = gsea_root,
    n_pathways = config$n_pathways,
    padj_cutoff = config$padj_cutoff
  )
}
message("\n MSigDB GSEA plots complete")

# -------------------------------------------------------------------- #
# 5.  Regenerate plots for SynGO GSEA                                   #
# -------------------------------------------------------------------- #
if (!is.null(syngo_gsea_results)) {
  message("\n", paste(rep("=", 70), collapse = ""))
  message("REGENERATING SynGO GSEA PLOTS")
  message(paste(rep("=", 70), collapse = ""))

  for (co in names(syngo_gsea_results)) {
    message("\n>>> Processing SynGO for contrast: ", co)
    syngo_res <- syngo_gsea_results[[co]]

    if (is.null(syngo_res)) {
      message("    Skipping ", co, " (no results)")
      next
    }

    # SynGO results are stored as single gseaResult objects, wrap in list
    syngo_list <- list(SynGO = syngo_res)

    plot_all_gsea_results(
      gsea_list = syngo_list,
      analysis_name = co,
      out_root = gsea_root,
      n_pathways = config$n_pathways,
      padj_cutoff = config$padj_cutoff
    )
  }
  message("\n SynGO GSEA plots complete")
}

# -------------------------------------------------------------------- #
# 6.  Regenerate plots for MitoCarta GSEA                               #
# -------------------------------------------------------------------- #
if (!is.null(mitocarta_gsea_results)) {
  message("\n", paste(rep("=", 70), collapse = ""))
  message("REGENERATING MitoCarta GSEA PLOTS")
  message(paste(rep("=", 70), collapse = ""))

  for (co in names(mitocarta_gsea_results)) {
    message("\n>>> Processing MitoCarta for contrast: ", co)
    mito_res <- mitocarta_gsea_results[[co]]

    if (is.null(mito_res)) {
      message("    Skipping ", co, " (no results)")
      next
    }

    # MitoCarta results are stored as single gseaResult objects, wrap in list
    mito_list <- list(MitoCarta = mito_res)

    plot_all_gsea_results(
      gsea_list = mito_list,
      analysis_name = co,
      out_root = gsea_root,
      n_pathways = config$n_pathways,
      padj_cutoff = config$padj_cutoff
    )
  }
  message("\n MitoCarta GSEA plots complete")
}

# -------------------------------------------------------------------- #
# 7.  Summary                                                           #
# -------------------------------------------------------------------- #
message("\n", paste(rep("=", 70), collapse = ""))
message("GSEA PLOT REGENERATION COMPLETE")
message(paste(rep("=", 70), collapse = ""))
message("\nOutput directory: ", gsea_root)
message("\nUsing unified color scheme:")
message("  Negative NES: ", DIVERGING_COLORS$negative, " (Blue)")
message("  Neutral:      ", DIVERGING_COLORS$neutral, " (White)")
message("  Positive NES: ", DIVERGING_COLORS$positive, " (Orange)")
message("\nRegenerated plots for:")
message("  - ", length(all_gsea_results), " contrasts (MSigDB)")
if (!is.null(syngo_gsea_results)) {
  message("  - ", length(syngo_gsea_results), " contrasts (SynGO)")
}
if (!is.null(mitocarta_gsea_results)) {
  message("  - ", length(mitocarta_gsea_results), " contrasts (MitoCarta)")
}
message(paste(rep("=", 70), collapse = ""))
