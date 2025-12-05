###############################################################################
##  MitoCarta GSEA Analysis – Separate Module                                ##
###############################################################################
##  Purpose: Run MitoCarta pathway enrichment analysis on differential
##           expression results from the main pipeline
##
##  Dependencies: Requires 1.Main_pipeline.R to have completed successfully
###############################################################################

# -------------------------------------------------------------------- #
# 0.  Configuration & Setup                                            #
# -------------------------------------------------------------------- #

# Load required packages
required_pkgs <- c("here", "clusterProfiler", "dplyr", "ggplot2")

for (p in required_pkgs){
  if (!requireNamespace(p, quietly = TRUE)){
    message(sprintf("• installing %s …", p))
    install.packages(p, repos = "https://cloud.r-project.org")
  }
  library(p, character.only = TRUE)
}

# -------------------------------------------------------------------- #
# 1.  Load Configuration & Checkpoint Data                             #
# -------------------------------------------------------------------- #

# Re-create config (must match main pipeline)
config <- list(
  out_root      = "03_Results/02_Analysis",
  helper_root   = "01_Scripts/RNAseq-toolkit",
  mitocarta_file = "00_Data/MitoCarta_3.0/MitoPathways3.0.gmx",
  force_recompute = FALSE   # Set TRUE to recompute MitoCarta GSEA
)

# Setup checkpoint directory
checkpoint_dir <- here::here(config$out_root, "checkpoints")
if (!dir.exists(checkpoint_dir)) {
  stop("❌ Checkpoint directory not found. Please run 1.Main_pipeline.R first.")
}

# Load checkpoint caching helper
load_or_compute <- function(checkpoint_file, compute_fn,
                            force_recompute = FALSE,
                            description = "computation") {
  if (!force_recompute && file.exists(checkpoint_file)) {
    message("📦 Loading cached ", description, " from ", basename(checkpoint_file))
    return(readRDS(checkpoint_file))
  }
  message("🔬 Computing ", description, "...")
  result <- compute_fn()
  message("💾 Saving ", description, " to ", basename(checkpoint_file))
  saveRDS(result, checkpoint_file)
  return(result)
}

# -------------------------------------------------------------------- #
# 2.  Load Required Data from Main Pipeline Checkpoints               #
# -------------------------------------------------------------------- #

message("\n📂 Loading data from main pipeline checkpoints...")

# Load contrast tables (contains DE results for all contrasts)
contrast_tables_file <- file.path(checkpoint_dir, "contrast_tables.rds")
if (!file.exists(contrast_tables_file)) {
  stop("❌ contrast_tables.rds not found. Please run 1.Main_pipeline.R first.")
}
contrast_tables <- readRDS(contrast_tables_file)
message("✓ Loaded contrast tables for ", length(contrast_tables), " contrasts")

# Load QC variables (contains sample annotation)
qc_vars_file <- file.path(checkpoint_dir, "qc_variables.rds")
if (!file.exists(qc_vars_file)) {
  stop("❌ qc_variables.rds not found. Please run 1.Main_pipeline.R first.")
}
qc_vars <- readRDS(qc_vars_file)
annot <- qc_vars$annot
message("✓ Loaded QC variables and sample annotations")

# -------------------------------------------------------------------- #
# 3.  Source Required Helper Functions                                 #
# -------------------------------------------------------------------- #

message("\n📚 Loading helper functions...")

source_if_present <- function(...) {
  path <- here::here(...)
  if (file.exists(path)) {
    source(path, echo = FALSE)
  } else {
    warning("helper not found → ", path)
  }
}

# Source MitoCarta-specific helpers
source_if_present("01_Scripts/R_scripts/parse_mitocarta_gmx.R")
source_if_present("01_Scripts/R_scripts/run_mitocarta_gsea.R")

# Source plotting utilities
source_if_present(config$helper_root, "scripts/custom_minimal_theme.R")
source_if_present(config$helper_root, "scripts/utils_plotting.R")
source_if_present(config$helper_root, "scripts/GSEA/GSEA_plotting/gsea_plotting_utils.R")
source_if_present(config$helper_root, "scripts/GSEA/GSEA_plotting/format_pathway_names.R")
source_if_present(config$helper_root, "scripts/GSEA/GSEA_plotting/gsea_dotplot.R")
source_if_present(config$helper_root, "scripts/GSEA/GSEA_plotting/gsea_barplot.R")
source_if_present(config$helper_root, "scripts/GSEA/GSEA_plotting/gsea_dotplot_facet.R")
source_if_present("01_Scripts/R_scripts/syngo_running_sum_plot.R")

message("✓ Helper functions loaded")

# -------------------------------------------------------------------- #
# 4.  Parse MitoCarta Pathways                                         #
# -------------------------------------------------------------------- #

message("\n🔬 Parsing MitoCarta pathways...")

# Check if MitoCarta file exists
mitocarta_path <- here::here(config$mitocarta_file)
if (!file.exists(mitocarta_path)) {
  stop("❌ MitoCarta GMX file not found at: ", mitocarta_path)
}

# Parse MitoCarta GMX file into TERM2GENE and TERM2NAME format
mitocarta_lists <- mitocarta_gmt(mitocarta_path)
message("✓ MitoCarta pathways parsed successfully")

# -------------------------------------------------------------------- #
# 5.  Run MitoCarta GSEA                                               #
# -------------------------------------------------------------------- #

message("\n🧬 Running MitoCarta GSEA analysis...")
message("   Contrasts to process: ", paste(names(contrast_tables), collapse = ", "))
message("   Force recompute: ", config$force_recompute)

# Use checkpoint caching for MitoCarta GSEA
mitocarta_gsea_results <- load_or_compute(
  checkpoint_file = file.path(checkpoint_dir, "mitocarta_gsea_results.rds"),
  force_recompute = config$force_recompute,
  description = "MitoCarta GSEA for all contrasts",
  compute_fn = function() {
    # Create result storage
    results <- list()

    # Run MitoCarta GSEA for each contrast
    message("Running MitoCarta GSEA for all contrasts...")
    for (co in names(contrast_tables)) {
      # Close any lingering graphic devices before starting new contrast
      close_all_devices()

      message("\n==== Working on contrast: ", co, " ====")
      message("   Genes in contrast table: ", nrow(contrast_tables[[co]]))

      # Capture the result with error handling
      gsea_out_dir <- here::here(config$out_root, "Plots/GSEA", co, "MitoCarta")

      safe_res <- tryCatch({
        run_mitocarta_gsea(
          contrast_tables[[co]],
          co,
          T2G = mitocarta_lists$T2G,
          T2N = mitocarta_lists$T2N,
          #nperm = 100000,  # default 100000
          output_dir = gsea_out_dir
        )
      }, error = function(e) {
        message("❗ ERROR in MitoCarta GSEA for contrast ", co, ": ", e$message)
        message("   Traceback:")
        print(traceback())
        NULL
      })

      if (!is.null(safe_res)) {
        message("✓ Completed ", co, " - ", nrow(safe_res@result), " pathways enriched")
      } else {
        message("⚠️  No results for ", co)
      }

      results[[co]] <- safe_res

      # Ensure all devices are closed after each iteration
      close_all_devices()
    }

    return(results)
  }
)

# -------------------------------------------------------------------- #
# 6.  Summary & Output                                                 #
# -------------------------------------------------------------------- #

message("\n", paste(rep("=", 78), collapse = ""))
message("✓ MitoCarta GSEA Analysis Complete")
message(paste(rep("=", 78), collapse = ""))

# Count successful results
null_checks <- sapply(mitocarta_gsea_results, function(x) !is.null(x))
n_success <- if(is.logical(null_checks)) sum(null_checks) else 0
message("\nResults summary:")
message("  Total contrasts processed: ", length(mitocarta_gsea_results))
message("  Successful analyses: ", n_success)
message("  Failed analyses: ", length(mitocarta_gsea_results) - n_success)

# Show pathway counts for each contrast
message("\nEnriched pathways per contrast:")
for (co in names(mitocarta_gsea_results)) {
  res <- mitocarta_gsea_results[[co]]
  if (!is.null(res)) {
    n_pathways <- nrow(res@result)
    n_sig <- sum(res@result$p.adjust < 0.05)
    message(sprintf("  %-30s: %3d total, %3d significant (padj < 0.05)",
                    co, n_pathways, n_sig))
  } else {
    message(sprintf("  %-30s: FAILED", co))
  }
}

message("\n📁 Output locations:")
message("  Checkpoint: ", file.path(checkpoint_dir, "mitocarta_gsea_results.rds"))
message("  Plots: ", here::here(config$out_root, "Plots/GSEA/[contrast]/MitoCarta/"))

message("\n✨ Done!")
