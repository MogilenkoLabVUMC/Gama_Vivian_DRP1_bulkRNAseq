###############################################################################
##  Standalone SynGO GSEA Analysis Runner
##  Loads checkpoints and runs SynGO analysis independently
###############################################################################

message(strrep("=", 60))
message("Starting standalone SynGO GSEA analysis...")
message(strrep("=", 60))

# -------------------------------------------------------------------- #
# 0.  Configuration (must match Analysis_pipeline.R)
# -------------------------------------------------------------------- #
config <- list(
  out_root      = "03_Results/02_Analysis",
  helper_root   = "01_Scripts/RNAseq-toolkit",
  syngo_dir     = "00_Data/SynGO_bulk_20231201",
  syngo_ns      = "CC",  # GO cellular-component ontology
  p_cutoff      = 0.05,
  fc_cutoff     = 2
)

checkpoint_dir <- here::here(config$out_root, "checkpoints")
gsea_root <- here::here(config$out_root, "Plots/GSEA")

# -------------------------------------------------------------------- #
# 1.  Load required packages
# -------------------------------------------------------------------- #
required_pkgs <- c(
  "limma", "dplyr", "ggplot2", "pheatmap", "RColorBrewer",
  "clusterProfiler", "fgsea", "here", "readxl"
)

for (p in required_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    message(sprintf("• installing %s …", p))
    install.packages(p, repos = "https://cloud.r-project.org")
  }
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

# -------------------------------------------------------------------- #
# 2.  Source helper functions
# -------------------------------------------------------------------- #
source_if_present <- function(...) {
  path <- here::here(...)
  if (file.exists(path)) {
    source(path, echo = FALSE)
    return(TRUE)
  } else {
    warning("helper not found → ", path)
    return(FALSE)
  }
}

# Source GSEA helpers
gsea_helpers <- c(
  "scripts/custom_minimal_theme.R",
  "scripts/GSEA/GSEA_plotting/gsea_plotting_utils.R",
  "scripts/GSEA/GSEA_plotting/format_pathway_names.R",
  "scripts/GSEA/GSEA_plotting/gsea_dotplot.R",
  "scripts/GSEA/GSEA_plotting/gsea_dotplot_facet.R",
  "scripts/GSEA/GSEA_plotting/gsea_barplot.R",
  "scripts/GSEA/GSEA_plotting/gsea_running_sum_plot.R",
  "scripts/GSEA/GSEA_plotting/gsea_heatmap.R",
  "scripts/GSEA/GSEA_processing/run_gsea.R"
)

message("\n📚 Sourcing helper functions...")
for (h in gsea_helpers) {
  source_if_present(config$helper_root, h)
}

# Source SynGO-specific functions
if (!source_if_present("01_Scripts/R_scripts/run_syngo_gsea.R")) {
  stop("ERROR: run_syngo_gsea.R not found!")
}

if (!source_if_present("01_Scripts/R_scripts/syngo_running_sum_plot.R")) {
  warning("syngo_running_sum_plot.R not found - running sum plots may fail")
}

message("✓ Helper functions loaded")

# -------------------------------------------------------------------- #
# 3.  Load checkpoints
# -------------------------------------------------------------------- #
message("\n📂 Loading checkpoints...")

# Check if checkpoint directory exists
if (!dir.exists(checkpoint_dir)) {
  stop("ERROR: Checkpoint directory not found: ", checkpoint_dir)
}

# Load fit object
fit_file <- file.path(checkpoint_dir, "fit_object.rds")
if (!file.exists(fit_file)) {
  stop("ERROR: fit_object.rds not found. Please run the main Analysis_pipeline.R first.")
}
fit <- readRDS(fit_file)
message("✓ Loaded fit_object.rds")

# Load contrasts
contrasts_file <- file.path(checkpoint_dir, "contrasts_matrix.rds")
if (!file.exists(contrasts_file)) {
  stop("ERROR: contrasts_matrix.rds not found. Please run the main Analysis_pipeline.R first.")
}
contrasts <- readRDS(contrasts_file)
message("✓ Loaded contrasts_matrix.rds (", ncol(contrasts), " contrasts)")

# Load QC variables (for annotation)
qc_file <- file.path(checkpoint_dir, "qc_variables.rds")
if (!file.exists(qc_file)) {
  warning("qc_variables.rds not found - sample annotation will not be available")
  annot <- NULL
} else {
  qc_vars <- readRDS(qc_file)
  annot <- qc_vars$annot
  message("✓ Loaded qc_variables.rds")
}

# -------------------------------------------------------------------- #
# 4.  Prepare SynGO GMT
# -------------------------------------------------------------------- #
message("\n🧬 Preparing SynGO gene sets...")

syngo_gmt <- function(syngo_dir, namespace = "CC") {
  requireNamespace("readxl", quietly = TRUE)

  # Check if files exist
  ann_file <- file.path(syngo_dir, "syngo_annotations.xlsx")
  ont_file <- file.path(syngo_dir, "syngo_ontologies.xlsx")

  if (!file.exists(ann_file)) {
    stop("ERROR: syngo_annotations.xlsx not found at: ", ann_file)
  }
  if (!file.exists(ont_file)) {
    stop("ERROR: syngo_ontologies.xlsx not found at: ", ont_file)
  }

  # Read xlsx files
  ann <- readxl::read_xlsx(ann_file)
  ont <- readxl::read_xlsx(ont_file)

  ann <- subset(ann, go_domain == namespace)

  if (nrow(ann) == 0) {
    stop("ERROR: No annotations found for namespace '", namespace, "'")
  }

  ## ── TERM2GENE (required) ───────────────────────────────────────────
  term2gene <- ann[, c("go_id", "hgnc_symbol")]
  colnames(term2gene) <- c("gs_name", "gene_symbol")

  ## ── TERM2NAME (optional but recommended) ───────────────────────────
  term2name <- unique(ont[, c("id", "name")])
  colnames(term2name) <- c("gs_name", "description")

  ## strip trailing " (GO:########)" → cleaner axis labels
  term2name$description <- sub(" \\(GO:[0-9]+\\)$", "", term2name$description)

  message("✓ Prepared ", nrow(term2gene), " gene-term mappings")
  message("✓ Prepared ", nrow(term2name), " term descriptions")

  list(T2G = term2gene, T2N = term2name)
}

# Execute GMT list preparation
tryCatch({
  syngo_lists <- syngo_gmt(here::here(config$syngo_dir), config$syngo_ns)
}, error = function(e) {
  stop("ERROR preparing SynGO GMT: ", e$message)
})

# -------------------------------------------------------------------- #
# 5.  Run SynGO GSEA for all contrasts
# -------------------------------------------------------------------- #
message("\n🔬 Running SynGO GSEA analysis...")
message("Will process ", ncol(contrasts), " contrasts")

# Create output directory
dir.create(gsea_root, recursive = TRUE, showWarnings = FALSE)

# Create a list to store all GSEA results
syngo_gsea_results <- list()

# Run GSEA for each contrast
for (co in colnames(contrasts)) {
  message("\n  → Processing: ", co)

  tryCatch({
    # Get DE results
    tbl <- topTable(fit, coef = co, number = Inf)
    message("    • Retrieved ", nrow(tbl), " genes from topTable")

    # Close any lingering graphic devices before starting new contrast
    while (dev.cur() > 1) dev.off()

    # Run SynGO GSEA
    result <- run_syngo_gsea(
      tbl,
      co,
      T2G = syngo_lists$T2G,
      T2N = syngo_lists$T2N,
      sample_annotation = annot,
      nperm = 100000,
      save_plots = TRUE
    )

    # Store result
    syngo_gsea_results[[co]] <- result

    # Report summary
    if (!is.null(result) && nrow(result@result) > 0) {
      n_sig <- sum(result@result$p.adjust < 0.05)
      message("    ✓ Found ", nrow(result@result), " pathways (",
              n_sig, " significant at FDR < 0.05)")
    } else {
      message("    ⚠ No pathways found")
    }

    # Ensure all devices are closed after each iteration
    while (dev.cur() > 1) dev.off()

  }, error = function(e) {
    message("    ✗ ERROR processing ", co, ": ", e$message)
    message("    Stack trace:")
    print(e)
    # Store NULL for failed contrasts
    syngo_gsea_results[[co]] <<- NULL
  })
}

# -------------------------------------------------------------------- #
# 6.  Save checkpoints
# -------------------------------------------------------------------- #
message("\n💾 Saving SynGO GSEA checkpoints...")

tryCatch({
  # Save SynGO GSEA results
  syngo_results_file <- file.path(checkpoint_dir, "syngo_gsea_results.rds")
  saveRDS(syngo_gsea_results, syngo_results_file)
  message("✓ Saved: syngo_gsea_results.rds (",
          format(file.size(syngo_results_file) / 1024^2, digits = 2), " MB)")

  # Save SynGO lists
  syngo_lists_file <- file.path(checkpoint_dir, "syngo_lists.rds")
  saveRDS(syngo_lists, syngo_lists_file)
  message("✓ Saved: syngo_lists.rds (",
          format(file.size(syngo_lists_file) / 1024^2, digits = 2), " MB)")

  message("\n", strrep("=", 60))
  message("✓ SynGO analysis complete!")
  message("  Checkpoints saved (", length(syngo_gsea_results), " contrasts)")
  message(strrep("=", 60))

}, error = function(e) {
  message("\n✗ ERROR saving checkpoints: ", e$message)
  stop(e)
})

# -------------------------------------------------------------------- #
# 7.  Summary report
# -------------------------------------------------------------------- #
message("\n📊 Summary Report:")
message(strrep("=", 60))

for (co in names(syngo_gsea_results)) {
  result <- syngo_gsea_results[[co]]
  if (!is.null(result) && nrow(result@result) > 0) {
    n_total <- nrow(result@result)
    n_sig_05 <- sum(result@result$p.adjust < 0.05)
    n_sig_01 <- sum(result@result$p.adjust < 0.01)
    message(sprintf("%-30s: %4d pathways (%d @ FDR<0.05, %d @ FDR<0.01)",
                    co, n_total, n_sig_05, n_sig_01))
  } else {
    message(sprintf("%-30s: FAILED or no results", co))
  }
}

message(strrep("=", 60))
message("\n✅ You can now run:")
message("  • viz_pooled_dotplots.R")
message("  • viz_ribosome_pathways.R")
message("  • Supp1.verify_enrichments.R")
