#' Run MitoCarta GSEA Analysis and Visualization
#'
#' @param tbl Differential expression results table
#' @param contrast Name of the contrast for plot labeling
#' @param T2G Term to gene mapping dataframe
#' @param T2N Term to name mapping dataframe
#' @param n_pathways Number of top pathways to display (default: 30)
#' @param padj_cutoff Adjusted p-value cutoff (default: 0.05)
#' @param save_plots Logical, save generated plots (default: TRUE)
#' @param sample_annotation Optional dataframe with sample annotations for heatmap
#' @param sample_order Optional vector specifying sample order for heatmap
#' @param nperm Number of permutations for GSEA (default: 100000)
#' @param output_dir Directory to save plots (default: file.path(gsea_root, contrast, "MitoCarta"))
#'
#' @return GSEA result object (invisible)

run_mitocarta_gsea <- function(tbl, contrast, T2G, T2N,
                               n_pathways = 30,
                               padj_cutoff = 0.05,
                               save_plots = TRUE,
                               sample_annotation = NULL,
                               sample_order = NULL,
                               nperm = 100000,
                               output_dir = NULL) {

  # Create ranked gene list for GSEA
  gene_vec <- sort(setNames(tbl$t, rownames(tbl)), decreasing = TRUE)

  # Set output directory
  if (is.null(output_dir)) {
    out_dir <- file.path(gsea_root, contrast, "MitoCarta")
  } else {
    out_dir <- output_dir
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # Run GSEA
  set.seed(123)
  res <- clusterProfiler::GSEA(
           geneList = gene_vec,
           TERM2GENE = T2G,
           TERM2NAME = T2N,
           pvalueCutoff = 1,
           pAdjustMethod = "fdr",
           eps = 0,
           by = "fgsea",
           nPermSimple = nperm)

  # Save results
  saveRDS(res, file = file.path(out_dir, "GSEA_MitoCarta_result.rds"))

  # Skip plotting if save_plots is FALSE
  if (!save_plots) {
    return(invisible(res))
  }

  # Get plot parameters for MitoCarta
  plot_par <- get_db_plot_params("mitocarta")  # falls back to default sizes

  # Generate and save all standard plots
  # 1. Combined dotplot
  save_gsea_plot(
      gsea_dotplot(res,
                  padj_cutoff = padj_cutoff,
                  showCategory = n_pathways,
                  title = paste(contrast, "MitoCarta")),
      filename = sprintf("%s_MitoCarta_dot.pdf", contrast),
      width = plot_par$width,
      height = plot_par$height,
      base_font_size = plot_par$font_size,
      dir = out_dir)

  # 2. NES barplot
  save_gsea_plot(
      gsea_barplot(res,
                  padj_cutoff = padj_cutoff,
                  top_n = n_pathways,
                  title = paste(contrast, "MitoCarta NES")),
      filename = sprintf("%s_MitoCarta_bar.pdf", contrast),
      width = plot_par$width,
      height = plot_par$height,
      base_font_size = plot_par$font_size,
      dir = out_dir)

  # 3. Up/down faceted dotplot
  save_gsea_plot(
      gsea_dotplot_facet(res,
                        padj_cutoff = padj_cutoff,
                        showCategory = n_pathways,
                        title = paste(contrast, "MitoCarta up / down")),
      filename = sprintf("%s_MitoCarta_facet.pdf", contrast),
      width = plot_par$width,
      height = plot_par$height * 1.4,
      base_font_size = plot_par$font_size,
      dir = out_dir)

  # 4. Running sum plot for top pathways
  tryCatch({
    # Source the specialized running sum plot function (same as SynGO uses)
    if (!exists("syngo_running_sum_plot")) {
      source(here::here("01_Scripts/R_scripts/syngo_running_sum_plot.R"))
    }

    # Find top pathways by absolute NES value
    if (nrow(res@result) > 0) {
      # Get top 5 pathways by absolute NES value
      top5 <- order(abs(res@result$NES), decreasing = TRUE)[1:min(5, nrow(res@result))]

      # Generate and save the plot using the specialized function
      save_gsea_plot(
        syngo_running_sum_plot(
          gsea_obj = res,
          gene_set_ids = top5,
          base_size = plot_par$font_size
        ),
        filename = sprintf("%s_MitoCarta_running_sum.pdf", contrast),
        width = plot_par$width,
        height = plot_par$height * 1.2,
        base_font_size = plot_par$font_size,
        dir = out_dir
      )
    } else {
      message("No pathways found in GSEA results for running sum plot")
    }
  }, error = function(e) {
    warning("Error generating MitoCarta running sum plot: ", e$message)
  })

  # 5. Split into up and down dotplots
  save_gsea_plot(
      gsea_dotplot(res,
                  filterBy = "NES_positive",
                  showCategory = n_pathways,
                  padj_cutoff = padj_cutoff,
                  title = sprintf("%s MitoCarta Up", contrast)),
      filename = sprintf("%s_MitoCarta_up_dot.pdf", contrast),
      width = plot_par$width,
      height = plot_par$height,
      base_font_size = plot_par$font_size,
      dir = out_dir)

  save_gsea_plot(
      gsea_dotplot(res,
                  filterBy = "NES_negative",
                  showCategory = n_pathways,
                  padj_cutoff = padj_cutoff,
                  title = sprintf("%s MitoCarta Down", contrast)),
      filename = sprintf("%s_MitoCarta_down_dot.pdf", contrast),
      width = plot_par$width,
      height = plot_par$height,
      base_font_size = plot_par$font_size,
      dir = out_dir)

  invisible(res)
}
