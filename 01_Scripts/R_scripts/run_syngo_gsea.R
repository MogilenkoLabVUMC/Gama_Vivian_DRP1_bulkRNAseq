run_syngo_gsea <- function(tbl, contrast, T2G, T2N,
                           n_pathways   = 30,
                           padj_cutoff  = 0.05) {

  gene_vec <- sort(setNames(tbl$t, rownames(tbl)), decreasing = TRUE)

  set.seed(123)
  res <- clusterProfiler::GSEA(
           geneList       = gene_vec,
           TERM2GENE      = T2G,
           TERM2NAME      = T2N,
           pvalueCutoff   = 1,
           pAdjustMethod  = "fdr",
           eps            = 0,
           by             = "fgsea",
           nPermSimple    = 100000)

  out_dir <- file.path(gsea_root, contrast, "SynGO")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(res, file = file.path(out_dir, "GSEA_SynGO_result.rds"))

  ## ---------- standard plots (same helpers you already have) ----------
  plot_par <- get_db_plot_params("syngo")   # falls back to default sizes

  save_gsea_plot(
      gsea_dotplot(res, padj_cutoff = padj_cutoff,
                   showCategory = n_pathways,
                   title = paste(contrast, "SynGO")),
      filename = sprintf("%s_SynGO_dot.pdf", contrast),
      width  = plot_par$width,
      height = plot_par$height,
      base_font_size = plot_par$font_size,
      dir = out_dir)

  save_gsea_plot(
      gsea_barplot(res, padj_cutoff = padj_cutoff,
                   top_n = n_pathways,
                   title = paste(contrast, "SynGO NES")),
      filename = sprintf("%s_SynGO_bar.pdf", contrast),
      width  = plot_par$width,
      height = plot_par$height,
      base_font_size = plot_par$font_size,
      dir = out_dir)

  save_gsea_plot(
      gsea_dotplot_facet(res, padj_cutoff = padj_cutoff,
                         showCategory = n_pathways,
                         title = paste(contrast, "SynGO up / down")),
      filename = sprintf("%s_SynGO_facet.pdf", contrast),
      width  = plot_par$width,
      height = plot_par$height * 1.4,
      base_font_size = plot_par$font_size,
      dir = out_dir)

  top5 <- head(order(abs(res@result$NES), decreasing = TRUE), 5)
  save_gsea_plot(
      gsea_running_sum_plot(res, gene_set_ids = top5),
      filename = sprintf("%s_SynGO_running_sum.pdf", contrast),
      width  = plot_par$width,
      height = plot_par$height * 1.2,
      base_font_size = plot_par$font_size,
      dir = out_dir)

  invisible(res)
}