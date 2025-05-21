#' Run SynGO GSEA Analysis and Visualization
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
#' @param output_dir Directory to save plots (default: file.path(gsea_root, contrast, "SynGO"))
#'
#' @return GSEA result object (invisible)

run_syngo_gsea <- function(tbl, contrast, T2G, T2N,
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
    out_dir <- file.path(gsea_root, contrast, "SynGO")
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
  saveRDS(res, file = file.path(out_dir, "GSEA_SynGO_result.rds"))
  
  # Skip plotting if save_plots is FALSE
  if (!save_plots) {
    return(invisible(res))
  }
  
  # Get plot parameters for SynGO
  plot_par <- get_db_plot_params("syngo")  # falls back to default sizes

  # Generate and save all standard plots
  # 1. Combined dotplot
  save_gsea_plot(
      gsea_dotplot(res, 
                  padj_cutoff = padj_cutoff,
                  showCategory = n_pathways,
                  title = paste(contrast, "SynGO")),
      filename = sprintf("%s_SynGO_dot.pdf", contrast),
      width = plot_par$width,
      height = plot_par$height,
      base_font_size = plot_par$font_size,
      dir = out_dir)

  # 2. NES barplot
  save_gsea_plot(
      gsea_barplot(res, 
                  padj_cutoff = padj_cutoff,
                  top_n = n_pathways,
                  title = paste(contrast, "SynGO NES")),
      filename = sprintf("%s_SynGO_bar.pdf", contrast),
      width = plot_par$width,
      height = plot_par$height,
      base_font_size = plot_par$font_size,
      dir = out_dir)

  # 3. Up/down faceted dotplot
  save_gsea_plot(
      gsea_dotplot_facet(res, 
                        padj_cutoff = padj_cutoff,
                        showCategory = n_pathways,
                        title = paste(contrast, "SynGO up / down")),
      filename = sprintf("%s_SynGO_facet.pdf", contrast),
      width = plot_par$width,
      height = plot_par$height * 1.4,
      base_font_size = plot_par$font_size,
      dir = out_dir)

 # Modify this section in run_syngo_gsea.R:

# In run_syngo_gsea.R, replace the running sum plot section with:

# 4. Running sum plot for top pathways
tryCatch({
  # Source the specialized SynGO running sum plot function
  if (!exists("syngo_running_sum_plot")) {
    source(file.path(dirname(getwd()), "01_Scripts/R_scripts/syngo_running_sum_plot.R"))
  }
  
  # Find top pathways by absolute NES value
  if (nrow(res@result) > 0) {
    # Get top 5 pathways by absolute NES value
    top5 <- order(abs(res@result$NES), decreasing = TRUE)[1:min(5, nrow(res@result))]
    
    # Generate and save the plot
    save_gsea_plot(
      syngo_running_sum_plot(
        gsea_obj = res, 
        gene_set_ids = top5,
        base_size = plot_par$font_size
      ),
      filename = sprintf("%s_SynGO_running_sum.pdf", contrast),
      width = plot_par$width,
      height = plot_par$height * 1.2,
      base_font_size = plot_par$font_size,
      dir = out_dir
    )
  } else {
    message("No pathways found in GSEA results for running sum plot")
  }
}, error = function(e) {
  warning("Error generating SynGO running sum plot: ", e$message)
})


  
  # 5. Split into up and down dotplots (like in the GSEA module)
  save_gsea_plot(
      gsea_dotplot(res,
                  filterBy = "NES_positive",
                  showCategory = n_pathways,
                  padj_cutoff = padj_cutoff,
                  title = sprintf("%s SynGO Up", contrast)),
      filename = sprintf("%s_SynGO_up_dot.pdf", contrast),
      width = plot_par$width,
      height = plot_par$height,
      base_font_size = plot_par$font_size,
      dir = out_dir)

  save_gsea_plot(
      gsea_dotplot(res,
                  filterBy = "NES_negative",
                  showCategory = n_pathways,
                  padj_cutoff = padj_cutoff,
                  title = sprintf("%s SynGO Down", contrast)),
      filename = sprintf("%s_SynGO_down_dot.pdf", contrast),
      width = plot_par$width,
      height = plot_par$height,
      base_font_size = plot_par$font_size,
      dir = out_dir)
  
#   # 6. Generate heatmap if sample annotations are provided
# if (!is.null(sample_annotation)) {
#   tryCatch({
#     top_tbl <- res@result |>
#       dplyr::filter(p.adjust < as.numeric(padj_cutoff))
    
#     # Check if we have any significant pathways
#     if (nrow(top_tbl) > 0) {
#       # Sort by p.adjust and take top n_pathways
#       top_tbl <- top_tbl |>
#         dplyr::arrange(p.adjust) |>
#         utils::head(n_pathways)
      
#       geneset <- top_tbl$ID
#       nes_vec <- top_tbl$NES
#       names(nes_vec) <- geneset
      
#       # Prepare matrix for heatmap
#       sample_ids <- rownames(sample_annotation)
      
#       # Ensure sample_ids is not NULL or empty
#       if (is.null(sample_ids) || length(sample_ids) == 0) {
#         warning("Sample annotation provided, but no row names found. Using numeric indices.")
#         sample_ids <- seq_len(nrow(sample_annotation))
#         rownames(sample_annotation) <- sample_ids
#       }
      
#       # Create matrix with NES values
#       mat <- matrix(rep(nes_vec, each = length(sample_ids)),
#                   nrow = length(sample_ids), byrow = TRUE,
#                   dimnames = list(sample_ids, geneset))
      
#       # Reorder samples if requested
#       if (!is.null(sample_order)) {
#         if (all(sample_order %in% rownames(mat))) {
#           mat <- mat[sample_order, , drop = FALSE]
#         } else {
#           warning("Sample order contains names not found in the sample annotation. Using original order.")
#         }
#       }
      
#       # Ensure annotation_col is properly formatted
#       ann_data <- sample_annotation[rownames(mat), , drop = FALSE]
      
#       # Set default colors for annotations if none provided
#       ann_col_list <- NULL
#       if (ncol(ann_data) > 0) {
#         ann_col_list <- list()
#         for (col_name in colnames(ann_data)) {
#           if (is.factor(ann_data[[col_name]]) || is.character(ann_data[[col_name]])) {
#             n_levels <- length(unique(ann_data[[col_name]]))
#             if (n_levels > 0) {
#               # Generate colors based on number of levels
#               col_palette <- RColorBrewer::brewer.pal(
#                 n = min(9, max(3, n_levels)),
#                 name = "Set2"
#               )
#               if (n_levels > 9) {
#                 col_palette <- colorRampPalette(col_palette)(n_levels)
#               }
              
#               # Create named vector
#               level_colors <- col_palette[1:n_levels]
#               names(level_colors) <- unique(ann_data[[col_name]])
#               ann_col_list[[col_name]] <- level_colors
#             }
#           }
#         }
#       }
      
#       # Generate and save heatmap
#       gsea_heatmap_save(mat,
#                        file = file.path(out_dir, sprintf("%s_SynGO_heatmap.pdf", contrast)),
#                        annotation_col = ann_data,
#                        ann_colors = ann_col_list,
#                        main = sprintf("%s - SynGO", contrast),
#                        gaps_col = NULL,
#                        cluster_cols = FALSE)
#     } else {
#       message("No significant pathways found for SynGO with padj_cutoff = ", padj_cutoff)
#     }
#   }, error = function(e) {
#     warning("Error generating heatmap for SynGO: ", e$message)
#   })
# }
  invisible(res)
}

# Helper function for SynGO GMT preparation
syngo_gmt <- function(syngo_dir, namespace = "CC") {
  requireNamespace("readxl", quietly = TRUE)
  
  # Read xlsx files
  ann <- readxl::read_xlsx(file.path(syngo_dir, "syngo_annotations.xlsx"))
  ont <- readxl::read_xlsx(file.path(syngo_dir, "syngo_ontologies.xlsx"))
  
  ann <- subset(ann, go_domain == namespace)
  
  ## ── TERM2GENE (required) ───────────────────────────────────────────
  term2gene <- ann[, c("go_id", "hgnc_symbol")]
  colnames(term2gene) <- c("gs_name", "gene_symbol")
  
  ## ── TERM2NAME (optional but recommended) ───────────────────────────
  term2name <- unique(ont[, c("id", "name")])
  colnames(term2name) <- c("gs_name", "description")
  
  ## strip trailing " (GO:########)" → cleaner axis labels
  term2name$description <- sub(" \\(GO:[0-9]+\\)$", "", term2name$description)
  
  list(T2G = term2gene, T2N = term2name)
}
