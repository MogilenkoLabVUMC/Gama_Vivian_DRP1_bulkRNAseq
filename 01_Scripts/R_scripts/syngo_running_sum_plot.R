#' SynGO Running Sum Plot (DEPRECATED)
#'
#' This function is deprecated and maintained only for backward compatibility.
#' Use `gsea_running_sum_plot()` from the RNAseq-toolkit instead, which now
#' works correctly for all databases including MSigDB, SynGO, and MitoCarta.
#'
#' @param gsea_obj SynGO gseaResult object
#' @param gene_set_ids numeric or character; 1-5 recommended
#' @param base_size Base font size for the theme
#' @param max_name_length Maximum length for pathway names in legend
#'
#' @return A patchwork object ready for saving
#'
#' @note DEPRECATED: This function now wraps gsea_running_sum_plot().
#'       The unified function was created on 2025-12-02 to consolidate
#'       all running sum plot implementations.
#'
#' @export
syngo_running_sum_plot <- function(gsea_obj,
                                   gene_set_ids = 1:5,
                                   base_size = 14,
                                   max_name_length = 40) {

  # Issue deprecation warning (once per session)
  .Deprecated(
    new = "gsea_running_sum_plot",
    package = "RNAseq-toolkit",
    msg = paste0(
      "syngo_running_sum_plot() is deprecated.\n",
      "Use gsea_running_sum_plot() which now works for all databases ",
      "(MSigDB, SynGO, MitoCarta).\n",
      "The function is located at: 01_Scripts/RNAseq-toolkit/scripts/GSEA/",
      "GSEA_plotting/gsea_running_sum_plot.R"
    )
  )

  # Check if the unified function is available
 if (!exists("gsea_running_sum_plot", mode = "function")) {
    # Try to source it
    toolkit_path <- here::here(
      "01_Scripts/RNAseq-toolkit/scripts/GSEA/GSEA_plotting/gsea_running_sum_plot.R"
    )
    if (file.exists(toolkit_path)) {
      source(toolkit_path)
    } else {
      stop("gsea_running_sum_plot() not found. Please source the RNAseq-toolkit.")
    }
  }

  # Delegate to the unified function
  gsea_running_sum_plot(
    gsea_obj = gsea_obj,
    gene_set_ids = gene_set_ids,
    base_size = base_size,
    max_name_length = max_name_length,
    title = "SynGO Gene Set Enrichment"
  )
}
