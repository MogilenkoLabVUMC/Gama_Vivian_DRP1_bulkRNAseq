#' SynGO-specific Running-sum plot
#'
#' A specialized version of running sum plot for SynGO GSEA results.
#' Handles both SYNGO: and GO: prefixed terms correctly.
#'
#' @param gsea_obj SynGO gseaResult object
#' @param gene_set_ids numeric or character; 1-5 recommended
#' @param base_size Base font size for the theme
#' @param max_name_length Maximum length for pathway names in legend
#' @return A patchwork object ready for saving
#' @export
syngo_running_sum_plot <- function(gsea_obj, 
                                  gene_set_ids = 1:5,
                                  base_size = 14,
                                  max_name_length = 40) {
  
  if (!requireNamespace("enrichplot", quietly = TRUE))
    stop("Package 'enrichplot' is required.")
  if (!requireNamespace("patchwork", quietly = TRUE))
    stop("Package 'patchwork' is required.")
  
  # Ensure gsea_obj is a valid object
  if (!methods::is(gsea_obj, "gseaResult") || is.null(gsea_obj@result))
    stop("Invalid GSEA result object")
  
  # Handle pathway IDs - either numeric indices or actual IDs
  if (is.numeric(gene_set_ids)) {
    # Check if we have enough results
    if (max(gene_set_ids) > nrow(gsea_obj@result))
      stop("Index out of bounds: requesting pathway index larger than available results")
    
    # Extract actual IDs from the result table
    path_ids <- gsea_obj@result$ID[gene_set_ids]
  } else {
    # User provided actual IDs
    path_ids <- gene_set_ids
    
    # Verify they exist in the results
    if (!all(path_ids %in% gsea_obj@result$ID))
      stop("Some specified pathway IDs not found in GSEA results")
  }
  
  # Get corresponding descriptions for better labels
  descriptions <- gsea_obj@result$Description[match(path_ids, gsea_obj@result$ID)]
  
  # Truncate long descriptions for readability in legend
  display_labels <- sapply(descriptions, function(x) {
    if (nchar(x) > max_name_length) {
      paste0(substr(x, 1, max_name_length-3), "...")
    } else {
      x
    }
  })
  names(display_labels) <- path_ids
  
  # Generate the running sum plots with specific parameters
  tryCatch({
    # Create a vibrant color palette that will be visible
    n_paths <- length(path_ids)
    colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999")
    
    # If we need more colors than the palette provides, interpolate
    if (n_paths > length(colors)) {
      colors <- grDevices::colorRampPalette(colors)(n_paths)
    } else {
      colors <- colors[1:n_paths]
    }
    
    # DO NOT assign names to the colors - this appears to be causing issues with enrichplot
    
    # Generate the plot using enrichplot directly with specified colors
    p_raw <- enrichplot::gseaplot2(gsea_obj,
                                  geneSetID = path_ids,
                                  title = "SynGO Gene Set Enrichment",
                                  subplots = c(1, 2, 3),
                                  pvalue_table = FALSE,
                                  rel_heights = c(1.5, 0.5, 0.5),
                                  color = colors)  # Pass colors directly
    
    # Define styling function without trying to override colors
    stylize <- function(p, show_x = TRUE, show_y = TRUE, show_legend = TRUE) {
      p + 
        # No scale_color_manual - this was causing the issue
        ggplot2::labs(color = NULL) +  # Remove legend title
        ggplot2::theme_classic(base_size = base_size) +
        ggplot2::theme(
          # Legend styling
          legend.position = if(show_legend) c(0.98, 0.98) else "none",
          legend.justification = c(1, 1),
          legend.background = ggplot2::element_rect(fill = "white", color = "grey90"),
          legend.margin = ggplot2::margin(5, 5, 5, 5),
          legend.key.size = ggplot2::unit(0.8, "lines"),
          
          # Panel and axis styling
          panel.background = ggplot2::element_rect(fill = "white", color = NA),
          axis.line = ggplot2::element_line(color = "black", size = 0.5),
          axis.ticks = ggplot2::element_line(color = "black", size = 0.5),
          
          # Text elements
          axis.title.x = if(show_x) ggplot2::element_text() else ggplot2::element_blank(),
          axis.title.y = if(show_y) ggplot2::element_text() else ggplot2::element_blank(),
          axis.text.x = if(show_x) ggplot2::element_text() else ggplot2::element_blank(),
          axis.text.y = if(show_y) ggplot2::element_text() else ggplot2::element_blank(),
          
          # Plot margins
          plot.margin = ggplot2::margin(5, 10, 5, 5)
        )
    }
    
    # Apply styling to each subplot without modifying the colors
    p1 <- stylize(p_raw[[1]], show_x = FALSE, show_y = TRUE, show_legend = TRUE)
    p2 <- stylize(p_raw[[2]], show_x = FALSE, show_y = FALSE, show_legend = FALSE)
    p3 <- stylize(p_raw[[3]], show_x = TRUE, show_y = TRUE, show_legend = FALSE)
    
    # Update the legend labels manually after styling
    if ("guides" %in% names(p1)) {
      p1 <- p1 + ggplot2::guides(
        color = ggplot2::guide_legend(
          override.aes = list(color = colors),
          labels = display_labels
        )
      )
    }
    
    # Combine plots with patchwork
    combined <- patchwork::wrap_plots(p1, p2, p3, 
                                    ncol = 1, 
                                    heights = c(2, 0.5, 0.5))
    
    return(combined)
    
  }, error = function(e) {
    stop("Error generating SynGO running sum plot: ", e$message)
  })
}
