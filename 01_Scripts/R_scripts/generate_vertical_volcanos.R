#' Generate Specific Sets of Vertical Volcano Plots
#'
#' Creates and saves predefined combinations of vertical volcano plots in both p-value and FDR modes,
#' with options to include or exclude calcium gene highlighting.
#'
#' @param contrast_tables List of data frames with DE results for each contrast.
#' @param config Configuration list with plotting parameters.
#' @param base_dir Base directory for saving output files.
#' @param highlight_calcium Logical. Whether to create a second set with calcium genes highlighted.
#'
#' @return NULL (saves plots to files)
#' @export
generate_vertical_volcano_sets <- function(contrast_tables, config, base_dir = NULL, highlight_calcium = TRUE) {
  if (is.null(base_dir)) {
    base_dir <- here::here(config$out_root, "Plots/Volcano")
  }
  
  # Create directories for both modes
  p_dir <- file.path(base_dir, "vertical_p")
  fdr_dir <- file.path(base_dir, "vertical_fdr")
  dir.create(p_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(fdr_dir, recursive = TRUE, showWarnings = FALSE)
  
  if (highlight_calcium) {
    # Create additional directories for calcium-highlighted versions
    p_dir_calcium <- file.path(base_dir, "vertical_p_calcium")
    fdr_dir_calcium <- file.path(base_dir, "vertical_fdr_calcium")
    dir.create(p_dir_calcium, recursive = TRUE, showWarnings = FALSE)
    dir.create(fdr_dir_calcium, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Define the contrast groups to generate
  contrast_groups <- list(
    group1 = c("G32A_vs_Ctrl_D35", "R403C_vs_Ctrl_D35"),
    group2 = c("G32A_vs_Ctrl_D65", "R403C_vs_Ctrl_D65"),
    group3 = c("G32A_vs_Ctrl_D35", "R403C_vs_Ctrl_D35", "G32A_vs_Ctrl_D65", "R403C_vs_Ctrl_D65"),
    group4 = c("Time_Ctrl", "Time_G32A", "Time_R403C"),
    group5 = c("Maturation_G32A_specific", "Maturation_R403C_specific")
  )
  
  # Names for the output files
  group_names <- c(
    "D35_comparisons",
    "D65_comparisons",
    "all_disease_vs_control",
    "time_effects",
    "maturation_effects"
  )
  
  # Function to create vertical volcanoes in specified mode, with or without calcium highlighting
  create_volcanos_by_mode <- function(mode, output_dir, highlight_genes = NULL) {
    # Set parameters based on mode
    decision_by <- if (mode == "p") "p" else "fdr"
    p_cutoff <- if (mode == "p") config$p_cutoff else 0.1
    
    # Process each group
    for (i in seq_along(contrast_groups)) {
      group <- contrast_groups[[i]]
      group_name <- group_names[i]
      
      # Check which contrasts are available
      available_contrasts <- names(contrast_tables)
      valid_contrasts <- group[group %in% available_contrasts]
      
      if (length(valid_contrasts) == 0) {
        warning("No matching contrasts found for group: ", group_name)
        next
      }
      
      # Create a list of volcano plots for this group
      volcano_list <- list()
      for (contrast in valid_contrasts) {
        volcano_list[[contrast]] <- create_vertical_volcano(
          contrast_tables[[contrast]],
          decision_by = decision_by,
          p_cutoff = p_cutoff,
          fc_cutoff = config$fc_cutoff,
          label_method = "top",
          highlight_gene = highlight_genes,
          title = contrast
        )
      }
      
      # Combine into one panel and save
      combined <- combine_volcano_row(volcano_list, keep_first_caption = FALSE)
      
      # Calculate appropriate width based on number of plots
      width <- 3 * length(volcano_list)
      height <- 6
      
      # Save the combined plot
      filename <- file.path(output_dir, paste0(group_name, "_vertical.pdf"))
      pdf(filename, width = width, height = height)
      print(combined)
      dev.off()
      
      # Save individual plots too
      for (contrast in names(volcano_list)) {
        pdf(file.path(output_dir, paste0(contrast, "_vertical.pdf")), 6, 7)
        print(volcano_list[[contrast]])
        dev.off()
      }
    }
  }
  
  # Version without calcium gene highlighting
  create_volcanos_by_mode("p", p_dir, highlight_genes = NULL)
  create_volcanos_by_mode("fdr", fdr_dir, highlight_genes = NULL)
  
  # Create "all contrasts" version without calcium highlighting
  create_all_contrasts <- function(mode, output_dir, highlight_genes = NULL) {
    decision_by <- if (mode == "p") "p" else "fdr"
    p_cutoff <- if (mode == "p") config$p_cutoff else 0.1
    
    volcano_list <- list()
    for (contrast in names(contrast_tables)) {
      volcano_list[[contrast]] <- create_vertical_volcano(
        contrast_tables[[contrast]],
        decision_by = decision_by,
        p_cutoff = p_cutoff,
        fc_cutoff = config$fc_cutoff,
        label_method = "top",
        highlight_gene = highlight_genes,
        title = contrast
      )
    }
    
    width <- 3 * length(volcano_list)
    height <- 6
    
    pdf(file.path(output_dir, "all_contrasts_vertical.pdf"), width = width, height = height)
    print(combine_volcano_row(volcano_list, keep_first_caption = FALSE))
    dev.off()
  }
  
  create_all_contrasts("p", p_dir, highlight_genes = NULL)
  create_all_contrasts("fdr", fdr_dir, highlight_genes = NULL)
  
  # Optional: Create version with calcium gene highlighting
  if (highlight_calcium) {
    message("Generating additional volcanoes with calcium genes highlighted...")
    create_volcanos_by_mode("p", p_dir_calcium, highlight_genes = config$calcium_genes)
    create_volcanos_by_mode("fdr", fdr_dir_calcium, highlight_genes = config$calcium_genes)
    
    create_all_contrasts("p", p_dir_calcium, highlight_genes = config$calcium_genes)
    create_all_contrasts("fdr", fdr_dir_calcium, highlight_genes = config$calcium_genes)
  }
  
  message("Vertical volcano plots generated in p-value and FDR modes.")
}
