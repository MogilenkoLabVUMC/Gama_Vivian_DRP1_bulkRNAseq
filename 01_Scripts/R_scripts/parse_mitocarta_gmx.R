###############################################################################
##  parse_mitocarta_gmx.R – Parse MitoCarta GMX file for GSEA                ##
###############################################################################
##  Purpose: Convert GMX (Gene Matrix Transposed) format to TERM2GENE and
##           TERM2NAME structures required by clusterProfiler/fgsea
##
##  GMX Format:
##    Row 1: Short pathway names
##    Row 2: Full hierarchical pathway names (pathway.subpathway.subsubpathway)
##    Row 3+: Gene symbols (one gene per row, distributed across pathway columns)
##
##  Output: List with T2G (TERM2GENE) and T2N (TERM2NAME) data frames
###############################################################################

#' Parse MitoCarta GMX file
#'
#' @param gmx_file Path to the GMX file
#' @return List with two data frames:
#'   - T2G: TERM2GENE mapping (gs_name, gene_symbol)
#'   - T2N: TERM2NAME mapping (gs_name, description)
#'
mitocarta_gmt <- function(gmx_file) {

  # Validate input
  if (!file.exists(gmx_file)) {
    stop("GMX file not found: ", gmx_file)
  }

  message("📖 Reading MitoCarta GMX file: ", gmx_file)

  # Read GMX file (tab-delimited, no header)
  gmx_data <- read.delim(gmx_file,
                         header = FALSE,
                         stringsAsFactors = FALSE,
                         check.names = FALSE,
                         fill = TRUE,  # Handle ragged arrays
                         comment.char = "",
                         quote = "")

  message("   Dimensions: ", nrow(gmx_data), " rows × ", ncol(gmx_data), " columns")

  # Extract the three components
  short_names <- as.character(gmx_data[1, ])
  full_names  <- as.character(gmx_data[2, ])
  gene_data   <- gmx_data[3:nrow(gmx_data), ]

  # Remove any completely empty columns
  non_empty_cols <- !is.na(full_names) & full_names != ""
  short_names <- short_names[non_empty_cols]
  full_names  <- full_names[non_empty_cols]
  gene_data   <- gene_data[, non_empty_cols, drop = FALSE]

  message("   Found ", length(full_names), " pathways")

  # Convert to TERM2GENE format (long format: one row per pathway-gene pair)
  term2gene_list <- vector("list", ncol(gene_data))

  for (i in seq_len(ncol(gene_data))) {
    # Extract genes for this pathway
    genes <- as.character(gene_data[, i])

    # Remove empty cells and NA values
    genes <- genes[!is.na(genes) & genes != "" & nchar(trimws(genes)) > 0]

    if (length(genes) > 0) {
      term2gene_list[[i]] <- data.frame(
        gs_name = rep(full_names[i], length(genes)),
        gene_symbol = genes,
        stringsAsFactors = FALSE
      )
    }
  }

  # Combine all pathway-gene pairs
  term2gene <- do.call(rbind, term2gene_list[!sapply(term2gene_list, is.null)])

  # Remove any duplicates (shouldn't be any, but just in case)
  term2gene <- unique(term2gene)

  message("   TERM2GENE: ", nrow(term2gene), " pathway-gene mappings")

  # Create TERM2NAME: map full hierarchical name to short name
  term2name <- data.frame(
    gs_name = full_names,
    description = short_names,
    stringsAsFactors = FALSE
  )

  # Remove empty rows
  term2name <- term2name[!is.na(term2name$gs_name) &
                         term2name$gs_name != "", ]

  # Remove duplicates
  term2name <- unique(term2name)

  message("   TERM2NAME: ", nrow(term2name), " pathway descriptions")

  # Validate output
  if (nrow(term2gene) == 0) {
    warning("No pathway-gene mappings found in GMX file!")
  }
  if (nrow(term2name) == 0) {
    warning("No pathway names found in GMX file!")
  }

  # Check for missing pathway names
  missing_names <- setdiff(term2gene$gs_name, term2name$gs_name)
  if (length(missing_names) > 0) {
    warning("Some pathways in TERM2GENE are missing from TERM2NAME: ",
            length(missing_names), " pathways")
  }

  message("✓ MitoCarta GMX parsing complete")

  return(list(T2G = term2gene, T2N = term2name))
}
