###############################################################################
##  GSEA Dotplot Helper Functions                                            ##
##                                                                             ##
##  Reusable functions for creating pooled GSEA dotplots with:               ##
##    - Gene Ratio calculation (leading edge / set size)                      ##
##    - Neuronal pathway filtering (iPSC cortical neuron-specific)            ##
##    - Significance outlines (black border for padj < 0.05)                  ##
##    - Colorblind-safe NES gradient (from color_config.R)                    ##
##                                                                             ##
##  Usage: source(here::here("01_Scripts/R_scripts/gsea_dotplot_helpers.R"))  ##
###############################################################################

library(ggplot2)
library(dplyr)
library(stringr)

# Ensure color config is loaded
if (!exists("DIVERGING_COLORS")) {
  source(here::here("01_Scripts/R_scripts/color_config.R"))
}

# ============================================================================
# GENE RATIO CALCULATION
# ============================================================================

#' Calculate Gene Ratio from Core Enrichment
#'
#' Computes the ratio of leading edge genes to total pathway genes.
#' Gene ratio = (number of core enrichment genes) / setSize
#'
#' @param core_enrichment Character vector of gene lists (genes separated by "/")
#' @param setSize Numeric vector of pathway set sizes
#' @return Numeric vector of gene ratios (0 to 1)
#'
#' @examples
#' gene_ratio <- calculate_gene_ratio("GENE1/GENE2/GENE3", 100)  # Returns 0.03
calculate_gene_ratio <- function(core_enrichment, setSize) {
  # Count genes in core_enrichment (genes separated by "/")
  count <- str_count(core_enrichment, "/") +
           ifelse(nchar(core_enrichment) > 0, 1, 0)
  return(count / setSize)
}

#' Add Gene Ratio Column to GSEA Results
#'
#' Wrapper function to add gene_ratio column to a GSEA results data frame.
#'
#' @param gsea_df Data frame with core_enrichment and setSize columns
#' @return Data frame with added GeneRatio column
add_gene_ratio <- function(gsea_df) {
  gsea_df$GeneRatio <- calculate_gene_ratio(
    gsea_df$core_enrichment,
    gsea_df$setSize
  )
  return(gsea_df)
}

# ============================================================================
# NEURONAL PATHWAY FILTERING
# ============================================================================

#' Filter Pathways for Neuronal Relevance
#'
#' Filters GSEA pathway results to keep only neuronal/brain-relevant pathways.
#' Removes tissue-specific pathways (cardiac, pancreatic, etc.) that are not
#' relevant to iPSC-derived cortical neurons.
#'
#' @param pathway_df Data frame with Description, ID, and Database columns
#' @return Filtered data frame with only neuronal-relevant pathways
#'
#' @details
#' Filtering rules:
#' - EXCLUDE: Pancreatic, cardiac, kidney, intestinal, lung, liver, breast,
#'            prostate, adipose, bile acid, xenobiotic, coagulation pathways
#' - INCLUDE (override): Neuronal, synaptic, axon, dendrite, brain, glial,
#'            ribosome, mitochondrial, cell cycle, apoptosis pathways
#' - ALWAYS KEEP: All SynGO pathways (100% neuronal-relevant)
filter_neuronal_pathways <- function(pathway_df) {

  if (nrow(pathway_df) == 0) return(pathway_df)

  ## EXCLUSION patterns (non-neuronal tissue-specific)
  exclude_patterns <- c(
    # Tissue-specific (non-neural)
    "PANCREA", "BETA_CELL", "ISLET",
    "CARDIAC", "HEART", "MYOCARDI", "CARDIOMYOCYTE",
    "KIDNEY", "RENAL", "NEPHRON",
    "INTESTIN", "COLON", "GUT", "DIGESTIVE",
    "LUNG", "ALVEOL", "BRONCH",
    "LIVER", "HEPAT",
    "BREAST", "MAMMARY",
    "PROSTAT",
    "MYOGENES",
    "ADIPOGEN", "ADIPOCYTE",
    # Metabolic (non-neuronal specific)
    "BILE_ACID",
    "XENOBIOTIC",
    "HEME_METABOL",
    "ANDROGEN_RESPONSE",
    # Blood/vascular
    "COAGULATION",
    "BLOOD_CLOT",
    # Other
    "UV_RESPONSE"
  )

  ## INCLUSION patterns (neuronal-relevant - override exclusions)
  include_patterns <- c(
    # Neuronal-specific
    "NEURO", "NEURAL", "NEURON", "NERVE",
    "SYNAP", "SYNAPT",
    "AXON", "DENDRIT",
    "BRAIN", "CORTEX", "HIPPOCAM", "CEREBR",
    "GLIAL", "ASTROCYT", "OLIGODENDRO",
    "MUSCLE.*NEURO", "NEUROMUSC",
    # Neuronal structures/functions
    "RIBOSOM", "TRANSLATION", "PROTEIN_SYNTHES",
    "MITOCHOND", "OXIDATIVE_PHOSPHORYL",
    # Relevant cell biology
    "E2F_TARGET", "G2M_CHECKPOINT",
    "MYC_TARGET",
    "APOPTOSIS", "P53_PATHWAY",
    "HYPOXIA",
    "DNA_REPAIR",
    "TNFA_SIGNALING",
    "INTERFERON", "INFLAMMAT",
    "EPITHELIAL_MESENCHYMAL",
    "MITOTIC_SPINDLE",
    "PROTEIN_SECRETION",
    "IL2_STAT5",
    "GLYCOLYSIS"
  )

  ## SynGO database: always keep (100% neuronal-relevant)
  syngo_mask <- grepl("syngo", pathway_df$Database, ignore.case = TRUE)

  ## Apply exclusions
  exclude_regex <- paste(exclude_patterns, collapse = "|")
  exclude_mask <- grepl(exclude_regex, pathway_df$Description, ignore.case = TRUE) |
                  grepl(exclude_regex, pathway_df$ID, ignore.case = TRUE)

  ## Apply inclusions (override exclusions)
  include_regex <- paste(include_patterns, collapse = "|")
  include_mask <- grepl(include_regex, pathway_df$Description, ignore.case = TRUE) |
                  grepl(include_regex, pathway_df$ID, ignore.case = TRUE)

  ## Final filter: (not excluded) OR (included) OR (syngo)
  keep_mask <- (!exclude_mask) | include_mask | syngo_mask

  filtered_df <- pathway_df[keep_mask, ]

  ## Report filtering
  n_removed <- nrow(pathway_df) - nrow(filtered_df)
  if (n_removed > 0) {
    message(sprintf("    Filtered: %d pathways removed (%d -> %d)",
                    n_removed, nrow(pathway_df), nrow(filtered_df)))
  }

  return(filtered_df)
}

# ============================================================================
# POOLED DOTPLOT CREATION
# ============================================================================

#' Create Pooled GSEA Dotplot with Significance Outlines
#'
#' Creates a publication-ready dotplot showing:
#' - X-axis: Gene Ratio (leading edge / set size)
#' - Y-axis: Pathway descriptions (sorted by gene ratio)
#' - Dot color: NES value (Blue-White-Orange gradient)
#' - Dot size: -log10(p.adjust)
#' - Black outline: Significant pathways (p.adjust < sig_cutoff)
#'
#' @param pathway_data Data frame with GSEA results
#' @param title Plot title
#' @param subtitle Plot subtitle
#' @param n_per_database Number of top pathways per database to show
#' @param highlight_sig Whether to add black outline for significant pathways
#' @param sig_cutoff P-value threshold for significance outline (default 0.05)
#' @param nes_limit Maximum |NES| for color scale (auto-calculated if NULL)
#' @return ggplot object
#'
#' @examples
#' p <- create_gsea_pooled_dotplot(
#'   filtered_data,
#'   title = "GSEA Pooled Dotplot",
#'   subtitle = "Top pathways by database",
#'   highlight_sig = TRUE
#' )
create_gsea_pooled_dotplot <- function(
    pathway_data,
    title = "GSEA Pooled Dotplot",
    subtitle = NULL,
    n_per_database = 5,
    highlight_sig = TRUE,
    sig_cutoff = 0.05,
    nes_limit = NULL
) {

  if (nrow(pathway_data) == 0) {
    return(ggplot() + labs(title = paste(title, "(No data)")))
  }

  # Calculate Gene Ratio if not present
  if (!"GeneRatio" %in% colnames(pathway_data)) {
    pathway_data <- add_gene_ratio(pathway_data)
  }

  # Select top N pathways per database
  top_pathways <- pathway_data %>%
    group_by(Database) %>%
    arrange(p.adjust) %>%
    slice_head(n = n_per_database) %>%
    ungroup()

  if (nrow(top_pathways) == 0) {
    return(ggplot() + labs(title = paste(title, "(No pathways after filtering)")))
  }

  # Calculate NES limits for color scale
  if (is.null(nes_limit)) {
    nes_max <- max(abs(top_pathways$NES), na.rm = TRUE)
    nes_limit <- ceiling(min(nes_max * 1.1, 4))  # Cap at 4
  }

  # Create base dotplot
  p <- ggplot(top_pathways,
              aes(x = GeneRatio, y = reorder(Description, GeneRatio))) +
    geom_point(
      aes(size = -log10(p.adjust), color = NES),
      alpha = 0.85
    ) +
    # Colorblind-safe NES gradient from color_config.R
    scale_color_gradient2(
      low = DIVERGING_COLORS$negative,   # Blue (#2166AC)
      mid = DIVERGING_COLORS$neutral,    # White (#F7F7F7)
      high = DIVERGING_COLORS$positive,  # Orange (#B35806)
      midpoint = 0,
      name = "NES",
      limits = c(-nes_limit, nes_limit),
      oob = scales::squish
    ) +
    scale_size_continuous(
      name = expression(-log[10](FDR)),
      range = c(3, 10)
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Gene Ratio (Leading Edge / Set Size)",
      y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.y = element_text(size = 9),
      legend.position = "right",
      panel.grid.major.y = element_line(color = "gray90"),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 9, color = "gray40")
    )

  # Add significance outlines
  if (highlight_sig) {
    sig_data <- top_pathways[top_pathways$p.adjust < sig_cutoff, ]

    if (nrow(sig_data) > 0) {
      p <- p + geom_point(
        data = sig_data,
        aes(x = GeneRatio, y = reorder(Description, GeneRatio),
            size = -log10(p.adjust)),
        shape = 21,
        color = "black",
        fill = NA,
        stroke = 1.2,
        show.legend = FALSE
      )
    }
  }

  return(p)
}

# ============================================================================
# CONFIRMATION MESSAGE
# ============================================================================
message("Loaded GSEA dotplot helpers from 01_Scripts/R_scripts/gsea_dotplot_helpers.R")
message("  Functions available: calculate_gene_ratio(), add_gene_ratio(),")
message("                       filter_neuronal_pathways(), create_gsea_pooled_dotplot()")
