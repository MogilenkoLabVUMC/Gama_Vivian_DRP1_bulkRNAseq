###############################################################################
##  Unified Color Configuration for DRP1 Bulk RNA-seq Analysis              ##
##                                                                           ##
##  This file provides a SINGLE SOURCE OF TRUTH for color schemes used      ##
##  across all R visualization scripts. It mirrors the Python color         ##
##  definitions in:                                                          ##
##    - 01_Scripts/Python/color_config.py (PRIMARY Python source)           ##
##    - 01_Scripts/RNAseq-toolkit/scripts/GSEA/GSEA_plotting_python/colormaps.py
##    - 01_Scripts/Python/pattern_definitions.py                             ##
##    - 01_Scripts/Python/semantic_categories.py                             ##
##                                                                           ##
##  Usage: source(here::here("01_Scripts/R_scripts/color_config.R"))        ##
##                                                                           ##
##  COLOR SYSTEM OVERVIEW:                                                   ##
##  ----------------------                                                   ##
##  1. DIVERGING_COLORS: Blue-White-Orange for NES/logFC heatmaps           ##
##     - AVOID using these colors for categorical annotations!              ##
##  2. HEATMAP_ANNOTATION_COLORS: Distinct colors for heatmap annotations   ##
##     - Use these for row/column annotations on diverging heatmaps         ##
##  3. MUTATION_COLORS: Standard mutation colors for general use            ##
##  4. Other categorical palettes for specific use cases                    ##
##                                                                           ##
##  COLORBLIND-SAFE: All palettes designed with deuteranopia/protanopia     ##
##  in mind. Avoid relying solely on red-green distinctions.                ##
###############################################################################

# =============================================================================
# DIVERGING COLOR PALETTE (Colorblind-safe Blue-White-Orange)
# =============================================================================
# This palette is used for NES values, logFC, and other diverging metrics
# where negative/positive values have distinct biological meaning.
#
# IMPORTANT: When using this gradient for heatmaps, use HEATMAP_ANNOTATION_COLORS
# for any row/column annotations to avoid color conflicts.

DIVERGING_COLORS <- list(
  negative = "#2166AC",  # Blue
  neutral  = "#F7F7F7",  # White
  positive = "#B35806"   # Orange/Brown
)

# Intermediate shades for 5-point gradients (logFC heatmaps)
DIVERGING_COLORS_5PT <- list(
  neg_strong = "#2166AC",
  neg_weak   = "#92c5de",
  neutral    = "#F7F7F7",
  pos_weak   = "#f4a582",
  pos_strong = "#B35806"
)

# =============================================================================
# NES/ENRICHMENT SCORE COLOR FUNCTIONS
# =============================================================================

#' Create NES color scale for ComplexHeatmap
#'
#' @param limits Numeric vector of length 2 for min/max NES values
#' @param n_colors Number of color steps (3 for simple, 5 for gradient)
#' @return colorRamp2 function for use in ComplexHeatmap
#' @examples
#' col_fun <- nes_color_scale(limits = c(-3.5, 3.5))
#' col_fun <- nes_color_scale(limits = c(-2, 2), n_colors = 5)
nes_color_scale <- function(limits = c(-3.5, 3.5), n_colors = 3) {
  if (!requireNamespace("circlize", quietly = TRUE)) {
    stop("Package 'circlize' required for colorRamp2. Install with: install.packages('circlize')")
  }

  if (n_colors == 3) {
    circlize::colorRamp2(
      c(limits[1], 0, limits[2]),
      c(DIVERGING_COLORS$negative, DIVERGING_COLORS$neutral, DIVERGING_COLORS$positive)
    )
  } else if (n_colors == 5) {
    mid <- (limits[1] + limits[2]) / 2
    q1 <- limits[1] / 2
    q3 <- limits[2] / 2
    circlize::colorRamp2(
      c(limits[1], q1, mid, q3, limits[2]),
      c(DIVERGING_COLORS_5PT$neg_strong, DIVERGING_COLORS_5PT$neg_weak,
        DIVERGING_COLORS_5PT$neutral,
        DIVERGING_COLORS_5PT$pos_weak, DIVERGING_COLORS_5PT$pos_strong)
    )
  } else {
    stop("n_colors must be 3 or 5")
  }
}

#' Create logFC color scale for ComplexHeatmap
#'
#' @param limits Numeric vector of length 2 for min/max logFC values
#' @return colorRamp2 function for use in ComplexHeatmap
#' @examples
#' col_fun <- logfc_color_scale(limits = c(-1.5, 1.5))
logfc_color_scale <- function(limits = c(-1.5, 1.5)) {
  if (!requireNamespace("circlize", quietly = TRUE)) {
    stop("Package 'circlize' required for colorRamp2")
  }

  # 5-point gradient for smoother logFC visualization
  q1 <- limits[1] / 2
  q3 <- limits[2] / 2
  circlize::colorRamp2(
    c(limits[1], q1, 0, q3, limits[2]),
    c(DIVERGING_COLORS_5PT$neg_strong, DIVERGING_COLORS_5PT$neg_weak,
      DIVERGING_COLORS_5PT$neutral,
      DIVERGING_COLORS_5PT$pos_weak, DIVERGING_COLORS_5PT$pos_strong)
  )
}

#' Get ggplot2 scale for NES values
#'
#' @param limits Numeric vector of length 2 for min/max values
#' @param name Legend title
#' @return ggplot2 scale_color_gradient2 object
#' @examples
#' p + nes_ggplot_scale(limits = c(-3.5, 3.5))
nes_ggplot_scale <- function(limits = c(-3.5, 3.5), name = "NES") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required")
  }

  ggplot2::scale_color_gradient2(
    low = DIVERGING_COLORS$negative,
    mid = DIVERGING_COLORS$neutral,
    high = DIVERGING_COLORS$positive,
    midpoint = 0,
    name = name,
    limits = limits,
    oob = scales::squish
  )
}

#' Get ggplot2 fill scale for NES values
#'
#' @param limits Numeric vector of length 2 for min/max values
#' @param name Legend title
#' @return ggplot2 scale_fill_gradient2 object
nes_ggplot_fill_scale <- function(limits = c(-3.5, 3.5), name = "NES") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required")
  }

  ggplot2::scale_fill_gradient2(
    low = DIVERGING_COLORS$negative,
    mid = DIVERGING_COLORS$neutral,
    high = DIVERGING_COLORS$positive,
    midpoint = 0,
    name = name,
    limits = limits,
    oob = scales::squish
  )
}

# =============================================================================
# HEATMAP ANNOTATION COLORS
# =============================================================================
# These colors are specifically designed for row/column annotations on heatmaps
# that use the Blue-White-Orange diverging gradient. They AVOID blue and orange
# to prevent visual confusion.
#
# Use these for: Mutation, Stage, Module annotations on ComplexHeatmap/pheatmap

HEATMAP_ANNOTATION_COLORS <- list(
  # Mutation colors for heatmap annotations (distinct from gradient blue/orange)
  mutation = c(
    "G32A"    = "#7B68EE",  # Medium Slate Blue (purple-ish, distinct from blue gradient)
    "R403C"   = "#DC143C",  # Crimson (red, distinct from orange gradient)
    "Ctrl"    = "#808080",  # Gray
    "Control" = "#808080"   # Gray (alias)
  ),

  # Stage/Trajectory colors (teal sequential - avoids red/pink confusion with mutations)
  stage = c(
    "Early"   = "#B2DFDB",  # Light teal
    "TrajDev" = "#4DB6AC",  # Medium teal
    "Late"    = "#00796B"   # Dark teal
  ),

  # Timepoint colors (for genotype x timepoint annotations)
  timepoint = c(
    "D35" = "#66C2A5",  # Teal/mint
    "D65" = "#8DA0CB"   # Periwinkle
  ),

  # Genotype colors (alias for mutation with full names)
  genotype = c(
    "Control" = "#808080",
    "G32A"    = "#7B68EE",
    "R403C"   = "#DC143C"
  )
)

# =============================================================================
# MUTATION COLORS (Standard - for bar charts, dotplots, non-heatmap figures)
# =============================================================================
# Colorblind-safe Okabe-Ito palette for mutation identity
# NOTE: For heatmap annotations, use HEATMAP_ANNOTATION_COLORS$mutation instead

MUTATION_COLORS <- c(
  "G32A"  = "#0072B2",  # Okabe-Ito Blue
  "R403C" = "#D55E00",  # Okabe-Ito Vermillion
  "Ctrl"  = "#999999"   # Gray
)

# =============================================================================
# PATTERN COLORS (Trajectory classification)
# =============================================================================
# These match the Python definitions in pattern_definitions.py

PATTERN_COLORS <- c(
  "Compensation"       = "#009E73",
  "Sign_reversal"      = "#9467BD",
  "Progressive"        = "#D55E00",
  "Natural_worsening"  = "#E69F00",
  "Natural_improvement"= "#56B4E9",
  "Late_onset"         = "#CC79A7",
  "Transient"          = "#0072B2",
  "Complex"            = "#F0E442",
  "Insufficient_data"  = "#DDDDDD"
)

# =============================================================================
# SUPER-CATEGORY COLORS (Simplified pattern groupings)
# =============================================================================

SUPER_CATEGORY_COLORS <- c(
  "Active_Compensation" = "#009E73",
  "Active_Reversal"     = "#9467BD",
  "Active_Progression"  = "#D55E00",
  "Passive"             = "#56B4E9",
  "Late_onset"          = "#CC79A7",
  "Other"               = "#999999",
  "Insufficient_data"   = "#DDDDDD"
)

# =============================================================================
# SEMANTIC CATEGORY COLORS (Biological pathway groupings)
# =============================================================================
# These match the Python definitions in semantic_categories.py

SEMANTIC_COLORS <- c(
  "Synapse"                    = "#8B4513",
  "Neuronal Development"       = "#A0522D",
  "Mitochondrial Dynamics"     = "#006400",
  "Electron Transport Chain"   = "#2E8B57",
  "ATP Synthase (Complex V)"   = "#DAA520",
  "Mitochondrial Metabolism"   = "#228B22",
  "Mitochondrial Function"     = "#32CD32",
  "Mitochondrial Ribosome"     = "#4169E1",
  "Mitochondrial Translation"  = "#1E90FF",
  "Ribosome Biogenesis"        = "#6A5ACD",
  "Cytoplasmic Ribosome"       = "#9932CC",
  "Cytoplasmic Translation"    = "#BA55D3",
  "Calcium Signaling"          = "#DC143C",
  "Other"                      = "#696969"
)

# =============================================================================
# TRAJECTORY STAGE COLORS (for standalone trajectory figures)
# =============================================================================
# NOTE: For heatmap annotations, use HEATMAP_ANNOTATION_COLORS$stage instead

TRAJECTORY_COLORS <- c(
  "Early"   = "#FEE0D2",  # Light pink/salmon
  "TrajDev" = "#FC9272",  # Medium salmon
  "Late"    = "#DE2D26"   # Dark red
)

# =============================================================================
# MODULE COLORS (for mechanistic cascade figures)
# =============================================================================
# Redesigned to avoid overlap with HEATMAP_ANNOTATION_COLORS
# Uses distinct hues: purple, gold, forest green, teal, maroon, coral

MODULE_COLORS <- c(
  "1. Mt Central Dogma"         = "#9467BD",  # Purple (distinct from teal stages)
  "2. Mt Ribosomes"             = "#DAA520",  # Goldenrod (distinct from orange gradient)
  "3. ATP Synthase (Complex V)" = "#2E8B57",  # Sea green (distinct)
  "4. Synaptic Ribo (Common)"   = "#17BECF",  # Cyan (distinct from mutations)
  "5. Postsynaptic Ribo (Only)" = "#8C564B",  # Brown (distinct)
  "6. Calcium Signaling"        = "#E377C2"   # Pink (softer than crimson mutation)
)

# =============================================================================
# RIBOSOME POOL COLORS (for ribosome paradox figures)
# =============================================================================

RIBOSOME_POOL_COLORS <- c(
  "1. Cytoplasmic Ribosome Biogenesis"   = "#B35806",
  "2. Synaptic Ribosomes (SynGO)"        = "#2166AC",
  "3. Mitochondrial Ribosomes (MitoCarta)"= "#D4A574"
)

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

#' Print color palette summary
#'
#' @return Invisible NULL, prints summary to console
print_color_summary <- function() {
  cat("\n=== DRP1 Analysis Color Palette Summary ===\n\n")

  cat("Diverging (NES/logFC):\n")
  cat(sprintf("  Negative: %s (Blue)\n", DIVERGING_COLORS$negative))
  cat(sprintf("  Neutral:  %s (White)\n", DIVERGING_COLORS$neutral))
  cat(sprintf("  Positive: %s (Orange)\n", DIVERGING_COLORS$positive))

  cat("\nMutation Colors:\n")
  for (name in names(MUTATION_COLORS)) {
    cat(sprintf("  %s: %s\n", name, MUTATION_COLORS[name]))
  }

  cat("\nPattern Colors:\n")
  for (name in names(PATTERN_COLORS)) {
    cat(sprintf("  %s: %s\n", name, PATTERN_COLORS[name]))
  }

  invisible(NULL)
}

# =============================================================================
# PHEATMAP HELPER FUNCTIONS
# =============================================================================

#' Get pheatmap-compatible diverging color palette
#'
#' @param n Number of colors to generate
#' @return Character vector of hex colors
#' @examples
#' pheatmap(..., color = get_diverging_palette(100))
get_diverging_palette <- function(n = 100) {
  if (!requireNamespace("grDevices", quietly = TRUE)) {
    stop("Package 'grDevices' required")
  }

  grDevices::colorRampPalette(c(
    DIVERGING_COLORS$negative,
    DIVERGING_COLORS$neutral,
    DIVERGING_COLORS$positive
  ))(n)
}

#' Get sample correlation heatmap color palette
#'
#' @param n Number of colors to generate
#' @return Character vector of hex colors for correlation values (0 to 1)
#' @examples
#' pheatmap(..., color = get_correlation_palette(50))
get_correlation_palette <- function(n = 50) {
  if (!requireNamespace("grDevices", quietly = TRUE)) {
    stop("Package 'grDevices' required")
  }

  # For correlation (typically 0.8-1.0), use a sequential palette
  # Light to dark: white -> teal for positive correlations
  grDevices::colorRampPalette(c(
    "#FFFFFF",   # White (low correlation)
    "#B2DFDB",   # Light teal
    "#4DB6AC",   # Medium teal
    "#00796B",   # Dark teal
    "#004D40"    # Very dark teal (high correlation)
  ))(n)
}

#' Get annotation colors for pheatmap/ComplexHeatmap
#'
#' Returns a list suitable for annotation_colors parameter in pheatmap
#' or col parameter in HeatmapAnnotation
#'
#' @return List of named color vectors
#' @examples
#' pheatmap(..., annotation_colors = get_heatmap_ann_colors())
get_heatmap_ann_colors <- function() {
  list(
    Genotype  = HEATMAP_ANNOTATION_COLORS$genotype,
    genotype  = HEATMAP_ANNOTATION_COLORS$genotype,
    Mutation  = HEATMAP_ANNOTATION_COLORS$mutation,
    mutation  = HEATMAP_ANNOTATION_COLORS$mutation,
    Stage     = HEATMAP_ANNOTATION_COLORS$stage,
    Trajectory = HEATMAP_ANNOTATION_COLORS$stage,
    Days      = HEATMAP_ANNOTATION_COLORS$timepoint,
    days      = HEATMAP_ANNOTATION_COLORS$timepoint,
    Timepoint = HEATMAP_ANNOTATION_COLORS$timepoint
  )
}

#' Create a color scale for z-scored expression heatmaps
#'
#' @param limits Numeric vector of length 2 for min/max z-score
#' @return colorRamp2 function for ComplexHeatmap or character vector for pheatmap
#' @examples
#' col_fun <- zscore_color_scale(limits = c(-2, 2))
zscore_color_scale <- function(limits = c(-2, 2), for_pheatmap = FALSE) {
  if (for_pheatmap) {
    return(get_diverging_palette(100))
  }

  if (!requireNamespace("circlize", quietly = TRUE)) {
    stop("Package 'circlize' required for colorRamp2")
  }

  circlize::colorRamp2(
    c(limits[1], 0, limits[2]),
    c(DIVERGING_COLORS$negative, DIVERGING_COLORS$neutral, DIVERGING_COLORS$positive)
  )
}

# =============================================================================
# BACKWARD COMPATIBILITY
# =============================================================================
# These aliases ensure old code continues to work

# Legacy alias for mutation colors
genotype_colors <- MUTATION_COLORS

# Legacy alias for ann_colors in main pipeline
ann_colors_heatmap <- get_heatmap_ann_colors()

message("Loaded unified color configuration from 01_Scripts/R_scripts/color_config.R")
