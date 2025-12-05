"""
Unified Color Configuration for DRP1 Bulk RNA-seq Analysis (Python)

This file provides a SINGLE SOURCE OF TRUTH for color schemes used across
all Python visualization scripts. It mirrors the R color definitions in:
  - 01_Scripts/R_scripts/color_config.R

COLOR SYSTEM OVERVIEW:
----------------------
1. DIVERGING_COLORS: Blue-White-Orange for NES/logFC heatmaps
   - AVOID using these colors for categorical annotations!
2. HEATMAP_ANNOTATION_COLORS: Distinct colors for heatmap annotations
   - Use these for row/column annotations on diverging heatmaps
3. MUTATION_COLORS: Standard mutation colors for general use
4. Other categorical palettes for specific use cases

COLORBLIND-SAFE: All palettes designed with deuteranopia/protanopia
in mind. Avoid relying solely on red-green distinctions.

Usage:
    from Python.color_config import DIVERGING_COLORS, MUTATION_COLORS
    from Python.color_config import create_diverging_cmap, HEATMAP_ANNOTATION_COLORS
"""

import numpy as np
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.colors as mcolors

# =============================================================================
# DIVERGING COLOR PALETTE (Colorblind-safe Blue-White-Orange)
# =============================================================================
# This palette is used for NES values, logFC, and other diverging metrics
# where negative/positive values have distinct biological meaning.
#
# IMPORTANT: When using this gradient for heatmaps, use HEATMAP_ANNOTATION_COLORS
# for any row/column annotations to avoid color conflicts.

DIVERGING_COLORS = {
    'negative': '#2166AC',  # Blue
    'neutral': '#F7F7F7',   # White
    'positive': '#B35806'   # Orange/Brown
}

# 5-point gradient for smoother diverging heatmaps
DIVERGING_COLORS_5PT = {
    'neg_strong': '#2166AC',
    'neg_weak': '#92c5de',
    'neutral': '#F7F7F7',
    'pos_weak': '#f4a582',
    'pos_strong': '#B35806'
}


def create_diverging_cmap(name='BlueWhiteOrange'):
    """
    Create a matplotlib colormap using the project's diverging colors.

    Returns
    -------
    LinearSegmentedColormap
        Colormap for NES/logFC heatmaps
    """
    colors = [DIVERGING_COLORS['negative'],
              DIVERGING_COLORS['neutral'],
              DIVERGING_COLORS['positive']]
    return LinearSegmentedColormap.from_list(name, colors)


# =============================================================================
# HEATMAP ANNOTATION COLORS
# =============================================================================
# These colors are specifically designed for row/column annotations on heatmaps
# that use the Blue-White-Orange diverging gradient. They AVOID blue and orange
# to prevent visual confusion.

HEATMAP_ANNOTATION_COLORS = {
    # Mutation colors for heatmap annotations (distinct from gradient blue/orange)
    'mutation': {
        'G32A': '#7B68EE',     # Medium Slate Blue (purple-ish, distinct from blue gradient)
        'R403C': '#DC143C',    # Crimson (red, distinct from orange gradient)
        'Ctrl': '#808080',     # Gray
        'Control': '#808080'   # Gray (alias)
    },

    # Stage/Trajectory colors (teal sequential - avoids red/pink confusion with mutations)
    'stage': {
        'Early': '#B2DFDB',    # Light teal
        'TrajDev': '#4DB6AC',  # Medium teal
        'Late': '#00796B'      # Dark teal
    },

    # Timepoint colors (for genotype x timepoint annotations)
    'timepoint': {
        'D35': '#66C2A5',      # Teal/mint
        'D65': '#8DA0CB'       # Periwinkle
    },

    # Genotype colors (alias for mutation with full names)
    'genotype': {
        'Control': '#808080',
        'G32A': '#7B68EE',
        'R403C': '#DC143C'
    }
}


# =============================================================================
# MUTATION COLORS (Standard - for bar charts, dotplots, non-heatmap figures)
# =============================================================================
# Colorblind-safe Okabe-Ito palette for mutation identity
# NOTE: For heatmap annotations, use HEATMAP_ANNOTATION_COLORS['mutation'] instead

MUTATION_COLORS = {
    'G32A': '#0072B2',   # Okabe-Ito Blue
    'R403C': '#D55E00',  # Okabe-Ito Vermillion
    'Ctrl': '#999999'    # Gray
}


# =============================================================================
# PATTERN COLORS (Trajectory classification)
# =============================================================================
# NOTE: These should match pattern_definitions.py - import from there for
# the canonical source. This is provided for convenience.

PATTERN_COLORS = {
    'Compensation': '#009E73',
    'Sign_reversal': '#9467BD',
    'Progressive': '#D55E00',
    'Natural_worsening': '#E69F00',
    'Natural_improvement': '#56B4E9',
    'Late_onset': '#CC79A7',
    'Transient': '#0072B2',
    'Complex': '#F0E442',
    'Insufficient_data': '#DDDDDD'
}


# =============================================================================
# SUPER-CATEGORY COLORS (Simplified pattern groupings)
# =============================================================================

SUPER_CATEGORY_COLORS = {
    'Active_Compensation': '#009E73',
    'Active_Reversal': '#9467BD',
    'Active_Progression': '#D55E00',
    'Passive': '#56B4E9',
    'Late_onset': '#CC79A7',
    'Other': '#999999',
    'Insufficient_data': '#DDDDDD'
}


# =============================================================================
# MODULE COLORS (for mechanistic cascade figures)
# =============================================================================
# Redesigned to avoid overlap with HEATMAP_ANNOTATION_COLORS

MODULE_COLORS = {
    '1. Mt Central Dogma': '#9467BD',         # Purple
    '2. Mt Ribosomes': '#DAA520',             # Goldenrod
    '3. ATP Synthase (Complex V)': '#2E8B57', # Sea green
    '4. Synaptic Ribo (Common)': '#17BECF',   # Cyan
    '5. Postsynaptic Ribo (Only)': '#8C564B', # Brown
    '6. Calcium Signaling': '#E377C2'         # Pink
}


# =============================================================================
# TRAJECTORY STAGE COLORS (for standalone trajectory figures)
# =============================================================================
# NOTE: For heatmap annotations, use HEATMAP_ANNOTATION_COLORS['stage'] instead

TRAJECTORY_COLORS = {
    'Early': '#FEE0D2',    # Light pink/salmon
    'TrajDev': '#FC9272',  # Medium salmon
    'Late': '#DE2D26'      # Dark red
}


# =============================================================================
# SEMANTIC CATEGORY COLORS (Biological pathway groupings)
# =============================================================================
# NOTE: This should match semantic_categories.py - consider importing from there

SEMANTIC_COLORS = {
    'Synapse': '#8B4513',
    'Neuronal Development': '#A0522D',
    'Mitochondrial Dynamics': '#006400',
    'Electron Transport Chain': '#2E8B57',
    'ATP Synthase (Complex V)': '#DAA520',
    'Mitochondrial Metabolism': '#228B22',
    'Mitochondrial Function': '#32CD32',
    'Mitochondrial Ribosome': '#4169E1',
    'Mitochondrial Translation': '#1E90FF',
    'Ribosome Biogenesis': '#6A5ACD',
    'Cytoplasmic Ribosome': '#9932CC',
    'Cytoplasmic Translation': '#BA55D3',
    'Calcium Signaling': '#DC143C',
    'Other': '#696969'
}


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def get_diverging_cmap():
    """Get the project's standard diverging colormap."""
    return create_diverging_cmap()


def get_correlation_palette(n=50):
    """
    Get a sequential palette for correlation heatmaps.

    Parameters
    ----------
    n : int
        Number of colors to generate

    Returns
    -------
    LinearSegmentedColormap
        Colormap for correlation values (typically 0.8-1.0)
    """
    colors = ['#FFFFFF', '#B2DFDB', '#4DB6AC', '#00796B', '#004D40']
    return LinearSegmentedColormap.from_list('correlation', colors, N=n)


def get_heatmap_ann_colors():
    """
    Get annotation colors for seaborn/matplotlib heatmaps.

    Returns
    -------
    dict
        Dictionary of color mappings for different annotation types
    """
    return HEATMAP_ANNOTATION_COLORS.copy()


def print_color_summary():
    """Print a summary of all color palettes."""
    print("\n=== DRP1 Analysis Color Palette Summary ===\n")

    print("Diverging (NES/logFC):")
    for name, color in DIVERGING_COLORS.items():
        print(f"  {name}: {color}")

    print("\nMutation Colors (standard):")
    for name, color in MUTATION_COLORS.items():
        print(f"  {name}: {color}")

    print("\nHeatmap Annotation - Mutation:")
    for name, color in HEATMAP_ANNOTATION_COLORS['mutation'].items():
        print(f"  {name}: {color}")

    print("\nHeatmap Annotation - Stage:")
    for name, color in HEATMAP_ANNOTATION_COLORS['stage'].items():
        print(f"  {name}: {color}")


if __name__ == '__main__':
    print_color_summary()
