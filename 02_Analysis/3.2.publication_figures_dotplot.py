#!/usr/bin/env python3
"""
Publication-Ready Trajectory Figures - DOTPLOT VERSION
=======================================================

Creates consolidated, colorblind-safe dotplot figures for manuscript:

1. Ribosome Paradox - Combined SynGO + MitoCarta showing opposite ribosome fates
2. Shared Compensation - MitoCarta pathways rescued in both mutations
3. Pattern Classification Summary - Quantitative overview across all databases
4. Semantic Pathway Overview - Comprehensive dotplot by biological category

Dotplot Features:
- Dot COLOR = NES value (Blue-White-Orange diverging palette)
- Dot SIZE = -log10(padj) - larger dots = more significant
- Dot EDGE = Black outline for significant (padj < 0.05), gray for non-significant
- CGP database excluded (cancer-focused, not relevant for neurobiology)
"""

import sys
from pathlib import Path

# Add module paths
sys.path.insert(0, str(Path(__file__).parent.parent / '01_Scripts'))
sys.path.insert(0, str(Path(__file__).parent.parent / '01_Scripts/RNAseq-toolkit/scripts/GSEA'))

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Patch
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import warnings
warnings.filterwarnings('ignore')

# Import from reusable toolkit
from GSEA_plotting_python import create_diverging_cmap, DotplotRenderer, DIVERGING_COLORS

# Import from project-specific modules
from Python.config import CONFIG, resolve_path, ensure_output_dir, EXCLUDE_PATHWAYS
from Python.semantic_categories import (
    SEMANTIC_CATEGORY_ORDER, SEMANTIC_COLORS, PATTERN_COLORS, MUTATION_COLORS,
    assign_semantic_category, categorize_pathways, sort_by_category, deduplicate_pathways
)
from Python.data_loader import load_classified_pathways, filter_pathways
from Python.patterns import is_compensation, add_pattern_classification
from Python.pattern_definitions import MEANINGFUL_PATTERNS, ACTIVE_PATTERNS

# Import unified color config for heatmap annotation colors
# Use these when mutations appear alongside NES gradient to avoid color confusion
from Python.color_config import HEATMAP_ANNOTATION_COLORS

# =============================================================================
# SETUP
# =============================================================================

# Create colormap
CMAP_DIVERGING = create_diverging_cmap()

# Ensure output directory exists with different name for dotplot version
OUTPUT_DIR = resolve_path('03_Results/02_Analysis/Plots/Publication_Figures_Dotplot')
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def add_size_legend_visual(ax, renderer, sig_cutoff=0.05):
    """
    Add a visual size legend showing the padj-to-size relationship.

    Creates example dots at different significance levels with labels.
    """
    ax.axis('off')

    # Define example padj values to show
    padj_values = [0.001, 0.01, 0.05, 0.1]
    labels = ['highly sig', 'very sig', 'threshold', 'not sig']

    # Calculate sizes using renderer's method (use actual sizes, not scaled)
    padj_array = np.array(padj_values).reshape(1, -1)
    dot_sizes = renderer._calculate_dot_sizes(padj_array)[0]

    # Position dots horizontally with more space
    x_positions = np.linspace(0.15, 0.85, len(padj_values))
    y_position = 0.5

    # Draw dots with ACTUAL sizes from the renderer
    for i, (x, size, padj, label) in enumerate(zip(x_positions, dot_sizes, padj_values, labels)):
        # Determine edge color
        edge_color = 'black' if padj < sig_cutoff else 'gray'
        edge_width = 2 if padj < sig_cutoff else 1

        # Draw dot using ACTUAL calculated size
        ax.scatter([x], [y_position], s=size, c='gray',
                  edgecolors=edge_color, linewidths=edge_width,
                  alpha=0.7, zorder=3)

        # Add label below with better spacing
        ax.text(x, 0.28, f'padj = {padj}',
               ha='center', va='top', fontsize=10, fontweight='bold')
        ax.text(x, 0.15, f'({label})',
               ha='center', va='top', fontsize=9, style='italic', color='#555')

    # Add title at top
    ax.text(0.5, 0.90, 'Dot Size Legend: Larger dots = more significant (lower padj)',
           ha='center', va='center', fontsize=11, fontweight='bold')

    # Add edge color explanation at bottom
    ax.text(0.5, 0.02, 'Black edge = padj < 0.05 (significant)  |  Gray edge = padj ≥ 0.05 (not significant)',
           ha='center', va='center', fontsize=9, style='italic')

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)


# =============================================================================
# DATA LOADING (with CGP exclusion)
# =============================================================================

def load_data():
    """Load classified pathway data with CGP excluded"""
    print("Loading data (excluding CGP database)...")
    df = load_classified_pathways()  # CGP excluded by default
    return df


# =============================================================================
# FIGURE 1: RIBOSOME PARADOX (DOTPLOT VERSION)
# =============================================================================

def create_ribosome_paradox_figure(df):
    """
    Create combined figure showing the Ribosome Paradox with 3 panels (dotplot version):
    - Panel A: Synaptic ribosomes (SynGO) - UP early, DOWN late [FAILURE]
    - Panel B: Mitochondrial ribosomes (MitoCarta) - DOWN early, UP recovery [PARTIAL RESCUE]
    - Panel C: Cytoplasmic Ribosome Biogenesis (GO:BP) - Shows compensatory context
    """
    print("\nCreating Ribosome Paradox figure (3 panels - DOTPLOT)...")

    # Extract ribosome-related pathways (same as original)
    syngo_ribosomes = df[
        (df['database'] == 'SynGO') &
        (df['Description'].str.contains('ribosome', case=False, na=False))
    ].copy()

    mito_ribosomes = df[
        (df['database'] == 'MitoCarta') &
        (df['Description'].str.contains('ribosome|translation|central_dogma', case=False, na=False))
    ].copy()

    gobp_ribosomes = df[
        (df['database'].isin(['gobp', 'gocc', 'gomf'])) &
        (df['Description'].str.contains('ribosom', case=False, na=False)) &
        (~df['Description'].str.contains('mitochondri', case=False, na=False))
    ].copy()

    # Filter Panel C: Keep only pathways with at least one significant result
    padj_cols = CONFIG['padj_cols']
    def has_any_significant(row):
        for col in padj_cols:
            if col in row and pd.notna(row[col]) and row[col] < CONFIG['padj_cutoff']:
                return True
        return False

    gobp_ribosomes['has_sig'] = gobp_ribosomes.apply(has_any_significant, axis=1)
    gobp_ribosomes = gobp_ribosomes[gobp_ribosomes['has_sig']].copy()

    print(f"  SynGO ribosomes: {len(syngo_ribosomes)}")
    print(f"  MitoCarta ribosomes: {len(mito_ribosomes)}")
    print(f"  GO ribosome biogenesis: {len(gobp_ribosomes)}")

    # Create figure with THREE panels + colorbar + size legend
    fig = plt.figure(figsize=(13, 15))
    gs = GridSpec(4, 2, width_ratios=[6, 0.3],
                  height_ratios=[1, 1, 1, 0.25],
                  wspace=0.05, hspace=0.35)

    ax_syngo = fig.add_subplot(gs[0, 0])
    ax_mito = fig.add_subplot(gs[1, 0])
    ax_cyto = fig.add_subplot(gs[2, 0])
    ax_cbar = fig.add_subplot(gs[0:3, 1])
    ax_legend = fig.add_subplot(gs[3, :])

    vmax = CONFIG['vmax']
    renderer = DotplotRenderer(cmap=CMAP_DIVERGING, vmax=vmax,
                               min_dot_size=20, max_dot_size=300)

    # Panel A: Synaptic Ribosomes
    _plot_trajectory_dotplot(
        ax_syngo, syngo_ribosomes, renderer,
        title='A. Synaptic Ribosomes (SynGO)',
        subtitle='Early compensatory upregulation fails during maturation'
    )

    # Panel B: Mitochondrial Ribosomes
    _plot_trajectory_dotplot(
        ax_mito, mito_ribosomes, renderer,
        title='B. Mitochondrial Ribosomes (MitoCarta)',
        subtitle='Early crisis with partial recovery during maturation'
    )

    # Panel C: Cytoplasmic Ribosome Biogenesis
    _plot_trajectory_dotplot(
        ax_cyto, gobp_ribosomes, renderer,
        title='C. Cytoplasmic Ribosome Biogenesis (GO)',
        subtitle='Broader ribosome biogenesis context'
    )

    # Shared colorbar
    norm = mcolors.Normalize(vmin=-vmax, vmax=vmax)
    sm = plt.cm.ScalarMappable(cmap=CMAP_DIVERGING, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, cax=ax_cbar)
    cbar.set_label('Normalized Enrichment Score (NES)', fontsize=11)
    cbar.ax.tick_params(labelsize=10)

    # Add size legend
    add_size_legend_visual(ax_legend, renderer)

    # Add annotation explaining the paradox
    fig.text(0.5, -0.01,
             'The Ribosome Paradox:\n'
             'A) Synaptic ribosomes show early UP (compensation attempt) → late DOWN (failure)\n'
             'B) Mitochondrial ribosomes show early DOWN (crisis) → late UP (recovery)\n'
             'C) Cytoplasmic ribosome biogenesis provides context for compensatory responses',
             ha='center', va='top', fontsize=9, style='italic',
             bbox=dict(boxstyle='round', facecolor='#F5F5F5', edgecolor='gray', alpha=0.8))

    # Main title
    fig.suptitle('The Ribosome Paradox: Opposite Fates of Translation Machineries (Dotplot)',
                fontsize=14, fontweight='bold', y=0.99)

    plt.tight_layout(rect=[0, 0.10, 0.95, 0.97])

    # Save PDF and PNG
    for ext, dpi in [('pdf', CONFIG['dpi']), ('png', 150)]:
        output_file = OUTPUT_DIR / f'Fig1_Ribosome_Paradox_dotplot.{ext}'
        fig.savefig(output_file, dpi=dpi, bbox_inches='tight', facecolor='white')
        print(f"  Saved: {output_file}")

    plt.close(fig)

    return syngo_ribosomes, mito_ribosomes, gobp_ribosomes


def _plot_trajectory_dotplot(ax, df_subset, renderer, title, subtitle):
    """
    Helper to plot a single trajectory dotplot panel.

    Uses scatter plot with:
    - Color = NES value
    - Size = -log10(padj)
    - Edge color = black for significant, gray for non-significant
    """
    if len(df_subset) == 0:
        ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
        ax.set_title(title, fontsize=12, fontweight='bold')
        return

    # Prepare data matrix
    nes_cols = CONFIG['nes_cols']
    padj_cols = CONFIG['padj_cols']

    # Sort by maximum absolute NES for visual clarity
    df_subset = df_subset.copy()
    df_subset['max_nes'] = df_subset[nes_cols].abs().max(axis=1)
    df_subset = df_subset.sort_values('max_nes', ascending=True)

    pathway_names = df_subset['Description'].values
    n_pathways = len(pathway_names)

    # Build combined matrix
    nes_matrix = df_subset[nes_cols].values
    padj_matrix = df_subset[padj_cols].values

    # Render dotplot
    renderer.render(ax, nes_matrix, padj_matrix)

    # Add thick vertical line separating mutations
    ax.axvline(2.5, color='white', linewidth=6, zorder=1)
    ax.axvline(2.5, color='black', linewidth=2, zorder=2)

    # Labels
    col_labels = ['Early\nG32A', 'TrajDev\nG32A', 'Late\nG32A',
                  'Early\nR403C', 'TrajDev\nR403C', 'Late\nR403C']
    ax.set_xticks(range(6))
    ax.set_xticklabels(col_labels, fontsize=9)
    ax.set_yticks(range(n_pathways))
    ax.set_yticklabels(pathway_names, fontsize=9)

    # Title
    ax.set_title(f'{title}\n{subtitle}', fontsize=11, fontweight='bold', pad=10)


# =============================================================================
# FIGURE 2: MitoCarta TRAJECTORY PATTERNS (DOTPLOT VERSION)
# =============================================================================

def create_mitocarta_trajectory_figure(df):
    """
    Create single-panel figure showing ALL MitoCarta pathways (dotplot version).
    """
    print("\nCreating MitoCarta Trajectory Patterns figure (all pathways - DOTPLOT)...")

    # Filter for MitoCarta
    df_mito = df[df['database'] == 'MitoCarta'].copy()
    print(f"  Total MitoCarta pathways: {len(df_mito)}")

    # Filter for pathways with at least one significant result
    padj_cols = CONFIG['padj_cols']

    def has_any_significant(row):
        for col in padj_cols:
            if col in row and pd.notna(row[col]) and row[col] < CONFIG['padj_cutoff']:
                return True
        return False

    df_mito['has_sig'] = df_mito.apply(has_any_significant, axis=1)
    df_mito = df_mito[df_mito['has_sig']].copy()
    print(f"  Pathways with ≥1 significant contrast: {len(df_mito)}")

    if len(df_mito) == 0:
        print("  No significant pathways found - skipping figure")
        return None

    # Sort by max |NES| for visual clarity
    nes_cols = CONFIG['nes_cols']
    df_mito['max_nes'] = df_mito[nes_cols].abs().max(axis=1)
    df_mito = df_mito.sort_values('max_nes', ascending=True)

    n_pathways = len(df_mito)

    # Create figure with size legend
    fig_height = max(10, n_pathways * 0.4 + 3)  # Extra space for legend
    fig = plt.figure(figsize=(14, fig_height))
    gs = GridSpec(2, 2, width_ratios=[6, 0.3],
                  height_ratios=[n_pathways * 0.4, 2.0],
                  hspace=0.15, wspace=0.05)

    ax_main = fig.add_subplot(gs[0, 0])
    ax_cbar = fig.add_subplot(gs[0, 1])
    ax_legend = fig.add_subplot(gs[1, :])

    vmax = CONFIG['vmax']
    renderer = DotplotRenderer(cmap=CMAP_DIVERGING, vmax=vmax,
                               min_dot_size=20, max_dot_size=300)

    # Plot all pathways
    _plot_trajectory_panel_dotplot(
        ax_main, df_mito, renderer,
        title='MitoCarta Mitochondrial Pathways',
        subtitle=f'All pathways with ≥1 significant contrast (n={n_pathways})',
        color=None
    )

    # Shared colorbar
    norm = mcolors.Normalize(vmin=-vmax, vmax=vmax)
    sm = plt.cm.ScalarMappable(cmap=CMAP_DIVERGING, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, cax=ax_cbar)
    cbar.set_label('Normalized Enrichment Score (NES)', fontsize=11)
    cbar.ax.tick_params(labelsize=10)

    # Add size legend
    add_size_legend_visual(ax_legend, renderer)

    # Main title
    fig.suptitle('MitoCarta Pathway Trajectories (Dotplot)\nG32A and R403C may show different trajectory patterns',
                fontsize=14, fontweight='bold', y=0.99)

    # Legend
    fig.text(0.5, 0.01,
             'Note: Each mutation follows its own trajectory - no forced classification',
             ha='center', va='bottom', fontsize=9,
             bbox=dict(boxstyle='round', facecolor='#F5F5F5', edgecolor='gray', alpha=0.8))

    plt.tight_layout(rect=[0, 0.02, 0.95, 0.97])

    # Save
    for ext, dpi in [('pdf', CONFIG['dpi']), ('png', 150)]:
        output_file = OUTPUT_DIR / f'Fig2_MitoCarta_Trajectory_Patterns_dotplot.{ext}'
        fig.savefig(output_file, dpi=dpi, bbox_inches='tight', facecolor='white')
        print(f"  Saved: {output_file}")

    plt.close(fig)

    return df_mito


def _plot_trajectory_panel_dotplot(ax, df_subset, renderer, title, subtitle, color=None):
    """
    Helper to plot a single trajectory pattern panel (dotplot version).
    """
    if len(df_subset) == 0:
        ax.text(0.5, 0.5, 'No pathways in this category',
                ha='center', va='center', transform=ax.transAxes, fontsize=10)
        ax.set_title(f'{title}\n{subtitle}', fontsize=11, fontweight='bold', pad=10)
        ax.axis('off')
        return

    nes_cols = CONFIG['nes_cols']
    padj_cols = CONFIG['padj_cols']

    # Sort by max |NES| for visual clarity
    df_subset = df_subset.copy()
    df_subset['max_nes'] = df_subset[nes_cols].abs().max(axis=1)
    df_subset = df_subset.sort_values('max_nes', ascending=True)

    pathway_names = df_subset['Description'].values
    n_pathways = len(pathway_names)

    nes_matrix = df_subset[nes_cols].values
    padj_matrix = df_subset[padj_cols].values

    # Render dotplot
    renderer.render(ax, nes_matrix, padj_matrix)

    # Thick vertical separator between mutations
    ax.axvline(2.5, color='white', linewidth=6, zorder=1)
    ax.axvline(2.5, color='black', linewidth=2, zorder=2)

    # Labels
    col_labels = ['Early\nG32A', 'TrajDev\nG32A', 'Late\nG32A',
                  'Early\nR403C', 'TrajDev\nR403C', 'Late\nR403C']
    ax.set_xticks(range(6))
    ax.set_xticklabels(col_labels, fontsize=9)
    ax.set_yticks(range(n_pathways))
    ax.set_yticklabels(pathway_names, fontsize=9)

    # Title (optionally with colored indicator)
    title_color = color if color else 'black'
    ax.set_title(f'{title}\n{subtitle}', fontsize=11, fontweight='bold', pad=10,
                color=title_color)


# =============================================================================
# FIGURE 3b: SynGO TRAJECTORY PATTERNS (DOTPLOT VERSION)
# =============================================================================

def create_syngo_trajectory_figure(df):
    """
    Create single-panel figure showing ALL SynGO pathways (dotplot version).
    """
    print("\nCreating SynGO Trajectory Patterns figure (all pathways - DOTPLOT)...")

    # Filter for SynGO
    df_syngo = df[df['database'] == 'SynGO'].copy()
    print(f"  Total SynGO pathways: {len(df_syngo)}")

    # Filter for pathways with at least one significant result
    padj_cols = CONFIG['padj_cols']

    def has_any_significant(row):
        for col in padj_cols:
            if col in row and pd.notna(row[col]) and row[col] < CONFIG['padj_cutoff']:
                return True
        return False

    df_syngo['has_sig'] = df_syngo.apply(has_any_significant, axis=1)
    df_syngo = df_syngo[df_syngo['has_sig']].copy()
    print(f"  Pathways with ≥1 significant contrast: {len(df_syngo)}")

    if len(df_syngo) == 0:
        print("  No significant pathways found - skipping figure")
        return None

    # Sort by max |NES| for visual clarity
    nes_cols = CONFIG['nes_cols']
    df_syngo['max_nes'] = df_syngo[nes_cols].abs().max(axis=1)
    df_syngo = df_syngo.sort_values('max_nes', ascending=True)

    n_pathways = len(df_syngo)

    # Create figure with size legend
    fig_height = max(10, n_pathways * 0.4 + 3)  # Extra space for legend
    fig = plt.figure(figsize=(14, fig_height))
    gs = GridSpec(2, 2, width_ratios=[6, 0.3],
                  height_ratios=[n_pathways * 0.4, 2.0],
                  hspace=0.15, wspace=0.05)

    ax_main = fig.add_subplot(gs[0, 0])
    ax_cbar = fig.add_subplot(gs[0, 1])
    ax_legend = fig.add_subplot(gs[1, :])

    vmax = CONFIG['vmax']
    renderer = DotplotRenderer(cmap=CMAP_DIVERGING, vmax=vmax,
                               min_dot_size=20, max_dot_size=300)

    # Plot all pathways
    _plot_trajectory_panel_dotplot(
        ax_main, df_syngo, renderer,
        title='SynGO Synaptic Pathways',
        subtitle=f'All pathways with ≥1 significant contrast (n={n_pathways})',
        color=None
    )

    # Shared colorbar
    norm = mcolors.Normalize(vmin=-vmax, vmax=vmax)
    sm = plt.cm.ScalarMappable(cmap=CMAP_DIVERGING, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, cax=ax_cbar)
    cbar.set_label('Normalized Enrichment Score (NES)', fontsize=11)
    cbar.ax.tick_params(labelsize=10)

    # Add size legend
    add_size_legend_visual(ax_legend, renderer)

    # Main title
    fig.suptitle('SynGO Synaptic Pathway Trajectories (Dotplot)\nG32A and R403C may show different trajectory patterns',
                fontsize=14, fontweight='bold', y=0.99)

    # Legend
    fig.text(0.5, 0.01,
             'Note: Each mutation follows its own trajectory - no forced classification',
             ha='center', va='bottom', fontsize=9,
             bbox=dict(boxstyle='round', facecolor='#F5F5F5', edgecolor='gray', alpha=0.8))

    plt.tight_layout(rect=[0, 0.02, 0.95, 0.97])

    # Save
    for ext, dpi in [('pdf', CONFIG['dpi']), ('png', 150)]:
        output_file = OUTPUT_DIR / f'Fig3b_SynGO_Trajectory_Patterns_dotplot.{ext}'
        fig.savefig(output_file, dpi=dpi, bbox_inches='tight', facecolor='white')
        print(f"  Saved: {output_file}")

    plt.close(fig)

    return df_syngo


# =============================================================================
# FIGURE 3: PATTERN CLASSIFICATION SUMMARY (same as original - no changes)
# =============================================================================

def create_pattern_summary_figure(df):
    """
    Create summary bar chart showing pattern classification across all databases.
    (This figure is unchanged from the original - it's not a heatmap/dotplot)
    """
    print("\nCreating Pattern Classification Summary...")

    # Ensure pattern columns exist
    if 'Pattern_G32A' not in df.columns or 'Pattern_R403C' not in df.columns:
        df = add_pattern_classification(df)

    # Count patterns per database and mutation
    pattern_counts = []

    for database in df['database'].unique():
        df_db = df[df['database'] == database]

        for mutation in ['G32A', 'R403C']:
            pattern_col = f'Pattern_{mutation}'
            if pattern_col not in df_db.columns:
                continue

            counts = df_db[pattern_col].value_counts()
            for pattern, count in counts.items():
                pattern_counts.append({
                    'Database': database,
                    'Mutation': mutation,
                    'Pattern': pattern,
                    'Count': count
                })

    df_counts = pd.DataFrame(pattern_counts)

    # Filter to meaningful patterns (using centralized constant)
    df_meaningful = df_counts[df_counts['Pattern'].isin(MEANINGFUL_PATTERNS)]

    # Collect all patterns present across BOTH mutations for legend
    all_patterns_present = set()
    for mut in ['G32A', 'R403C']:
        df_mut = df_meaningful[df_meaningful['Mutation'] == mut]
        all_patterns_present.update(df_mut['Pattern'].unique())

    # Create figure
    fig, axes = plt.subplots(1, 2, figsize=(14, 8), sharey=False)

    # Use heatmap annotation colors for mutations (distinct from NES gradient blue/orange)
    # This prevents visual confusion when figures are viewed alongside NES heatmaps
    mutation_bar_colors = HEATMAP_ANNOTATION_COLORS['mutation']

    # Panel A: Counts by database
    ax1 = axes[0]

    # Get consistent database order
    databases = df_meaningful['Database'].unique()
    databases = sorted(databases)
    x = np.arange(len(databases))
    width = 0.35

    for mut_idx, mutation in enumerate(['G32A', 'R403C']):
        df_mut = df_meaningful[df_meaningful['Mutation'] == mutation]
        pivot = df_mut.pivot_table(index='Database', columns='Pattern',
                                   values='Count', fill_value=0)

        # Reindex to ensure all databases are present
        pivot = pivot.reindex(databases, fill_value=0)

        pattern_order = [p for p in MEANINGFUL_PATTERNS if p in pivot.columns]
        pivot = pivot[pattern_order]

        bottom = np.zeros(len(pivot))
        offset = (mut_idx - 0.5) * width

        for pattern in pattern_order:
            if pattern in pivot.columns:
                color = PATTERN_COLORS.get(pattern, '#999999')
                ax1.barh(x + offset, pivot[pattern], width, left=bottom,
                        color=color, alpha=0.85)
                bottom += pivot[pattern].values

        # Add mutation labels next to each bar
        # Use heatmap annotation colors for consistency with Panel B
        for i, db in enumerate(databases):
            bar_total = pivot.loc[db].sum() if db in pivot.index else 0
            if bar_total > 0:
                ax1.text(bar_total + 1, i + offset, mutation[:4],
                        fontsize=8, fontweight='bold', va='center',
                        color=mutation_bar_colors[mutation])

    ax1.set_yticks(x)
    ax1.set_yticklabels(databases, fontsize=10)
    ax1.set_xlabel('Number of Pathways', fontsize=11)
    ax1.set_title('A. Pattern Distribution by Database', fontsize=12, fontweight='bold')

    # Create legend handles manually to include ALL patterns present
    legend_handles = [mpatches.Patch(color=PATTERN_COLORS[p], label=p, alpha=0.85)
                     for p in MEANINGFUL_PATTERNS if p in all_patterns_present]
    ax1.legend(handles=legend_handles, loc='center right', fontsize=9, framealpha=0.9,
               bbox_to_anchor=(1.0, 0.5))

    # Panel B: Focus on compensation vs progressive
    ax2 = axes[1]

    summary_data = []
    for mutation in ['G32A', 'R403C']:
        df_mut = df_counts[df_counts['Mutation'] == mutation]

        comp_count = df_mut[df_mut['Pattern'] == 'Compensation']['Count'].sum()
        rev_count = df_mut[df_mut['Pattern'] == 'Sign_reversal']['Count'].sum()
        worsening_count = df_mut[df_mut['Pattern'].isin(['Progressive', 'Natural_worsening'])]['Count'].sum()
        other_count = df_mut[~df_mut['Pattern'].isin(['Compensation', 'Sign_reversal', 'Progressive',
                                                       'Natural_worsening', 'Insufficient_data'])]['Count'].sum()

        summary_data.append({'Mutation': mutation, 'Pattern': 'Compensation', 'Count': comp_count})
        summary_data.append({'Mutation': mutation, 'Pattern': 'Sign_reversal', 'Count': rev_count})
        summary_data.append({'Mutation': mutation, 'Pattern': 'Worsening', 'Count': worsening_count})
        summary_data.append({'Mutation': mutation, 'Pattern': 'Other', 'Count': other_count})

    df_summary = pd.DataFrame(summary_data)

    x = np.arange(4)
    width = 0.35

    for mut_idx, mutation in enumerate(['G32A', 'R403C']):
        df_mut = df_summary[df_summary['Mutation'] == mutation]
        counts = [
            df_mut[df_mut['Pattern'] == 'Compensation']['Count'].values[0],
            df_mut[df_mut['Pattern'] == 'Sign_reversal']['Count'].values[0],
            df_mut[df_mut['Pattern'] == 'Worsening']['Count'].values[0],
            df_mut[df_mut['Pattern'] == 'Other']['Count'].values[0]
        ]

        bars = ax2.bar(x + (mut_idx - 0.5) * width, counts, width,
                      label=mutation, color=mutation_bar_colors[mutation], alpha=0.85)

        for bar, count in zip(bars, counts):
            ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 2,
                    str(count), ha='center', va='bottom', fontsize=10, fontweight='bold')

    ax2.set_xticks(x)
    ax2.set_xticklabels(['Compensation\n(rescue)', 'Sign reversal\n(trajectory flip)',
                         'Worsening\n(progressive)', 'Other\n(complex)'],
                        fontsize=11)
    ax2.set_ylabel('Number of Pathways', fontsize=11)
    ax2.set_title('B. Summary: Compensation vs Worsening', fontsize=12, fontweight='bold')
    ax2.legend(loc='upper right', fontsize=10)

    # Add interpretation
    comp_g32a = df_summary[(df_summary['Mutation']=='G32A') & (df_summary['Pattern']=='Compensation')]['Count'].values[0]
    comp_r403c = df_summary[(df_summary['Mutation']=='R403C') & (df_summary['Pattern']=='Compensation')]['Count'].values[0]

    interpretation = (
        f'Key Finding: R403C shows {comp_r403c} compensation pathways vs {comp_g32a} for G32A\n'
        f'suggesting stronger adaptive responses in R403C mutation.'
    )
    fig.text(0.5, 0.02, interpretation, ha='center', fontsize=10, style='italic',
             bbox=dict(boxstyle='round', facecolor='#E8F4E8', edgecolor='green', alpha=0.8))

    fig.suptitle('Trajectory Pattern Classification Across All Databases',
                fontsize=14, fontweight='bold', y=0.98)

    plt.tight_layout(rect=[0, 0.06, 1, 0.95])

    # Save (note: using same filename as original since this figure is unchanged)
    for ext, dpi in [('pdf', CONFIG['dpi']), ('png', 150)]:
        output_file = OUTPUT_DIR / f'Fig3_Pattern_Classification_Summary.{ext}'
        fig.savefig(output_file, dpi=dpi, bbox_inches='tight', facecolor='white')
        print(f"  Saved: {output_file}")

    plt.close(fig)

    return df_counts


# =============================================================================
# FIGURE 4: SEMANTIC PATHWAY OVERVIEW (DOTPLOT VERSION)
# =============================================================================

def create_semantic_grouping_figure(df):
    """
    Create comprehensive semantic-grouped pathway dotplot.
    """
    print("\nCreating Semantic Pathway Overview figure (DOTPLOT)...")

    # Prepare pathway data with NES columns
    nes_cols = CONFIG['nes_cols']
    padj_cols = CONFIG['padj_cols']

    # Filter pathways by data points
    df_filtered = filter_pathways(df, min_data_points=CONFIG['min_data_points'])
    print(f"  Pathways with ≥{CONFIG['min_data_points']} data points: {len(df_filtered)}")

    # Assign semantic categories WITH exclusion list
    df_filtered = categorize_pathways(df_filtered, exclude_pathways=EXCLUDE_PATHWAYS)
    print(f"  After excluding irrelevant pathways: {len(df_filtered)}")

    # Deduplicate pathways across databases
    df_filtered = deduplicate_pathways(df_filtered, nes_cols=nes_cols)
    print(f"  After deduplication: {len(df_filtered)}")

    # Further filter: only pathways with at least one significant result
    def has_significant(row):
        for col in padj_cols:
            if col in row and pd.notna(row[col]) and row[col] < CONFIG['padj_cutoff']:
                return True
        return False

    df_filtered['has_sig'] = df_filtered.apply(has_significant, axis=1)
    df_sig = df_filtered[df_filtered['has_sig']].copy()
    print(f"  Pathways with significant results: {len(df_sig)}")

    if len(df_sig) == 0:
        print("  No significant pathways - skipping figure")
        return None

    # Calculate max |NES| for sorting within category
    df_sig['max_nes'] = df_sig[nes_cols].abs().max(axis=1)

    # Sort by semantic category and NES
    df_sig = sort_by_category(df_sig, nes_col='max_nes')

    # Limit to top pathways per category
    max_per_category = CONFIG['max_per_category']
    df_plot = df_sig.groupby('Semantic_Category').head(max_per_category).reset_index(drop=True)
    print(f"  Pathways to plot (max {max_per_category} per category): {len(df_plot)}")

    # Group pathways by category for separate dotplot blocks
    category_data = []

    for cat in SEMANTIC_CATEGORY_ORDER:
        cat_df = df_plot[df_plot['Semantic_Category'] == cat]
        if len(cat_df) > 0:
            color = SEMANTIC_COLORS.get(cat, '#333333')
            category_data.append((cat, cat_df, color))

    n_categories = len(category_data)
    print(f"  Categories with data: {n_categories}")

    if n_categories == 0:
        print("  No categories with data - skipping figure")
        return None

    # Calculate heights for each category block
    GAP_FRACTION = 0.8
    ROW_HEIGHT = 0.35

    category_heights = []
    total_pathways = 0
    for cat, cat_df, color in category_data:
        n_rows = len(cat_df)
        category_heights.append(n_rows)
        total_pathways += n_rows

    # Total figure height
    total_gap_height = (n_categories - 1) * GAP_FRACTION * ROW_HEIGHT * 3
    fig_height = max(10, total_pathways * ROW_HEIGHT + total_gap_height + 6)  # +6 for legend and spacing

    # Create figure with GridSpec for true separation
    height_ratios = []
    for i, h in enumerate(category_heights):
        height_ratios.append(h)
        if i < len(category_heights) - 1:
            height_ratios.append(GAP_FRACTION * 3)

    # Add legend row at the end
    height_ratios.append(5)  # Legend height (increased to prevent overlap)

    n_grid_rows = len(height_ratios)

    fig = plt.figure(figsize=(14, fig_height))
    gs = GridSpec(n_grid_rows, 2, width_ratios=[10, 0.4], height_ratios=height_ratios,
                  hspace=0, wspace=0.02)

    vmax = CONFIG['vmax']
    renderer = DotplotRenderer(cmap=CMAP_DIVERGING, vmax=vmax,
                               min_dot_size=15, max_dot_size=200)

    col_labels = ['Early\nG32A', 'TrajDev\nG32A', 'Late\nG32A',
                  'Early\nR403C', 'TrajDev\nR403C', 'Late\nR403C']

    axes = []
    grid_row = 0

    for cat_idx, (cat, cat_df, color) in enumerate(category_data):
        # Create axis for this category block
        ax = fig.add_subplot(gs[grid_row, 0])
        axes.append(ax)

        # Prepare data for this category
        nes_matrix = cat_df[nes_cols].values
        padj_matrix = cat_df[padj_cols].values
        pathway_names = cat_df['Description'].values
        n_rows = len(cat_df)

        # Render dotplot
        renderer.render(ax, nes_matrix, padj_matrix)

        # Thick vertical separator between mutations
        ax.axvline(2.5, color='white', linewidth=4, zorder=1)
        ax.axvline(2.5, color='black', linewidth=1.5, zorder=2)

        # Y-axis: pathway names
        ax.set_yticks(range(n_rows))
        ax.set_yticklabels(pathway_names, fontsize=8)

        # X-axis: only show labels on the last category
        ax.set_xticks(range(6))
        if cat_idx == len(category_data) - 1:
            ax.set_xticklabels(col_labels, fontsize=9)
        else:
            ax.set_xticklabels([])

        # Category label at top-left
        ax.text(-0.02, 1.1, cat, fontsize=11, fontweight='bold',
                va='bottom', ha='left', color=color,
                transform=ax.transAxes)

        # Add box around this category block
        for spine in ax.spines.values():
            spine.set_linewidth(1.5)
            spine.set_color('black')

        # Move to next grid position
        grid_row += 1
        if cat_idx < len(category_data) - 1:
            # Create invisible gap axis
            ax_gap = fig.add_subplot(gs[grid_row, 0])
            ax_gap.axis('off')
            grid_row += 1

    # Colorbar (spans all rows except legend)
    ax_cbar = fig.add_subplot(gs[:-1, 1])
    norm = mcolors.Normalize(vmin=-vmax, vmax=vmax)
    sm = plt.cm.ScalarMappable(cmap=CMAP_DIVERGING, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, cax=ax_cbar)
    cbar.set_label('NES', fontsize=11)

    # Size legend (bottom row, spanning both columns)
    ax_legend = fig.add_subplot(gs[-1, :])
    add_size_legend_visual(ax_legend, renderer)

    # Title and legend
    fig.suptitle('Comprehensive Pathway Overview by Semantic Category (Dotplot)\n'
                 f'Filtered: ≥{CONFIG["min_data_points"]} trajectory points, at least one significant result',
                 fontsize=12, fontweight='bold', y=0.995)

    plt.tight_layout(rect=[0, 0.05, 1, 0.97])

    # Save
    for ext, dpi in [('pdf', CONFIG['dpi']), ('png', 150)]:
        output_file = OUTPUT_DIR / f'Fig4_Semantic_Pathway_Overview_dotplot.{ext}'
        fig.savefig(output_file, dpi=dpi, bbox_inches='tight', facecolor='white')
        print(f"  Saved: {output_file}")

    plt.close(fig)

    return df_plot


# =============================================================================
# MAIN EXECUTION
# =============================================================================

def main():
    """Generate all publication figures (dotplot version)"""
    print("="*80)
    print("GENERATING PUBLICATION-READY FIGURES (DOTPLOT VERSION)")
    print("="*80)
    print(f"Output directory: {OUTPUT_DIR}")
    print(f"CGP database: EXCLUDED")
    print(f"Visualization: Dotplot (color=NES, size=padj, edge=significance)")

    # Load data
    df = load_data()

    # Create figures
    print("\n" + "-"*80)
    syngo_ribo, mito_ribo, gobp_ribo = create_ribosome_paradox_figure(df)

    print("\n" + "-"*80)
    df_mito = create_mitocarta_trajectory_figure(df)

    print("\n" + "-"*80)
    df_syngo = create_syngo_trajectory_figure(df)

    print("\n" + "-"*80)
    df_patterns = create_pattern_summary_figure(df)

    print("\n" + "-"*80)
    df_semantic = create_semantic_grouping_figure(df)

    # Summary
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    print(f"Output directory: {OUTPUT_DIR}")
    print("\nGenerated files:")
    for f in sorted(OUTPUT_DIR.glob('*')):
        print(f"  - {f.name}")

    print("\n" + "="*80)
    print("COMPLETE!")
    print("="*80)


if __name__ == '__main__':
    main()
