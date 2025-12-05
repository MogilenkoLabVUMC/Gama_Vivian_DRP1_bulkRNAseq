#!/usr/bin/env python3
"""
Publication-Ready Trajectory Figures
=====================================

Creates consolidated, colorblind-safe figures for manuscript:

1. Ribosome Paradox - Combined SynGO + MitoCarta showing opposite ribosome fates
2. Shared Compensation - MitoCarta pathways rescued in both mutations
3. Pattern Classification Summary - Quantitative overview across all databases
4. Semantic Pathway Overview - Comprehensive heatmap by biological category

Features:
- Colorblind-safe Blue-White-Orange diverging palette
- All cells show NES values; Bold* for significant, normal text for non-significant
  (color naturally fades as NES approaches 0, preserving gradient legend validity)
- CGP database excluded (cancer-focused, not relevant for neurobiology)
- New semantic categories: Mitochondrial Dynamics, Neuronal Development
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
from matplotlib.patches import Patch, Rectangle
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import warnings
warnings.filterwarnings('ignore')

# Import from reusable toolkit
from GSEA_plotting_python import create_diverging_cmap, HeatmapRenderer, DIVERGING_COLORS

# Import from project-specific modules
from Python.config import CONFIG, resolve_path, ensure_output_dir, EXCLUDE_PATHWAYS
from Python.semantic_categories import (
    SEMANTIC_CATEGORY_ORDER, SEMANTIC_COLORS, PATTERN_COLORS, MUTATION_COLORS,
    assign_semantic_category, categorize_pathways, sort_by_category, deduplicate_pathways
)
from Python.data_loader import load_classified_pathways, filter_pathways
from Python.patterns import is_compensation, add_pattern_classification
from Python.pattern_definitions import MEANINGFUL_PATTERNS, ACTIVE_PATTERNS

# =============================================================================
# SETUP
# =============================================================================

# Create colormap
CMAP_DIVERGING = create_diverging_cmap()

# Ensure output directory exists
OUTPUT_DIR = ensure_output_dir()


# =============================================================================
# DATA LOADING (with CGP exclusion)
# =============================================================================

def load_data():
    """Load classified pathway data with CGP excluded"""
    print("Loading data (excluding CGP database)...")
    df = load_classified_pathways()  # CGP excluded by default
    return df


# =============================================================================
# FIGURE 1: RIBOSOME PARADOX
# =============================================================================

def create_ribosome_paradox_figure(df):
    """
    Create combined figure showing the Ribosome Paradox with 3 panels:
    - Panel A: Synaptic ribosomes (SynGO) - UP early, DOWN late [FAILURE]
    - Panel B: Mitochondrial ribosomes (MitoCarta) - DOWN early, UP recovery [PARTIAL RESCUE]
    - Panel C: Cytoplasmic Ribosome Biogenesis (GO:BP) - Shows compensatory context
    """
    print("\nCreating Ribosome Paradox figure (3 panels)...")

    # Extract ribosome-related pathways
    syngo_ribosomes = df[
        (df['database'] == 'SynGO') &
        (df['Description'].str.contains('ribosome', case=False, na=False))
    ].copy()

    mito_ribosomes = df[
        (df['database'] == 'MitoCarta') &
        (df['Description'].str.contains('ribosome|translation|central_dogma', case=False, na=False))
    ].copy()

    # Panel C: Cytoplasmic ribosome biogenesis from GO:BP (not mitochondrial)
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

    # Create figure with THREE panels
    fig = plt.figure(figsize=(12, 14))
    gs = GridSpec(3, 2, width_ratios=[6, 0.3], height_ratios=[1, 1, 1],
                  wspace=0.05, hspace=0.35)

    ax_syngo = fig.add_subplot(gs[0, 0])
    ax_mito = fig.add_subplot(gs[1, 0])
    ax_cyto = fig.add_subplot(gs[2, 0])
    ax_cbar = fig.add_subplot(gs[:, 1])

    vmax = CONFIG['vmax']
    renderer = HeatmapRenderer(cmap=CMAP_DIVERGING, vmax=vmax,
                               nonsig_style='show_values')

    # Panel A: Synaptic Ribosomes
    _plot_trajectory_heatmap(
        ax_syngo, syngo_ribosomes, renderer,
        title='A. Synaptic Ribosomes (SynGO)',
        subtitle='Early compensatory upregulation fails during maturation'
    )

    # Panel B: Mitochondrial Ribosomes
    _plot_trajectory_heatmap(
        ax_mito, mito_ribosomes, renderer,
        title='B. Mitochondrial Ribosomes (MitoCarta)',
        subtitle='Early crisis with partial recovery during maturation'
    )

    # Panel C: Cytoplasmic Ribosome Biogenesis
    _plot_trajectory_heatmap(
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

    # Add annotation explaining the paradox (positioned lower to avoid overlap)
    fig.text(0.5, -0.01,
             'The Ribosome Paradox:\n'
             'A) Synaptic ribosomes show early UP (compensation attempt) → late DOWN (failure)\n'
             'B) Mitochondrial ribosomes show early DOWN (crisis) → late UP (recovery)\n'
             'C) Cytoplasmic ribosome biogenesis provides context for compensatory responses\n\n'
             'Bold* = Significant (padj < 0.05)  |  Normal = Not significant',
             ha='center', va='top', fontsize=9, style='italic',
             bbox=dict(boxstyle='round', facecolor='#F5F5F5', edgecolor='gray', alpha=0.8))

    # Main title
    fig.suptitle('The Ribosome Paradox: Opposite Fates of Translation Machineries',
                fontsize=14, fontweight='bold', y=0.99)

    plt.tight_layout(rect=[0, 0.10, 0.95, 0.97])

    # Save PDF and PNG
    for ext, dpi in [('pdf', CONFIG['dpi']), ('png', 150)]:
        output_file = OUTPUT_DIR / f'Fig1_Ribosome_Paradox.{ext}'
        fig.savefig(output_file, dpi=dpi, bbox_inches='tight', facecolor='white')
        print(f"  Saved: {output_file}")

    plt.close(fig)

    return syngo_ribosomes, mito_ribosomes, gobp_ribosomes


def _plot_trajectory_heatmap(ax, df_subset, renderer, title, subtitle):
    """
    Helper to plot a single trajectory heatmap panel.

    Shows all NES values with consistent coloring. Significant cells marked with *.
    No hatching - non-significant values shown with normal text.
    """
    if len(df_subset) == 0:
        ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
        ax.set_title(title, fontsize=12, fontweight='bold')
        return

    # Prepare data matrix: rows = pathways, cols = [Early_G32A, TrajDev_G32A, Late_G32A, Early_R403C, ...]
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
    sig_matrix = (df_subset[padj_cols].values < CONFIG['padj_cutoff'])

    # Plot heatmap
    im = ax.imshow(nes_matrix, cmap=CMAP_DIVERGING, aspect='auto',
                   vmin=-renderer.vmax, vmax=renderer.vmax, interpolation='nearest')

    # Add thick vertical line separating mutations with small gap effect
    ax.axvline(2.5, color='white', linewidth=6)  # White gap
    ax.axvline(2.5, color='black', linewidth=2)  # Black line on top

    # NO hatching - show values for all cells with consistent coloring
    # Add NES values with significance markers (bold* for significant, normal for non-sig)
    renderer.annotate_cells(ax, nes_matrix, sig_matrix)

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
# FIGURE 2: MitoCarta TRAJECTORY PATTERNS (all pathways, no subclassification)
# =============================================================================

def create_mitocarta_trajectory_figure(df):
    """
    Create single-panel figure showing ALL MitoCarta pathways with any significance.

    No subclassification - each mutation may follow different trajectories.
    Shows all pathways that have significance in at least one contrast for at least one mutant.
    """
    print("\nCreating MitoCarta Trajectory Patterns figure (all pathways)...")

    # Filter for MitoCarta
    df_mito = df[df['database'] == 'MitoCarta'].copy()
    print(f"  Total MitoCarta pathways: {len(df_mito)}")

    # Filter for pathways with at least one significant result in any contrast
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

    # Create figure with single panel
    fig_height = max(10, n_pathways * 0.4)
    fig = plt.figure(figsize=(14, fig_height))
    gs = GridSpec(1, 2, width_ratios=[6, 0.3], wspace=0.05)

    ax_main = fig.add_subplot(gs[0, 0])
    ax_cbar = fig.add_subplot(gs[0, 1])

    vmax = CONFIG['vmax']
    renderer = HeatmapRenderer(cmap=CMAP_DIVERGING, vmax=vmax, font_size=8,
                               nonsig_style='show_values')

    # Plot all pathways
    _plot_trajectory_panel(
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

    # Main title
    fig.suptitle('MitoCarta Pathway Trajectories\nG32A and R403C may show different trajectory patterns',
                fontsize=14, fontweight='bold', y=0.99)

    # Legend
    fig.text(0.5, 0.01,
             'Bold* = Significant (padj < 0.05)  |  Normal = Not significant\n'
             'Note: Each mutation follows its own trajectory - no forced classification',
             ha='center', va='bottom', fontsize=9,
             bbox=dict(boxstyle='round', facecolor='#F5F5F5', edgecolor='gray', alpha=0.8))

    plt.tight_layout(rect=[0, 0.05, 0.95, 0.97])

    # Save
    for ext, dpi in [('pdf', CONFIG['dpi']), ('png', 150)]:
        output_file = OUTPUT_DIR / f'Fig2_MitoCarta_Trajectory_Patterns.{ext}'
        fig.savefig(output_file, dpi=dpi, bbox_inches='tight', facecolor='white')
        print(f"  Saved: {output_file}")

    plt.close(fig)

    return df_mito


def _plot_trajectory_panel(ax, df_subset, renderer, title, subtitle, color=None):
    """
    Helper to plot a single trajectory pattern panel.

    Shows all NES values with consistent coloring. No hatching.
    Significant cells marked with bold*.
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
    sig_matrix = (df_subset[padj_cols].values < CONFIG['padj_cutoff'])

    # Plot heatmap
    im = ax.imshow(nes_matrix, cmap=CMAP_DIVERGING, aspect='auto',
                   vmin=-renderer.vmax, vmax=renderer.vmax, interpolation='nearest')

    # Thick vertical separator between mutations with gap effect
    ax.axvline(2.5, color='white', linewidth=6)  # White gap
    ax.axvline(2.5, color='black', linewidth=2)  # Black line on top

    # NO hatching - show values for all cells with consistent coloring
    # Add NES values with significance markers (bold* for significant, normal for non-sig)
    renderer.annotate_cells(ax, nes_matrix, sig_matrix)

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
# FIGURE 3b: SynGO TRAJECTORY PATTERNS (all pathways, no subclassification)
# =============================================================================

def create_syngo_trajectory_figure(df):
    """
    Create single-panel figure showing ALL SynGO pathways with any significance.

    No subclassification - each mutation may follow different trajectories.
    Shows all pathways that have significance in at least one contrast for at least one mutant.
    """
    print("\nCreating SynGO Trajectory Patterns figure (all pathways)...")

    # Filter for SynGO
    df_syngo = df[df['database'] == 'SynGO'].copy()
    print(f"  Total SynGO pathways: {len(df_syngo)}")

    # Filter for pathways with at least one significant result in any contrast
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

    # Create figure with single panel
    fig_height = max(10, n_pathways * 0.4)
    fig = plt.figure(figsize=(14, fig_height))
    gs = GridSpec(1, 2, width_ratios=[6, 0.3], wspace=0.05)

    ax_main = fig.add_subplot(gs[0, 0])
    ax_cbar = fig.add_subplot(gs[0, 1])

    vmax = CONFIG['vmax']
    renderer = HeatmapRenderer(cmap=CMAP_DIVERGING, vmax=vmax, font_size=8,
                               nonsig_style='show_values')

    # Plot all pathways
    _plot_trajectory_panel(
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

    # Main title
    fig.suptitle('SynGO Synaptic Pathway Trajectories\nG32A and R403C may show different trajectory patterns',
                fontsize=14, fontweight='bold', y=0.99)

    # Legend
    fig.text(0.5, 0.01,
             'Bold* = Significant (padj < 0.05)  |  Normal = Not significant\n'
             'Note: Each mutation follows its own trajectory - no forced classification',
             ha='center', va='bottom', fontsize=9,
             bbox=dict(boxstyle='round', facecolor='#F5F5F5', edgecolor='gray', alpha=0.8))

    plt.tight_layout(rect=[0, 0.05, 0.95, 0.97])

    # Save
    for ext, dpi in [('pdf', CONFIG['dpi']), ('png', 150)]:
        output_file = OUTPUT_DIR / f'Fig3b_SynGO_Trajectory_Patterns.{ext}'
        fig.savefig(output_file, dpi=dpi, bbox_inches='tight', facecolor='white')
        print(f"  Saved: {output_file}")

    plt.close(fig)

    return df_syngo


# =============================================================================
# FIGURE 3: PATTERN CLASSIFICATION SUMMARY
# =============================================================================

def create_pattern_summary_figure(df):
    """
    Create summary bar chart showing pattern classification across all databases.

    Fixed issues:
    - Mutation labels positioned next to bars, not far right
    - Legend includes ALL patterns present (not just G32A patterns)
    - Clear visual distinction between G32A and R403C bars
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
                # Don't add label here - we'll create legend manually
                ax1.barh(x + offset, pivot[pattern], width, left=bottom,
                        color=color, alpha=0.85)
                bottom += pivot[pattern].values

        # FIX: Add mutation labels next to each bar (at the end of each database's bar)
        for i, db in enumerate(databases):
            bar_total = pivot.loc[db].sum() if db in pivot.index else 0
            if bar_total > 0:
                ax1.text(bar_total + 1, i + offset, mutation[:4],
                        fontsize=8, fontweight='bold', va='center',
                        color=MUTATION_COLORS[mutation])

    ax1.set_yticks(x)
    ax1.set_yticklabels(databases, fontsize=10)
    ax1.set_xlabel('Number of Pathways', fontsize=11)
    ax1.set_title('A. Pattern Distribution by Database', fontsize=12, fontweight='bold')

    # FIX: Create legend handles manually to include ALL patterns present
    # Position legend at center-right of the plot area
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
                      label=mutation, color=MUTATION_COLORS[mutation], alpha=0.85)

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

    # Save
    for ext, dpi in [('pdf', CONFIG['dpi']), ('png', 150)]:
        output_file = OUTPUT_DIR / f'Fig3_Pattern_Classification_Summary.{ext}'
        fig.savefig(output_file, dpi=dpi, bbox_inches='tight', facecolor='white')
        print(f"  Saved: {output_file}")

    plt.close(fig)

    return df_counts


# =============================================================================
# FIGURE 4: SEMANTIC PATHWAY OVERVIEW
# =============================================================================

def create_semantic_grouping_figure(df):
    """
    Create comprehensive semantic-grouped pathway heatmap.

    Features:
    - Semantic categories with TRUE whitespace breaks between groups
    - Category labels placed at top-left of each group (above first pathway)
    - Filter: Include pathways with data in ≥3 of 6 columns
    - Excludes irrelevant pathways (cilium, sperm, cancer, viral, cytoskeletal)
    - Deduplicates pathways across databases
    """
    print("\nCreating Semantic Pathway Overview figure...")

    # Prepare pathway data with NES columns
    nes_cols = CONFIG['nes_cols']
    padj_cols = CONFIG['padj_cols']

    # Filter pathways by data points
    df_filtered = filter_pathways(df, min_data_points=CONFIG['min_data_points'])
    print(f"  Pathways with ≥{CONFIG['min_data_points']} data points: {len(df_filtered)}")

    # Assign semantic categories WITH exclusion list
    # This removes irrelevant pathways (cilium, sperm, cancer, viral, cytoskeletal)
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

    # =========================================================================
    # Group pathways by category for separate heatmap blocks
    # =========================================================================
    category_data = []  # List of (category_name, df_subset, color)

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

    # =========================================================================
    # Calculate heights for each category block
    # =========================================================================
    GAP_FRACTION = 0.8   # Increased for true whitespace separation between categories
    ROW_HEIGHT = 0.35    # Height per pathway row in inches

    category_heights = []
    total_pathways = 0
    for cat, cat_df, color in category_data:
        n_rows = len(cat_df)
        category_heights.append(n_rows)
        total_pathways += n_rows

    # Total figure height
    total_gap_height = (n_categories - 1) * GAP_FRACTION * ROW_HEIGHT * 3  # Gap space
    fig_height = max(10, total_pathways * ROW_HEIGHT + total_gap_height + 2)  # +2 for title/legend

    # =========================================================================
    # Create figure with GridSpec for true separation
    # =========================================================================
    # Height ratios: category heights with gaps between them
    height_ratios = []
    for i, h in enumerate(category_heights):
        height_ratios.append(h)
        if i < len(category_heights) - 1:
            height_ratios.append(GAP_FRACTION * 3)  # Gap row

    n_grid_rows = len(height_ratios)

    fig = plt.figure(figsize=(14, fig_height))
    gs = GridSpec(n_grid_rows, 2, width_ratios=[10, 0.4], height_ratios=height_ratios,
                  hspace=0, wspace=0.02)

    vmax = CONFIG['vmax']
    renderer = HeatmapRenderer(cmap=CMAP_DIVERGING, vmax=vmax, font_size=7,
                               nonsig_style='show_values')

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
        sig_matrix = (cat_df[padj_cols].values < CONFIG['padj_cutoff'])
        pathway_names = cat_df['Description'].values
        n_rows = len(cat_df)

        # Plot heatmap
        masked_nes = np.ma.masked_invalid(nes_matrix)
        im = ax.imshow(masked_nes, cmap=CMAP_DIVERGING, aspect='auto',
                       vmin=-vmax, vmax=vmax, interpolation='nearest')

        # Add NES values with significance markers
        renderer.annotate_cells(ax, nes_matrix, sig_matrix)

        # Thick vertical separator between mutations
        ax.axvline(2.5, color='white', linewidth=4, zorder=3)
        ax.axvline(2.5, color='black', linewidth=1.5, zorder=4)

        # Y-axis: pathway names
        ax.set_yticks(range(n_rows))
        ax.set_yticklabels(pathway_names, fontsize=8)

        # X-axis: only show labels on the last category
        ax.set_xticks(range(6))
        if cat_idx == len(category_data) - 1:
            ax.set_xticklabels(col_labels, fontsize=9)
        else:
            ax.set_xticklabels([])

        # Category label at top-left (above first pathway)
        # Positioned to the left (-0.2) and above (1.1) to clear pathway names
        ax.text(-0.02, 1.1, cat, fontsize=11, fontweight='bold',
                va='bottom', ha='left', color=color,
                transform=ax.transAxes)

        # Add box around this category block
        for spine in ax.spines.values():
            spine.set_linewidth(1.5)
            spine.set_color('black')

        # Move to next grid position (skip gap row)
        grid_row += 1
        if cat_idx < len(category_data) - 1:
            # Create invisible gap axis
            ax_gap = fig.add_subplot(gs[grid_row, 0])
            ax_gap.axis('off')
            grid_row += 1

    # =========================================================================
    # Colorbar (spans all rows)
    # =========================================================================
    ax_cbar = fig.add_subplot(gs[:, 1])
    norm = mcolors.Normalize(vmin=-vmax, vmax=vmax)
    sm = plt.cm.ScalarMappable(cmap=CMAP_DIVERGING, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, cax=ax_cbar)
    cbar.set_label('NES', fontsize=11)

    # =========================================================================
    # Title and legend
    # =========================================================================
    fig.suptitle('Comprehensive Pathway Overview by Semantic Category\n'
                 f'Filtered: ≥{CONFIG["min_data_points"]} trajectory points, at least one significant result',
                 fontsize=12, fontweight='bold', y=0.995)

    # Legend at bottom
    fig.text(0.5, 0.002, 'Bold* = Significant (padj < 0.05)  |  Normal = Not significant',
             ha='center', va='bottom', fontsize=9,
             bbox=dict(boxstyle='round', facecolor='#F5F5F5', edgecolor='gray', alpha=0.8))

    plt.tight_layout(rect=[0, 0.02, 1, 0.97])

    # Save
    for ext, dpi in [('pdf', CONFIG['dpi']), ('png', 150)]:
        output_file = OUTPUT_DIR / f'Fig4_Semantic_Pathway_Overview.{ext}'
        fig.savefig(output_file, dpi=dpi, bbox_inches='tight', facecolor='white')
        print(f"  Saved: {output_file}")

    plt.close(fig)

    return df_plot


# =============================================================================
# SUPPLEMENTARY: COLOR PALETTE REFERENCE
# =============================================================================

def create_palette_reference():
    """Create a reference figure showing the colorblind-safe palettes used"""
    print("\nCreating color palette reference...")

    fig, axes = plt.subplots(1, 3, figsize=(15, 4))

    # Panel 1: Diverging colormap
    ax1 = axes[0]
    gradient = np.linspace(-3, 3, 256).reshape(1, -1)
    ax1.imshow(gradient, cmap=CMAP_DIVERGING, aspect='auto', extent=[-3, 3, 0, 1])
    ax1.set_yticks([])
    ax1.set_xlabel('NES Value', fontsize=11)
    ax1.set_title('Diverging Palette (NES)\nBlue (down) - White - Orange (up)',
                  fontsize=11, fontweight='bold')

    ax1.text(-3, -0.3, f'Down: {DIVERGING_COLORS["negative"]}', fontsize=9, ha='left')
    ax1.text(0, -0.3, f'Neutral: {DIVERGING_COLORS["neutral"]}', fontsize=9, ha='center')
    ax1.text(3, -0.3, f'Up: {DIVERGING_COLORS["positive"]}', fontsize=9, ha='right')

    # Panel 2: Pattern colors
    ax2 = axes[1]
    patterns = list(PATTERN_COLORS.keys())
    colors = list(PATTERN_COLORS.values())
    y_pos = np.arange(len(patterns))

    bars = ax2.barh(y_pos, [1]*len(patterns), color=colors, edgecolor='black', linewidth=0.5)
    ax2.set_yticks(y_pos)
    ax2.set_yticklabels(patterns, fontsize=10)
    ax2.set_xlim(0, 1.5)
    ax2.set_xticks([])
    ax2.set_title('Pattern Categories\n(Colorblind-safe)', fontsize=11, fontweight='bold')

    for i, (pattern, color) in enumerate(PATTERN_COLORS.items()):
        ax2.text(1.05, i, color, fontsize=9, va='center', family='monospace')

    # Panel 3: Mutation colors
    ax3 = axes[2]
    mutations = list(MUTATION_COLORS.keys())
    mut_colors = list(MUTATION_COLORS.values())

    bars = ax3.barh([0, 1], [1, 1], color=mut_colors, edgecolor='black', linewidth=0.5)
    ax3.set_yticks([0, 1])
    ax3.set_yticklabels(mutations, fontsize=12, fontweight='bold')
    ax3.set_xlim(0, 1.5)
    ax3.set_xticks([])
    ax3.set_title('Mutation Colors\n(Colorblind-safe)', fontsize=11, fontweight='bold')

    for i, (mut, color) in enumerate(MUTATION_COLORS.items()):
        ax3.text(1.05, i, color, fontsize=10, va='center', family='monospace')

    fig.suptitle('Colorblind-Safe Palette Reference', fontsize=14, fontweight='bold')
    plt.tight_layout()

    for ext, dpi in [('pdf', CONFIG['dpi']), ('png', 150)]:
        output_file = OUTPUT_DIR / f'Supplementary_Color_Palette_Reference.{ext}'
        fig.savefig(output_file, dpi=dpi, bbox_inches='tight', facecolor='white')
        print(f"  Saved: {output_file}")

    plt.close(fig)


# =============================================================================
# MAIN EXECUTION
# =============================================================================

def main():
    """Generate all publication figures"""
    print("="*80)
    print("GENERATING PUBLICATION-READY FIGURES")
    print("="*80)
    print(f"Output directory: {OUTPUT_DIR}")
    print(f"CGP database: EXCLUDED")
    print(f"Significance style: Bold* for padj < 0.05, normal text for non-significant")

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

    print("\n" + "-"*80)
    create_palette_reference()

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
