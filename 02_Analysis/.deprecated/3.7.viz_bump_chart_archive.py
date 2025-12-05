#!/usr/bin/env python3
"""
Bump Chart / Slope Graph Visualization
======================================

Creates bump charts (slope graphs) showing individual pathway trajectories
through the Early → TrajDev → Late stages.

Unlike alluvial diagrams that show categorical flows, bump charts show:
- Y-axis: NES value (continuous) or Rank (ordinal)
- X-axis: Stage (Early, TrajDev, Late)
- Each line: One pathway's trajectory
- Color: Pattern classification

Generates 6 variants to empirically determine best visualization:
- 2 Y-axis types (NES, Rank) × 3 scopes (Focused, Significant, All)

Author: Claude Code
Date: 2025-11-27
"""

import sys
from pathlib import Path

# Add module paths
sys.path.insert(0, str(Path(__file__).parent.parent / '01_Scripts'))

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import warnings
warnings.filterwarnings('ignore')

# Import project modules
from Python.config import resolve_path
from Python.pattern_definitions import get_pattern_colors, PADJ_SIGNIFICANT, NES_EFFECT
from Python.semantic_categories import (
    is_relevant_for_highlight,
    get_highlight_priority,
    assign_semantic_category,
    PRIORITY_CATEGORIES_FOR_HIGHLIGHT,
    MAX_OTHER_CATEGORY_HIGHLIGHTS,
)

# =============================================================================
# CONFIGURATION
# =============================================================================

OUTPUT_DIR = resolve_path('03_Results/02_Analysis/Plots/Trajectory_Flow')
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Pattern colors (colorblind-safe)
PATTERN_COLORS = get_pattern_colors()

# Mutation colors
MUTATION_COLORS = {
    'G32A': '#0072B2',    # Blue
    'R403C': '#D55E00',   # Vermillion
}

# Define scope filters
SCOPE_CONFIGS = {
    'focused': {
        'name': 'Focused (MitoCarta + SynGO)',
        'databases': ['MitoCarta', 'SynGO'],
        'expected_n': '~100',
    },
    'significant': {
        'name': 'Significant (p<0.05 any stage)',
        'databases': None,  # All databases
        'expected_n': '~2-4k',
    },
    'all': {
        'name': 'All Pathways',
        'databases': None,
        'expected_n': '~12k',
    },
}

# Alpha values for different scopes (less alpha for denser plots)
ALPHA_BY_SCOPE = {
    'focused': 0.7,
    'significant': 0.15,
    'all': 0.05,
}

# Line width by scope
LINEWIDTH_BY_SCOPE = {
    'focused': 1.5,
    'significant': 0.5,
    'all': 0.3,
}


# =============================================================================
# DATA LOADING
# =============================================================================

def load_data():
    """Load pathway data from master GSEA table."""
    print("Loading pathway data...")

    master_file = resolve_path('03_Results/02_Analysis/master_gsea_table.csv')
    df = pd.read_csv(master_file)

    # Get required columns
    id_cols = ['pathway_id', 'database', 'Description']
    nes_cols = [
        'NES_Early_G32A', 'NES_TrajDev_G32A', 'NES_Late_G32A',
        'NES_Early_R403C', 'NES_TrajDev_R403C', 'NES_Late_R403C'
    ]
    padj_cols = [
        'p.adjust_Early_G32A', 'p.adjust_TrajDev_G32A', 'p.adjust_Late_G32A',
        'p.adjust_Early_R403C', 'p.adjust_TrajDev_R403C', 'p.adjust_Late_R403C'
    ]
    pattern_cols = ['Pattern_G32A', 'Pattern_R403C', 'Confidence_G32A', 'Confidence_R403C']

    keep_cols = id_cols + nes_cols + padj_cols + pattern_cols
    available_cols = [c for c in keep_cols if c in df.columns]

    # Deduplicate to unique pathways
    df_unique = df[available_cols].drop_duplicates(subset=['pathway_id', 'database']).copy()

    print(f"  Loaded {len(df_unique)} unique pathways")
    print(f"  Databases: {df_unique['database'].nunique()}")

    return df_unique


def filter_by_scope(df, scope):
    """Filter dataframe based on scope configuration."""
    config = SCOPE_CONFIGS[scope]

    df_filtered = df.copy()

    # Filter by database if specified
    if config['databases']:
        df_filtered = df_filtered[df_filtered['database'].isin(config['databases'])]

    # For 'significant' scope, filter to pathways with meaningful patterns
    # (patterns other than Complex/Insufficient_data require significance)
    if scope == 'significant':
        # Significant = has a meaningful pattern in either mutation
        meaningful_patterns = ['Compensation', 'Progressive', 'Natural_improvement',
                               'Natural_worsening', 'Late_onset', 'Transient']

        sig_mask = pd.Series(False, index=df_filtered.index)

        for mutation in ['G32A', 'R403C']:
            pattern_col = f'Pattern_{mutation}'
            if pattern_col in df_filtered.columns:
                sig_mask |= df_filtered[pattern_col].isin(meaningful_patterns)

        df_filtered = df_filtered[sig_mask]

    # Remove rows with all-NA NES values (for either mutation)
    nes_cols_g32a = ['NES_Early_G32A', 'NES_TrajDev_G32A', 'NES_Late_G32A']
    nes_cols_r403c = ['NES_Early_R403C', 'NES_TrajDev_R403C', 'NES_Late_R403C']

    has_g32a = df_filtered[nes_cols_g32a].notna().any(axis=1)
    has_r403c = df_filtered[nes_cols_r403c].notna().any(axis=1)

    df_filtered = df_filtered[has_g32a | has_r403c]

    print(f"  Scope '{scope}': {len(df_filtered)} pathways")

    return df_filtered


# =============================================================================
# BUMP CHART FUNCTIONS
# =============================================================================

def create_bump_chart_nes(df, mutation, ax, scope='focused'):
    """
    Create bump chart using NES values on Y-axis.

    Parameters
    ----------
    df : pd.DataFrame
        Filtered pathway data
    mutation : str
        'G32A' or 'R403C'
    ax : matplotlib.axes.Axes
        Axes to plot on
    scope : str
        'focused', 'significant', or 'all' (affects visual parameters)
    """
    stages = ['Early', 'TrajDev', 'Late']
    nes_cols = [f'NES_{s}_{mutation}' for s in stages]
    pattern_col = f'Pattern_{mutation}'

    # Check columns exist
    if pattern_col not in df.columns:
        print(f"  Warning: {pattern_col} not found")
        return

    # Filter to rows with complete NES data for this mutation
    df_plot = df[df[nes_cols].notna().all(axis=1)].copy()

    if len(df_plot) == 0:
        ax.text(0.5, 0.5, f'No complete data for {mutation}',
                ha='center', va='center', transform=ax.transAxes)
        return

    # Get visual parameters for scope
    alpha = ALPHA_BY_SCOPE.get(scope, 0.5)
    linewidth = LINEWIDTH_BY_SCOPE.get(scope, 1.0)

    # Plot each pathway trajectory
    x_positions = [0, 1, 2]

    for _, row in df_plot.iterrows():
        nes_values = [row[col] for col in nes_cols]
        pattern = row[pattern_col]
        color = PATTERN_COLORS.get(pattern, '#999999')

        ax.plot(x_positions, nes_values,
                color=color, alpha=alpha, linewidth=linewidth,
                solid_capstyle='round')

    # Add zero line
    ax.axhline(y=0, color='black', linewidth=0.5, linestyle='-', alpha=0.3)

    # Add threshold lines
    ax.axhline(y=NES_EFFECT, color='gray', linewidth=0.5, linestyle='--', alpha=0.3)
    ax.axhline(y=-NES_EFFECT, color='gray', linewidth=0.5, linestyle='--', alpha=0.3)

    # Formatting
    ax.set_xticks(x_positions)
    ax.set_xticklabels(['Early\n(D35)', 'TrajDev\n(Maturation)', 'Late\n(D65)'], fontsize=10)
    ax.set_xlim(-0.3, 2.3)

    # Y-axis limits based on data
    y_max = max(4, df_plot[nes_cols].max().max() * 1.1)
    y_min = min(-4, df_plot[nes_cols].min().min() * 1.1)
    ax.set_ylim(y_min, y_max)

    ax.set_title(f'{mutation}\n(n={len(df_plot)})',
                 fontsize=12, fontweight='bold',
                 color=MUTATION_COLORS.get(mutation, 'black'))

    # Grid
    ax.grid(axis='y', alpha=0.3, linestyle='-', linewidth=0.5)
    ax.set_axisbelow(True)


def create_bump_chart_rank(df, mutation, ax, scope='focused'):
    """
    Create bump chart using rank on Y-axis.

    Each stage has its own ranking (1 = most negative NES, N = most positive).
    Lines connect ranks across stages.

    Parameters
    ----------
    df : pd.DataFrame
        Filtered pathway data
    mutation : str
        'G32A' or 'R403C'
    ax : matplotlib.axes.Axes
        Axes to plot on
    scope : str
        'focused', 'significant', or 'all' (affects visual parameters)
    """
    stages = ['Early', 'TrajDev', 'Late']
    nes_cols = [f'NES_{s}_{mutation}' for s in stages]
    pattern_col = f'Pattern_{mutation}'

    # Check columns exist
    if pattern_col not in df.columns:
        print(f"  Warning: {pattern_col} not found")
        return

    # Filter to rows with complete NES data for this mutation
    df_plot = df[df[nes_cols].notna().all(axis=1)].copy()

    if len(df_plot) == 0:
        ax.text(0.5, 0.5, f'No complete data for {mutation}',
                ha='center', va='center', transform=ax.transAxes)
        return

    # Compute ranks for each stage
    for col in nes_cols:
        rank_col = col.replace('NES_', 'Rank_')
        df_plot[rank_col] = df_plot[col].rank(ascending=True, method='first')

    rank_cols = [col.replace('NES_', 'Rank_') for col in nes_cols]

    # Get visual parameters for scope
    alpha = ALPHA_BY_SCOPE.get(scope, 0.5)
    linewidth = LINEWIDTH_BY_SCOPE.get(scope, 1.0)

    # Plot each pathway trajectory
    x_positions = [0, 1, 2]
    n_pathways = len(df_plot)

    for _, row in df_plot.iterrows():
        rank_values = [row[col] for col in rank_cols]
        pattern = row[pattern_col]
        color = PATTERN_COLORS.get(pattern, '#999999')

        ax.plot(x_positions, rank_values,
                color=color, alpha=alpha, linewidth=linewidth,
                solid_capstyle='round')

    # Add middle line
    mid_rank = n_pathways / 2
    ax.axhline(y=mid_rank, color='black', linewidth=0.5, linestyle='-', alpha=0.3)

    # Formatting
    ax.set_xticks(x_positions)
    ax.set_xticklabels(['Early\n(D35)', 'TrajDev\n(Maturation)', 'Late\n(D65)'], fontsize=10)
    ax.set_xlim(-0.3, 2.3)

    # Y-axis: ranks
    ax.set_ylim(0, n_pathways + 1)
    ax.invert_yaxis()  # Rank 1 at top

    ax.set_title(f'{mutation}\n(n={len(df_plot)})',
                 fontsize=12, fontweight='bold',
                 color=MUTATION_COLORS.get(mutation, 'black'))

    # Grid
    ax.grid(axis='y', alpha=0.3, linestyle='-', linewidth=0.5)
    ax.set_axisbelow(True)


def create_legend(ax):
    """Create legend showing pattern colors."""
    # Define patterns to show in legend (excluding Complex and Insufficient_data for clarity)
    legend_patterns = [
        ('Compensation', '#009E73', 'Early defect + TrajDev opposes + Late improved'),
        ('Natural_improvement', '#56B4E9', 'Early defect + passive recovery'),
        ('Progressive', '#D55E00', 'Early defect + TrajDev amplifies + Late worsened'),
        ('Natural_worsening', '#E69F00', 'Early defect + passive worsening'),
        ('Late_onset', '#CC79A7', 'No early defect + Late emerges'),
        ('Transient', '#0072B2', 'Strong early + fully resolved'),
        ('Complex', '#F0E442', 'Multiphasic or inconsistent'),
    ]

    handles = []
    for pattern, color, description in legend_patterns:
        patch = mpatches.Patch(color=color, label=f'{pattern}')
        handles.append(patch)

    ax.legend(handles=handles, loc='center left', fontsize=9,
              title='Pattern', title_fontsize=10, framealpha=0.9)
    ax.axis('off')


def create_faceted_bump_chart(df, y_type='nes', scope='focused'):
    """
    Create side-by-side bump charts for G32A and R403C.

    Parameters
    ----------
    df : pd.DataFrame
        Pathway data
    y_type : str
        'nes' for NES values, 'rank' for rank positions
    scope : str
        'focused', 'significant', or 'all'

    Returns
    -------
    matplotlib.figure.Figure
    """
    # Filter data by scope
    df_scope = filter_by_scope(df, scope)

    if len(df_scope) == 0:
        print(f"  Warning: No data for scope '{scope}'")
        return None

    # Create figure with 3 columns: G32A, R403C, Legend
    fig, axes = plt.subplots(1, 3, figsize=(14, 8),
                              gridspec_kw={'width_ratios': [3, 3, 1.5]})

    # Select plotting function
    plot_fn = create_bump_chart_nes if y_type == 'nes' else create_bump_chart_rank

    # Plot G32A
    plot_fn(df_scope, 'G32A', axes[0], scope=scope)
    axes[0].set_ylabel('NES' if y_type == 'nes' else 'Rank', fontsize=11)

    # Plot R403C
    plot_fn(df_scope, 'R403C', axes[1], scope=scope)
    # Share y-axis limits for NES
    if y_type == 'nes':
        y_min = min(axes[0].get_ylim()[0], axes[1].get_ylim()[0])
        y_max = max(axes[0].get_ylim()[1], axes[1].get_ylim()[1])
        axes[0].set_ylim(y_min, y_max)
        axes[1].set_ylim(y_min, y_max)

    # Legend
    create_legend(axes[2])

    # Main title
    scope_name = SCOPE_CONFIGS[scope]['name']
    y_label = 'NES Value' if y_type == 'nes' else 'Rank (by NES)'
    fig.suptitle(f'Pathway Trajectory Bump Chart\n{scope_name} | Y-axis: {y_label}',
                 fontsize=14, fontweight='bold', y=0.98)

    plt.tight_layout(rect=[0, 0, 1, 0.95])

    return fig


# =============================================================================
# WEIGHTED BUMP CHART FUNCTIONS (Pass 2 and Pass 3)
# =============================================================================

def compute_line_weights(df, mutation):
    """
    Compute line thickness weights based on pattern frequency.
    Rarer patterns get thicker lines and are rendered on top.

    Parameters
    ----------
    df : pd.DataFrame
        Filtered pathway data
    mutation : str
        'G32A' or 'R403C'

    Returns
    -------
    dict
        {pattern: {'linewidth': float, 'alpha': float, 'zorder': int}}
    """
    pattern_col = f'Pattern_{mutation}'

    if pattern_col not in df.columns:
        return {}

    counts = df[pattern_col].value_counts()
    total = len(df)

    weights = {}
    for pattern, count in counts.items():
        frequency = count / total

        if frequency > 0.3:  # Dominant (Complex, sometimes Compensation)
            # Background: thin, transparent, low z-order
            weights[pattern] = {'linewidth': 0.3, 'alpha': 0.15, 'zorder': 1}
        elif frequency > 0.1:  # Common (Natural_improvement, Compensation)
            weights[pattern] = {'linewidth': 0.8, 'alpha': 0.4, 'zorder': 2}
        elif frequency > 0.01:  # Uncommon (Natural_worsening)
            weights[pattern] = {'linewidth': 1.5, 'alpha': 0.7, 'zorder': 3}
        else:  # Rare (Late_onset, Transient, Progressive)
            # Foreground: thick, opaque, high z-order
            weights[pattern] = {'linewidth': 2.5, 'alpha': 0.9, 'zorder': 4}

    return weights


def add_pathway_labels(ax, df, mutation, weights, y_type='nes',
                       include_trajdev=True, max_labels=5):
    """
    Add pathway name labels for rare patterns at line endpoints.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to add labels to
    df : pd.DataFrame
        Pathway data (already filtered)
    mutation : str
        'G32A' or 'R403C'
    weights : dict
        Line weight dictionary from compute_line_weights
    y_type : str
        'nes' or 'rank'
    include_trajdev : bool
        Whether TrajDev is included (affects x position of labels)
    max_labels : int
        Maximum labels per pattern
    """
    pattern_col = f'Pattern_{mutation}'
    early_col = f'NES_Early_{mutation}'
    late_col = f'NES_Late_{mutation}'

    # Only label rare patterns (zorder >= 4)
    rare_patterns = [p for p, w in weights.items() if w.get('zorder', 0) >= 4]

    if not rare_patterns:
        return

    # X position for labels
    x_end = 2 if include_trajdev else 1

    for pattern in rare_patterns:
        df_pattern = df[df[pattern_col] == pattern].copy()

        if len(df_pattern) == 0:
            continue

        # Sort by |NES change| to get most dramatic trajectories
        df_pattern['abs_change'] = (df_pattern[late_col] - df_pattern[early_col]).abs()
        df_pattern = df_pattern.nlargest(min(max_labels, len(df_pattern)), 'abs_change')

        color = PATTERN_COLORS.get(pattern, '#999999')

        for i, (_, row) in enumerate(df_pattern.iterrows()):
            if y_type == 'nes':
                y_end = row[late_col]
            else:
                rank_col = f'Rank_Late_{mutation}'
                y_end = row.get(rank_col, row[late_col])

            # Short pathway name (truncate if needed)
            desc = str(row.get('Description', row.get('pathway_id', '')))
            label = desc[:35] + '...' if len(desc) > 35 else desc

            # Offset labels vertically to reduce overlap
            y_offset = (i - len(df_pattern)/2) * 0.15

            ax.annotate(
                label,
                xy=(x_end, y_end),
                xytext=(x_end + 0.12, y_end + y_offset),
                fontsize=6,
                color=color,
                ha='left', va='center',
                arrowprops=dict(arrowstyle='-', color=color, lw=0.5, alpha=0.5),
                zorder=10
            )


def create_bump_chart_weighted(df, mutation, ax, y_type='nes',
                                include_trajdev=True, scope='focused'):
    """
    Create weighted bump chart with inverse-frequency line weighting.

    Parameters
    ----------
    df : pd.DataFrame
        Filtered pathway data
    mutation : str
        'G32A' or 'R403C'
    ax : matplotlib.axes.Axes
        Axes to plot on
    y_type : str
        'nes' for NES values, 'rank' for rank positions
    include_trajdev : bool
        If True, 3-point (Early→TrajDev→Late)
        If False, 2-point (Early→Late direct)
    scope : str
        'focused', 'significant', or 'all'
    """
    pattern_col = f'Pattern_{mutation}'

    # Set up stages based on include_trajdev
    if include_trajdev:
        stages = ['Early', 'TrajDev', 'Late']
        x_positions = [0, 1, 2]
        x_labels = ['Early\n(D35)', 'TrajDev\n(Maturation)', 'Late\n(D65)']
    else:
        stages = ['Early', 'Late']
        x_positions = [0, 1]
        x_labels = ['Early\n(D35)', 'Late\n(D65)']

    nes_cols = [f'NES_{s}_{mutation}' for s in stages]

    # Check columns exist
    if pattern_col not in df.columns:
        print(f"  Warning: {pattern_col} not found")
        return

    # Filter to rows with complete NES data for this mutation
    df_plot = df[df[nes_cols].notna().all(axis=1)].copy()

    if len(df_plot) == 0:
        ax.text(0.5, 0.5, f'No complete data for {mutation}',
                ha='center', va='center', transform=ax.transAxes)
        return

    # Compute ranks if needed
    if y_type == 'rank':
        for col in nes_cols:
            rank_col = col.replace('NES_', 'Rank_')
            df_plot[rank_col] = df_plot[col].rank(ascending=True, method='first')

    # Compute line weights based on pattern frequency
    weights = compute_line_weights(df_plot, mutation)

    # Get patterns sorted by z-order (background first, rare last)
    pattern_order = sorted(
        df_plot[pattern_col].unique(),
        key=lambda p: weights.get(p, {}).get('zorder', 0)
    )

    # Plot each pattern group in z-order
    for pattern in pattern_order:
        df_pattern = df_plot[df_plot[pattern_col] == pattern]
        style = weights.get(pattern, {'linewidth': 1, 'alpha': 0.5, 'zorder': 1})
        color = PATTERN_COLORS.get(pattern, '#999999')

        for _, row in df_pattern.iterrows():
            if y_type == 'nes':
                y_values = [row[col] for col in nes_cols]
            else:
                rank_cols = [col.replace('NES_', 'Rank_') for col in nes_cols]
                y_values = [row[col] for col in rank_cols]

            ax.plot(x_positions, y_values,
                    color=color,
                    linewidth=style['linewidth'],
                    alpha=style['alpha'],
                    zorder=style['zorder'],
                    solid_capstyle='round')

    # Add pathway labels for rare patterns
    add_pathway_labels(ax, df_plot, mutation, weights, y_type=y_type,
                       include_trajdev=include_trajdev, max_labels=5)

    # Add reference lines
    if y_type == 'nes':
        ax.axhline(y=0, color='black', linewidth=0.5, linestyle='-', alpha=0.3, zorder=0)
        ax.axhline(y=NES_EFFECT, color='gray', linewidth=0.5, linestyle='--', alpha=0.3, zorder=0)
        ax.axhline(y=-NES_EFFECT, color='gray', linewidth=0.5, linestyle='--', alpha=0.3, zorder=0)

        # Y-axis limits
        y_max = max(4, df_plot[nes_cols].max().max() * 1.1)
        y_min = min(-4, df_plot[nes_cols].min().min() * 1.1)
        ax.set_ylim(y_min, y_max)
    else:
        n_pathways = len(df_plot)
        ax.axhline(y=n_pathways/2, color='black', linewidth=0.5, linestyle='-', alpha=0.3, zorder=0)
        ax.set_ylim(0, n_pathways + 1)
        ax.invert_yaxis()

    # Formatting
    ax.set_xticks(x_positions)
    ax.set_xticklabels(x_labels, fontsize=10)
    ax.set_xlim(-0.3, x_positions[-1] + 0.8)  # Extra space for labels

    ax.set_title(f'{mutation}\n(n={len(df_plot)})',
                 fontsize=12, fontweight='bold',
                 color=MUTATION_COLORS.get(mutation, 'black'))

    ax.grid(axis='y', alpha=0.3, linestyle='-', linewidth=0.5)
    ax.set_axisbelow(True)


def create_legend_weighted(ax, df, mutation):
    """Create legend showing pattern colors and line weights."""
    pattern_col = f'Pattern_{mutation}'

    if pattern_col not in df.columns:
        ax.axis('off')
        return

    weights = compute_line_weights(df, mutation)

    # Sort patterns by zorder for legend ordering
    patterns_sorted = sorted(
        weights.keys(),
        key=lambda p: weights.get(p, {}).get('zorder', 0),
        reverse=True  # Rare first in legend
    )

    handles = []
    labels_list = []

    for pattern in patterns_sorted:
        style = weights.get(pattern, {})
        color = PATTERN_COLORS.get(pattern, '#999999')
        count = len(df[df[pattern_col] == pattern])

        # Create line with appropriate width
        from matplotlib.lines import Line2D
        line = Line2D([0], [0], color=color,
                      linewidth=min(style.get('linewidth', 1) * 1.5, 4),
                      alpha=min(style.get('alpha', 0.5) * 1.5, 1.0))
        handles.append(line)

        # Label includes count
        zorder = style.get('zorder', 0)
        tier = {4: 'RARE', 3: 'uncommon', 2: 'common', 1: 'background'}.get(zorder, '')
        labels_list.append(f'{pattern} (n={count}) [{tier}]')

    ax.legend(handles, labels_list, loc='center left', fontsize=8,
              title='Pattern (line weight)', title_fontsize=9, framealpha=0.9)
    ax.axis('off')


def create_faceted_bump_weighted(df, y_type='nes', scope='focused', include_trajdev=True):
    """
    Create side-by-side weighted bump charts for G32A and R403C.

    Parameters
    ----------
    df : pd.DataFrame
        Pathway data
    y_type : str
        'nes' for NES values, 'rank' for rank positions
    scope : str
        'focused', 'significant', or 'all'
    include_trajdev : bool
        If True, 3-point trajectory; if False, 2-point (direct)

    Returns
    -------
    matplotlib.figure.Figure
    """
    # Filter data by scope
    df_scope = filter_by_scope(df, scope)

    if len(df_scope) == 0:
        print(f"  Warning: No data for scope '{scope}'")
        return None

    # Create figure with 3 columns: G32A, R403C, Legend
    fig, axes = plt.subplots(1, 3, figsize=(16, 8),
                              gridspec_kw={'width_ratios': [3.5, 3.5, 1.5]})

    # Plot G32A
    create_bump_chart_weighted(df_scope, 'G32A', axes[0], y_type=y_type,
                               include_trajdev=include_trajdev, scope=scope)
    axes[0].set_ylabel('NES' if y_type == 'nes' else 'Rank', fontsize=11)

    # Plot R403C
    create_bump_chart_weighted(df_scope, 'R403C', axes[1], y_type=y_type,
                               include_trajdev=include_trajdev, scope=scope)

    # Share y-axis limits for NES
    if y_type == 'nes':
        y_min = min(axes[0].get_ylim()[0], axes[1].get_ylim()[0])
        y_max = max(axes[0].get_ylim()[1], axes[1].get_ylim()[1])
        axes[0].set_ylim(y_min, y_max)
        axes[1].set_ylim(y_min, y_max)

    # Legend
    create_legend_weighted(axes[2], df_scope, 'G32A')

    # Main title
    scope_name = SCOPE_CONFIGS[scope]['name']
    y_label = 'NES Value' if y_type == 'nes' else 'Rank (by NES)'
    traj_type = '3-Point (Early→TrajDev→Late)' if include_trajdev else '2-Point Direct (Early→Late)'
    fig.suptitle(f'Weighted Bump Chart | {scope_name}\nY-axis: {y_label} | {traj_type}',
                 fontsize=14, fontweight='bold', y=0.98)

    plt.tight_layout(rect=[0, 0, 1, 0.94])

    return fig


# =============================================================================
# SUMMARY RIBBON VERSION (for dense plots)
# =============================================================================

def create_summary_bump_chart(df, mutation, ax, scope='all', include_trajdev=True):
    """
    Create summary bump chart with median trajectory and IQR ribbon per pattern.

    For dense plots (all pathways), individual lines become unreadable.
    This version shows:
    - Solid line: Median NES per pattern
    - Shaded ribbon: IQR (25th-75th percentile)

    Parameters
    ----------
    df : pd.DataFrame
        Pathway data
    mutation : str
        'G32A' or 'R403C'
    ax : matplotlib.axes.Axes
        Axes to plot on
    scope : str
        'focused', 'significant', or 'all'
    include_trajdev : bool
        If True, 3-point (Early→TrajDev→Late)
        If False, 2-point (Early→Late direct)
    """
    if include_trajdev:
        stages = ['Early', 'TrajDev', 'Late']
        x_positions = [0, 1, 2]
        x_labels = ['Early\n(D35)', 'TrajDev\n(Maturation)', 'Late\n(D65)']
    else:
        stages = ['Early', 'Late']
        x_positions = [0, 1]
        x_labels = ['Early\n(D35)', 'Late\n(D65)']

    nes_cols = [f'NES_{s}_{mutation}' for s in stages]
    pattern_col = f'Pattern_{mutation}'

    # Filter to rows with complete NES data
    df_plot = df[df[nes_cols].notna().all(axis=1)].copy()

    if len(df_plot) == 0:
        ax.text(0.5, 0.5, f'No complete data for {mutation}',
                ha='center', va='center', transform=ax.transAxes)
        return

    # Patterns to summarize (in order of importance)
    patterns_to_show = ['Compensation', 'Natural_improvement', 'Progressive',
                        'Natural_worsening', 'Late_onset', 'Transient']

    for pattern in patterns_to_show:
        df_pattern = df_plot[df_plot[pattern_col] == pattern]

        if len(df_pattern) < 3:  # Need at least 3 for meaningful summary
            continue

        color = PATTERN_COLORS.get(pattern, '#999999')

        # Compute statistics
        medians = [df_pattern[col].median() for col in nes_cols]
        q25 = [df_pattern[col].quantile(0.25) for col in nes_cols]
        q75 = [df_pattern[col].quantile(0.75) for col in nes_cols]

        # Plot ribbon (IQR)
        ax.fill_between(x_positions, q25, q75, color=color, alpha=0.2)

        # Plot median line
        ax.plot(x_positions, medians, color=color, linewidth=2.5,
                label=f'{pattern} (n={len(df_pattern)})', marker='o', markersize=5)

    # Add zero line
    ax.axhline(y=0, color='black', linewidth=0.5, linestyle='-', alpha=0.5)

    # Formatting
    ax.set_xticks(x_positions)
    ax.set_xticklabels(x_labels, fontsize=10)
    ax.set_xlim(-0.3, x_positions[-1] + 0.3)
    ax.set_ylim(-3, 3)

    ax.set_title(f'{mutation}\n(n={len(df_plot)} total)',
                 fontsize=12, fontweight='bold',
                 color=MUTATION_COLORS.get(mutation, 'black'))

    ax.legend(loc='upper right', fontsize=8, framealpha=0.9)
    ax.grid(axis='y', alpha=0.3, linestyle='-', linewidth=0.5)


def create_summary_faceted(df, scope='all', include_trajdev=True):
    """Create summary bump chart (median + IQR) for dense data."""
    df_scope = filter_by_scope(df, scope)

    if len(df_scope) == 0:
        return None

    fig, axes = plt.subplots(1, 2, figsize=(12, 6), sharey=True)

    create_summary_bump_chart(df_scope, 'G32A', axes[0], scope=scope,
                               include_trajdev=include_trajdev)
    axes[0].set_ylabel('NES', fontsize=11)

    create_summary_bump_chart(df_scope, 'R403C', axes[1], scope=scope,
                               include_trajdev=include_trajdev)

    scope_name = SCOPE_CONFIGS[scope]['name']
    traj_type = '3-Point (Early→TrajDev→Late)' if include_trajdev else '2-Point Direct (Early→Late)'
    fig.suptitle(f'Pathway Trajectory Summary\n{scope_name} | Median + IQR by Pattern\n{traj_type}',
                 fontsize=14, fontweight='bold', y=0.98)

    plt.tight_layout(rect=[0, 0, 1, 0.93])

    return fig


# =============================================================================
# PASS 5: REFINED WEIGHTING WITH INVERSE ALPHA
# =============================================================================

def compute_line_weights_refined(df, mutation):
    """
    Refined weighting with inverse alpha (thicker = more opaque).

    Changes from original compute_line_weights:
    - Reduced max linewidth: 2.5 → 2.0 (less dominating)
    - Increased min linewidth: 0.3 → 0.5 (more visible, especially in legend)
    - Inverse alpha: thicker = more opaque (60% → 90%)

    Parameters
    ----------
    df : pd.DataFrame
        Filtered pathway data
    mutation : str
        'G32A' or 'R403C'

    Returns
    -------
    dict
        {pattern: {'linewidth': float, 'alpha': float, 'zorder': int}}
    """
    pattern_col = f'Pattern_{mutation}'

    if pattern_col not in df.columns:
        return {}

    counts = df[pattern_col].value_counts()
    total = len(df)

    weights = {}
    for pattern, count in counts.items():
        frequency = count / total

        if frequency > 0.3:  # Dominant (Complex)
            weights[pattern] = {'linewidth': 0.5, 'alpha': 0.60, 'zorder': 1}
        elif frequency > 0.1:  # Common (Compensation, Natural_improvement)
            weights[pattern] = {'linewidth': 0.8, 'alpha': 0.70, 'zorder': 2}
        elif frequency > 0.01:  # Uncommon (Natural_worsening)
            weights[pattern] = {'linewidth': 1.2, 'alpha': 0.80, 'zorder': 3}
        else:  # Rare (Late_onset, Transient, Progressive)
            weights[pattern] = {'linewidth': 2.0, 'alpha': 0.90, 'zorder': 4}

    return weights


def compute_global_line_weights(df, mutations=['G32A', 'R403C']):
    """
    Compute line weights using COMBINED pattern frequencies across all mutations.

    This ensures consistent z-order/layering across both mutation panels,
    so the same pattern appears at the same layer in both panels.

    Parameters
    ----------
    df : pd.DataFrame
        Filtered pathway data
    mutations : list
        List of mutations to combine

    Returns
    -------
    dict
        {pattern: {'linewidth': float, 'alpha': float, 'zorder': int}}
    """
    # Collect pattern counts from all mutations
    pattern_counts = {}

    for mutation in mutations:
        pattern_col = f'Pattern_{mutation}'
        if pattern_col not in df.columns:
            continue

        # Get NES columns to filter to complete data
        nes_cols = [f'NES_Early_{mutation}', f'NES_Late_{mutation}']
        df_mut = df[df[nes_cols].notna().all(axis=1)]

        counts = df_mut[pattern_col].value_counts().to_dict()
        for pattern, count in counts.items():
            pattern_counts[pattern] = pattern_counts.get(pattern, 0) + count

    total = sum(pattern_counts.values())
    if total == 0:
        return {}

    weights = {}
    for pattern, count in pattern_counts.items():
        frequency = count / total

        if frequency > 0.3:  # Dominant (Complex)
            weights[pattern] = {'linewidth': 0.5, 'alpha': 0.60, 'zorder': 1}
        elif frequency > 0.1:  # Common (Compensation, Natural_improvement)
            weights[pattern] = {'linewidth': 0.8, 'alpha': 0.70, 'zorder': 2}
        elif frequency > 0.01:  # Uncommon
            weights[pattern] = {'linewidth': 1.2, 'alpha': 0.80, 'zorder': 3}
        else:  # Rare (Late_onset, Transient, Progressive)
            weights[pattern] = {'linewidth': 2.0, 'alpha': 0.90, 'zorder': 4}

    return weights


def create_bump_chart_refined(df, mutation, ax, y_type='nes', scope='focused'):
    """
    Create refined weighted bump chart (2-point only) with inverse alpha.

    Parameters
    ----------
    df : pd.DataFrame
        Filtered pathway data
    mutation : str
        'G32A' or 'R403C'
    ax : matplotlib.axes.Axes
        Axes to plot on
    y_type : str
        'nes' for NES values, 'rank' for rank positions
    scope : str
        'focused', 'significant', or 'all'
    """
    pattern_col = f'Pattern_{mutation}'

    # 2-point only (Early → Late)
    stages = ['Early', 'Late']
    x_positions = [0, 1]
    x_labels = ['Early\n(D35)', 'Late\n(D65)']

    nes_cols = [f'NES_{s}_{mutation}' for s in stages]

    # Check columns exist
    if pattern_col not in df.columns:
        print(f"  Warning: {pattern_col} not found")
        return

    # Filter to rows with complete NES data for this mutation
    df_plot = df[df[nes_cols].notna().all(axis=1)].copy()

    if len(df_plot) == 0:
        ax.text(0.5, 0.5, f'No complete data for {mutation}',
                ha='center', va='center', transform=ax.transAxes)
        return

    # Compute ranks if needed
    if y_type == 'rank':
        for col in nes_cols:
            rank_col = col.replace('NES_', 'Rank_')
            df_plot[rank_col] = df_plot[col].rank(ascending=True, method='first')

    # Compute REFINED line weights (inverse alpha)
    weights = compute_line_weights_refined(df_plot, mutation)

    # Get patterns sorted by z-order (background first, rare last)
    pattern_order = sorted(
        df_plot[pattern_col].unique(),
        key=lambda p: weights.get(p, {}).get('zorder', 0)
    )

    # Plot each pattern group in z-order
    for pattern in pattern_order:
        df_pattern = df_plot[df_plot[pattern_col] == pattern]
        style = weights.get(pattern, {'linewidth': 1, 'alpha': 0.5, 'zorder': 1})
        color = PATTERN_COLORS.get(pattern, '#999999')

        for _, row in df_pattern.iterrows():
            if y_type == 'nes':
                y_values = [row[col] for col in nes_cols]
            else:
                rank_cols = [col.replace('NES_', 'Rank_') for col in nes_cols]
                y_values = [row[col] for col in rank_cols]

            ax.plot(x_positions, y_values,
                    color=color,
                    linewidth=style['linewidth'],
                    alpha=style['alpha'],
                    zorder=style['zorder'],
                    solid_capstyle='round')

    # Add reference lines
    if y_type == 'nes':
        ax.axhline(y=0, color='black', linewidth=0.5, linestyle='-', alpha=0.3, zorder=0)
        ax.axhline(y=NES_EFFECT, color='gray', linewidth=0.5, linestyle='--', alpha=0.3, zorder=0)
        ax.axhline(y=-NES_EFFECT, color='gray', linewidth=0.5, linestyle='--', alpha=0.3, zorder=0)

        # Y-axis limits
        y_max = max(4, df_plot[nes_cols].max().max() * 1.1)
        y_min = min(-4, df_plot[nes_cols].min().min() * 1.1)
        ax.set_ylim(y_min, y_max)
    else:
        n_pathways = len(df_plot)
        ax.axhline(y=n_pathways/2, color='black', linewidth=0.5, linestyle='-', alpha=0.3, zorder=0)
        ax.set_ylim(0, n_pathways + 1)
        ax.invert_yaxis()

    # Formatting
    ax.set_xticks(x_positions)
    ax.set_xticklabels(x_labels, fontsize=10)
    ax.set_xlim(-0.3, 1.3)

    ax.set_title(f'{mutation}\n(n={len(df_plot)})',
                 fontsize=12, fontweight='bold',
                 color=MUTATION_COLORS.get(mutation, 'black'))

    ax.grid(axis='y', alpha=0.3, linestyle='-', linewidth=0.5)
    ax.set_axisbelow(True)


def create_legend_refined(ax, df, mutation):
    """Create legend showing pattern colors and refined line weights."""
    pattern_col = f'Pattern_{mutation}'

    if pattern_col not in df.columns:
        ax.axis('off')
        return

    weights = compute_line_weights_refined(df, mutation)

    # Sort patterns by zorder for legend ordering
    patterns_sorted = sorted(
        weights.keys(),
        key=lambda p: weights.get(p, {}).get('zorder', 0),
        reverse=True  # Rare first in legend
    )

    handles = []
    labels_list = []

    for pattern in patterns_sorted:
        style = weights.get(pattern, {})
        color = PATTERN_COLORS.get(pattern, '#999999')
        count = len(df[df[pattern_col] == pattern])

        # Create line with appropriate width (scaled up for visibility)
        from matplotlib.lines import Line2D
        line = Line2D([0], [0], color=color,
                      linewidth=min(style.get('linewidth', 1) * 2, 5),
                      alpha=min(style.get('alpha', 0.5), 1.0))
        handles.append(line)

        # Label includes count and tier
        zorder = style.get('zorder', 0)
        tier = {4: 'RARE', 3: 'uncommon', 2: 'common', 1: 'background'}.get(zorder, '')
        labels_list.append(f'{pattern} (n={count}) [{tier}]')

    ax.legend(handles, labels_list, loc='center left', fontsize=8,
              title='Pattern (refined weight)', title_fontsize=9, framealpha=0.9)
    ax.axis('off')


def create_legend_unified(ax, df, mutations=['G32A', 'R403C']):
    """
    Create unified legend showing pattern colors with counts from BOTH mutations.

    This ensures the legend reflects all patterns present in either panel,
    not just the left panel.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes for the legend
    df : pd.DataFrame
        Filtered pathway data
    mutations : list
        List of mutations to include (default: ['G32A', 'R403C'])
    """
    from matplotlib.lines import Line2D

    # Collect pattern counts from all mutations
    pattern_counts = {}  # {pattern: {mutation: count}}
    all_patterns = set()

    for mutation in mutations:
        pattern_col = f'Pattern_{mutation}'
        if pattern_col not in df.columns:
            continue

        # Get NES columns to filter to complete data
        nes_cols = [f'NES_Early_{mutation}', f'NES_Late_{mutation}']
        df_mut = df[df[nes_cols].notna().all(axis=1)]

        counts = df_mut[pattern_col].value_counts().to_dict()

        for pattern, count in counts.items():
            if pattern not in pattern_counts:
                pattern_counts[pattern] = {}
            pattern_counts[pattern][mutation] = count
            all_patterns.add(pattern)

    if not all_patterns:
        ax.axis('off')
        return

    # Compute weights using combined frequencies for z-order
    # Use total count across mutations for frequency calculation
    total_pathways = sum(
        sum(pattern_counts.get(p, {}).values())
        for p in all_patterns
    )

    pattern_weights = {}
    for pattern in all_patterns:
        total_count = sum(pattern_counts.get(pattern, {}).values())
        frequency = total_count / total_pathways if total_pathways > 0 else 0

        if frequency > 0.3:
            pattern_weights[pattern] = {'linewidth': 0.5, 'alpha': 0.60, 'zorder': 1, 'tier': 'background'}
        elif frequency > 0.1:
            pattern_weights[pattern] = {'linewidth': 0.8, 'alpha': 0.70, 'zorder': 2, 'tier': 'common'}
        elif frequency > 0.01:
            pattern_weights[pattern] = {'linewidth': 1.2, 'alpha': 0.80, 'zorder': 3, 'tier': 'uncommon'}
        else:
            pattern_weights[pattern] = {'linewidth': 2.0, 'alpha': 0.90, 'zorder': 4, 'tier': 'RARE'}

    # Sort patterns by zorder (rare first in legend)
    patterns_sorted = sorted(
        all_patterns,
        key=lambda p: pattern_weights.get(p, {}).get('zorder', 0),
        reverse=True
    )

    handles = []
    labels_list = []

    for pattern in patterns_sorted:
        style = pattern_weights.get(pattern, {})
        color = PATTERN_COLORS.get(pattern, '#999999')

        # Build count string: "G32A: X, R403C: Y"
        counts_str_parts = []
        for mutation in mutations:
            count = pattern_counts.get(pattern, {}).get(mutation, 0)
            counts_str_parts.append(f'{mutation}: {count}')
        counts_str = ', '.join(counts_str_parts)

        # Create line
        line = Line2D([0], [0], color=color,
                      linewidth=min(style.get('linewidth', 1) * 2, 5),
                      alpha=min(style.get('alpha', 0.5), 1.0))
        handles.append(line)

        # Label format: "Pattern (G32A: X, R403C: Y) [tier]"
        tier = style.get('tier', '')
        labels_list.append(f'{pattern}\n  ({counts_str}) [{tier}]')

    ax.legend(handles, labels_list, loc='center left', fontsize=7,
              title='Pattern (both mutations)', title_fontsize=9, framealpha=0.9,
              labelspacing=1.2)
    ax.axis('off')


def create_faceted_bump_refined(df, y_type='nes', scope='focused'):
    """
    Create side-by-side refined bump charts for G32A and R403C.
    Uses refined weighting (inverse alpha, adjusted linewidths).
    2-point only (Early → Late).
    """
    # Filter data by scope
    df_scope = filter_by_scope(df, scope)

    if len(df_scope) == 0:
        print(f"  Warning: No data for scope '{scope}'")
        return None

    # Create figure with 3 columns: G32A, R403C, Legend
    fig, axes = plt.subplots(1, 3, figsize=(14, 8),
                              gridspec_kw={'width_ratios': [3, 3, 1.5]})

    # Plot G32A
    create_bump_chart_refined(df_scope, 'G32A', axes[0], y_type=y_type, scope=scope)
    axes[0].set_ylabel('NES' if y_type == 'nes' else 'Rank', fontsize=11)

    # Plot R403C
    create_bump_chart_refined(df_scope, 'R403C', axes[1], y_type=y_type, scope=scope)

    # Share y-axis limits for NES
    if y_type == 'nes':
        y_min = min(axes[0].get_ylim()[0], axes[1].get_ylim()[0])
        y_max = max(axes[0].get_ylim()[1], axes[1].get_ylim()[1])
        axes[0].set_ylim(y_min, y_max)
        axes[1].set_ylim(y_min, y_max)

    # Legend
    create_legend_refined(axes[2], df_scope, 'G32A')

    # Main title
    scope_name = SCOPE_CONFIGS[scope]['name']
    y_label = 'NES Value' if y_type == 'nes' else 'Rank (by NES)'
    fig.suptitle(f'Refined Bump Chart (Inverse Alpha) | {scope_name}\nY-axis: {y_label} | 2-Point Direct (Early→Late)',
                 fontsize=14, fontweight='bold', y=0.98)

    plt.tight_layout(rect=[0, 0, 1, 0.94])

    return fig


# =============================================================================
# PASS 6: HIGHLIGHTED SIGNIFICANT PATHWAYS
# =============================================================================

def identify_highlight_pathways(df, mutation, max_per_pattern=5):
    """
    Identify pathways to highlight based on:
    1. Biological relevance (exclude irrelevant tissue/disease contexts)
    2. Significant in at least one contrast (Early, TrajDev, Late)
    3. Priority by semantic category (mito, synapse, translation first)
    4. Top N by |NES change| (Early→Late) within each Pattern

    Parameters
    ----------
    df : pd.DataFrame
        Filtered pathway data
    mutation : str
        'G32A' or 'R403C'
    max_per_pattern : int
        Maximum pathways to highlight per pattern

    Returns
    -------
    set
        Set of pathway_ids to highlight
    """
    pattern_col = f'Pattern_{mutation}'
    early_col = f'NES_Early_{mutation}'
    late_col = f'NES_Late_{mutation}'

    highlight_ids = set()

    # Get unique patterns (excluding Complex and Insufficient_data)
    meaningful_patterns = ['Compensation', 'Progressive', 'Natural_improvement',
                           'Natural_worsening', 'Late_onset', 'Transient']

    for pattern in meaningful_patterns:
        df_pattern = df[df[pattern_col] == pattern].copy()

        if len(df_pattern) == 0:
            continue

        # STEP 1: Filter out biologically irrelevant pathways
        df_pattern = df_pattern[
            df_pattern['Description'].apply(is_relevant_for_highlight)
        ]

        if len(df_pattern) == 0:
            continue

        # STEP 2: Filter to significant pathways
        # Use ever_significant_trajectory column if available (pre-computed)
        # Otherwise, meaningful patterns already require significance
        if 'ever_significant_trajectory' in df_pattern.columns:
            df_pattern = df_pattern[df_pattern['ever_significant_trajectory'] == True]
        elif 'ever_significant' in df_pattern.columns:
            df_pattern = df_pattern[df_pattern['ever_significant'] == True]
        # Note: If neither column exists, meaningful patterns already imply significance

        if len(df_pattern) == 0:
            continue

        # STEP 3: Assign semantic category and priority
        df_pattern['_semantic_cat'] = df_pattern.apply(
            lambda row: assign_semantic_category(row), axis=1
        )
        df_pattern['_priority'] = df_pattern['_semantic_cat'].apply(get_highlight_priority)

        # Calculate |NES change| from Early to Late
        if early_col in df_pattern.columns and late_col in df_pattern.columns:
            df_pattern['abs_nes_change'] = (df_pattern[late_col] - df_pattern[early_col]).abs()

            # STEP 4: Select top pathways with priority weighting
            # First, get priority pathways (semantic category in priority list)
            priority_mask = df_pattern['_priority'] < 99
            df_priority = df_pattern[priority_mask].copy()
            df_other = df_pattern[~priority_mask].copy()

            selected_ids = []

            # From priority categories: sort by priority (lower=better), then by |NES change|
            if len(df_priority) > 0:
                df_priority = df_priority.sort_values(
                    ['_priority', 'abs_nes_change'],
                    ascending=[True, False]
                )
                n_priority = min(max_per_pattern, len(df_priority))
                selected_ids.extend(df_priority.head(n_priority)['pathway_id'].tolist())

            # From 'Other' category: only if we haven't filled quota
            remaining_slots = max_per_pattern - len(selected_ids)
            if remaining_slots > 0 and len(df_other) > 0:
                # Limit 'Other' category to MAX_OTHER_CATEGORY_HIGHLIGHTS
                n_other = min(remaining_slots, MAX_OTHER_CATEGORY_HIGHLIGHTS, len(df_other))
                df_other = df_other.nlargest(n_other, 'abs_nes_change')
                selected_ids.extend(df_other['pathway_id'].tolist())

            highlight_ids.update(selected_ids)

    return highlight_ids


def create_bump_chart_highlight(df, mutation, ax, y_type='nes', scope='focused',
                                 global_weights=None):
    """
    Create highlighted bump chart with significant pathways at full opacity and labeled.

    Builds on refined weighting but:
    - Highlighted pathways: alpha=1.0, labeled
    - Non-highlighted: use refined alpha from compute_line_weights_refined

    Parameters
    ----------
    df : pd.DataFrame
        Filtered pathway data
    mutation : str
        'G32A' or 'R403C'
    ax : matplotlib.axes.Axes
        Axes to plot on
    y_type : str
        'nes' for NES values, 'rank' for rank positions
    scope : str
        'focused', 'significant', or 'all'
    global_weights : dict, optional
        Pre-computed global weights for consistent z-order across panels.
        If None, computes per-panel weights.
    """
    pattern_col = f'Pattern_{mutation}'

    # 2-point only (Early → Late)
    stages = ['Early', 'Late']
    x_positions = [0, 1]
    x_labels = ['Early\n(D35)', 'Late\n(D65)']

    nes_cols = [f'NES_{s}_{mutation}' for s in stages]

    # Check columns exist
    if pattern_col not in df.columns:
        print(f"  Warning: {pattern_col} not found")
        return

    # Filter to rows with complete NES data for this mutation
    df_plot = df[df[nes_cols].notna().all(axis=1)].copy()

    if len(df_plot) == 0:
        ax.text(0.5, 0.5, f'No complete data for {mutation}',
                ha='center', va='center', transform=ax.transAxes)
        return

    # Compute ranks if needed
    if y_type == 'rank':
        for col in nes_cols:
            rank_col = col.replace('NES_', 'Rank_')
            df_plot[rank_col] = df_plot[col].rank(ascending=True, method='first')

    # Identify pathways to highlight
    highlight_ids = identify_highlight_pathways(df_plot, mutation, max_per_pattern=5)

    # Use global weights if provided (for consistent z-order across panels)
    # Otherwise compute per-panel weights
    if global_weights is not None:
        weights = global_weights
    else:
        weights = compute_line_weights_refined(df_plot, mutation)

    # Get patterns sorted by z-order (background first, rare last)
    pattern_order = sorted(
        df_plot[pattern_col].unique(),
        key=lambda p: weights.get(p, {}).get('zorder', 0)
    )

    # First pass: Plot all non-highlighted pathways
    for pattern in pattern_order:
        df_pattern = df_plot[df_plot[pattern_col] == pattern]
        df_non_highlight = df_pattern[~df_pattern['pathway_id'].isin(highlight_ids)]

        style = weights.get(pattern, {'linewidth': 1, 'alpha': 0.5, 'zorder': 1})
        color = PATTERN_COLORS.get(pattern, '#999999')

        for _, row in df_non_highlight.iterrows():
            if y_type == 'nes':
                y_values = [row[col] for col in nes_cols]
            else:
                rank_cols = [col.replace('NES_', 'Rank_') for col in nes_cols]
                y_values = [row[col] for col in rank_cols]

            ax.plot(x_positions, y_values,
                    color=color,
                    linewidth=style['linewidth'],
                    alpha=style['alpha'],
                    zorder=style['zorder'],
                    solid_capstyle='round')

    # Second pass: Plot highlighted pathways on top with full opacity + collect labels
    # Use adjustText for collision avoidance
    try:
        from adjustText import adjust_text
        use_adjust_text = True
    except ImportError:
        use_adjust_text = False

    label_texts = []  # Collect text objects for adjustText
    label_points_x = []  # X coordinates of line endpoints
    label_points_y = []  # Y coordinates of line endpoints

    for pattern in pattern_order:
        df_pattern = df_plot[df_plot[pattern_col] == pattern]
        df_highlight = df_pattern[df_pattern['pathway_id'].isin(highlight_ids)]

        if len(df_highlight) == 0:
            continue

        style = weights.get(pattern, {'linewidth': 1, 'alpha': 0.5, 'zorder': 1})
        color = PATTERN_COLORS.get(pattern, '#999999')

        for _, row in df_highlight.iterrows():
            if y_type == 'nes':
                y_values = [row[col] for col in nes_cols]
            else:
                rank_cols = [col.replace('NES_', 'Rank_') for col in nes_cols]
                y_values = [row[col] for col in rank_cols]

            # Full opacity for highlighted, keep same linewidth
            ax.plot(x_positions, y_values,
                    color=color,
                    linewidth=style['linewidth'] * 1.2,  # Slightly thicker
                    alpha=1.0,  # Full opacity
                    zorder=10,  # On top
                    solid_capstyle='round')

            # Collect label info for later adjustment
            y_end = y_values[-1]
            desc = str(row.get('Description', row.get('pathway_id', '')))
            label = desc[:30] + '...' if len(desc) > 30 else desc

            # Create text object at initial position
            txt = ax.text(1.08, y_end, label,
                         fontsize=6, color=color, ha='left', va='center',
                         fontweight='bold', zorder=11)
            label_texts.append(txt)
            label_points_x.append(1)  # Line endpoint x
            label_points_y.append(y_end)  # Line endpoint y

    # Apply adjustText for label collision avoidance
    if use_adjust_text and label_texts:
        adjust_text(
            label_texts,
            x=label_points_x,
            y=label_points_y,
            ax=ax,
            arrowprops=dict(arrowstyle='-', color='gray', lw=0.5, alpha=0.5),
            force_text=(0.5, 1.0),  # Force to push text away
            force_points=(0.2, 0.5),  # Force from line endpoints
            expand_text=(1.2, 1.5),  # Expand bounding boxes
            only_move={'points': 'y', 'text': 'xy'},  # Allow vertical movement
            lim=500,  # Max iterations
        )

    # Add reference lines
    if y_type == 'nes':
        ax.axhline(y=0, color='black', linewidth=0.5, linestyle='-', alpha=0.3, zorder=0)
        ax.axhline(y=NES_EFFECT, color='gray', linewidth=0.5, linestyle='--', alpha=0.3, zorder=0)
        ax.axhline(y=-NES_EFFECT, color='gray', linewidth=0.5, linestyle='--', alpha=0.3, zorder=0)

        # Y-axis limits
        y_max = max(4, df_plot[nes_cols].max().max() * 1.1)
        y_min = min(-4, df_plot[nes_cols].min().min() * 1.1)
        ax.set_ylim(y_min, y_max)
    else:
        n_pathways = len(df_plot)
        ax.axhline(y=n_pathways/2, color='black', linewidth=0.5, linestyle='-', alpha=0.3, zorder=0)
        ax.set_ylim(0, n_pathways + 1)
        ax.invert_yaxis()

    # Formatting
    ax.set_xticks(x_positions)
    ax.set_xticklabels(x_labels, fontsize=10)
    ax.set_xlim(-0.3, 1.6)  # Extra space for labels

    n_highlighted = len(df_plot[df_plot['pathway_id'].isin(highlight_ids)])
    ax.set_title(f'{mutation}\n(n={len(df_plot)}, {n_highlighted} highlighted)',
                 fontsize=12, fontweight='bold',
                 color=MUTATION_COLORS.get(mutation, 'black'))

    ax.grid(axis='y', alpha=0.3, linestyle='-', linewidth=0.5)
    ax.set_axisbelow(True)


def create_bump_chart_highlight_v8(df, mutation, ax, y_type='nes', scope='focused',
                                    global_weights=None):
    """
    Create highlighted bump chart with IMPROVED LAYOUT (Pass 8).

    Improvements over v7:
    - Lines start closer to y-axis (clear NES correspondence)
    - Labels positioned further right (no overlap with lines)
    - Stronger adjustText parameters for better separation

    Parameters
    ----------
    df : pd.DataFrame
        Filtered pathway data
    mutation : str
        'G32A' or 'R403C'
    ax : matplotlib.axes.Axes
        Axes to plot on
    y_type : str
        'nes' for NES values, 'rank' for rank positions
    scope : str
        'focused', 'significant', or 'all'
    global_weights : dict, optional
        Pre-computed global weights for consistent z-order across panels.
    """
    pattern_col = f'Pattern_{mutation}'

    # 2-point only (Early → Late)
    stages = ['Early', 'Late']
    x_positions = [0, 1]
    x_labels = ['Early\n(D35)', 'Late\n(D65)']

    nes_cols = [f'NES_{s}_{mutation}' for s in stages]

    # Check columns exist
    if pattern_col not in df.columns:
        print(f"  Warning: {pattern_col} not found")
        return

    # Filter to rows with complete NES data for this mutation
    df_plot = df[df[nes_cols].notna().all(axis=1)].copy()

    if len(df_plot) == 0:
        ax.text(0.5, 0.5, f'No complete data for {mutation}',
                ha='center', va='center', transform=ax.transAxes)
        return

    # Compute ranks if needed
    if y_type == 'rank':
        for col in nes_cols:
            rank_col = col.replace('NES_', 'Rank_')
            df_plot[rank_col] = df_plot[col].rank(ascending=True, method='first')

    # Identify pathways to highlight
    highlight_ids = identify_highlight_pathways(df_plot, mutation, max_per_pattern=5)

    # Use global weights if provided
    if global_weights is not None:
        weights = global_weights
    else:
        weights = compute_line_weights_refined(df_plot, mutation)

    # Get patterns sorted by z-order (background first, rare last)
    pattern_order = sorted(
        df_plot[pattern_col].unique(),
        key=lambda p: weights.get(p, {}).get('zorder', 0)
    )

    # First pass: Plot all non-highlighted pathways
    for pattern in pattern_order:
        df_pattern = df_plot[df_plot[pattern_col] == pattern]
        df_non_highlight = df_pattern[~df_pattern['pathway_id'].isin(highlight_ids)]

        style = weights.get(pattern, {'linewidth': 1, 'alpha': 0.5, 'zorder': 1})
        color = PATTERN_COLORS.get(pattern, '#999999')

        for _, row in df_non_highlight.iterrows():
            if y_type == 'nes':
                y_values = [row[col] for col in nes_cols]
            else:
                rank_cols = [col.replace('NES_', 'Rank_') for col in nes_cols]
                y_values = [row[col] for col in rank_cols]

            ax.plot(x_positions, y_values,
                    color=color,
                    linewidth=style['linewidth'],
                    alpha=style['alpha'],
                    zorder=style['zorder'],
                    solid_capstyle='round')

    # Second pass: Plot highlighted pathways + collect labels
    try:
        from adjustText import adjust_text
        use_adjust_text = True
    except ImportError:
        use_adjust_text = False

    label_texts = []
    label_points_x = []
    label_points_y = []

    for pattern in pattern_order:
        df_pattern = df_plot[df_plot[pattern_col] == pattern]
        df_highlight = df_pattern[df_pattern['pathway_id'].isin(highlight_ids)]

        if len(df_highlight) == 0:
            continue

        style = weights.get(pattern, {'linewidth': 1, 'alpha': 0.5, 'zorder': 1})
        color = PATTERN_COLORS.get(pattern, '#999999')

        for _, row in df_highlight.iterrows():
            if y_type == 'nes':
                y_values = [row[col] for col in nes_cols]
            else:
                rank_cols = [col.replace('NES_', 'Rank_') for col in nes_cols]
                y_values = [row[col] for col in rank_cols]

            # Plot highlighted line
            ax.plot(x_positions, y_values,
                    color=color,
                    linewidth=style['linewidth'] * 1.2,
                    alpha=1.0,
                    zorder=10,
                    solid_capstyle='round')

            # Collect label info - MOVED FURTHER RIGHT (1.15 instead of 1.08)
            y_end = y_values[-1]
            desc = str(row.get('Description', row.get('pathway_id', '')))
            label = desc[:30] + '...' if len(desc) > 30 else desc

            txt = ax.text(1.15, y_end, label,
                         fontsize=6, color=color, ha='left', va='center',
                         fontweight='bold', zorder=11)
            label_texts.append(txt)
            label_points_x.append(1)
            label_points_y.append(y_end)

    # Apply adjustText with STRONGER parameters
    if use_adjust_text and label_texts:
        adjust_text(
            label_texts,
            x=label_points_x,
            y=label_points_y,
            ax=ax,
            arrowprops=dict(arrowstyle='-', color='gray', lw=0.5, alpha=0.5),
            force_text=(0.8, 1.2),      # Stronger push
            force_points=(0.3, 0.8),    # Stronger repulsion
            expand_text=(1.3, 1.6),     # More expansion
            only_move={'points': 'y', 'text': 'xy'},
            lim=500,
        )

    # Add reference lines
    if y_type == 'nes':
        ax.axhline(y=0, color='black', linewidth=0.5, linestyle='-', alpha=0.3, zorder=0)
        ax.axhline(y=NES_EFFECT, color='gray', linewidth=0.5, linestyle='--', alpha=0.3, zorder=0)
        ax.axhline(y=-NES_EFFECT, color='gray', linewidth=0.5, linestyle='--', alpha=0.3, zorder=0)

        y_max = max(4, df_plot[nes_cols].max().max() * 1.1)
        y_min = min(-4, df_plot[nes_cols].min().min() * 1.1)
        ax.set_ylim(y_min, y_max)
    else:
        n_pathways = len(df_plot)
        ax.axhline(y=n_pathways/2, color='black', linewidth=0.5, linestyle='-', alpha=0.3, zorder=0)
        ax.set_ylim(0, n_pathways + 1)
        ax.invert_yaxis()

    # Formatting - TIGHTER LEFT MARGIN, MORE SPACE RIGHT
    ax.set_xticks(x_positions)
    ax.set_xticklabels(x_labels, fontsize=10)
    ax.set_xlim(-0.05, 1.8)  # Lines start close to y-axis, more label space

    n_highlighted = len(df_plot[df_plot['pathway_id'].isin(highlight_ids)])
    ax.set_title(f'{mutation}\n(n={len(df_plot)}, {n_highlighted} highlighted)',
                 fontsize=12, fontweight='bold',
                 color=MUTATION_COLORS.get(mutation, 'black'))

    ax.grid(axis='y', alpha=0.3, linestyle='-', linewidth=0.5)
    ax.set_axisbelow(True)


def create_bump_chart_highlight_v9(df, mutation, ax, y_type='nes', scope='focused',
                                    global_weights=None):
    """
    Create highlighted bump chart with FIXED CONNECTORS (Pass 9).

    Improvements over v8:
    - Labels moved further right (x=1.25) to avoid overlapping lines
    - Manual horizontal connectors (fixed position, don't follow label movement)
    - adjustText without arrow connectors

    Parameters
    ----------
    df : pd.DataFrame
        Filtered pathway data
    mutation : str
        'G32A' or 'R403C'
    ax : matplotlib.axes.Axes
        Axes to plot on
    y_type : str
        'nes' for NES values, 'rank' for rank positions
    scope : str
        'focused', 'significant', or 'all'
    global_weights : dict, optional
        Pre-computed global weights for consistent z-order across panels.
    """
    pattern_col = f'Pattern_{mutation}'

    # 2-point only (Early → Late)
    stages = ['Early', 'Late']
    x_positions = [0, 1]
    x_labels = ['Early\n(D35)', 'Late\n(D65)']

    nes_cols = [f'NES_{s}_{mutation}' for s in stages]

    # Check columns exist
    if pattern_col not in df.columns:
        print(f"  Warning: {pattern_col} not found")
        return

    # Filter to rows with complete NES data for this mutation
    df_plot = df[df[nes_cols].notna().all(axis=1)].copy()

    if len(df_plot) == 0:
        ax.text(0.5, 0.5, f'No complete data for {mutation}',
                ha='center', va='center', transform=ax.transAxes)
        return

    # Compute ranks if needed
    if y_type == 'rank':
        for col in nes_cols:
            rank_col = col.replace('NES_', 'Rank_')
            df_plot[rank_col] = df_plot[col].rank(ascending=True, method='first')

    # Identify pathways to highlight
    highlight_ids = identify_highlight_pathways(df_plot, mutation, max_per_pattern=5)

    # Use global weights if provided
    if global_weights is not None:
        weights = global_weights
    else:
        weights = compute_line_weights_refined(df_plot, mutation)

    # Get patterns sorted by z-order (background first, rare last)
    pattern_order = sorted(
        df_plot[pattern_col].unique(),
        key=lambda p: weights.get(p, {}).get('zorder', 0)
    )

    # First pass: Plot all non-highlighted pathways
    for pattern in pattern_order:
        df_pattern = df_plot[df_plot[pattern_col] == pattern]
        df_non_highlight = df_pattern[~df_pattern['pathway_id'].isin(highlight_ids)]

        style = weights.get(pattern, {'linewidth': 1, 'alpha': 0.5, 'zorder': 1})
        color = PATTERN_COLORS.get(pattern, '#999999')

        for _, row in df_non_highlight.iterrows():
            if y_type == 'nes':
                y_values = [row[col] for col in nes_cols]
            else:
                rank_cols = [col.replace('NES_', 'Rank_') for col in nes_cols]
                y_values = [row[col] for col in rank_cols]

            ax.plot(x_positions, y_values,
                    color=color,
                    linewidth=style['linewidth'],
                    alpha=style['alpha'],
                    zorder=style['zorder'],
                    solid_capstyle='round')

    # Second pass: Plot highlighted pathways + connectors + labels
    try:
        from adjustText import adjust_text
        use_adjust_text = True
    except ImportError:
        use_adjust_text = False

    label_texts = []
    label_points_x = []
    label_points_y = []

    # Connector parameters
    connector_x_start = 1.02   # Just after line endpoint
    connector_x_end = 1.23     # Just before label
    label_x = 1.25             # Label start position

    for pattern in pattern_order:
        df_pattern = df_plot[df_plot[pattern_col] == pattern]
        df_highlight = df_pattern[df_pattern['pathway_id'].isin(highlight_ids)]

        if len(df_highlight) == 0:
            continue

        style = weights.get(pattern, {'linewidth': 1, 'alpha': 0.5, 'zorder': 1})
        color = PATTERN_COLORS.get(pattern, '#999999')

        for _, row in df_highlight.iterrows():
            if y_type == 'nes':
                y_values = [row[col] for col in nes_cols]
            else:
                rank_cols = [col.replace('NES_', 'Rank_') for col in nes_cols]
                y_values = [row[col] for col in rank_cols]

            # Plot highlighted line
            ax.plot(x_positions, y_values,
                    color=color,
                    linewidth=style['linewidth'] * 1.2,
                    alpha=1.0,
                    zorder=10,
                    solid_capstyle='round')

            # Draw MANUAL horizontal connector (fixed position)
            y_end = y_values[-1]
            ax.plot([connector_x_start, connector_x_end], [y_end, y_end],
                    color='gray', lw=0.5, alpha=0.5, zorder=9)

            # Create label FURTHER RIGHT
            desc = str(row.get('Description', row.get('pathway_id', '')))
            label = desc[:30] + '...' if len(desc) > 30 else desc

            txt = ax.text(label_x, y_end, label,
                         fontsize=6, color=color, ha='left', va='center',
                         fontweight='bold', zorder=11)
            label_texts.append(txt)
            label_points_x.append(connector_x_end)  # Connector end, not line endpoint
            label_points_y.append(y_end)

    # Apply adjustText WITHOUT arrow connectors (we drew our own)
    if use_adjust_text and label_texts:
        adjust_text(
            label_texts,
            x=label_points_x,
            y=label_points_y,
            ax=ax,
            arrowprops=None,  # No automatic connectors - we drew manual ones
            force_text=(0.5, 0.8),
            force_points=(0.2, 0.5),
            expand_text=(1.2, 1.4),
            only_move={'points': 'y', 'text': 'xy'},
            lim=500,
        )

    # Add reference lines
    if y_type == 'nes':
        ax.axhline(y=0, color='black', linewidth=0.5, linestyle='-', alpha=0.3, zorder=0)
        ax.axhline(y=NES_EFFECT, color='gray', linewidth=0.5, linestyle='--', alpha=0.3, zorder=0)
        ax.axhline(y=-NES_EFFECT, color='gray', linewidth=0.5, linestyle='--', alpha=0.3, zorder=0)

        y_max = max(4, df_plot[nes_cols].max().max() * 1.1)
        y_min = min(-4, df_plot[nes_cols].min().min() * 1.1)
        ax.set_ylim(y_min, y_max)
    else:
        n_pathways = len(df_plot)
        ax.axhline(y=n_pathways/2, color='black', linewidth=0.5, linestyle='-', alpha=0.3, zorder=0)
        ax.set_ylim(0, n_pathways + 1)
        ax.invert_yaxis()

    # Formatting - EVEN MORE SPACE for labels
    ax.set_xticks(x_positions)
    ax.set_xticklabels(x_labels, fontsize=10)
    ax.set_xlim(-0.05, 2.0)  # More label space

    n_highlighted = len(df_plot[df_plot['pathway_id'].isin(highlight_ids)])
    ax.set_title(f'{mutation}\n(n={len(df_plot)}, {n_highlighted} highlighted)',
                 fontsize=12, fontweight='bold',
                 color=MUTATION_COLORS.get(mutation, 'black'))

    ax.grid(axis='y', alpha=0.3, linestyle='-', linewidth=0.5)
    ax.set_axisbelow(True)


def create_bump_chart_highlight_v10(df, mutation, ax, y_type='nes', scope='focused',
                                    global_weights=None):
    """
    Create highlighted bump chart with DIAGONAL CONNECTORS (Pass 10).

    Improvements over v9:
    - Labels moved far right (x=2.0) to clearly separate from lines
    - Diagonal connectors via adjustText (follow label movement)
    - No manual horizontal connectors
    - Left-aligned labels (no ladder pattern)

    Parameters
    ----------
    df : pd.DataFrame
        Filtered pathway data
    mutation : str
        'G32A' or 'R403C'
    ax : matplotlib.axes.Axes
        Axes to plot on
    y_type : str
        'nes' for NES values, 'rank' for rank positions
    scope : str
        'focused', 'significant', or 'all'
    global_weights : dict, optional
        Pre-computed global weights for consistent z-order across panels.
    """
    pattern_col = f'Pattern_{mutation}'

    # 2-point only (Early → Late)
    stages = ['Early', 'Late']
    x_positions = [0, 1]
    x_labels = ['Early\n(D35)', 'Late\n(D65)']

    nes_cols = [f'NES_{s}_{mutation}' for s in stages]

    # Check columns exist
    if pattern_col not in df.columns:
        print(f"  Warning: {pattern_col} not found")
        return

    # Filter to rows with complete NES data for this mutation
    df_plot = df[df[nes_cols].notna().all(axis=1)].copy()

    if len(df_plot) == 0:
        ax.text(0.5, 0.5, f'No complete data for {mutation}',
                ha='center', va='center', transform=ax.transAxes)
        return

    # Compute ranks if needed
    if y_type == 'rank':
        for col in nes_cols:
            rank_col = col.replace('NES_', 'Rank_')
            df_plot[rank_col] = df_plot[col].rank(ascending=True, method='first')

    # Identify pathways to highlight
    highlight_ids = identify_highlight_pathways(df_plot, mutation, max_per_pattern=5)

    # Use global weights if provided
    if global_weights is not None:
        weights = global_weights
    else:
        weights = compute_line_weights_refined(df_plot, mutation)

    # Get patterns sorted by z-order (background first, rare last)
    pattern_order = sorted(
        df_plot[pattern_col].unique(),
        key=lambda p: weights.get(p, {}).get('zorder', 0)
    )

    # First pass: Plot all non-highlighted pathways
    for pattern in pattern_order:
        df_pattern = df_plot[df_plot[pattern_col] == pattern]
        df_non_highlight = df_pattern[~df_pattern['pathway_id'].isin(highlight_ids)]

        style = weights.get(pattern, {'linewidth': 1, 'alpha': 0.5, 'zorder': 1})
        color = PATTERN_COLORS.get(pattern, '#999999')

        for _, row in df_non_highlight.iterrows():
            if y_type == 'nes':
                y_values = [row[col] for col in nes_cols]
            else:
                rank_cols = [col.replace('NES_', 'Rank_') for col in nes_cols]
                y_values = [row[col] for col in rank_cols]

            ax.plot(x_positions, y_values,
                    color=color,
                    linewidth=style['linewidth'],
                    alpha=style['alpha'],
                    zorder=style['zorder'],
                    solid_capstyle='round')

    # Second pass: Plot highlighted pathways + collect labels
    try:
        from adjustText import adjust_text
        use_adjust_text = True
    except ImportError:
        use_adjust_text = False

    label_texts = []
    label_points_x = []
    label_points_y = []

    # Label position - FAR RIGHT
    label_x = 2.0

    for pattern in pattern_order:
        df_pattern = df_plot[df_plot[pattern_col] == pattern]
        df_highlight = df_pattern[df_pattern['pathway_id'].isin(highlight_ids)]

        if len(df_highlight) == 0:
            continue

        style = weights.get(pattern, {'linewidth': 1, 'alpha': 0.5, 'zorder': 1})
        color = PATTERN_COLORS.get(pattern, '#999999')

        for _, row in df_highlight.iterrows():
            if y_type == 'nes':
                y_values = [row[col] for col in nes_cols]
            else:
                rank_cols = [col.replace('NES_', 'Rank_') for col in nes_cols]
                y_values = [row[col] for col in rank_cols]

            # Plot highlighted line
            ax.plot(x_positions, y_values,
                    color=color,
                    linewidth=style['linewidth'] * 1.2,
                    alpha=1.0,
                    zorder=10,
                    solid_capstyle='round')

            # Create label FAR RIGHT (no manual connectors - adjustText will draw them)
            y_end = y_values[-1]
            desc = str(row.get('Description', row.get('pathway_id', '')))
            label = desc[:30] + '...' if len(desc) > 30 else desc

            txt = ax.text(label_x, y_end, label,
                         fontsize=6, color=color, ha='left', va='center',
                         fontweight='bold', zorder=11)
            label_texts.append(txt)
            label_points_x.append(1)  # Line endpoint x position (anchor for connector)
            label_points_y.append(y_end)

    # Apply adjustText WITH diagonal connectors
    if use_adjust_text and label_texts:
        adjust_text(
            label_texts,
            x=label_points_x,
            y=label_points_y,
            ax=ax,
            arrowprops=dict(arrowstyle='-', color='gray', lw=0.5, alpha=0.5),
            force_text=(0.5, 1.0),
            force_points=(0.3, 0.8),
            expand_text=(1.3, 1.5),
            only_move={'points': 'y', 'text': 'xy'},
            lim=500,
        )

    # Add reference lines
    if y_type == 'nes':
        ax.axhline(y=0, color='black', linewidth=0.5, linestyle='-', alpha=0.3, zorder=0)
        ax.axhline(y=NES_EFFECT, color='gray', linewidth=0.5, linestyle='--', alpha=0.3, zorder=0)
        ax.axhline(y=-NES_EFFECT, color='gray', linewidth=0.5, linestyle='--', alpha=0.3, zorder=0)

        y_max = max(4, df_plot[nes_cols].max().max() * 1.1)
        y_min = min(-4, df_plot[nes_cols].min().min() * 1.1)
        ax.set_ylim(y_min, y_max)
    else:
        n_pathways = len(df_plot)
        ax.axhline(y=n_pathways/2, color='black', linewidth=0.5, linestyle='-', alpha=0.3, zorder=0)
        ax.set_ylim(0, n_pathways + 1)
        ax.invert_yaxis()

    # Formatting - MUCH MORE SPACE for labels at x=2
    ax.set_xticks(x_positions)
    ax.set_xticklabels(x_labels, fontsize=10)
    ax.set_xlim(-0.05, 3.5)  # Wide for labels at x=2

    n_highlighted = len(df_plot[df_plot['pathway_id'].isin(highlight_ids)])
    ax.set_title(f'{mutation}\n(n={len(df_plot)}, {n_highlighted} highlighted)',
                 fontsize=12, fontweight='bold',
                 color=MUTATION_COLORS.get(mutation, 'black'))

    ax.grid(axis='y', alpha=0.3, linestyle='-', linewidth=0.5)
    ax.set_axisbelow(True)


def create_bump_chart_highlight_v11(df, mutation, ax, y_type='nes', scope='focused',
                                    global_weights=None):
    """
    Create highlighted bump chart with MANUAL CONNECTORS after adjustText (Pass 11).

    Improvements over v10:
    - Labels at fixed x=2.0 (no horizontal movement)
    - adjustText for vertical-only collision avoidance
    - Manual connectors drawn AFTER adjustText, from line endpoints to labels
    - Diagonal connectors that follow label vertical shifts

    Parameters
    ----------
    df : pd.DataFrame
        Filtered pathway data
    mutation : str
        'G32A' or 'R403C'
    ax : matplotlib.axes.Axes
        Axes to plot on
    y_type : str
        'nes' for NES values, 'rank' for rank positions
    scope : str
        'focused', 'significant', or 'all'
    global_weights : dict, optional
        Pre-computed global weights for consistent z-order across panels.
    """
    pattern_col = f'Pattern_{mutation}'

    # 2-point only (Early → Late)
    stages = ['Early', 'Late']
    x_positions = [0, 1]
    x_labels = ['Early\n(D35)', 'Late\n(D65)']

    nes_cols = [f'NES_{s}_{mutation}' for s in stages]

    # Check columns exist
    if pattern_col not in df.columns:
        print(f"  Warning: {pattern_col} not found")
        return

    # Filter to rows with complete NES data for this mutation
    df_plot = df[df[nes_cols].notna().all(axis=1)].copy()

    if len(df_plot) == 0:
        ax.text(0.5, 0.5, f'No complete data for {mutation}',
                ha='center', va='center', transform=ax.transAxes)
        return

    # Compute ranks if needed
    if y_type == 'rank':
        for col in nes_cols:
            rank_col = col.replace('NES_', 'Rank_')
            df_plot[rank_col] = df_plot[col].rank(ascending=True, method='first')

    # Identify pathways to highlight
    highlight_ids = identify_highlight_pathways(df_plot, mutation, max_per_pattern=5)

    # Use global weights if provided
    if global_weights is not None:
        weights = global_weights
    else:
        weights = compute_line_weights_refined(df_plot, mutation)

    # Get patterns sorted by z-order (background first, rare last)
    pattern_order = sorted(
        df_plot[pattern_col].unique(),
        key=lambda p: weights.get(p, {}).get('zorder', 0)
    )

    # First pass: Plot all non-highlighted pathways
    for pattern in pattern_order:
        df_pattern = df_plot[df_plot[pattern_col] == pattern]
        df_non_highlight = df_pattern[~df_pattern['pathway_id'].isin(highlight_ids)]

        style = weights.get(pattern, {'linewidth': 1, 'alpha': 0.5, 'zorder': 1})
        color = PATTERN_COLORS.get(pattern, '#999999')

        for _, row in df_non_highlight.iterrows():
            if y_type == 'nes':
                y_values = [row[col] for col in nes_cols]
            else:
                rank_cols = [col.replace('NES_', 'Rank_') for col in nes_cols]
                y_values = [row[col] for col in rank_cols]

            ax.plot(x_positions, y_values,
                    color=color,
                    linewidth=style['linewidth'],
                    alpha=style['alpha'],
                    zorder=style['zorder'],
                    solid_capstyle='round')

    # Second pass: Plot highlighted pathways + collect labels
    try:
        from adjustText import adjust_text
        use_adjust_text = True
    except ImportError:
        use_adjust_text = False

    label_texts = []
    label_points_y = []  # Original y positions (line endpoints)

    # Label position - FIXED at x=2.0
    label_x = 2.0

    for pattern in pattern_order:
        df_pattern = df_plot[df_plot[pattern_col] == pattern]
        df_highlight = df_pattern[df_pattern['pathway_id'].isin(highlight_ids)]

        if len(df_highlight) == 0:
            continue

        style = weights.get(pattern, {'linewidth': 1, 'alpha': 0.5, 'zorder': 1})
        color = PATTERN_COLORS.get(pattern, '#999999')

        for _, row in df_highlight.iterrows():
            if y_type == 'nes':
                y_values = [row[col] for col in nes_cols]
            else:
                rank_cols = [col.replace('NES_', 'Rank_') for col in nes_cols]
                y_values = [row[col] for col in rank_cols]

            # Plot highlighted line
            ax.plot(x_positions, y_values,
                    color=color,
                    linewidth=style['linewidth'] * 1.2,
                    alpha=1.0,
                    zorder=10,
                    solid_capstyle='round')

            # Create label at FIXED x position
            y_end = y_values[-1]
            desc = str(row.get('Description', row.get('pathway_id', '')))
            label = desc[:30] + '...' if len(desc) > 30 else desc

            txt = ax.text(label_x, y_end, label,
                         fontsize=6, color=color, ha='left', va='center',
                         fontweight='bold', zorder=11)
            label_texts.append(txt)
            label_points_y.append(y_end)  # Store original y for connector

    # Apply adjustText for VERTICAL-ONLY collision avoidance (no arrows)
    if use_adjust_text and label_texts:
        adjust_text(
            label_texts,
            x=[label_x] * len(label_texts),  # Original x = label x
            y=label_points_y,
            ax=ax,
            arrowprops=None,  # NO arrows - we'll draw manually
            force_text=(0.0, 1.0),  # Only vertical force
            force_points=(0.0, 0.5),
            expand_text=(1.0, 1.4),
            only_move={'points': 'y', 'text': 'y'},  # ONLY vertical movement
            lim=500,
        )

    # Draw MANUAL connectors AFTER adjustText
    for txt, y_original in zip(label_texts, label_points_y):
        # Get adjusted label position
        label_pos = txt.get_position()
        label_y_adjusted = label_pos[1]

        # Get label color from text
        label_color = txt.get_color()

        # Draw connector from line endpoint to label
        ax.plot([1.02, label_x - 0.02], [y_original, label_y_adjusted],
                color='gray', lw=0.5, alpha=0.5, zorder=9)

    # Add reference lines
    if y_type == 'nes':
        ax.axhline(y=0, color='black', linewidth=0.5, linestyle='-', alpha=0.3, zorder=0)
        ax.axhline(y=NES_EFFECT, color='gray', linewidth=0.5, linestyle='--', alpha=0.3, zorder=0)
        ax.axhline(y=-NES_EFFECT, color='gray', linewidth=0.5, linestyle='--', alpha=0.3, zorder=0)

        y_max = max(4, df_plot[nes_cols].max().max() * 1.1)
        y_min = min(-4, df_plot[nes_cols].min().min() * 1.1)
        ax.set_ylim(y_min, y_max)
    else:
        n_pathways = len(df_plot)
        ax.axhline(y=n_pathways/2, color='black', linewidth=0.5, linestyle='-', alpha=0.3, zorder=0)
        ax.set_ylim(0, n_pathways + 1)
        ax.invert_yaxis()

    # Formatting - Wide for labels at x=2
    ax.set_xticks(x_positions)
    ax.set_xticklabels(x_labels, fontsize=10)
    ax.set_xlim(-0.05, 3.5)

    n_highlighted = len(df_plot[df_plot['pathway_id'].isin(highlight_ids)])
    ax.set_title(f'{mutation}\n(n={len(df_plot)}, {n_highlighted} highlighted)',
                 fontsize=12, fontweight='bold',
                 color=MUTATION_COLORS.get(mutation, 'black'))

    ax.grid(axis='y', alpha=0.3, linestyle='-', linewidth=0.5)
    ax.set_axisbelow(True)


def create_faceted_bump_highlight_v9(df, y_type='nes', scope='focused'):
    """
    Create side-by-side highlighted bump charts with FIXED CONNECTORS (Pass 9).

    Improvements:
    - Labels further right (no overlap with lines)
    - Manual horizontal connectors (fixed position)
    - No adjustText arrows (cleaner appearance)
    """
    # Filter data by scope
    df_scope = filter_by_scope(df, scope)

    if len(df_scope) == 0:
        print(f"  Warning: No data for scope '{scope}'")
        return None

    # Compute global weights
    global_weights = compute_global_line_weights(df_scope, mutations=['G32A', 'R403C'])

    # WIDER figure for more label space
    fig, axes = plt.subplots(1, 3, figsize=(20, 8),
                              gridspec_kw={'width_ratios': [4.5, 4.5, 1.5]})

    # Plot G32A
    create_bump_chart_highlight_v9(df_scope, 'G32A', axes[0], y_type=y_type, scope=scope,
                                   global_weights=global_weights)
    axes[0].set_ylabel('NES' if y_type == 'nes' else 'Rank', fontsize=11)

    # Plot R403C
    create_bump_chart_highlight_v9(df_scope, 'R403C', axes[1], y_type=y_type, scope=scope,
                                   global_weights=global_weights)

    # Share y-axis limits for NES
    if y_type == 'nes':
        y_min = min(axes[0].get_ylim()[0], axes[1].get_ylim()[0])
        y_max = max(axes[0].get_ylim()[1], axes[1].get_ylim()[1])
        axes[0].set_ylim(y_min, y_max)
        axes[1].set_ylim(y_min, y_max)

    # Unified legend
    create_legend_unified(axes[2], df_scope, mutations=['G32A', 'R403C'])

    # Main title
    scope_name = SCOPE_CONFIGS[scope]['name']
    y_label = 'NES Value' if y_type == 'nes' else 'Rank (by NES)'
    fig.suptitle(f'Bump Chart (Fixed Connectors) | {scope_name}\nY-axis: {y_label} | Labels right of lines',
                 fontsize=14, fontweight='bold', y=0.98)

    plt.tight_layout(rect=[0, 0, 1, 0.94])

    return fig


def create_faceted_bump_highlight_v10(df, y_type='nes', scope='focused'):
    """
    Create side-by-side highlighted bump charts with DIAGONAL CONNECTORS (Pass 10).

    Improvements over v9:
    - Labels at x=2 (far right, clear separation from lines)
    - Diagonal connectors via adjustText (follow label movement)
    - No manual horizontal connectors
    - Left-aligned labels (no ladder pattern)
    """
    # Filter data by scope
    df_scope = filter_by_scope(df, scope)

    if len(df_scope) == 0:
        print(f"  Warning: No data for scope '{scope}'")
        return None

    # Compute global weights
    global_weights = compute_global_line_weights(df_scope, mutations=['G32A', 'R403C'])

    # EVEN WIDER figure for labels at x=2
    fig, axes = plt.subplots(1, 3, figsize=(24, 8),
                              gridspec_kw={'width_ratios': [5.5, 5.5, 1.5]})

    # Plot G32A
    create_bump_chart_highlight_v10(df_scope, 'G32A', axes[0], y_type=y_type, scope=scope,
                                   global_weights=global_weights)
    axes[0].set_ylabel('NES' if y_type == 'nes' else 'Rank', fontsize=11)

    # Plot R403C
    create_bump_chart_highlight_v10(df_scope, 'R403C', axes[1], y_type=y_type, scope=scope,
                                   global_weights=global_weights)

    # Share y-axis limits for NES
    if y_type == 'nes':
        y_min = min(axes[0].get_ylim()[0], axes[1].get_ylim()[0])
        y_max = max(axes[0].get_ylim()[1], axes[1].get_ylim()[1])
        axes[0].set_ylim(y_min, y_max)
        axes[1].set_ylim(y_min, y_max)

    # Unified legend
    create_legend_unified(axes[2], df_scope, mutations=['G32A', 'R403C'])

    # Main title
    scope_name = SCOPE_CONFIGS[scope]['name']
    y_label = 'NES Value' if y_type == 'nes' else 'Rank (by NES)'
    fig.suptitle(f'Bump Chart (Diagonal Connectors) | {scope_name}\nY-axis: {y_label} | Labels far right',
                 fontsize=14, fontweight='bold', y=0.98)

    plt.tight_layout(rect=[0, 0, 1, 0.94])

    return fig


def create_faceted_bump_highlight_v11(df, y_type='nes', scope='focused'):
    """
    Create side-by-side highlighted bump charts with MANUAL CONNECTORS (Pass 11).

    Improvements over v10:
    - Labels at fixed x=2.0 (no horizontal movement)
    - adjustText for vertical-only collision avoidance
    - Manual connectors drawn AFTER adjustText
    - Diagonal connectors from line endpoints to adjusted label positions
    """
    # Filter data by scope
    df_scope = filter_by_scope(df, scope)

    if len(df_scope) == 0:
        print(f"  Warning: No data for scope '{scope}'")
        return None

    # Compute global weights
    global_weights = compute_global_line_weights(df_scope, mutations=['G32A', 'R403C'])

    # Wide figure for labels at x=2
    fig, axes = plt.subplots(1, 3, figsize=(24, 8),
                              gridspec_kw={'width_ratios': [5.5, 5.5, 1.5]})

    # Plot G32A
    create_bump_chart_highlight_v11(df_scope, 'G32A', axes[0], y_type=y_type, scope=scope,
                                   global_weights=global_weights)
    axes[0].set_ylabel('NES' if y_type == 'nes' else 'Rank', fontsize=11)

    # Plot R403C
    create_bump_chart_highlight_v11(df_scope, 'R403C', axes[1], y_type=y_type, scope=scope,
                                   global_weights=global_weights)

    # Share y-axis limits for NES
    if y_type == 'nes':
        y_min = min(axes[0].get_ylim()[0], axes[1].get_ylim()[0])
        y_max = max(axes[0].get_ylim()[1], axes[1].get_ylim()[1])
        axes[0].set_ylim(y_min, y_max)
        axes[1].set_ylim(y_min, y_max)

    # Unified legend
    create_legend_unified(axes[2], df_scope, mutations=['G32A', 'R403C'])

    # Main title
    scope_name = SCOPE_CONFIGS[scope]['name']
    y_label = 'NES Value' if y_type == 'nes' else 'Rank (by NES)'
    fig.suptitle(f'Bump Chart (Manual Connectors) | {scope_name}\nY-axis: {y_label} | Labels at x=2, vertical adjust',
                 fontsize=14, fontweight='bold', y=0.98)

    plt.tight_layout(rect=[0, 0, 1, 0.94])

    return fig


def create_faceted_bump_highlight_v8(df, y_type='nes', scope='focused'):
    """
    Create side-by-side highlighted bump charts with IMPROVED LAYOUT (Pass 8).

    Improvements:
    - Lines start closer to y-axis
    - Labels positioned further right with better separation
    - Wider figure for more label space
    """
    # Filter data by scope
    df_scope = filter_by_scope(df, scope)

    if len(df_scope) == 0:
        print(f"  Warning: No data for scope '{scope}'")
        return None

    # Compute global weights
    global_weights = compute_global_line_weights(df_scope, mutations=['G32A', 'R403C'])

    # WIDER figure for label space
    fig, axes = plt.subplots(1, 3, figsize=(18, 8),
                              gridspec_kw={'width_ratios': [4, 4, 1.5]})

    # Plot G32A
    create_bump_chart_highlight_v8(df_scope, 'G32A', axes[0], y_type=y_type, scope=scope,
                                   global_weights=global_weights)
    axes[0].set_ylabel('NES' if y_type == 'nes' else 'Rank', fontsize=11)

    # Plot R403C
    create_bump_chart_highlight_v8(df_scope, 'R403C', axes[1], y_type=y_type, scope=scope,
                                   global_weights=global_weights)

    # Share y-axis limits for NES
    if y_type == 'nes':
        y_min = min(axes[0].get_ylim()[0], axes[1].get_ylim()[0])
        y_max = max(axes[0].get_ylim()[1], axes[1].get_ylim()[1])
        axes[0].set_ylim(y_min, y_max)
        axes[1].set_ylim(y_min, y_max)

    # Unified legend
    create_legend_unified(axes[2], df_scope, mutations=['G32A', 'R403C'])

    # Main title
    scope_name = SCOPE_CONFIGS[scope]['name']
    y_label = 'NES Value' if y_type == 'nes' else 'Rank (by NES)'
    fig.suptitle(f'Bump Chart (Improved Layout) | {scope_name}\nY-axis: {y_label} | Lines left, labels right',
                 fontsize=14, fontweight='bold', y=0.98)

    plt.tight_layout(rect=[0, 0, 1, 0.94])

    return fig


def create_faceted_bump_highlight_original(df, y_type='nes', scope='focused'):
    """
    Create side-by-side highlighted bump charts for G32A and R403C.
    ORIGINAL version (Pass 6) - uses per-panel weights and simple labels.

    This is preserved as p6 for comparison with the refined p7 version.
    """
    # Filter data by scope
    df_scope = filter_by_scope(df, scope)

    if len(df_scope) == 0:
        print(f"  Warning: No data for scope '{scope}'")
        return None

    # Create figure with 3 columns: G32A, R403C, Legend
    fig, axes = plt.subplots(1, 3, figsize=(16, 8),
                              gridspec_kw={'width_ratios': [3.5, 3.5, 1.5]})

    # Plot G32A with per-panel weights (original behavior)
    create_bump_chart_highlight(df_scope, 'G32A', axes[0], y_type=y_type, scope=scope,
                                global_weights=None)
    axes[0].set_ylabel('NES' if y_type == 'nes' else 'Rank', fontsize=11)

    # Plot R403C with per-panel weights
    create_bump_chart_highlight(df_scope, 'R403C', axes[1], y_type=y_type, scope=scope,
                                global_weights=None)

    # Share y-axis limits for NES
    if y_type == 'nes':
        y_min = min(axes[0].get_ylim()[0], axes[1].get_ylim()[0])
        y_max = max(axes[0].get_ylim()[1], axes[1].get_ylim()[1])
        axes[0].set_ylim(y_min, y_max)
        axes[1].set_ylim(y_min, y_max)

    # Legend - use original single-mutation legend (G32A only)
    create_legend_refined(axes[2], df_scope, 'G32A')

    # Main title
    scope_name = SCOPE_CONFIGS[scope]['name']
    y_label = 'NES Value' if y_type == 'nes' else 'Rank (by NES)'
    fig.suptitle(f'Highlighted Bump Chart (Original) | {scope_name}\nY-axis: {y_label} | Top 5 per pattern labeled',
                 fontsize=14, fontweight='bold', y=0.98)

    plt.tight_layout(rect=[0, 0, 1, 0.94])

    return fig


def create_faceted_bump_highlight(df, y_type='nes', scope='focused'):
    """
    Create side-by-side highlighted bump charts for G32A and R403C.
    REFINED version (Pass 7) with all improvements:
    - Biological relevance filtering
    - Unified legend (both mutations)
    - Global z-order (consistent layering)
    - Label collision avoidance (adjustText)

    Uses GLOBAL line weights for consistent z-order across both panels.
    """
    # Filter data by scope
    df_scope = filter_by_scope(df, scope)

    if len(df_scope) == 0:
        print(f"  Warning: No data for scope '{scope}'")
        return None

    # Compute GLOBAL weights for consistent layering across both panels
    global_weights = compute_global_line_weights(df_scope, mutations=['G32A', 'R403C'])

    # Create figure with 3 columns: G32A, R403C, Legend
    fig, axes = plt.subplots(1, 3, figsize=(16, 8),
                              gridspec_kw={'width_ratios': [3.5, 3.5, 1.5]})

    # Plot G32A with global weights
    create_bump_chart_highlight(df_scope, 'G32A', axes[0], y_type=y_type, scope=scope,
                                global_weights=global_weights)
    axes[0].set_ylabel('NES' if y_type == 'nes' else 'Rank', fontsize=11)

    # Plot R403C with global weights
    create_bump_chart_highlight(df_scope, 'R403C', axes[1], y_type=y_type, scope=scope,
                                global_weights=global_weights)

    # Share y-axis limits for NES
    if y_type == 'nes':
        y_min = min(axes[0].get_ylim()[0], axes[1].get_ylim()[0])
        y_max = max(axes[0].get_ylim()[1], axes[1].get_ylim()[1])
        axes[0].set_ylim(y_min, y_max)
        axes[1].set_ylim(y_min, y_max)

    # Legend - use unified legend showing both mutations
    create_legend_unified(axes[2], df_scope, mutations=['G32A', 'R403C'])

    # Main title
    scope_name = SCOPE_CONFIGS[scope]['name']
    y_label = 'NES Value' if y_type == 'nes' else 'Rank (by NES)'
    fig.suptitle(f'Highlighted Bump Chart | {scope_name}\nY-axis: {y_label} | Top 5 significant per pattern labeled',
                 fontsize=14, fontweight='bold', y=0.98)

    plt.tight_layout(rect=[0, 0, 1, 0.94])

    return fig


# =============================================================================
# MAIN EXECUTION
# =============================================================================

def main():
    """Generate all bump chart variants."""
    print("=" * 80)
    print("BUMP CHART VISUALIZATION")
    print("=" * 80)
    print(f"Output directory: {OUTPUT_DIR}")

    # Load data
    df = load_data()

    # =========================================================================
    # PASS 1: Original (uniform line thickness)
    # Naming: bump_p1_uniform_{ytype}_{scope}
    # Key feature: Uniform line thickness, no weighting
    # =========================================================================
    print("\n" + "=" * 60)
    print("PASS 1: Original bump charts (uniform lines)")
    print("  Naming: bump_p1_uniform_{ytype}_{scope}")
    print("=" * 60)

    variants = [
        ('nes', 'focused'), ('rank', 'focused'),
        ('nes', 'significant'), ('rank', 'significant'),
        ('nes', 'all'), ('rank', 'all'),
    ]
    for y_type, scope in variants:
        print(f"\n  Creating: p1_uniform_{y_type}_{scope}")
        fig = create_faceted_bump_chart(df, y_type=y_type, scope=scope)
        if fig is not None:
            filename = f'bump_p1_uniform_{y_type}_{scope}'
            for ext, dpi in [('pdf', 300), ('png', 150)]:
                output_file = OUTPUT_DIR / f'{filename}.{ext}'
                fig.savefig(output_file, dpi=dpi, bbox_inches='tight', facecolor='white')
                print(f"    Saved: {output_file.name}")
            plt.close(fig)

    # =========================================================================
    # PASS 2: Weighted + 3-Point Trajectory (Early → TrajDev → Late)
    # Naming: bump_p2_weighted3pt_{ytype}_{scope}
    # Key feature: Frequency-based line weighting, includes TrajDev middle point
    # =========================================================================
    print("\n" + "=" * 60)
    print("PASS 2: Weighted bump charts (3-point, with TrajDev)")
    print("  Naming: bump_p2_weighted3pt_{ytype}_{scope}")
    print("=" * 60)

    for scope in ['focused', 'significant', 'all']:
        for y_type in ['nes', 'rank']:
            print(f"\n  Creating: p2_weighted3pt_{y_type}_{scope}")

            fig = create_faceted_bump_weighted(df, y_type=y_type, scope=scope,
                                               include_trajdev=True)

            if fig is not None:
                filename = f'bump_p2_weighted3pt_{y_type}_{scope}'

                for ext, dpi in [('pdf', 300), ('png', 150)]:
                    output_file = OUTPUT_DIR / f'{filename}.{ext}'
                    fig.savefig(output_file, dpi=dpi, bbox_inches='tight', facecolor='white')
                    print(f"    Saved: {output_file.name}")

                plt.close(fig)

    # =========================================================================
    # PASS 3: Weighted + 2-Point Direct Trajectory (Early → Late)
    # Naming: bump_p3_direct2pt_{ytype}_{scope}
    # Key feature: 2-point direct trajectory (no TrajDev), cleaner visualization
    # =========================================================================
    print("\n" + "=" * 60)
    print("PASS 3: Weighted bump charts (2-point direct, no TrajDev)")
    print("  Naming: bump_p3_direct2pt_{ytype}_{scope}")
    print("=" * 60)

    for scope in ['focused', 'significant', 'all']:
        for y_type in ['nes', 'rank']:
            print(f"\n  Creating: p3_direct2pt_{y_type}_{scope}")

            fig = create_faceted_bump_weighted(df, y_type=y_type, scope=scope,
                                               include_trajdev=False)

            if fig is not None:
                filename = f'bump_p3_direct2pt_{y_type}_{scope}'

                for ext, dpi in [('pdf', 300), ('png', 150)]:
                    output_file = OUTPUT_DIR / f'{filename}.{ext}'
                    fig.savefig(output_file, dpi=dpi, bbox_inches='tight', facecolor='white')
                    print(f"    Saved: {output_file.name}")

                plt.close(fig)

    # =========================================================================
    # PASS 4: Summary Direct (2-point, no TrajDev)
    # Naming: bump_p4_summary_{scope}
    # Key feature: Median + IQR ribbons per pattern (for dense data)
    # =========================================================================
    print("\n" + "=" * 60)
    print("PASS 4: Summary Direct (2-point, no TrajDev)")
    print("  Naming: bump_p4_summary_{scope}")
    print("=" * 60)

    for scope in ['significant', 'all']:
        print(f"\n  Creating: p4_summary_{scope}")

        fig = create_summary_faceted(df, scope=scope, include_trajdev=False)

        if fig is not None:
            filename = f'bump_p4_summary_{scope}'

            for ext, dpi in [('pdf', 300), ('png', 150)]:
                output_file = OUTPUT_DIR / f'{filename}.{ext}'
                fig.savefig(output_file, dpi=dpi, bbox_inches='tight', facecolor='white')
                print(f"    Saved: {output_file.name}")

            plt.close(fig)

    # =========================================================================
    # PASS 5: Refined Weighting (2-point, inverse alpha)
    # Naming: bump_p5_inversealpha_{ytype}_{scope}
    # Key feature: Inverse alpha (thicker=more opaque), adjusted linewidths
    # =========================================================================
    print("\n" + "=" * 60)
    print("PASS 5: Refined bump charts (2-point, inverse alpha)")
    print("  Naming: bump_p5_inversealpha_{ytype}_{scope}")
    print("=" * 60)

    for scope in ['focused', 'significant', 'all']:
        for y_type in ['nes', 'rank']:
            print(f"\n  Creating: p5_inversealpha_{y_type}_{scope}")

            fig = create_faceted_bump_refined(df, y_type=y_type, scope=scope)

            if fig is not None:
                filename = f'bump_p5_inversealpha_{y_type}_{scope}'

                for ext, dpi in [('pdf', 300), ('png', 150)]:
                    output_file = OUTPUT_DIR / f'{filename}.{ext}'
                    fig.savefig(output_file, dpi=dpi, bbox_inches='tight', facecolor='white')
                    print(f"    Saved: {output_file.name}")

                plt.close(fig)

    # =========================================================================
    # PASS 6: Highlighted Significant (2-point, labels) - ORIGINAL
    # Naming: bump_p6_labeled_{ytype}_{scope}
    # Key feature: Top 5 significant pathways per pattern highlighted + labeled
    # Note: This is the baseline version before Pass 7 refinements
    # =========================================================================
    print("\n" + "=" * 60)
    print("PASS 6: Highlighted bump charts (significant pathways labeled)")
    print("  Naming: bump_p6_labeled_{ytype}_{scope}")
    print("  Note: Original version - see Pass 7 for refined version")
    print("=" * 60)

    for scope in ['focused', 'significant', 'all']:
        for y_type in ['nes', 'rank']:
            print(f"\n  Creating: p6_labeled_{y_type}_{scope}")

            # Use per-panel weights (original behavior) for p6
            fig = create_faceted_bump_highlight_original(df, y_type=y_type, scope=scope)

            if fig is not None:
                filename = f'bump_p6_labeled_{y_type}_{scope}'

                for ext, dpi in [('pdf', 300), ('png', 150)]:
                    output_file = OUTPUT_DIR / f'{filename}.{ext}'
                    fig.savefig(output_file, dpi=dpi, bbox_inches='tight', facecolor='white')
                    print(f"    Saved: {output_file.name}")

                plt.close(fig)

    # =========================================================================
    # PASS 7: Refined Highlighted (2-point, all fixes applied)
    # Naming: bump_p7_refined_{ytype}_{scope}
    # Key features:
    #   - 7D: Biological relevance filtering (exclude irrelevant pathways)
    #   - 7A: Unified legend (counts from both mutations)
    #   - 7B: Global z-order (consistent layering across panels)
    #   - 7C: Label collision avoidance (adjustText)
    # =========================================================================
    print("\n" + "=" * 60)
    print("PASS 7: Refined bump charts (all improvements)")
    print("  Naming: bump_p7_refined_{ytype}_{scope}")
    print("  Fixes: biological filtering, unified legend, global z-order, label collision")
    print("=" * 60)

    for scope in ['focused', 'significant', 'all']:
        for y_type in ['nes', 'rank']:
            print(f"\n  Creating: p7_refined_{y_type}_{scope}")

            fig = create_faceted_bump_highlight(df, y_type=y_type, scope=scope)

            if fig is not None:
                filename = f'bump_p7_refined_{y_type}_{scope}'

                for ext, dpi in [('pdf', 300), ('png', 150)]:
                    output_file = OUTPUT_DIR / f'{filename}.{ext}'
                    fig.savefig(output_file, dpi=dpi, bbox_inches='tight', facecolor='white')
                    print(f"    Saved: {output_file.name}")

                plt.close(fig)

    # =========================================================================
    # PASS 8: Improved Layout (lines left, labels right)
    # Naming: bump_p8_layout_{ytype}_{scope}
    # Key features:
    #   - Lines start closer to y-axis (clear NES correspondence)
    #   - Labels positioned further right (no overlap with lines)
    #   - Wider figure for better label spacing
    # =========================================================================
    print("\n" + "=" * 60)
    print("PASS 8: Improved layout (lines left, labels right)")
    print("  Naming: bump_p8_layout_{ytype}_{scope}")
    print("  Key: tighter left margin, labels further right")
    print("=" * 60)

    for scope in ['focused', 'significant', 'all']:
        for y_type in ['nes', 'rank']:
            print(f"\n  Creating: p8_layout_{y_type}_{scope}")

            fig = create_faceted_bump_highlight_v8(df, y_type=y_type, scope=scope)

            if fig is not None:
                filename = f'bump_p8_layout_{y_type}_{scope}'

                for ext, dpi in [('pdf', 300), ('png', 150)]:
                    output_file = OUTPUT_DIR / f'{filename}.{ext}'
                    fig.savefig(output_file, dpi=dpi, bbox_inches='tight', facecolor='white')
                    print(f"    Saved: {output_file.name}")

                plt.close(fig)

    # =========================================================================
    # PASS 9: Fixed Connectors (labels right, manual connectors)
    # Naming: bump_p9_fixconn_{ytype}_{scope}
    # Key features:
    #   - Labels moved further right (x=1.25) - no overlap with lines
    #   - Manual horizontal connectors (fixed position)
    #   - No adjustText arrows (cleaner appearance)
    # =========================================================================
    print("\n" + "=" * 60)
    print("PASS 9: Fixed connectors (labels right, manual connectors)")
    print("  Naming: bump_p9_fixconn_{ytype}_{scope}")
    print("  Key: labels at x=1.25, horizontal connectors 1.02->1.23")
    print("=" * 60)

    for scope in ['focused', 'significant', 'all']:
        for y_type in ['nes', 'rank']:
            print(f"\n  Creating: p9_fixconn_{y_type}_{scope}")

            fig = create_faceted_bump_highlight_v9(df, y_type=y_type, scope=scope)

            if fig is not None:
                filename = f'bump_p9_fixconn_{y_type}_{scope}'

                for ext, dpi in [('pdf', 300), ('png', 150)]:
                    output_file = OUTPUT_DIR / f'{filename}.{ext}'
                    fig.savefig(output_file, dpi=dpi, bbox_inches='tight', facecolor='white')
                    print(f"    Saved: {output_file.name}")

                plt.close(fig)

    # =========================================================================
    # PASS 10: Diagonal Connectors (labels far right, adjustText arrows)
    # Naming: bump_p10_diagconn_{ytype}_{scope}
    # Key features:
    #   - Labels at x=2 (far right, clear separation)
    #   - Diagonal connectors via adjustText (follow label movement)
    #   - No manual horizontal connectors
    #   - Left-aligned labels (no ladder pattern)
    # =========================================================================
    print("\n" + "=" * 60)
    print("PASS 10: Diagonal connectors (labels far right, adjustText arrows)")
    print("  Naming: bump_p10_diagconn_{ytype}_{scope}")
    print("  Key: labels at x=2, diagonal connectors from adjustText")
    print("=" * 60)

    for scope in ['focused', 'significant', 'all']:
        for y_type in ['nes', 'rank']:
            print(f"\n  Creating: p10_diagconn_{y_type}_{scope}")

            fig = create_faceted_bump_highlight_v10(df, y_type=y_type, scope=scope)

            if fig is not None:
                filename = f'bump_p10_diagconn_{y_type}_{scope}'

                for ext, dpi in [('pdf', 300), ('png', 150)]:
                    output_file = OUTPUT_DIR / f'{filename}.{ext}'
                    fig.savefig(output_file, dpi=dpi, bbox_inches='tight', facecolor='white')
                    print(f"    Saved: {output_file.name}")

                plt.close(fig)

    # =========================================================================
    # PASS 11: Manual Connectors (labels fixed at x=2, connectors after adjustText)
    # Naming: bump_p11_manual_{ytype}_{scope}
    # Key features:
    #   - Labels at fixed x=2.0 (no horizontal movement)
    #   - adjustText for vertical-only collision avoidance
    #   - Manual connectors drawn AFTER adjustText
    #   - Diagonal connectors from line endpoints to adjusted label positions
    # =========================================================================
    print("\n" + "=" * 60)
    print("PASS 11: Manual connectors (labels fixed, connectors after adjustText)")
    print("  Naming: bump_p11_manual_{ytype}_{scope}")
    print("  Key: labels at x=2, manual diagonal connectors after vertical adjust")
    print("=" * 60)

    for scope in ['focused', 'significant', 'all']:
        for y_type in ['nes', 'rank']:
            print(f"\n  Creating: p11_manual_{y_type}_{scope}")

            fig = create_faceted_bump_highlight_v11(df, y_type=y_type, scope=scope)

            if fig is not None:
                filename = f'bump_p11_manual_{y_type}_{scope}'

                for ext, dpi in [('pdf', 300), ('png', 150)]:
                    output_file = OUTPUT_DIR / f'{filename}.{ext}'
                    fig.savefig(output_file, dpi=dpi, bbox_inches='tight', facecolor='white')
                    print(f"    Saved: {output_file.name}")

                plt.close(fig)

    # Generate README
    generate_readme()

    print("\n" + "=" * 80)
    print("COMPLETE!")
    print("=" * 80)
    print(f"\nOutput directory: {OUTPUT_DIR}")
    print("\nGenerated bump chart files:")
    for f in sorted(OUTPUT_DIR.glob('bump_*')):
        print(f"  - {f.name}")


def generate_readme():
    """Generate README documentation for bump charts."""
    readme_content = """# Bump Chart Visualizations

## Purpose

Bump charts (slope graphs) show **individual pathway trajectories** through the
Early → TrajDev → Late stages. Unlike alluvial diagrams that show categorical flows,
bump charts reveal:

1. **Actual NES values** or ranks at each stage
2. **Pattern-specific trajectory shapes** (e.g., Compensation vs Natural_improvement)
3. **Individual pathway journeys** not just aggregate counts

## Naming Convention

Files are named with chronological pass numbers and metadata:
`bump_p{N}_{metadata}_{ytype}_{scope}.{ext}`

| Pass | Prefix | Key Change |
|------|--------|------------|
| 1 | `p1_uniform` | Uniform line thickness (baseline) |
| 2 | `p2_weighted3pt` | Frequency-weighted lines, 3-point (with TrajDev) |
| 3 | `p3_direct2pt` | 2-point direct trajectory (Early→Late only) |
| 4 | `p4_summary` | Median + IQR ribbons (for dense data) |
| 5 | `p5_inversealpha` | Inverse alpha (thicker = more opaque) |
| 6 | `p6_labeled` | Top 5 significant per pattern labeled (original) |
| 7 | `p7_refined` | All refinements (filtering, legend, z-order) |
| 8 | `p8_layout` | Improved layout (lines left, labels right) |
| 9 | `p9_fixconn` | Fixed connectors (no label-line overlap) |
| 10 | `p10_diagconn` | Diagonal connectors, labels far right |
| 11 | `p11_manual` | **RECOMMENDED** - Manual connectors after adjustText |

## Y-axis Types

| Type | Description | Best For |
|------|-------------|----------|
| **nes** | Continuous enrichment score | Seeing actual magnitude changes |
| **rank** | Ordinal position (1=most negative) | Seeing relative position changes |

## Scopes

| Scope | Pathways | Description |
|-------|----------|-------------|
| **focused** | ~100 | MitoCarta + SynGO only (domain-specific) |
| **significant** | ~2-4k | At least one stage p<0.05 |
| **all** | ~12k | All tested pathways |

## File Reference

### Pass 1: Uniform (baseline)
| File | Description |
|------|-------------|
| `bump_p1_uniform_nes_{scope}` | Uniform lines, NES on Y-axis |
| `bump_p1_uniform_rank_{scope}` | Uniform lines, Rank on Y-axis |

### Pass 2: Weighted 3-Point (with TrajDev)
| File | Description |
|------|-------------|
| `bump_p2_weighted3pt_nes_{scope}` | Weighted lines, 3-point trajectory, NES |
| `bump_p2_weighted3pt_rank_{scope}` | Weighted lines, 3-point trajectory, Rank |

### Pass 3: Direct 2-Point (Early→Late)
| File | Description |
|------|-------------|
| `bump_p3_direct2pt_nes_{scope}` | Weighted lines, 2-point direct, NES |
| `bump_p3_direct2pt_rank_{scope}` | Weighted lines, 2-point direct, Rank |

### Pass 4: Summary Ribbons
| File | Description |
|------|-------------|
| `bump_p4_summary_{scope}` | Median line + IQR ribbon per pattern |

### Pass 5: Inverse Alpha
| File | Description |
|------|-------------|
| `bump_p5_inversealpha_nes_{scope}` | Inverse alpha weighting, NES |
| `bump_p5_inversealpha_rank_{scope}` | Inverse alpha weighting, Rank |

### Pass 6: Labeled Highlights (Original)
| File | Description |
|------|-------------|
| `bump_p6_labeled_nes_{scope}` | Top 5 significant labeled, NES |
| `bump_p6_labeled_rank_{scope}` | Top 5 significant labeled, Rank |

### Pass 7: Refined Highlights (RECOMMENDED)
| File | Description |
|------|-------------|
| `bump_p7_refined_nes_{scope}` | All refinements applied, NES |
| `bump_p7_refined_rank_{scope}` | All refinements applied, Rank |

**Pass 7 Improvements:**
- **Biological filtering**: Excludes irrelevant pathways (cardiac, muscle, cancer, etc.)
- **Priority categories**: Prioritizes mito, synapse, translation pathways
- **Unified legend**: Shows counts from both G32A and R403C mutations
- **Global z-order**: Consistent layering across both mutation panels
- **Label collision**: Uses adjustText for non-overlapping labels

### Pass 8: Improved Layout
| File | Description |
|------|-------------|
| `bump_p8_layout_nes_{scope}` | Lines left, labels right, NES |
| `bump_p8_layout_rank_{scope}` | Lines left, labels right, Rank |

**Pass 8 Improvements (over Pass 7):**
- **Lines closer to y-axis**: Clear correspondence between line start and NES value
- **Labels further right**: No overlap between labels and trajectory lines
- **Wider figure**: More space for labels, cleaner appearance

### Pass 9: Fixed Connectors
| File | Description |
|------|-------------|
| `bump_p9_fixconn_nes_{scope}` | Fixed connectors, no label-line overlap, NES |
| `bump_p9_fixconn_rank_{scope}` | Fixed connectors, no label-line overlap, Rank |

**Pass 9 Improvements (over Pass 8):**
- **Labels further right (x=1.25)**: Clear separation from trajectory line endpoints
- **Manual fixed connectors**: Horizontal lines from x=1.02 to x=1.23, stay fixed even when labels shift
- **No whitespace gaps**: Connectors bridge the gap between line endpoints and labels
- **adjustText without arrows**: Labels still avoid collision but don't drag connectors around
- **Wider figure (20in)**: More space for label region

### Pass 10: Diagonal Connectors
| File | Description |
|------|-------------|
| `bump_p10_diagconn_nes_{scope}` | Diagonal connectors, labels far right, NES |
| `bump_p10_diagconn_rank_{scope}` | Diagonal connectors, labels far right, Rank |

**Pass 10 Improvements (over Pass 9):**
- **Labels far right (x=2)**: Completely clear of trajectory lines
- **Diagonal connectors**: adjustText draws connectors that follow label movement
- **Left-aligned labels**: All labels start at same x position (no ladder pattern)
- **Wider figure (24in)**: Accommodates labels at x=2

### Pass 11: Manual Connectors (RECOMMENDED)
| File | Description |
|------|-------------|
| `bump_p11_manual_nes_{scope}` | Manual diagonal connectors, labels at x=2, NES |
| `bump_p11_manual_rank_{scope}` | Manual diagonal connectors, labels at x=2, Rank |

**Pass 11 Improvements (over Pass 10):**
- **Labels at fixed x=2**: No horizontal movement (left-aligned, no ladder)
- **Vertical-only adjustText**: Labels can shift up/down to avoid collision
- **Manual connectors AFTER adjustText**: Drawn from line endpoints to adjusted label positions
- **Proper diagonal connectors**: From (x=1.02, y_original) to (x=1.98, y_adjusted)

## Interpretation Guide

### Reading the Plots

- **Left panel**: G32A mutation (GTPase domain)
- **Right panel**: R403C mutation (Stalk domain)
- **Each line**: One pathway's trajectory
- **Color**: Pattern classification

### Pattern Trajectory Shapes

| Pattern | Expected Shape |
|---------|----------------|
| **Compensation** | Starts away from 0 → TrajDev moves toward 0 → Ends near 0 |
| **Natural_improvement** | Starts away from 0 → Flat TrajDev → Ends nearer 0 |
| **Progressive** | Starts away from 0 → TrajDev moves further → Ends further from 0 |
| **Late_onset** | Starts near 0 → Ends away from 0 |

### NES vs Rank

- **NES plots**: Show actual magnitude. Compensation should show large |NES| at Early,
  opposing TrajDev, small |NES| at Late.
- **Rank plots**: Show relative ordering. Pathways that improve move toward middle ranks
  at Late stage.

## Color Legend

- Compensation: Bluish green (#009E73)
- Progressive: Vermillion (#D55E00)
- Natural_improvement: Sky blue (#56B4E9)
- Natural_worsening: Orange (#E69F00)
- Late_onset: Reddish purple (#CC79A7)
- Transient: Blue (#0072B2)
- Complex: Yellow (#F0E442)

## See Also

- `docs/PATTERN_CLASSIFICATION.md` - Pattern definitions
- `01_Scripts/Python/pattern_definitions.py` - Classification code
- `README_alluvial_goals.md` - Motivation for bump chart approach

---
Generated: """ + pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')

    with open(OUTPUT_DIR / 'README_bump_charts.md', 'w') as f:
        f.write(readme_content)
    print("\n  Saved: README_bump_charts.md")


if __name__ == '__main__':
    main()
