#!/usr/bin/env python3
"""
Normalized Pattern Summary Visualization
=========================================

Creates normalized (100% stacked) pattern summary figures to address the issue
that large databases (gobp with ~5000 pathways) mask the relative pattern
distributions in smaller databases.

Key improvements over original:
1. Each database bar normalized to 100% for proportion comparison
2. Absolute counts shown as text annotations (e.g., "Compensation: 45 (60%)")
3. Optional dual-panel view: normalized + absolute reference

This script does NOT replace the original pattern_summary figures - it provides
an alternative perspective focused on relative patterns.
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

# Import from project-specific modules
from Python.config import CONFIG, resolve_path
from Python.semantic_categories import PATTERN_COLORS, MUTATION_COLORS
from Python.patterns import add_pattern_classification
from Python.pattern_definitions import MEANINGFUL_PATTERNS

# =============================================================================
# SETUP
# =============================================================================

OUTPUT_DIR = resolve_path('03_Results/02_Analysis/Plots/Pattern_Summary_Normalized')
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# =============================================================================
# DATA LOADING
# =============================================================================

def load_classified_data():
    """Load pre-classified pathway data from authoritative master table"""
    print("Loading classified pathway data from master_gsea_table.csv...")

    # Use authoritative master table with significance-based patterns
    data_file = resolve_path('03_Results/02_Analysis/master_gsea_table.csv')
    df_long = pd.read_csv(data_file)

    # Pivot to wide format (one row per pathway)
    # Pattern columns are consistent across contrasts, so take first occurrence
    id_cols = ['pathway_id', 'database', 'Description']
    pattern_cols = ['Pattern_G32A', 'Confidence_G32A', 'Super_Category_G32A',
                    'Pattern_R403C', 'Confidence_R403C', 'Super_Category_R403C']
    nes_cols = ['NES_Early_G32A', 'NES_Early_R403C', 'NES_TrajDev_G32A',
                'NES_TrajDev_R403C', 'NES_Late_G32A', 'NES_Late_R403C']

    keep_cols = id_cols + pattern_cols + nes_cols
    available_cols = [c for c in keep_cols if c in df_long.columns]

    df = df_long[available_cols].drop_duplicates(subset=['pathway_id', 'database']).copy()

    print(f"  Loaded {len(df)} unique pathways from {df['database'].nunique()} databases")
    print(f"  Using significance-based pattern classifications")

    return df


def prepare_pattern_counts(df):
    """
    Prepare pattern count data for both mutations.

    Returns DataFrame with columns: Database, Mutation, Pattern, Count
    """
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

    return pd.DataFrame(pattern_counts)


# =============================================================================
# OPTION 1: NORMALIZED 100% STACKED BARS
# =============================================================================

def create_normalized_pattern_summary(df_counts, show_annotations=True):
    """
    Create normalized (100% stacked) pattern summary with absolute count annotations.

    Parameters
    ----------
    df_counts : pd.DataFrame
        Pattern counts with columns: Database, Mutation, Pattern, Count
    show_annotations : bool
        If True, add text annotations showing absolute counts and percentages
    """
    print("\nCreating normalized pattern summary (100% stacked bars)...")

    # Filter to meaningful patterns (exclude Complex for cleaner visualization)
    patterns_to_plot = [p for p in MEANINGFUL_PATTERNS if p != 'Complex']
    df_meaningful = df_counts[df_counts['Pattern'].isin(patterns_to_plot)].copy()

    if len(df_meaningful) == 0:
        print("  No meaningful patterns found")
        return None

    # Create figure
    fig, axes = plt.subplots(1, 2, figsize=(16, 8), sharey=True)

    # Get databases
    databases = sorted(df_meaningful['Database'].unique())
    n_databases = len(databases)

    for idx, (mutation, ax) in enumerate(zip(['G32A', 'R403C'], axes)):
        df_mut = df_meaningful[df_meaningful['Mutation'] == mutation].copy()

        # Pivot to get pattern counts per database
        pivot = df_mut.pivot_table(
            index='Database',
            columns='Pattern',
            values='Count',
            fill_value=0
        )

        # Reindex to ensure all databases present
        pivot = pivot.reindex(databases, fill_value=0)

        # Calculate totals and proportions
        totals = pivot.sum(axis=1)
        pivot_pct = pivot.div(totals, axis=0) * 100  # Convert to percentage

        # Order patterns
        pattern_order = [p for p in MEANINGFUL_PATTERNS if p in pivot.columns]
        pivot = pivot[pattern_order]
        pivot_pct = pivot_pct[pattern_order]

        # Plot normalized stacked bars
        x = np.arange(n_databases)
        left = np.zeros(n_databases)

        for pattern in pattern_order:
            if pattern not in pivot_pct.columns:
                continue

            color = PATTERN_COLORS.get(pattern, '#999999')
            widths = pivot_pct[pattern].values

            bars = ax.barh(x, widths, left=left,
                          color=color, alpha=0.85,
                          edgecolor='white', linewidth=0.5)

            # Add annotations if requested
            if show_annotations:
                for i, (width, count) in enumerate(zip(widths, pivot[pattern].values)):
                    if width > 5:  # Only annotate if segment is large enough
                        x_pos = left[i] + width / 2

                        # Format: "Pattern: N (X%)"
                        if count > 0:
                            label = f'{int(count)}\n({width:.0f}%)'
                            ax.text(x_pos, i, label,
                                   ha='center', va='center',
                                   fontsize=7, fontweight='bold',
                                   color='white' if width > 15 else 'black')

            left += widths

        # Add total counts at end of bars
        for i, (db, total) in enumerate(zip(databases, totals)):
            ax.text(101, i, f'n={int(total)}',
                   ha='left', va='center',
                   fontsize=8, fontweight='bold',
                   color=MUTATION_COLORS[mutation])

        # Formatting
        ax.set_yticks(x)
        ax.set_yticklabels(databases, fontsize=10)
        ax.set_xlabel('Proportion of Pathways (%)', fontsize=11)
        ax.set_xlim(0, 115)  # Extra space for total annotations

        ax.set_title(f'{mutation}',
                    fontsize=14, fontweight='bold',
                    color=MUTATION_COLORS[mutation])

        # Add vertical line at 50%
        ax.axvline(50, color='gray', linestyle='--', alpha=0.3, linewidth=1)

        if idx == 0:
            ax.set_ylabel('Database', fontsize=11)

    # Create shared legend
    legend_handles = [
        mpatches.Patch(color=PATTERN_COLORS[p], label=p, alpha=0.85)
        for p in MEANINGFUL_PATTERNS if p in df_meaningful['Pattern'].unique()
    ]
    fig.legend(handles=legend_handles,
              loc='lower center', ncol=len(legend_handles),
              fontsize=10, frameon=True, framealpha=0.9,
              bbox_to_anchor=(0.5, -0.10))

    # Title
    fig.suptitle('Normalized Pattern Distribution (100% Stacked)\n'
                'Each database scaled to 100% to show relative pattern proportions',
                fontsize=14, fontweight='bold', y=0.98)

    # Interpretation text
    interpretation = (
        'Note: Bars normalized to 100% per database. Absolute counts shown as "N (%)" inside bars.\n'
        'Total pathways per database shown as "n=X" at right. Enables comparison of pattern proportions\n'
        'across databases of different sizes (e.g., gobp with 5000 vs hallmark with 10 pathways).'
    )
    fig.text(0.5, 0.01, interpretation,
            ha='center', va='bottom', fontsize=9, style='italic',
            bbox=dict(boxstyle='round', facecolor='#FFF9E6',
                     edgecolor='orange', alpha=0.8))

    plt.tight_layout(rect=[0, 0.13, 1, 0.96])

    return fig


# =============================================================================
# OPTION 2: DUAL-PANEL COMPARISON
# =============================================================================

def create_dual_panel_comparison(df_counts, mutation='G32A'):
    """
    Create dual-panel figure showing both normalized and absolute counts.

    Parameters
    ----------
    df_counts : pd.DataFrame
        Pattern counts
    mutation : str
        Which mutation to show
    """
    print(f"\nCreating dual-panel comparison for {mutation}...")

    # Filter data (exclude Complex for cleaner visualization)
    patterns_to_plot = [p for p in MEANINGFUL_PATTERNS if p != 'Complex']
    df_mut = df_counts[
        (df_counts['Mutation'] == mutation) &
        (df_counts['Pattern'].isin(patterns_to_plot))
    ].copy()

    if len(df_mut) == 0:
        print(f"  No data for {mutation}")
        return None

    # Prepare data
    databases = sorted(df_mut['Database'].unique())
    pivot = df_mut.pivot_table(
        index='Database',
        columns='Pattern',
        values='Count',
        fill_value=0
    ).reindex(databases, fill_value=0)

    pattern_order = [p for p in MEANINGFUL_PATTERNS if p in pivot.columns]
    pivot = pivot[pattern_order]

    # Calculate proportions
    totals = pivot.sum(axis=1)
    pivot_pct = pivot.div(totals, axis=0) * 100

    # Create figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8), sharey=True)

    x = np.arange(len(databases))

    # Panel 1: Normalized (100% stacked)
    left = np.zeros(len(databases))
    for pattern in pattern_order:
        color = PATTERN_COLORS.get(pattern, '#999999')
        widths = pivot_pct[pattern].values
        ax1.barh(x, widths, left=left, color=color, alpha=0.85,
                edgecolor='white', linewidth=0.5, label=pattern)
        left += widths

    ax1.set_xlabel('Proportion (%)', fontsize=11)
    ax1.set_xlim(0, 100)
    ax1.set_title('A. Normalized Pattern Distribution',
                 fontsize=12, fontweight='bold')
    ax1.axvline(50, color='gray', linestyle='--', alpha=0.3)

    # Panel 2: Absolute counts (traditional stacked)
    left = np.zeros(len(databases))
    for pattern in pattern_order:
        color = PATTERN_COLORS.get(pattern, '#999999')
        counts = pivot[pattern].values
        ax2.barh(x, counts, left=left, color=color, alpha=0.85,
                edgecolor='white', linewidth=0.5)
        left += counts

    ax2.set_xlabel('Absolute Count', fontsize=11)
    ax2.set_title('B. Absolute Pathway Counts (Reference)',
                 fontsize=12, fontweight='bold')

    # Shared y-axis
    ax1.set_yticks(x)
    ax1.set_yticklabels(databases, fontsize=10)
    ax1.set_ylabel('Database', fontsize=11)

    # Legend
    ax1.legend(loc='center left', bbox_to_anchor=(1.02, 0.5),
              fontsize=9, frameon=True)

    fig.suptitle(f'{mutation} - Normalized vs Absolute Pattern Comparison',
                fontsize=14, fontweight='bold',
                color=MUTATION_COLORS[mutation], y=0.98)

    plt.tight_layout(rect=[0, 0.02, 1, 0.96])

    return fig


# =============================================================================
# OPTION 3: PERCENTAGE TABLE WITH HEATMAP
# =============================================================================

def create_percentage_heatmap(df_counts):
    """
    Create heatmap showing percentage of each pattern per database.

    Alternative visualization that shows all data in a compact format.
    """
    print("\nCreating percentage heatmap...")

    # Filter to meaningful patterns (exclude Complex for cleaner visualization)
    patterns_to_plot = [p for p in MEANINGFUL_PATTERNS if p != 'Complex']
    df_meaningful = df_counts[df_counts['Pattern'].isin(patterns_to_plot)].copy()

    fig, axes = plt.subplots(1, 2, figsize=(14, 8), sharey=True)

    for idx, (mutation, ax) in enumerate(zip(['G32A', 'R403C'], axes)):
        df_mut = df_meaningful[df_meaningful['Mutation'] == mutation]

        # Pivot and calculate percentages
        pivot = df_mut.pivot_table(
            index='Database',
            columns='Pattern',
            values='Count',
            fill_value=0
        )

        databases = sorted(pivot.index)
        pivot = pivot.reindex(databases, fill_value=0)
        pattern_order = [p for p in MEANINGFUL_PATTERNS if p in pivot.columns]
        pivot = pivot[pattern_order]

        # Calculate percentages
        totals = pivot.sum(axis=1)
        pivot_pct = pivot.div(totals, axis=0) * 100

        # Plot heatmap
        im = ax.imshow(pivot_pct.values, cmap='YlOrRd', aspect='auto',
                      vmin=0, vmax=100)

        # Add percentage text
        for i in range(len(databases)):
            for j in range(len(pattern_order)):
                pct = pivot_pct.iloc[i, j]
                count = pivot.iloc[i, j]
                if pct > 0:
                    text = ax.text(j, i, f'{pct:.0f}%\n(n={int(count)})',
                                 ha='center', va='center',
                                 fontsize=8, fontweight='bold',
                                 color='white' if pct > 50 else 'black')

        # Labels
        ax.set_xticks(range(len(pattern_order)))
        ax.set_xticklabels(pattern_order, rotation=45, ha='right', fontsize=9)
        ax.set_yticks(range(len(databases)))
        ax.set_yticklabels(databases, fontsize=10)

        ax.set_title(f'{mutation}', fontsize=12, fontweight='bold',
                    color=MUTATION_COLORS[mutation])

        if idx == 0:
            ax.set_ylabel('Database', fontsize=11)

    # Add colorbar with much more spacing to avoid overlap
    cbar = fig.colorbar(im, ax=axes, orientation='horizontal',
                       pad=1.10, shrink=0.8)
    cbar.set_label('Percentage of Pathways (%)', fontsize=10)

    fig.suptitle('Pattern Distribution Heatmap (% of total per database)',
                fontsize=14, fontweight='bold', y=0.97)

    plt.tight_layout(rect=[0, 0.36, 1, 0.78])

    return fig


# =============================================================================
# MAIN EXECUTION
# =============================================================================

def main():
    """Generate all normalized pattern summary figures"""
    print("="*80)
    print("NORMALIZED PATTERN SUMMARY VISUALIZATION")
    print("="*80)
    print(f"Output directory: {OUTPUT_DIR}")
    print("\nPurpose: Address database size bias in pattern visualization")
    print("  - Large databases (gobp ~5000) mask smaller databases")
    print("  - Normalized views enable proportion comparison")
    print("  - Absolute counts preserved as annotations")

    # Load data
    df = load_classified_data()
    df_counts = prepare_pattern_counts(df)

    print(f"\nPattern counts prepared:")
    print(f"  Databases: {df_counts['Database'].nunique()}")
    print(f"  Mutations: {df_counts['Mutation'].nunique()}")
    print(f"  Patterns: {df_counts['Pattern'].nunique()}")

    # Create figures
    print("\n" + "-"*80)

    # Option 1: Normalized 100% stacked bars (RECOMMENDED)
    fig1 = create_normalized_pattern_summary(df_counts, show_annotations=True)
    if fig1:
        for ext, dpi in [('pdf', 300), ('png', 150)]:
            output_file = OUTPUT_DIR / f'pattern_summary_normalized.{ext}'
            fig1.savefig(output_file, dpi=dpi, bbox_inches='tight', facecolor='white')
            print(f"  Saved: {output_file.name}")
        plt.close(fig1)

    # Option 2: Dual-panel comparison (for each mutation)
    print("\n" + "-"*80)
    for mutation in ['G32A', 'R403C']:
        fig2 = create_dual_panel_comparison(df_counts, mutation=mutation)
        if fig2:
            for ext, dpi in [('pdf', 300), ('png', 150)]:
                output_file = OUTPUT_DIR / f'pattern_comparison_dual_{mutation}.{ext}'
                fig2.savefig(output_file, dpi=dpi, bbox_inches='tight', facecolor='white')
                print(f"  Saved: {output_file.name}")
            plt.close(fig2)

    # Option 3: Percentage heatmap
    print("\n" + "-"*80)
    fig3 = create_percentage_heatmap(df_counts)
    if fig3:
        for ext, dpi in [('pdf', 300), ('png', 150)]:
            output_file = OUTPUT_DIR / f'pattern_summary_heatmap.{ext}'
            fig3.savefig(output_file, dpi=dpi, bbox_inches='tight', facecolor='white')
            print(f"  Saved: {output_file.name}")
        plt.close(fig3)

    # Summary
    print("\n" + "="*80)
    print("COMPLETE!")
    print("="*80)
    print(f"Output directory: {OUTPUT_DIR}")
    print("\nGenerated files:")
    for f in sorted(OUTPUT_DIR.glob('*')):
        print(f"  - {f.name}")

    print("\nComparison with original figures:")
    print("  Original: Absolute counts (database size dominates)")
    print("  New:      Normalized proportions (equal comparison)")
    print("\nRecommendation: Use normalized version for pattern comparison,")
    print("                reference original for absolute context.")
    print("="*80)


if __name__ == '__main__':
    main()
