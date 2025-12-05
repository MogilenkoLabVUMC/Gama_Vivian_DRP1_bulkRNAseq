"""
Data loading utilities with database filtering.

Provides functions to load GSEA results and pathway classifications
with automatic exclusion of irrelevant databases (CGP, canon).
"""

import pandas as pd
import numpy as np
from pathlib import Path

from .config import CONFIG, resolve_path


def load_classified_pathways(exclude_databases=None, min_data_points=None):
    """
    Load classified pathway data with database filtering.

    Automatically excludes CGP and canon databases by default.
    Uses authoritative master_gsea_table.csv with significance-based pattern classifications.

    Parameters
    ----------
    exclude_databases : list, optional
        Databases to exclude. Default: CONFIG['excluded_databases']
    min_data_points : int, optional
        Minimum non-NA values required across trajectory columns.
        Default: None (no filtering)

    Returns
    -------
    pd.DataFrame
        Filtered pathway data with pattern classifications
    """
    # Load from authoritative master table
    data_path = resolve_path(CONFIG['master_gsea_table'])

    if not data_path.exists():
        raise FileNotFoundError(f"Master GSEA table not found: {data_path}")

    print(f"Loading classified pathways from: master_gsea_table.csv")
    df = pd.read_csv(data_path)

    # Pivot to wide format (one row per pathway) to match expected structure
    # Group by pathway_id, database, Description and keep pattern columns
    id_cols = ['pathway_id', 'database', 'Description']
    pattern_cols = ['Pattern_G32A', 'Confidence_G32A', 'Super_Category_G32A',
                    'Pattern_R403C', 'Confidence_R403C', 'Super_Category_R403C']
    nes_cols = ['NES_Early_G32A', 'NES_Early_R403C', 'NES_TrajDev_G32A',
                'NES_TrajDev_R403C', 'NES_Late_G32A', 'NES_Late_R403C']

    # Take first occurrence of each pathway (patterns are consistent across contrasts)
    keep_cols = id_cols + pattern_cols + nes_cols
    available_cols = [c for c in keep_cols if c in df.columns]

    df = df[available_cols].drop_duplicates(subset=['pathway_id', 'database']).copy()

    # Merge p.adjust columns from wide format file (needed for some visualizations)
    wide_file = resolve_path('03_Results/02_Analysis/Python_exports/gsea_results_wide.csv')
    if wide_file.exists():
        df_wide = pd.read_csv(wide_file)

        # Get all p.adjust columns and rename them to trajectory format
        contrast_to_traj = {
            'p.adjust_G32A_vs_Ctrl_D35': 'p.adjust_Early_G32A',
            'p.adjust_Maturation_G32A_specific': 'p.adjust_TrajDev_G32A',
            'p.adjust_G32A_vs_Ctrl_D65': 'p.adjust_Late_G32A',
            'p.adjust_R403C_vs_Ctrl_D35': 'p.adjust_Early_R403C',
            'p.adjust_Maturation_R403C_specific': 'p.adjust_TrajDev_R403C',
            'p.adjust_R403C_vs_Ctrl_D65': 'p.adjust_Late_R403C',
        }

        # Rename columns in wide file
        df_wide = df_wide.rename(columns=contrast_to_traj)

        # Get renamed p.adjust columns that exist
        padj_cols_to_merge = [v for k, v in contrast_to_traj.items() if v in df_wide.columns]
        merge_cols = ['pathway_id', 'database'] + padj_cols_to_merge
        available_merge_cols = [c for c in merge_cols if c in df_wide.columns]

        if len(available_merge_cols) > 2:  # More than just ID cols
            df = df.merge(
                df_wide[available_merge_cols].drop_duplicates(),
                on=['pathway_id', 'database'],
                how='left'
            )
            print(f"  Merged p.adjust columns: {len(padj_cols_to_merge)} columns")

    print(f"  Loaded {len(df)} unique pathways across {df['database'].nunique()} databases")

    # Apply database exclusions
    if exclude_databases is None:
        exclude_databases = CONFIG['excluded_databases']

    if exclude_databases:
        n_before = len(df)
        df = df[~df['database'].isin(exclude_databases)].copy()
        n_excluded = n_before - len(df)
        if n_excluded > 0:
            print(f"  Excluded {n_excluded} pathways from: {exclude_databases}")

    # Apply data point filtering if requested
    if min_data_points is not None:
        df = filter_pathways(df, min_data_points=min_data_points)

    print(f"  Final: {len(df)} pathways from {df['database'].nunique()} databases")

    return df


def filter_pathways(df, min_data_points=3):
    """
    Filter pathways by data availability.

    Keeps pathways that have at least `min_data_points` non-NA values
    across the 6 trajectory columns (Early/TrajDev/Late × G32A/R403C).

    Parameters
    ----------
    df : pd.DataFrame
        Pathway data with NES columns
    min_data_points : int
        Minimum non-NA values required (default: 3)

    Returns
    -------
    pd.DataFrame
        Filtered DataFrame
    """
    nes_cols = CONFIG['nes_cols']

    # Ensure columns exist
    available_cols = [c for c in nes_cols if c in df.columns]

    if not available_cols:
        print("  WARNING: No NES columns found for filtering")
        return df

    df = df.copy()
    df['_data_points'] = df[available_cols].notna().sum(axis=1)

    n_before = len(df)
    df_filtered = df[df['_data_points'] >= min_data_points].copy()
    df_filtered = df_filtered.drop(columns=['_data_points'])

    n_filtered = n_before - len(df_filtered)
    if n_filtered > 0:
        print(f"  Filtered {n_filtered} pathways with <{min_data_points} data points")

    return df_filtered


def get_database_summary(df):
    """
    Get summary statistics by database.

    Parameters
    ----------
    df : pd.DataFrame
        Pathway data

    Returns
    -------
    pd.DataFrame
        Summary with counts and percentages
    """
    summary = df.groupby('database').agg(
        n_pathways=('pathway_id', 'count'),
    ).reset_index()

    summary['pct'] = (summary['n_pathways'] / summary['n_pathways'].sum() * 100).round(1)
    summary = summary.sort_values('n_pathways', ascending=False)

    return summary
