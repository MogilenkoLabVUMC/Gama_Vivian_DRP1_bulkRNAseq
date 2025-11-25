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
        Filtered pathway data
    """
    data_path = resolve_path(CONFIG['classified_data'])

    if not data_path.exists():
        raise FileNotFoundError(f"Classified data not found: {data_path}")

    print(f"Loading classified pathways from: {data_path}")
    df = pd.read_csv(data_path)
    print(f"  Loaded {len(df)} pathways across {df['database'].nunique()} databases")

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


def load_gsea_trajectory_data(contrast_file=None, include_nonsignificant=False):
    """
    Load GSEA trajectory data in long format.

    Parameters
    ----------
    contrast_file : str or Path, optional
        Path to GSEA trajectory CSV. Default: CONFIG['data_dir'] / 'gsea_trajectory.csv'
        (or 'gsea_trajectory_all.csv' if include_nonsignificant=True)
    include_nonsignificant : bool, optional
        If True, load ALL tested pathways including those never significant.
        If False (default), load only pathways significant in at least one contrast.
        This allows distinguishing "tested but not significant" from "not tested".

    Returns
    -------
    tuple
        (df_long, contrast_mapping)
        df_long: Long-format DataFrame with pathway GSEA results
        contrast_mapping: Dict mapping original contrasts to trajectory names

    Notes
    -----
    When include_nonsignificant=True, the returned DataFrame includes columns:
    - ever_significant: True if pathway is significant (padj < 0.05) in any contrast
    - ever_significant_trajectory: True if significant in any trajectory contrast
    """
    if contrast_file is None:
        data_dir = resolve_path(CONFIG['data_dir'])
        if include_nonsignificant:
            contrast_file = data_dir / 'gsea_trajectory_all.csv'
        else:
            contrast_file = data_dir / 'gsea_trajectory.csv'
    else:
        contrast_file = Path(contrast_file)

    if not contrast_file.exists():
        raise FileNotFoundError(f"GSEA trajectory data not found: {contrast_file}")

    file_type = "all pathways (incl. non-significant)" if include_nonsignificant else "significant pathways only"
    print(f"Loading GSEA trajectory data ({file_type})")
    print(f"  Source: {contrast_file}")
    df = pd.read_csv(contrast_file)

    # Apply contrast name mapping
    contrast_mapping = CONFIG['contrast_mapping']
    df['trajectory_stage'] = df['contrast'].map(contrast_mapping)

    # Filter out unmapped contrasts
    df = df[df['trajectory_stage'].notna()].copy()

    print(f"  Loaded {len(df)} GSEA results for {df['trajectory_stage'].nunique()} stages")

    # Report significance breakdown if available
    if 'ever_significant_trajectory' in df.columns:
        n_sig = df['ever_significant_trajectory'].sum()
        n_total = len(df)
        print(f"  Significance: {n_sig}/{n_total} records from pathways with ≥1 significant result")

    return df, contrast_mapping


def pivot_to_wide(df_long, id_cols=None):
    """
    Pivot long-format GSEA data to wide format.

    Parameters
    ----------
    df_long : pd.DataFrame
        Long-format DataFrame from load_gsea_trajectory_data
    id_cols : list, optional
        Columns to use as index. Default: ['pathway_id', 'database', 'Description']

    Returns
    -------
    pd.DataFrame
        Wide-format DataFrame with NES_* and p.adjust_* columns
    """
    if id_cols is None:
        id_cols = ['pathway_id', 'database', 'Description']

    # Ensure required columns exist
    available_id_cols = [c for c in id_cols if c in df_long.columns]

    df_pivot = df_long.pivot_table(
        index=available_id_cols,
        columns='trajectory_stage',
        values=['NES', 'p.adjust'],
        aggfunc='first'
    ).reset_index()

    # Flatten column names
    df_pivot.columns = [
        '_'.join(col).strip('_') if col[1] else col[0]
        for col in df_pivot.columns.values
    ]

    return df_pivot


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
