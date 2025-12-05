"""
Trajectory pattern classification for DRP1 analysis.

This module provides pattern classification functions that use the canonical
definitions from pattern_definitions.py. It serves as the interface for
existing code while ensuring consistency with the new significance-based
classification system.

Pattern Types (see pattern_definitions.py for full documentation):
- Compensation: Active adaptive response (TrajDev significant + opposes Early)
- Progressive: Active worsening (TrajDev significant + amplifies Early)
- Natural_worsening: Passive deterioration (TrajDev not significant + worsened)
- Natural_improvement: Passive recovery (TrajDev not significant + improved)
- Late_onset: Maturation-dependent (no Early defect, Late defect emerges)
- Transient: Developmental delay (strong Early defect, fully resolved)
- Complex: Multiphasic or inconsistent patterns

IMPORTANT: The new classification requires BOTH significance (p.adjust < 0.05)
AND magnitude (|NES| > 0.5) for defect detection. This is more stringent than
the previous magnitude-only approach.
"""

import pandas as pd
import numpy as np
import warnings
from typing import Optional, Tuple

from .config import CONFIG

# Import canonical definitions
from .pattern_definitions import (
    # Thresholds
    PADJ_SIGNIFICANT,
    PADJ_TRENDING,
    NES_EFFECT,
    NES_STRONG,
    IMPROVEMENT_RATIO,
    WORSENING_RATIO,
    # Functions
    classify_pattern,
    add_pattern_classification as _add_pattern_classification_new,
    get_pattern_colors,
    get_pattern_summary as _get_pattern_summary_new,
    get_strict_counts,
    get_potential_counts,
    # Definitions
    PATTERN_DEFINITIONS
)


# =============================================================================
# BACKWARD COMPATIBLE FUNCTIONS
# =============================================================================

def classify_trajectory_pattern(
    row,
    mutation: Optional[str] = None,
    early_col: Optional[str] = None,
    trajdev_col: Optional[str] = None,
    late_col: Optional[str] = None,
    use_significance: bool = True
) -> str:
    """
    Classify trajectory pattern based on Early, TrajDev, and Late values.

    This function wraps the canonical classify_pattern() function from
    pattern_definitions.py. It can operate in two modes:

    1. Significance-based (default, recommended): Uses both p.adjust and NES
    2. Magnitude-only (legacy): Uses only NES values (for backward compatibility)

    Parameters
    ----------
    row : pd.Series
        Row with NES and p.adjust columns
    mutation : str, optional
        Mutation name ('G32A' or 'R403C') to auto-construct column names
    early_col, trajdev_col, late_col : str, optional
        Explicit column names for NES values
    use_significance : bool, default True
        If True, use significance-based classification (recommended).
        If False, use legacy magnitude-only classification.

    Returns
    -------
    str
        Pattern classification (without confidence level)
    """
    # Get column names
    if mutation:
        early_nes_col = f'NES_Early_{mutation}'
        trajdev_nes_col = f'NES_TrajDev_{mutation}'
        late_nes_col = f'NES_Late_{mutation}'
        early_padj_col = f'p.adjust_Early_{mutation}'
        trajdev_padj_col = f'p.adjust_TrajDev_{mutation}'
        late_padj_col = f'p.adjust_Late_{mutation}'
    elif early_col and trajdev_col and late_col:
        early_nes_col = early_col
        trajdev_nes_col = trajdev_col
        late_nes_col = late_col
        # Derive p.adjust column names
        early_padj_col = early_col.replace('NES_', 'p.adjust_')
        trajdev_padj_col = trajdev_col.replace('NES_', 'p.adjust_')
        late_padj_col = late_col.replace('NES_', 'p.adjust_')
    else:
        early_nes_col = 'Early_NES'
        trajdev_nes_col = 'TrajDev_NES'
        late_nes_col = 'Late_NES'
        early_padj_col = 'Early_padj'
        trajdev_padj_col = 'TrajDev_padj'
        late_padj_col = 'Late_padj'

    # Get values
    early_nes = row.get(early_nes_col, np.nan)
    trajdev_nes = row.get(trajdev_nes_col, np.nan)
    late_nes = row.get(late_nes_col, np.nan)

    if use_significance:
        # Use canonical significance-based classification
        early_padj = row.get(early_padj_col, 1.0)
        trajdev_padj = row.get(trajdev_padj_col, 1.0)
        late_padj = row.get(late_padj_col, 1.0)

        # Handle missing p.adjust values
        if pd.isna(early_padj):
            early_padj = 1.0
        if pd.isna(trajdev_padj):
            trajdev_padj = 1.0
        if pd.isna(late_padj):
            late_padj = 1.0

        pattern, _ = classify_pattern(
            early_nes, early_padj,
            trajdev_nes, trajdev_padj,
            late_nes, late_padj
        )
        return pattern

    else:
        # Legacy magnitude-only classification (backward compatibility)
        return _classify_trajectory_pattern_legacy(early_nes, trajdev_nes, late_nes)


def _classify_trajectory_pattern_legacy(
    early: float,
    trajdev: float,
    late: float
) -> str:
    """
    Legacy magnitude-only classification (for backward compatibility).

    DEPRECATED: This function does not use significance testing. Use the
    new significance-based classify_pattern() from pattern_definitions.py
    for rigorous analysis.
    """
    # Handle missing values
    if pd.isna([early, trajdev, late]).any():
        return 'Insufficient_data'

    early_abs = abs(early)
    late_abs = abs(late)
    trajdev_abs = abs(trajdev)

    # Thresholds (legacy)
    low_threshold = CONFIG.get('low_threshold', 0.5)
    defect_threshold = 1.0

    # Late-onset: no early defect, but late defect appears
    if early_abs < low_threshold and late_abs > defect_threshold:
        return 'Late_onset'

    # Transient: early defect resolves
    if early_abs > defect_threshold and late_abs < low_threshold:
        return 'Transient'

    # Compensation: defect gets better with active trajectory change
    if early_abs > low_threshold and late_abs < early_abs:
        if np.sign(early) != np.sign(trajdev) and trajdev_abs > low_threshold:
            return 'Compensation'
        else:
            return 'Natural_improvement'

    # Progressive: defect gets worse with active trajectory change
    if early_abs > low_threshold and late_abs > early_abs:
        if np.sign(early) == np.sign(trajdev) and trajdev_abs > low_threshold:
            return 'Progressive'
        else:
            return 'Natural_worsening'

    return 'Complex'


def is_compensation(
    row,
    mutation: str,
    nes_thresh: Optional[float] = None,
    padj_cutoff: Optional[float] = None
) -> bool:
    """
    Check if a pathway shows compensation pattern for given mutation.

    This function uses the NEW significance-based criteria:
    - Early defect: p.adjust < 0.05 AND |NES| > threshold
    - TrajDev significant: p.adjust < 0.05 AND |NES| > 0.5
    - TrajDev opposes Early direction
    - Late outcome improved: |Late| < |Early|

    Parameters
    ----------
    row : pd.Series
        Row with NES and p.adjust columns
    mutation : str
        Mutation name ('G32A' or 'R403C')
    nes_thresh : float, optional
        Minimum |NES| for early defect. Default: NES_EFFECT (0.5)
    padj_cutoff : float, optional
        Significance cutoff. Default: PADJ_SIGNIFICANT (0.05)

    Returns
    -------
    bool
        True if pathway shows compensation pattern with High confidence
    """
    if nes_thresh is None:
        nes_thresh = NES_EFFECT
    if padj_cutoff is None:
        padj_cutoff = PADJ_SIGNIFICANT

    early_nes = row.get(f'NES_Early_{mutation}', np.nan)
    trajdev_nes = row.get(f'NES_TrajDev_{mutation}', np.nan)
    late_nes = row.get(f'NES_Late_{mutation}', np.nan)
    early_padj = row.get(f'p.adjust_Early_{mutation}', 1.0)
    trajdev_padj = row.get(f'p.adjust_TrajDev_{mutation}', 1.0)
    late_padj = row.get(f'p.adjust_Late_{mutation}', 1.0)

    # Handle missing data
    if pd.isna([early_nes, trajdev_nes, late_nes]).any():
        return False

    # Use canonical classification
    pattern, confidence = classify_pattern(
        early_nes, early_padj,
        trajdev_nes, trajdev_padj,
        late_nes, late_padj
    )

    return pattern == 'Compensation' and confidence == 'High'


def add_pattern_classification(
    df: pd.DataFrame,
    mutations: Optional[list] = None,
    use_significance: bool = True
) -> pd.DataFrame:
    """
    Add pattern classification columns to dataframe.

    This function wraps the canonical add_pattern_classification() from
    pattern_definitions.py. It adds both Pattern and Confidence columns.

    Parameters
    ----------
    df : pd.DataFrame
        Wide-format DataFrame with NES and p.adjust columns
    mutations : list, optional
        Mutations to classify. Default: ['G32A', 'R403C']
    use_significance : bool, default True
        If True, use significance-based classification (adds Confidence columns).
        If False, use legacy magnitude-only classification.

    Returns
    -------
    pd.DataFrame
        DataFrame with Pattern_{mutation} and Confidence_{mutation} columns
    """
    if use_significance:
        # Use canonical significance-based classification
        return _add_pattern_classification_new(df, mutations)
    else:
        # Legacy magnitude-only classification
        warnings.warn(
            "Magnitude-only classification is deprecated. "
            "Use use_significance=True for rigorous analysis.",
            DeprecationWarning
        )
        return _add_pattern_classification_legacy(df, mutations)


def _add_pattern_classification_legacy(
    df: pd.DataFrame,
    mutations: Optional[list] = None
) -> pd.DataFrame:
    """Legacy magnitude-only pattern classification."""
    if mutations is None:
        mutations = ['G32A', 'R403C']

    df = df.copy()

    for mutation in mutations:
        early_col = f'NES_Early_{mutation}'
        trajdev_col = f'NES_TrajDev_{mutation}'
        late_col = f'NES_Late_{mutation}'

        if all(col in df.columns for col in [early_col, trajdev_col, late_col]):
            df[f'Pattern_{mutation}'] = df.apply(
                lambda row: classify_trajectory_pattern(
                    row, mutation=mutation, use_significance=False
                ),
                axis=1
            )

            pattern_counts = df[f'Pattern_{mutation}'].value_counts()
            print(f"\n{mutation} pattern distribution (legacy magnitude-only):")
            for pattern, count in pattern_counts.items():
                print(f"  {pattern}: {count}")

    return df


def get_pattern_summary(
    df: pd.DataFrame,
    mutations: Optional[list] = None
) -> pd.DataFrame:
    """
    Get summary of pattern classifications across databases.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with Pattern_* columns
    mutations : list, optional
        Mutations to summarize. Default: ['G32A', 'R403C']

    Returns
    -------
    pd.DataFrame
        Summary DataFrame with counts per database and mutation
    """
    if mutations is None:
        mutations = ['G32A', 'R403C']

    summary_rows = []

    for mutation in mutations:
        pattern_col = f'Pattern_{mutation}'
        confidence_col = f'Confidence_{mutation}'

        if pattern_col not in df.columns:
            continue

        has_confidence = confidence_col in df.columns

        for database in df['database'].unique():
            df_db = df[df['database'] == database]
            counts = df_db[pattern_col].value_counts()

            for pattern, count in counts.items():
                row_data = {
                    'Database': database,
                    'Mutation': mutation,
                    'Pattern': pattern,
                    'Count': count
                }

                # Add confidence breakdown if available
                if has_confidence:
                    subset = df_db[df_db[pattern_col] == pattern]
                    row_data['High_Confidence'] = (subset[confidence_col] == 'High').sum()
                    row_data['Medium_Confidence'] = (subset[confidence_col] == 'Medium').sum()

                summary_rows.append(row_data)

    return pd.DataFrame(summary_rows)


def identify_shared_compensation(
    df: pd.DataFrame,
    nes_thresh: Optional[float] = None,
    padj_cutoff: Optional[float] = None,
    high_confidence_only: bool = True
) -> pd.DataFrame:
    """
    Identify pathways showing compensation in BOTH G32A and R403C.

    Parameters
    ----------
    df : pd.DataFrame
        Wide-format DataFrame with NES and p.adjust columns
    nes_thresh : float, optional
        Minimum |NES| for early defect. Default: NES_EFFECT
    padj_cutoff : float, optional
        Significance cutoff. Default: PADJ_SIGNIFICANT
    high_confidence_only : bool, default True
        If True, only include High confidence compensation.
        If False, include High and Medium confidence.

    Returns
    -------
    pd.DataFrame
        Pathways showing shared compensation
    """
    df = df.copy()

    # Check compensation for each mutation
    df['comp_G32A'] = df.apply(
        lambda r: is_compensation(r, 'G32A', nes_thresh, padj_cutoff),
        axis=1
    )
    df['comp_R403C'] = df.apply(
        lambda r: is_compensation(r, 'R403C', nes_thresh, padj_cutoff),
        axis=1
    )

    # Filter to shared compensation
    df_shared = df[df['comp_G32A'] & df['comp_R403C']].copy()

    print(f"Found {len(df_shared)} pathways with shared compensation (High confidence)")
    print(f"  G32A only: {(df['comp_G32A'] & ~df['comp_R403C']).sum()}")
    print(f"  R403C only: {(~df['comp_G32A'] & df['comp_R403C']).sum()}")

    # Clean up
    df_shared = df_shared.drop(columns=['comp_G32A', 'comp_R403C'])

    return df_shared


# =============================================================================
# EXPORTS
# =============================================================================

# Re-export key items from pattern_definitions for convenience
__all__ = [
    # Functions
    'classify_trajectory_pattern',
    'classify_pattern',
    'add_pattern_classification',
    'get_pattern_summary',
    'is_compensation',
    'identify_shared_compensation',
    'get_pattern_colors',
    'get_strict_counts',
    'get_potential_counts',
    # Thresholds
    'PADJ_SIGNIFICANT',
    'PADJ_TRENDING',
    'NES_EFFECT',
    'NES_STRONG',
    'IMPROVEMENT_RATIO',
    'WORSENING_RATIO',
    # Definitions
    'PATTERN_DEFINITIONS'
]
