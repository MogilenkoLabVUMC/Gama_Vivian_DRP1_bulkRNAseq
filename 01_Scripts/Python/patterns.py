"""
Trajectory pattern classification for DRP1 analysis.

Classifies pathways into trajectory patterns based on their
Early → TrajDev → Late NES dynamics:

- Compensation: Early defect that improves through active correction
- Progressive: Early defect that worsens over time
- Natural_worsening: Passive worsening without active trajectory change
- Natural_improvement: Passive improvement without active trajectory change
- Late_onset: No early defect, but defect appears later
- Transient: Early defect that resolves
- Persistent: Stable defect throughout
- Complex: Doesn't fit other patterns
- Insufficient_data: Missing trajectory data
"""

import pandas as pd
import numpy as np

from .config import CONFIG


def classify_trajectory_pattern(row, mutation=None,
                                early_col=None, trajdev_col=None, late_col=None):
    """
    Classify trajectory pattern based on Early, TrajDev, and Late NES values.

    Parameters
    ----------
    row : pd.Series
        Row with NES columns
    mutation : str, optional
        Mutation name ('G32A' or 'R403C') to auto-construct column names
    early_col, trajdev_col, late_col : str, optional
        Explicit column names for NES values

    Returns
    -------
    str
        Pattern classification
    """
    # Get NES values
    if mutation:
        early = row.get(f'NES_Early_{mutation}', np.nan)
        trajdev = row.get(f'NES_TrajDev_{mutation}', np.nan)
        late = row.get(f'NES_Late_{mutation}', np.nan)
    elif early_col and trajdev_col and late_col:
        early = row.get(early_col, np.nan)
        trajdev = row.get(trajdev_col, np.nan)
        late = row.get(late_col, np.nan)
    else:
        # Try generic column names
        early = row.get('Early_NES', np.nan)
        trajdev = row.get('TrajDev_NES', np.nan)
        late = row.get('Late_NES', np.nan)

    # Handle missing values
    if pd.isna([early, trajdev, late]).any():
        return 'Insufficient_data'

    early_abs = abs(early)
    late_abs = abs(late)
    trajdev_abs = abs(trajdev)

    # Thresholds
    low_threshold = CONFIG.get('low_threshold', 0.5)
    defect_threshold = 1.0  # |NES| > 1.0 indicates meaningful defect

    # Late-onset: no early defect, but late defect appears
    if early_abs < low_threshold and late_abs > defect_threshold:
        return 'Late_onset'

    # Transient: early defect resolves
    if early_abs > defect_threshold and late_abs < low_threshold:
        return 'Transient'

    # Compensation: defect gets better with active trajectory change
    if early_abs > low_threshold and late_abs < early_abs:
        # Check if trajectory deviation opposes early defect
        if np.sign(early) != np.sign(trajdev) and trajdev_abs > low_threshold:
            return 'Compensation'
        else:
            return 'Natural_improvement'

    # Progressive: defect gets worse with active trajectory change
    if early_abs > low_threshold and late_abs > early_abs:
        # Check if trajectory deviation amplifies early defect
        if np.sign(early) == np.sign(trajdev) and trajdev_abs > low_threshold:
            return 'Progressive'
        else:
            return 'Natural_worsening'

    # Persistent: stable defect with no trajectory change
    if (early_abs > low_threshold and
        abs(late_abs - early_abs) < low_threshold and
        trajdev_abs < low_threshold):
        return 'Persistent'

    return 'Complex'


def is_compensation(row, mutation, nes_thresh=None, padj_cutoff=None):
    """
    Check if a pathway shows compensation pattern for given mutation.

    Compensation criteria:
    - Early defect: |NES| > threshold, padj < cutoff
    - TrajDev opposes Early direction
    - Late outcome improved: |Late| < |Early|

    Parameters
    ----------
    row : pd.Series
        Row with NES and p.adjust columns
    mutation : str
        Mutation name ('G32A' or 'R403C')
    nes_thresh : float, optional
        Minimum |NES| for early defect. Default: CONFIG['nes_threshold']
    padj_cutoff : float, optional
        Significance cutoff. Default: CONFIG['padj_cutoff']

    Returns
    -------
    bool
        True if pathway shows compensation pattern
    """
    if nes_thresh is None:
        nes_thresh = CONFIG['nes_threshold']
    if padj_cutoff is None:
        padj_cutoff = CONFIG['padj_cutoff']

    early = row.get(f'NES_Early_{mutation}', np.nan)
    trajdev = row.get(f'NES_TrajDev_{mutation}', np.nan)
    late = row.get(f'NES_Late_{mutation}', np.nan)
    early_padj = row.get(f'p.adjust_Early_{mutation}', 1.0)

    # Check for missing data
    if pd.isna([early, trajdev, late]).any():
        return False

    # Criteria
    has_early_defect = abs(early) > nes_thresh and early_padj < padj_cutoff
    trajdev_opposes = np.sign(early) != np.sign(trajdev) and abs(trajdev) > 0.5
    late_improved = abs(late) < abs(early)

    return has_early_defect and trajdev_opposes and late_improved


def add_pattern_classification(df, mutations=None):
    """
    Add pattern classification columns to dataframe.

    Parameters
    ----------
    df : pd.DataFrame
        Wide-format DataFrame with NES columns
    mutations : list, optional
        Mutations to classify. Default: ['G32A', 'R403C']

    Returns
    -------
    pd.DataFrame
        DataFrame with Pattern_G32A and Pattern_R403C columns
    """
    if mutations is None:
        mutations = ['G32A', 'R403C']

    df = df.copy()

    for mutation in mutations:
        early_col = f'NES_Early_{mutation}'
        trajdev_col = f'NES_TrajDev_{mutation}'
        late_col = f'NES_Late_{mutation}'

        # Check if all required columns exist
        if all(col in df.columns for col in [early_col, trajdev_col, late_col]):
            df[f'Pattern_{mutation}'] = df.apply(
                lambda row: classify_trajectory_pattern(row, mutation=mutation),
                axis=1
            )

            # Count patterns
            pattern_counts = df[f'Pattern_{mutation}'].value_counts()
            print(f"\n{mutation} pattern distribution:")
            for pattern, count in pattern_counts.items():
                print(f"  {pattern}: {count}")

    return df


def get_pattern_summary(df, mutations=None):
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
        if pattern_col not in df.columns:
            continue

        for database in df['database'].unique():
            df_db = df[df['database'] == database]
            counts = df_db[pattern_col].value_counts()

            for pattern, count in counts.items():
                summary_rows.append({
                    'Database': database,
                    'Mutation': mutation,
                    'Pattern': pattern,
                    'Count': count
                })

    return pd.DataFrame(summary_rows)


def identify_shared_compensation(df, nes_thresh=None, padj_cutoff=None):
    """
    Identify pathways showing compensation in BOTH G32A and R403C.

    Parameters
    ----------
    df : pd.DataFrame
        Wide-format DataFrame with NES and p.adjust columns
    nes_thresh, padj_cutoff : float, optional
        Thresholds for compensation criteria

    Returns
    -------
    pd.DataFrame
        Pathways showing shared compensation
    """
    df = df.copy()

    # Check compensation for each mutation
    df['comp_G32A'] = df.apply(lambda r: is_compensation(r, 'G32A', nes_thresh, padj_cutoff), axis=1)
    df['comp_R403C'] = df.apply(lambda r: is_compensation(r, 'R403C', nes_thresh, padj_cutoff), axis=1)

    # Filter to shared compensation
    df_shared = df[df['comp_G32A'] & df['comp_R403C']].copy()

    print(f"Found {len(df_shared)} pathways with shared compensation")
    print(f"  G32A only: {(df['comp_G32A'] & ~df['comp_R403C']).sum()}")
    print(f"  R403C only: {(~df['comp_G32A'] & df['comp_R403C']).sum()}")

    # Clean up
    df_shared = df_shared.drop(columns=['comp_G32A', 'comp_R403C'])

    return df_shared
