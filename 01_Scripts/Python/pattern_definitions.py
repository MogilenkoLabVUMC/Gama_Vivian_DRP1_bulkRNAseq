"""
Canonical pattern definitions for DRP1 trajectory analysis.

This module serves as the SINGLE SOURCE OF TRUTH for pattern classification
across both GSEA (NES-based) and GSVA (enrichment score-based) analyses.

Biological Framework
--------------------
- Early: Mutation vs Control at D35 (initial state)
- TrajDev: Mutation-specific maturation trajectory (deviation from control maturation)
- Late: Mutation vs Control at D65 (mature state)

Key Insight
-----------
TrajDev = (D65_mut - D35_mut) - (D65_ctrl - D35_ctrl)

This represents how the mutant's maturation DEVIATES from normal development.
A significant TrajDev indicates active transcriptional plasticity, not just
passive developmental changes.

Pattern Classification
---------------------
1. Compensation: Active adaptive response (TrajDev significant + opposes Early)
2. Progressive: Active worsening (TrajDev significant + amplifies Early)
3. Natural_improvement: Passive recovery (TrajDev not significant + improved)
4. Natural_worsening: Passive deterioration (TrajDev not significant + worsened)
5. Late_onset: Maturation-dependent (no Early defect, Late defect emerges)
6. Transient: Developmental delay (strong Early defect, fully resolved)
7. Complex: Multiphasic or inconsistent patterns

Author: Claude Code
Date: 2025-11-26
"""

import numpy as np
import pandas as pd
from typing import Tuple, Optional, Dict, Any


# =============================================================================
# THRESHOLDS
# =============================================================================

# Significance thresholds
PADJ_SIGNIFICANT = 0.05      # p.adjust threshold for "significant"
PADJ_TRENDING = 0.10         # p.adjust threshold for "trending" (medium confidence)

# Effect size thresholds (GSEA NES)
NES_EFFECT = 0.5             # Minimum |NES| for any effect
NES_STRONG = 1.0             # Minimum |NES| for "strong" defect (used for Transient, Late_onset)

# Improvement/worsening ratios
IMPROVEMENT_RATIO = 0.7      # |Late|/|Early| < 0.7 = "improved" (30% reduction)
WORSENING_RATIO = 1.3        # |Late|/|Early| > 1.3 = "worsened" (30% increase)

# GSVA-specific thresholds (scaled for typical GSVA score range of -1 to 1)
GSVA_EFFECT = 0.15           # Equivalent to NES 0.5
GSVA_STRONG = 0.30           # Equivalent to NES 1.0


# =============================================================================
# SUPER-CATEGORY MAPPING (Simplified for main figures/text)
# =============================================================================
#
# The 7-pattern system provides detailed trajectory classification, but for
# main text interpretation, a simpler 5-category system is often more useful.
#
# Super-categories:
# - Active_Compensation: Patterns with significant TrajDev opposing early defects
# - Active_Progression: Patterns with significant TrajDev amplifying defects (rare)
# - Passive: Recovery/worsening without active TrajDev involvement
# - Late_onset: Maturation-dependent dysfunction (biologically distinct)
# - Other: Transient, Complex patterns
#
# Usage: Main figures use super-categories; detailed 7-pattern taxonomy reserved
# for Methods section and Supplementary materials.

SUPER_CATEGORY_MAP = {
    'Compensation': 'Active_Compensation',
    'Progressive': 'Active_Progression',
    'Natural_improvement': 'Passive',
    'Natural_worsening': 'Passive',
    'Late_onset': 'Late_onset',
    'Transient': 'Other',
    'Complex': 'Other',
    'Insufficient_data': 'Insufficient_data'
}

# Colors for super-categories (colorblind-safe)
SUPER_CATEGORY_COLORS = {
    'Active_Compensation': '#009E73',  # Bluish green - good outcome
    'Active_Progression': '#D55E00',   # Vermillion - worsening
    'Passive': '#56B4E9',              # Sky blue - neutral
    'Late_onset': '#CC79A7',           # Reddish purple - maturation-dependent
    'Other': '#999999',                # Gray
    'Insufficient_data': '#DDDDDD'     # Light gray
}

# Order for plotting
SUPER_CATEGORY_ORDER = [
    'Active_Compensation',
    'Active_Progression',
    'Passive',
    'Late_onset',
    'Other',
    'Insufficient_data'
]


# =============================================================================
# PATTERN DEFINITIONS (for documentation and export)
# =============================================================================

PATTERN_DEFINITIONS: Dict[str, Dict[str, Any]] = {
    'Compensation': {
        'criteria': 'Early sig defect + TrajDev sig opposing + Late improved/resolved',
        'interpretation': (
            'Active adaptive response; mutant significantly deviates from normal '
            'maturation in a direction that opposes the initial defect'
        ),
        'biological_significance': (
            'Reveals neuronal plasticity and active compensation capacity. '
            'These pathways are targets for understanding adaptive mechanisms.'
        ),
        'confidence_levels': ['High (p<0.05)', 'Medium (0.05<=p<0.10)'],
        'color': '#009E73'  # Bluish green - good outcome
    },
    'Progressive': {
        'criteria': 'Early sig defect + TrajDev sig amplifying + Late worsened',
        'interpretation': (
            'Cumulative damage; mutant significantly deviates from normal '
            'maturation in a direction that amplifies the initial defect'
        ),
        'biological_significance': (
            'Indicates failed compensation and progressive pathophysiology. '
            'These pathways represent therapeutic targets to halt progression.'
        ),
        'confidence_levels': ['High (p<0.05)', 'Medium (0.05<=p<0.10)'],
        'color': '#D55E00'  # Vermillion - worsening
    },
    'Natural_worsening': {
        'criteria': 'Early sig defect + TrajDev NOT significant + Late worsened',
        'interpretation': (
            'Passive deterioration; defect worsens during development without '
            'active trajectory change from normal maturation'
        ),
        'biological_significance': (
            'Lacks adaptive plasticity; follows normal maturation trajectory '
            'but cannot buffer the initial defect.'
        ),
        'confidence_levels': ['High'],
        'color': '#E69F00'  # Orange - passive worsening
    },
    'Natural_improvement': {
        'criteria': 'Early sig defect + TrajDev NOT significant + Late improved/resolved',
        'interpretation': (
            'Passive recovery; defect improves through normal developmental '
            'buffering without active compensatory transcriptional programs'
        ),
        'biological_significance': (
            'Developmental robustness handles the defect without requiring '
            'active compensation. Indicates inherent system resilience.'
        ),
        'confidence_levels': ['High'],
        'color': '#56B4E9'  # Sky blue - passive improvement
    },
    'Late_onset': {
        'criteria': 'Early NO defect (p>=0.10 or |NES|<=0.5) + Late sig strong defect (|NES|>1.0, p<0.05)',
        'interpretation': (
            'Maturation-dependent dysfunction; pathway appears normal at D35 '
            'but becomes significantly dysregulated at D65'
        ),
        'biological_significance': (
            'Indicates maturation-dependent vulnerability. These pathways may '
            'be critical for mature neuronal function.'
        ),
        'confidence_levels': ['High'],
        'color': '#CC79A7'  # Reddish purple
    },
    'Transient': {
        'criteria': 'Early sig strong defect (|NES|>1.0, p<0.05) + Late resolved (|NES|<0.5)',
        'interpretation': (
            'Developmental delay that fully recovers by maturation. The strong '
            'early disruption does not persist to the mature state.'
        ),
        'biological_significance': (
            'Temporary disruption during critical developmental window. '
            'May indicate timing-specific vulnerability.'
        ),
        'confidence_levels': ['High'],
        'color': '#0072B2'  # Blue
    },
    'Complex': {
        'criteria': 'Does not fit other patterns',
        'interpretation': 'Non-linear or multiphasic dynamics across trajectory',
        'biological_significance': (
            'Requires individual pathway inspection. May represent pathways '
            'with multiple regulatory mechanisms.'
        ),
        'confidence_levels': ['-'],
        'color': '#F0E442'  # Yellow
    },
    'Insufficient_data': {
        'criteria': 'Missing NES or p.adjust values for required stages',
        'interpretation': 'Cannot classify due to missing data',
        'biological_significance': 'N/A',
        'confidence_levels': ['-'],
        'color': '#DDDDDD'  # Light gray
    }
}


# =============================================================================
# CLASSIFICATION FUNCTION
# =============================================================================

def classify_pattern(
    early_nes: float,
    early_padj: float,
    trajdev_nes: float,
    trajdev_padj: float,
    late_nes: float,
    late_padj: float
) -> Tuple[str, Optional[str]]:
    """
    Classify trajectory pattern based on NES and p.adjust values.

    This function implements the canonical pattern classification algorithm
    that distinguishes between ACTIVE compensation/progression (significant
    TrajDev) and PASSIVE improvement/worsening (non-significant TrajDev).

    Parameters
    ----------
    early_nes : float
        Normalized enrichment score for Early stage (Mutation_vs_Ctrl_D35)
    early_padj : float
        Adjusted p-value for Early stage
    trajdev_nes : float
        Normalized enrichment score for TrajDev (Maturation_Mutation_specific)
    trajdev_padj : float
        Adjusted p-value for TrajDev
    late_nes : float
        Normalized enrichment score for Late stage (Mutation_vs_Ctrl_D65)
    late_padj : float
        Adjusted p-value for Late stage

    Returns
    -------
    tuple
        (pattern_name: str, confidence_level: str or None)
        Confidence is 'High' for p<0.05, 'Medium' for 0.05<=p<0.10, None for Complex

    Examples
    --------
    >>> classify_pattern(-1.5, 0.01, 1.2, 0.02, -0.3, 0.5)
    ('Compensation', 'High')

    >>> classify_pattern(1.2, 0.03, 0.8, 0.2, 1.8, 0.01)
    ('Natural_worsening', 'High')
    """
    # Handle missing values
    if pd.isna([early_nes, trajdev_nes, late_nes, early_padj]).any():
        return ('Insufficient_data', None)

    early_abs = abs(early_nes)
    late_abs = abs(late_nes)
    trajdev_abs = abs(trajdev_nes)

    # -------------------------------------------------------------------------
    # Step 1: Early defect assessment
    # -------------------------------------------------------------------------
    early_sig_defect = (early_padj < PADJ_SIGNIFICANT) and (early_abs > NES_EFFECT)
    early_trending = (early_padj < PADJ_TRENDING) and (early_abs > NES_EFFECT)
    early_strong = (early_padj < PADJ_SIGNIFICANT) and (early_abs > NES_STRONG)
    early_no_defect = (early_padj >= PADJ_TRENDING) or (early_abs <= NES_EFFECT)

    # -------------------------------------------------------------------------
    # Step 2: Late outcome assessment
    # -------------------------------------------------------------------------
    late_sig_defect = (late_padj < PADJ_SIGNIFICANT) and (late_abs > NES_STRONG)
    late_resolved = late_abs < NES_EFFECT

    # Improvement/worsening ratios (protect against division by zero)
    if early_abs > 0.1:
        ratio = late_abs / early_abs
        improved = (ratio < IMPROVEMENT_RATIO) or late_resolved
        worsened = ratio > WORSENING_RATIO
    else:
        improved = late_resolved
        worsened = late_abs > NES_STRONG

    # -------------------------------------------------------------------------
    # Step 3: TrajDev assessment (Active vs Passive)
    # -------------------------------------------------------------------------
    trajdev_sig = (trajdev_padj < PADJ_SIGNIFICANT) and (trajdev_abs > NES_EFFECT)

    # Direction assessment (only meaningful if early has a direction)
    if early_abs > 0.1:
        trajdev_opposes = np.sign(trajdev_nes) != np.sign(early_nes)
        trajdev_amplifies = np.sign(trajdev_nes) == np.sign(early_nes)
    else:
        trajdev_opposes = False
        trajdev_amplifies = False

    # -------------------------------------------------------------------------
    # Step 4: Pattern classification
    # -------------------------------------------------------------------------

    # Late_onset: no early defect, significant late defect
    if early_no_defect and late_sig_defect:
        return ('Late_onset', 'High')

    # Patterns requiring early defect
    if early_sig_defect or early_trending:
        confidence = 'High' if early_sig_defect else 'Medium'

        # Active patterns (TrajDev significant) - check BEFORE Transient
        # because significant opposing TrajDev indicates active compensation,
        # not passive transient recovery
        if trajdev_sig and trajdev_opposes and improved:
            return ('Compensation', confidence)

        if trajdev_sig and trajdev_amplifies and worsened:
            return ('Progressive', confidence)

        # Transient: strong early defect, fully resolved, but NOT actively compensated
        # (Must come AFTER active pattern checks)
        if early_strong and late_resolved:
            return ('Transient', 'High')

        # Passive patterns (TrajDev not significant)
        if not trajdev_sig:
            if improved:
                return ('Natural_improvement', 'High' if early_sig_defect else 'Medium')
            if worsened:
                return ('Natural_worsening', 'High' if early_sig_defect else 'Medium')

    # Transient for cases without early_sig_defect but with strong early
    # (edge case: early_trending but strong)
    if early_strong and late_resolved:
        return ('Transient', 'High')

    return ('Complex', None)


def add_pattern_classification(
    df: pd.DataFrame,
    mutations: Optional[list] = None
) -> pd.DataFrame:
    """
    Add pattern classification columns to dataframe.

    Parameters
    ----------
    df : pd.DataFrame
        Wide-format DataFrame with NES and p.adjust columns for each stage.
        Expected columns: NES_Early_{mutation}, p.adjust_Early_{mutation},
        NES_TrajDev_{mutation}, p.adjust_TrajDev_{mutation},
        NES_Late_{mutation}, p.adjust_Late_{mutation}
    mutations : list, optional
        Mutations to classify. Default: ['G32A', 'R403C']

    Returns
    -------
    pd.DataFrame
        DataFrame with Pattern_{mutation} and Confidence_{mutation} columns
    """
    if mutations is None:
        mutations = ['G32A', 'R403C']

    df = df.copy()

    for mutation in mutations:
        early_nes_col = f'NES_Early_{mutation}'
        early_padj_col = f'p.adjust_Early_{mutation}'
        trajdev_nes_col = f'NES_TrajDev_{mutation}'
        trajdev_padj_col = f'p.adjust_TrajDev_{mutation}'
        late_nes_col = f'NES_Late_{mutation}'
        late_padj_col = f'p.adjust_Late_{mutation}'

        required_cols = [early_nes_col, early_padj_col, trajdev_nes_col,
                        trajdev_padj_col, late_nes_col, late_padj_col]

        # Check if all required columns exist
        if not all(col in df.columns for col in required_cols):
            missing = [c for c in required_cols if c not in df.columns]
            print(f"  Warning: Missing columns for {mutation}: {missing}")
            continue

        # Apply classification to each row
        results = df.apply(
            lambda row: classify_pattern(
                row[early_nes_col], row[early_padj_col],
                row[trajdev_nes_col], row[trajdev_padj_col],
                row[late_nes_col], row[late_padj_col]
            ),
            axis=1
        )

        # Unpack results into separate columns
        df[f'Pattern_{mutation}'] = results.apply(lambda x: x[0])
        df[f'Confidence_{mutation}'] = results.apply(lambda x: x[1])

        # Print summary
        pattern_counts = df[f'Pattern_{mutation}'].value_counts()
        print(f"\n{mutation} pattern distribution:")
        for pattern, count in pattern_counts.items():
            pct = count / len(df) * 100
            print(f"  {pattern}: {count} ({pct:.1f}%)")

        # Print confidence breakdown for active patterns
        for pattern in ['Compensation', 'Progressive']:
            subset = df[df[f'Pattern_{mutation}'] == pattern]
            if len(subset) > 0:
                conf_counts = subset[f'Confidence_{mutation}'].value_counts()
                high = conf_counts.get('High', 0)
                med = conf_counts.get('Medium', 0)
                print(f"    {pattern}: {high} High, {med} Medium confidence")

    return df


def get_pattern_colors() -> Dict[str, str]:
    """
    Get colorblind-safe colors for each pattern.

    Returns
    -------
    dict
        Mapping of pattern name to hex color code
    """
    return {name: info['color'] for name, info in PATTERN_DEFINITIONS.items()}


def get_pattern_summary() -> pd.DataFrame:
    """
    Get summary table of all pattern definitions.

    Returns
    -------
    pd.DataFrame
        Summary table with Pattern, Criteria, Interpretation columns
    """
    rows = []
    for name, info in PATTERN_DEFINITIONS.items():
        if name != 'Insufficient_data':
            rows.append({
                'Pattern': name,
                'Criteria': info['criteria'],
                'Interpretation': info['interpretation'],
                'Biological_Significance': info['biological_significance']
            })
    return pd.DataFrame(rows)


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def is_high_confidence(pattern: str, confidence: str) -> bool:
    """Check if classification is high confidence."""
    return confidence == 'High' and pattern not in ['Complex', 'Insufficient_data']


def get_strict_counts(df: pd.DataFrame, mutation: str) -> pd.Series:
    """
    Get pattern counts using strict (High confidence only) criteria.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with Pattern and Confidence columns
    mutation : str
        Mutation name (G32A or R403C)

    Returns
    -------
    pd.Series
        Counts per pattern (High confidence only)
    """
    mask = df[f'Confidence_{mutation}'] == 'High'
    return df.loc[mask, f'Pattern_{mutation}'].value_counts()


def get_potential_counts(df: pd.DataFrame, mutation: str) -> pd.Series:
    """
    Get pattern counts using potential (High + Medium confidence) criteria.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with Pattern and Confidence columns
    mutation : str
        Mutation name (G32A or R403C)

    Returns
    -------
    pd.Series
        Counts per pattern (High + Medium confidence)
    """
    mask = df[f'Confidence_{mutation}'].isin(['High', 'Medium'])
    return df.loc[mask, f'Pattern_{mutation}'].value_counts()


# =============================================================================
# SUPER-CATEGORY FUNCTIONS
# =============================================================================

def get_super_category(pattern: str) -> str:
    """
    Map a detailed pattern to its super-category.

    Parameters
    ----------
    pattern : str
        One of the 7 pattern names (Compensation, Progressive, etc.)

    Returns
    -------
    str
        Super-category name (Active_Compensation, Active_Progression,
        Passive, Late_onset, Other, or Insufficient_data)
    """
    return SUPER_CATEGORY_MAP.get(pattern, 'Other')


def add_super_category_columns(
    df: pd.DataFrame,
    mutations: Optional[list] = None
) -> pd.DataFrame:
    """
    Add Super_Category columns based on existing Pattern columns.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with Pattern_{mutation} columns
    mutations : list, optional
        Mutations to process. Default: ['G32A', 'R403C']

    Returns
    -------
    pd.DataFrame
        DataFrame with added Super_Category_{mutation} columns
    """
    if mutations is None:
        mutations = ['G32A', 'R403C']

    df = df.copy()

    for mutation in mutations:
        pattern_col = f'Pattern_{mutation}'
        super_col = f'Super_Category_{mutation}'

        if pattern_col in df.columns:
            df[super_col] = df[pattern_col].map(SUPER_CATEGORY_MAP).fillna('Other')

            # Print summary
            super_counts = df[super_col].value_counts()
            print(f"\n{mutation} super-category distribution:")
            for cat in SUPER_CATEGORY_ORDER:
                if cat in super_counts.index:
                    count = super_counts[cat]
                    pct = count / len(df) * 100
                    print(f"  {cat}: {count} ({pct:.1f}%)")

    return df


def get_super_category_colors() -> Dict[str, str]:
    """
    Get colorblind-safe colors for each super-category.

    Returns
    -------
    dict
        Mapping of super-category name to hex color code
    """
    return SUPER_CATEGORY_COLORS.copy()


def get_super_category_order() -> list:
    """
    Get the standard plotting order for super-categories.

    Returns
    -------
    list
        Ordered list of super-category names
    """
    return SUPER_CATEGORY_ORDER.copy()


def get_super_category_summary(df: pd.DataFrame, mutation: str) -> pd.DataFrame:
    """
    Get summary statistics for super-categories.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with Pattern and Super_Category columns
    mutation : str
        Mutation name (G32A or R403C)

    Returns
    -------
    pd.DataFrame
        Summary with counts and percentages per super-category
    """
    super_col = f'Super_Category_{mutation}'

    if super_col not in df.columns:
        df = add_super_category_columns(df, mutations=[mutation])

    counts = df[super_col].value_counts()
    total = len(df)

    summary = pd.DataFrame({
        'Super_Category': SUPER_CATEGORY_ORDER,
        'Count': [counts.get(cat, 0) for cat in SUPER_CATEGORY_ORDER],
        'Percentage': [counts.get(cat, 0) / total * 100 for cat in SUPER_CATEGORY_ORDER]
    })

    return summary[summary['Count'] > 0]
