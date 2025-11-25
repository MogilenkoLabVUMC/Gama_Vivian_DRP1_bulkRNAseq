#!/usr/bin/env python3
"""
Master Pathway Table Generator
===============================

Creates a comprehensive master table of all pathways tested in GSEA analysis with:
1. Pathway name and database
2. Contrast information
3. All GSEA statistics (NES, pvalue, p.adjust, qvalue, enrichmentScore, setSize)
4. Pattern classification for G32A and R403C
5. Change consistency between mutants

This table provides complete information for all pathways across all contrasts,
enabling comprehensive downstream analysis and interpretation.

Output: 03_Results/02_Analysis/master_pathway_table.csv
"""

import sys
from pathlib import Path

# Add module paths
sys.path.insert(0, str(Path(__file__).parent.parent / '01_Scripts'))

import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings('ignore')

# Import project modules
from Python.config import CONFIG, resolve_path

# =============================================================================
# SETUP
# =============================================================================

OUTPUT_DIR = resolve_path('03_Results/02_Analysis')
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# =============================================================================
# DATA LOADING
# =============================================================================

def load_gsea_results():
    """Load complete GSEA results from long format CSV"""
    print("Loading GSEA results...")

    data_dir = resolve_path(CONFIG['data_dir'])
    df_long = pd.read_csv(data_dir / 'gsea_results_long.csv')

    print(f"  Loaded {len(df_long)} GSEA results")
    print(f"  Contrasts: {df_long['contrast'].nunique()}")
    print(f"  Databases: {df_long['database'].nunique()}")
    print(f"  Unique pathways: {df_long['pathway_id'].nunique()}")

    return df_long


def load_pattern_classifications():
    """Load pattern classifications from pathways_classified.csv"""
    print("\nLoading pattern classifications...")

    pattern_file = resolve_path('03_Results/02_Analysis/Plots/Cross_database_validation/pathways_classified.csv')
    df_patterns = pd.read_csv(pattern_file)

    print(f"  Loaded {len(df_patterns)} classified pathways")

    # Select relevant columns for merging
    pattern_cols = ['pathway_id', 'database', 'Description',
                   'Pattern_G32A', 'Pattern_R403C',
                   'NES_Early_G32A', 'NES_Early_R403C',
                   'NES_Late_G32A', 'NES_Late_R403C',
                   'NES_TrajDev_G32A', 'NES_TrajDev_R403C']

    # Only keep columns that exist
    existing_cols = [col for col in pattern_cols if col in df_patterns.columns]
    df_patterns = df_patterns[existing_cols]

    return df_patterns


# =============================================================================
# CONSISTENCY CLASSIFICATION
# =============================================================================

def classify_change_consistency(row):
    """
    Classify consistency of change between G32A and R403C mutants.

    Categories:
    - Co-upregulated: Both mutants show positive NES
    - Co-downregulated: Both mutants show negative NES
    - Opposite: One up, one down
    - G32A-specific: Only G32A significant
    - R403C-specific: Only R403C significant
    - Both-nonsignificant: Neither mutant significant
    - Mixed: Complex pattern across stages

    Parameters
    ----------
    row : pd.Series
        Row with Pattern_G32A and Pattern_R403C columns

    Returns
    -------
    str
        Consistency classification
    """
    pattern_g32a = row.get('Pattern_G32A', 'Unknown')
    pattern_r403c = row.get('Pattern_R403C', 'Unknown')

    # Handle missing patterns
    if pd.isna(pattern_g32a) or pd.isna(pattern_r403c):
        return 'Insufficient_data'

    if pattern_g32a == 'Insufficient_data' or pattern_r403c == 'Insufficient_data':
        return 'Insufficient_data'

    # Get NES values for consistency check (Early stage as primary indicator)
    nes_early_g32a = row.get('NES_Early_G32A', np.nan)
    nes_early_r403c = row.get('NES_Early_R403C', np.nan)

    # Both show same pattern type
    if pattern_g32a == pattern_r403c:
        # Check directionality
        if not pd.isna(nes_early_g32a) and not pd.isna(nes_early_r403c):
            if nes_early_g32a > 0 and nes_early_r403c > 0:
                return f'Consistent_{pattern_g32a}_co-upregulated'
            elif nes_early_g32a < 0 and nes_early_r403c < 0:
                return f'Consistent_{pattern_g32a}_co-downregulated'
            else:
                return f'Consistent_{pattern_g32a}_opposite-direction'
        else:
            return f'Consistent_{pattern_g32a}'

    # Different patterns
    else:
        # Check if one is more severe
        severity_order = ['Insufficient_data', 'Complex', 'Persistent',
                         'Natural_improvement', 'Natural_worsening',
                         'Transient', 'Late_onset', 'Progressive', 'Compensation']

        try:
            g32a_severity = severity_order.index(pattern_g32a)
            r403c_severity = severity_order.index(pattern_r403c)

            if g32a_severity > r403c_severity:
                return f'Inconsistent_G32A-more-severe'
            elif r403c_severity > g32a_severity:
                return f'Inconsistent_R403C-more-severe'
            else:
                return f'Inconsistent_{pattern_g32a}_vs_{pattern_r403c}'
        except ValueError:
            return f'Inconsistent_{pattern_g32a}_vs_{pattern_r403c}'


def add_consistency_classification(df):
    """Add change_consistency column to dataframe"""
    print("\nClassifying change consistency between mutants...")

    df['change_consistency'] = df.apply(classify_change_consistency, axis=1)

    # Print summary
    consistency_counts = df['change_consistency'].value_counts()
    print("\nConsistency classification summary:")
    for consistency, count in consistency_counts.head(15).items():
        print(f"  {consistency}: {count}")

    return df


# =============================================================================
# MASTER TABLE GENERATION
# =============================================================================

def create_master_table(df_long, df_patterns):
    """
    Create master table by merging GSEA results with pattern classifications.

    The resulting table has one row per pathway per contrast, with pattern
    classifications and consistency annotations.
    """
    print("\nCreating master table...")

    # Merge GSEA results with pattern classifications
    df_master = df_long.merge(
        df_patterns,
        on=['pathway_id', 'database', 'Description'],
        how='left'
    )

    print(f"  Master table: {len(df_master)} rows")

    # Add consistency classification
    # For each unique pathway, calculate consistency once
    pathway_consistency = df_patterns.copy()
    pathway_consistency = add_consistency_classification(pathway_consistency)

    # Merge consistency back to master table
    df_master = df_master.merge(
        pathway_consistency[['pathway_id', 'change_consistency']],
        on='pathway_id',
        how='left'
    )

    # Reorder columns for clarity
    column_order = [
        # Pathway identification
        'pathway_id',
        'database',
        'ID',
        'Description',

        # Contrast information
        'contrast',
        'new_name',
        'category',
        'mutation',

        # GSEA statistics
        'NES',
        'pvalue',
        'p.adjust',
        'qvalue',
        'enrichmentScore',
        'setSize',

        # Pattern classifications
        'Pattern_G32A',
        'Pattern_R403C',
        'change_consistency',

        # Trajectory NES values (for reference)
        'NES_Early_G32A',
        'NES_Early_R403C',
        'NES_TrajDev_G32A',
        'NES_TrajDev_R403C',
        'NES_Late_G32A',
        'NES_Late_R403C',

        # Significance flags
        'ever_significant',
        'ever_significant_trajectory'
    ]

    # Only keep columns that exist
    existing_columns = [col for col in column_order if col in df_master.columns]
    df_master = df_master[existing_columns]

    # Sort by database, pathway, and contrast
    df_master = df_master.sort_values(['database', 'Description', 'contrast'])

    return df_master


# =============================================================================
# SUMMARY STATISTICS
# =============================================================================

def print_summary_statistics(df_master):
    """Print summary statistics about the master table"""
    print("\n" + "="*80)
    print("MASTER TABLE SUMMARY")
    print("="*80)

    print(f"\nTotal rows: {len(df_master):,}")
    print(f"Unique pathways: {df_master['pathway_id'].nunique():,}")
    print(f"Databases: {df_master['database'].nunique()}")
    print(f"Contrasts: {df_master['contrast'].nunique()}")

    print("\n--- Pathways per database ---")
    db_counts = df_master.groupby('database')['pathway_id'].nunique().sort_values(ascending=False)
    for db, count in db_counts.items():
        print(f"  {db}: {count:,} pathways")

    print("\n--- Results per contrast ---")
    contrast_counts = df_master.groupby('contrast').size().sort_values(ascending=False)
    for contrast, count in contrast_counts.items():
        print(f"  {contrast}: {count:,} results")

    print("\n--- Pattern distribution ---")
    if 'Pattern_G32A' in df_master.columns:
        g32a_patterns = df_master.drop_duplicates('pathway_id')['Pattern_G32A'].value_counts()
        print("\nG32A patterns (unique pathways):")
        for pattern, count in g32a_patterns.head(10).items():
            print(f"  {pattern}: {count}")

    if 'Pattern_R403C' in df_master.columns:
        r403c_patterns = df_master.drop_duplicates('pathway_id')['Pattern_R403C'].value_counts()
        print("\nR403C patterns (unique pathways):")
        for pattern, count in r403c_patterns.head(10).items():
            print(f"  {pattern}: {count}")

    print("\n--- Change consistency ---")
    if 'change_consistency' in df_master.columns:
        consistency_counts = df_master.drop_duplicates('pathway_id')['change_consistency'].value_counts()
        print("\nTop consistency categories:")
        for consistency, count in consistency_counts.head(15).items():
            print(f"  {consistency}: {count}")

    print("\n--- Significance statistics ---")
    if 'ever_significant' in df_master.columns:
        n_sig = df_master['ever_significant'].sum()
        n_total = len(df_master)
        print(f"  Results with significance in at least one contrast: {n_sig:,} ({100*n_sig/n_total:.1f}%)")

    if 'ever_significant_trajectory' in df_master.columns:
        n_sig_traj = df_master['ever_significant_trajectory'].sum()
        print(f"  Results with trajectory significance: {n_sig_traj:,} ({100*n_sig_traj/n_total:.1f}%)")


# =============================================================================
# MAIN EXECUTION
# =============================================================================

def main():
    """Main execution function"""
    print("="*80)
    print("MASTER PATHWAY TABLE GENERATOR")
    print("="*80)
    print(f"\nOutput directory: {OUTPUT_DIR}")

    # Load data
    df_long = load_gsea_results()
    df_patterns = load_pattern_classifications()

    # Create master table
    df_master = create_master_table(df_long, df_patterns)

    # Print summary statistics
    print_summary_statistics(df_master)

    # Save master table
    output_file = OUTPUT_DIR / 'master_pathway_table.csv'
    df_master.to_csv(output_file, index=False)

    print("\n" + "="*80)
    print("COMPLETE!")
    print("="*80)
    print(f"\nSaved master table: {output_file}")
    print(f"  Rows: {len(df_master):,}")
    print(f"  Columns: {len(df_master.columns)}")
    print(f"\nColumns in master table:")
    for i, col in enumerate(df_master.columns, 1):
        print(f"  {i:2d}. {col}")

    print("\n" + "="*80)
    print("Use this table for comprehensive pathway analysis!")
    print("="*80)


if __name__ == '__main__':
    main()
