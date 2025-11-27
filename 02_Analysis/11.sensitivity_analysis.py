#!/usr/bin/env python3
"""
Sensitivity Analysis for Pattern Classification Thresholds

Tests robustness of pattern classification across reasonable threshold variations.
Generates supplementary figure/table showing qualitative conclusions remain stable.

Thresholds tested:
- NES_EFFECT: 0.4, 0.5 (default), 0.6
- NES_STRONG: 0.8, 1.0 (default), 1.2
- IMPROVEMENT_RATIO: 0.6, 0.7 (default), 0.8
- WORSENING_RATIO: 1.25, 1.3 (default), 1.4

Author: Claude Code
Date: 2025-11-26
"""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import Tuple, Optional, Dict, List
from itertools import product
import matplotlib.pyplot as plt
import seaborn as sns

# =============================================================================
# PATHS
# =============================================================================

PROJECT_ROOT = Path(__file__).parent.parent
DATA_DIR = PROJECT_ROOT / "03_Results" / "02_Analysis"
OUTPUT_DIR = DATA_DIR / "Sensitivity_Analysis"
OUTPUT_DIR.mkdir(exist_ok=True)

# =============================================================================
# PARAMETERIZED CLASSIFICATION FUNCTION
# =============================================================================

def classify_pattern_parameterized(
    early_nes: float,
    early_padj: float,
    trajdev_nes: float,
    trajdev_padj: float,
    late_nes: float,
    late_padj: float,
    # Thresholds as parameters
    padj_significant: float = 0.05,
    padj_trending: float = 0.10,
    nes_effect: float = 0.5,
    nes_strong: float = 1.0,
    improvement_ratio: float = 0.7,
    worsening_ratio: float = 1.3
) -> Tuple[str, Optional[str]]:
    """
    Classify trajectory pattern with parameterized thresholds.

    This allows testing different threshold combinations for sensitivity analysis.
    """
    # Handle missing values
    if pd.isna([early_nes, trajdev_nes, late_nes, early_padj]).any():
        return ('Insufficient_data', None)

    early_abs = abs(early_nes)
    late_abs = abs(late_nes)
    trajdev_abs = abs(trajdev_nes)

    # Step 1: Early defect assessment
    early_sig_defect = (early_padj < padj_significant) and (early_abs > nes_effect)
    early_trending = (early_padj < padj_trending) and (early_abs > nes_effect)
    early_strong = (early_padj < padj_significant) and (early_abs > nes_strong)
    early_no_defect = (early_padj >= padj_trending) or (early_abs <= nes_effect)

    # Step 2: Late outcome assessment
    late_sig_defect = (late_padj < padj_significant) and (late_abs > nes_strong)
    late_resolved = late_abs < nes_effect

    # Improvement/worsening ratios
    if early_abs > 0.1:
        ratio = late_abs / early_abs
        improved = (ratio < improvement_ratio) or late_resolved
        worsened = ratio > worsening_ratio
    else:
        improved = late_resolved
        worsened = late_abs > nes_strong

    # Step 3: TrajDev assessment
    trajdev_sig = (trajdev_padj < padj_significant) and (trajdev_abs > nes_effect)

    # Direction assessment
    if early_abs > 0.1:
        trajdev_opposes = np.sign(trajdev_nes) != np.sign(early_nes)
        trajdev_amplifies = np.sign(trajdev_nes) == np.sign(early_nes)
    else:
        trajdev_opposes = False
        trajdev_amplifies = False

    # Step 4: Pattern classification
    # Late_onset: no early defect, significant late defect
    if early_no_defect and late_sig_defect:
        return ('Late_onset', 'High')

    # Patterns requiring early defect
    if early_sig_defect or early_trending:
        confidence = 'High' if early_sig_defect else 'Medium'

        # Active patterns - check BEFORE Transient
        if trajdev_sig and trajdev_opposes and improved:
            return ('Compensation', confidence)

        if trajdev_sig and trajdev_amplifies and worsened:
            return ('Progressive', confidence)

        # Transient: strong early, fully resolved
        if early_strong and late_resolved:
            return ('Transient', 'High')

        # Passive patterns
        if not trajdev_sig:
            if improved:
                return ('Natural_improvement', 'High' if early_sig_defect else 'Medium')
            if worsened:
                return ('Natural_worsening', 'High' if early_sig_defect else 'Medium')

    # Edge case transient
    if early_strong and late_resolved:
        return ('Transient', 'High')

    return ('Complex', None)


# =============================================================================
# SENSITIVITY ANALYSIS
# =============================================================================

def run_sensitivity_analysis(df: pd.DataFrame, mutations: List[str] = ['G32A', 'R403C']) -> pd.DataFrame:
    """
    Run classification with all threshold combinations and compare results.

    Parameters
    ----------
    df : pd.DataFrame
        Wide-format GSEA results with NES and p.adjust columns
    mutations : list
        Mutations to analyze

    Returns
    -------
    pd.DataFrame
        Summary of pattern distributions across threshold combinations
    """
    # Threshold combinations to test
    nes_effect_values = [0.4, 0.5, 0.6]
    nes_strong_values = [0.8, 1.0, 1.2]
    improvement_ratio_values = [0.6, 0.7, 0.8]
    worsening_ratio_values = [1.25, 1.3, 1.4]

    results = []

    # Column mappings
    col_map = {
        'G32A': {
            'early_nes': 'NES_G32A_vs_Ctrl_D35',
            'early_padj': 'p.adjust_G32A_vs_Ctrl_D35',
            'trajdev_nes': 'NES_Maturation_G32A_specific',
            'trajdev_padj': 'p.adjust_Maturation_G32A_specific',
            'late_nes': 'NES_G32A_vs_Ctrl_D65',
            'late_padj': 'p.adjust_G32A_vs_Ctrl_D65'
        },
        'R403C': {
            'early_nes': 'NES_R403C_vs_Ctrl_D35',
            'early_padj': 'p.adjust_R403C_vs_Ctrl_D35',
            'trajdev_nes': 'NES_Maturation_R403C_specific',
            'trajdev_padj': 'p.adjust_Maturation_R403C_specific',
            'late_nes': 'NES_R403C_vs_Ctrl_D65',
            'late_padj': 'p.adjust_R403C_vs_Ctrl_D65'
        }
    }

    # Generate all combinations
    combinations = list(product(
        nes_effect_values,
        nes_strong_values,
        improvement_ratio_values,
        worsening_ratio_values
    ))

    print(f"Testing {len(combinations)} threshold combinations...")

    for nes_effect, nes_strong, imp_ratio, wors_ratio in combinations:
        # Skip invalid combinations (nes_effect should be < nes_strong)
        if nes_effect >= nes_strong:
            continue

        combo_id = f"NES_eff={nes_effect}_strong={nes_strong}_imp={imp_ratio}_wors={wors_ratio}"

        for mutation in mutations:
            cols = col_map[mutation]

            # Classify each pathway
            patterns = []
            for _, row in df.iterrows():
                pattern, conf = classify_pattern_parameterized(
                    early_nes=row[cols['early_nes']],
                    early_padj=row[cols['early_padj']],
                    trajdev_nes=row[cols['trajdev_nes']],
                    trajdev_padj=row[cols['trajdev_padj']],
                    late_nes=row[cols['late_nes']],
                    late_padj=row[cols['late_padj']],
                    nes_effect=nes_effect,
                    nes_strong=nes_strong,
                    improvement_ratio=imp_ratio,
                    worsening_ratio=wors_ratio
                )
                patterns.append(pattern)

            # Count patterns
            pattern_counts = pd.Series(patterns).value_counts()
            total = len(patterns)

            # Calculate super-category percentages
            active_comp = pattern_counts.get('Compensation', 0)
            active_prog = pattern_counts.get('Progressive', 0)
            passive = (pattern_counts.get('Natural_improvement', 0) +
                      pattern_counts.get('Natural_worsening', 0))
            late_onset = pattern_counts.get('Late_onset', 0)
            transient = pattern_counts.get('Transient', 0)
            complex_pat = pattern_counts.get('Complex', 0)
            insufficient = pattern_counts.get('Insufficient_data', 0)

            results.append({
                'combination_id': combo_id,
                'NES_EFFECT': nes_effect,
                'NES_STRONG': nes_strong,
                'IMPROVEMENT_RATIO': imp_ratio,
                'WORSENING_RATIO': wors_ratio,
                'mutation': mutation,
                'n_pathways': total,
                'Compensation': active_comp,
                'Compensation_pct': active_comp / total * 100,
                'Progressive': active_prog,
                'Progressive_pct': active_prog / total * 100,
                'Natural_improvement': pattern_counts.get('Natural_improvement', 0),
                'Natural_worsening': pattern_counts.get('Natural_worsening', 0),
                'Passive_total': passive,
                'Passive_pct': passive / total * 100,
                'Late_onset': late_onset,
                'Late_onset_pct': late_onset / total * 100,
                'Transient': transient,
                'Transient_pct': transient / total * 100,
                'Complex': complex_pat,
                'Complex_pct': complex_pat / total * 100,
                'Insufficient_data': insufficient,
                'is_default': (nes_effect == 0.5 and nes_strong == 1.0 and
                              imp_ratio == 0.7 and wors_ratio == 1.3)
            })

    return pd.DataFrame(results)


def analyze_key_claims(sensitivity_df: pd.DataFrame) -> pd.DataFrame:
    """
    Check whether key biological claims remain stable across thresholds.

    Key claims to verify:
    1. Compensation is dominant classifiable pattern (among non-Complex)
    2. R403C shows more compensation than G32A
    3. Progressive patterns are rare (<5%)
    4. Passive patterns exist but are secondary to Active compensation
    """
    claims = []

    for combo_id in sensitivity_df['combination_id'].unique():
        combo_data = sensitivity_df[sensitivity_df['combination_id'] == combo_id]
        g32a = combo_data[combo_data['mutation'] == 'G32A'].iloc[0]
        r403c = combo_data[combo_data['mutation'] == 'R403C'].iloc[0]

        # Calculate non-Complex totals for more meaningful "dominance"
        g32a_classifiable = g32a['n_pathways'] - g32a['Complex'] - g32a.get('Insufficient_data', 0)
        r403c_classifiable = r403c['n_pathways'] - r403c['Complex'] - r403c.get('Insufficient_data', 0)

        # Compensation fraction of classifiable pathways
        g32a_comp_of_classifiable = g32a['Compensation'] / g32a_classifiable * 100 if g32a_classifiable > 0 else 0
        r403c_comp_of_classifiable = r403c['Compensation'] / r403c_classifiable * 100 if r403c_classifiable > 0 else 0

        claims.append({
            'combination_id': combo_id,
            'NES_EFFECT': g32a['NES_EFFECT'],
            'NES_STRONG': g32a['NES_STRONG'],
            'IMPROVEMENT_RATIO': g32a['IMPROVEMENT_RATIO'],
            'WORSENING_RATIO': g32a['WORSENING_RATIO'],
            'is_default': g32a['is_default'],
            # Actual counts and percentages
            'comp_count_G32A': g32a['Compensation'],
            'comp_count_R403C': r403c['Compensation'],
            'comp_pct_G32A': g32a['Compensation_pct'],
            'comp_pct_R403C': r403c['Compensation_pct'],
            'classifiable_G32A': g32a_classifiable,
            'classifiable_R403C': r403c_classifiable,
            'comp_of_classifiable_G32A': g32a_comp_of_classifiable,
            'comp_of_classifiable_R403C': r403c_comp_of_classifiable,
            # Claim 1: Compensation is largest classifiable pattern (> 50% of non-Complex)
            'comp_dominates_classifiable_G32A': g32a_comp_of_classifiable > 50,
            'comp_dominates_classifiable_R403C': r403c_comp_of_classifiable > 50,
            # Claim 2: Which mutation shows more compensation
            'R403C_more_compensation': r403c['Compensation'] > g32a['Compensation'],
            'compensation_diff': r403c['Compensation'] - g32a['Compensation'],
            # Claim 3: Progressive is rare (<5% of total)
            'progressive_rare_G32A': g32a['Progressive_pct'] < 5,
            'progressive_rare_R403C': r403c['Progressive_pct'] < 5,
            'prog_pct_G32A': g32a['Progressive_pct'],
            'prog_pct_R403C': r403c['Progressive_pct'],
            # Claim 4: Compensation > Passive
            'comp_exceeds_passive_G32A': g32a['Compensation'] > g32a['Passive_total'],
            'comp_exceeds_passive_R403C': r403c['Compensation'] > r403c['Passive_total'],
        })

    return pd.DataFrame(claims)


def create_sensitivity_heatmap(sensitivity_df: pd.DataFrame, mutation: str, output_dir: Path):
    """
    Create heatmap showing pattern percentages across threshold combinations.
    """
    # Filter for this mutation
    df = sensitivity_df[sensitivity_df['mutation'] == mutation].copy()

    # Focus on NES_EFFECT vs IMPROVEMENT_RATIO (fixing others at default)
    df_subset = df[
        (df['NES_STRONG'] == 1.0) &
        (df['WORSENING_RATIO'] == 1.3)
    ].copy()

    if len(df_subset) == 0:
        print(f"No data for heatmap subset for {mutation}")
        return

    # Create pivot table for Compensation percentage
    pivot = df_subset.pivot_table(
        values='Compensation_pct',
        index='IMPROVEMENT_RATIO',
        columns='NES_EFFECT',
        aggfunc='mean'
    )

    fig, ax = plt.subplots(figsize=(8, 6))
    sns.heatmap(
        pivot,
        annot=True,
        fmt='.1f',
        cmap='YlGnBu',
        ax=ax,
        cbar_kws={'label': 'Compensation %'}
    )
    ax.set_title(f'{mutation}: Compensation % by Threshold\n(NES_STRONG=1.0, WORSENING_RATIO=1.3)')
    ax.set_xlabel('NES_EFFECT threshold')
    ax.set_ylabel('IMPROVEMENT_RATIO threshold')

    # Mark default
    default_row = list(pivot.index).index(0.7)
    default_col = list(pivot.columns).index(0.5)
    ax.add_patch(plt.Rectangle(
        (default_col, default_row), 1, 1,
        fill=False, edgecolor='red', linewidth=3
    ))

    plt.tight_layout()
    plt.savefig(output_dir / f'sensitivity_heatmap_{mutation}.pdf', dpi=150)
    plt.savefig(output_dir / f'sensitivity_heatmap_{mutation}.png', dpi=150)
    plt.close()
    print(f"  Saved heatmap for {mutation}")


def create_summary_figure(sensitivity_df: pd.DataFrame, claims_df: pd.DataFrame, output_dir: Path):
    """
    Create comprehensive summary figure showing threshold robustness.
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    # Get default combination for reference
    default_data = sensitivity_df[sensitivity_df['is_default']]

    # 1. Compensation % across all combinations (boxplot)
    ax = axes[0, 0]
    g32a_comp = sensitivity_df[sensitivity_df['mutation'] == 'G32A']['Compensation_pct']
    r403c_comp = sensitivity_df[sensitivity_df['mutation'] == 'R403C']['Compensation_pct']

    bp = ax.boxplot([g32a_comp, r403c_comp], labels=['G32A', 'R403C'], patch_artist=True)
    bp['boxes'][0].set_facecolor('#1f77b4')
    bp['boxes'][1].set_facecolor('#ff7f0e')

    # Add default values as stars
    for i, (mut, col) in enumerate([(g32a_comp, '#1f77b4'), (r403c_comp, '#ff7f0e')]):
        default_val = default_data[default_data['mutation'] == ['G32A', 'R403C'][i]]['Compensation_pct'].values[0]
        ax.scatter([i+1], [default_val], marker='*', s=200, c='red', zorder=5, label='Default' if i==0 else '')

    ax.set_ylabel('Compensation %')
    ax.set_title('A. Compensation % Across All Threshold Combinations\n(Red star = default thresholds)')
    ax.legend()

    # 2. Pattern distribution comparison (default thresholds)
    ax = axes[0, 1]
    patterns = ['Compensation', 'Progressive', 'Natural_improvement',
                'Natural_worsening', 'Late_onset', 'Transient', 'Complex']

    x = np.arange(len(patterns))
    width = 0.35

    g32a_default = default_data[default_data['mutation'] == 'G32A'].iloc[0]
    r403c_default = default_data[default_data['mutation'] == 'R403C'].iloc[0]

    g32a_vals = [g32a_default[p] for p in patterns]
    r403c_vals = [r403c_default[p] for p in patterns]

    ax.bar(x - width/2, g32a_vals, width, label='G32A', color='#1f77b4')
    ax.bar(x + width/2, r403c_vals, width, label='R403C', color='#ff7f0e')
    ax.set_xticks(x)
    ax.set_xticklabels(patterns, rotation=45, ha='right')
    ax.set_ylabel('Number of Pathways')
    ax.set_title('B. Pattern Distribution (Default Thresholds)')
    ax.legend()

    # 3. Claim stability: Multiple claims
    ax = axes[1, 0]
    total_combos = len(claims_df)

    # Key stable claims
    claims_to_check = {
        'R403C > G32A\nCompensation': claims_df['R403C_more_compensation'].sum(),
        'Progressive\nRare (G32A)': claims_df['progressive_rare_G32A'].sum(),
        'Progressive\nRare (R403C)': claims_df['progressive_rare_R403C'].sum(),
        'Comp > Passive\n(G32A)': claims_df['comp_exceeds_passive_G32A'].sum(),
        'Comp > Passive\n(R403C)': claims_df['comp_exceeds_passive_R403C'].sum(),
    }

    claim_names = list(claims_to_check.keys())
    claim_pcts = [v / total_combos * 100 for v in claims_to_check.values()]

    bars = ax.bar(range(len(claim_names)), claim_pcts, color=['#2ecc71', '#3498db', '#e74c3c', '#9b59b6', '#f39c12'])
    ax.axhline(y=100, color='green', linestyle='--', alpha=0.5)
    ax.set_xticks(range(len(claim_names)))
    ax.set_xticklabels(claim_names, fontsize=8)
    ax.set_ylabel('% of Combinations Where Claim Holds')
    ax.set_title(f'C. Claim Stability Across {total_combos} Threshold Combinations')
    ax.set_ylim(0, 115)

    for bar, pct in zip(bars, claim_pcts):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 2,
                f'{pct:.0f}%', ha='center', fontsize=10, fontweight='bold')

    # 4. Range of compensation percentages
    ax = axes[1, 1]
    for i, mutation in enumerate(['G32A', 'R403C']):
        mut_data = sensitivity_df[sensitivity_df['mutation'] == mutation]
        comp_vals = mut_data['Compensation_pct'].values

        ax.scatter(
            [i] * len(comp_vals),
            comp_vals,
            alpha=0.3,
            s=50,
            c=['#1f77b4', '#ff7f0e'][i]
        )

        # Add range annotation
        min_val, max_val = comp_vals.min(), comp_vals.max()
        median_val = np.median(comp_vals)
        ax.hlines(median_val, i-0.2, i+0.2, colors='black', linewidths=2)
        ax.text(i+0.3, median_val, f'Median: {median_val:.1f}%\nRange: {min_val:.1f}-{max_val:.1f}%',
                fontsize=9, va='center')

    ax.set_xticks([0, 1])
    ax.set_xticklabels(['G32A', 'R403C'])
    ax.set_ylabel('Compensation %')
    ax.set_title('D. Compensation % Range Across Thresholds')

    plt.tight_layout()
    plt.savefig(output_dir / 'sensitivity_analysis_summary.pdf', dpi=150)
    plt.savefig(output_dir / 'sensitivity_analysis_summary.png', dpi=150)
    plt.close()
    print("  Saved summary figure")


# =============================================================================
# MAIN
# =============================================================================

def main():
    print("=" * 70)
    print("SENSITIVITY ANALYSIS FOR PATTERN CLASSIFICATION")
    print("=" * 70)

    # Load data
    print("\n1. Loading GSEA results (wide format)...")
    wide_file = DATA_DIR / "Python_exports" / "gsea_results_wide.csv"
    df = pd.read_csv(wide_file)
    print(f"   Loaded {len(df)} unique pathways")

    # Run sensitivity analysis
    print("\n2. Running sensitivity analysis...")
    sensitivity_df = run_sensitivity_analysis(df)
    print(f"   Generated {len(sensitivity_df)} classification results")

    # Analyze key claims
    print("\n3. Analyzing claim stability...")
    claims_df = analyze_key_claims(sensitivity_df)

    # Save results
    print("\n4. Saving results...")
    sensitivity_df.to_csv(OUTPUT_DIR / 'sensitivity_results.csv', index=False)
    claims_df.to_csv(OUTPUT_DIR / 'claim_stability.csv', index=False)
    print(f"   Saved to {OUTPUT_DIR}")

    # Create visualizations
    print("\n5. Creating visualizations...")
    for mutation in ['G32A', 'R403C']:
        create_sensitivity_heatmap(sensitivity_df, mutation, OUTPUT_DIR)
    create_summary_figure(sensitivity_df, claims_df, OUTPUT_DIR)

    # Print summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    default_data = sensitivity_df[sensitivity_df['is_default']]
    print("\nDefault threshold results:")
    for mutation in ['G32A', 'R403C']:
        row = default_data[default_data['mutation'] == mutation].iloc[0]
        print(f"\n  {mutation}:")
        print(f"    Compensation: {row['Compensation']:.0f} ({row['Compensation_pct']:.1f}%)")
        print(f"    Progressive:  {row['Progressive']:.0f} ({row['Progressive_pct']:.1f}%)")
        print(f"    Passive:      {row['Passive_total']:.0f} ({row['Passive_pct']:.1f}%)")
        print(f"    Late_onset:   {row['Late_onset']:.0f} ({row['Late_onset_pct']:.1f}%)")
        print(f"    Transient:    {row['Transient']:.0f} ({row['Transient_pct']:.1f}%)")
        print(f"    Complex:      {row['Complex']:.0f} ({row['Complex_pct']:.1f}%)")

    print("\nClaim stability across threshold variations:")
    n_combos = len(claims_df)

    # Key claims
    claims_summary = {
        'R403C shows more compensation than G32A': claims_df['R403C_more_compensation'].sum(),
        'Progressive patterns rare for G32A (<5%)': claims_df['progressive_rare_G32A'].sum(),
        'Progressive patterns rare for R403C (<5%)': claims_df['progressive_rare_R403C'].sum(),
        'Compensation > Passive for G32A': claims_df['comp_exceeds_passive_G32A'].sum(),
        'Compensation > Passive for R403C': claims_df['comp_exceeds_passive_R403C'].sum(),
    }

    for claim, count in claims_summary.items():
        pct = count / n_combos * 100
        status = "✓ STABLE" if pct == 100 else ("⚠ MOSTLY STABLE" if pct > 80 else "✗ UNSTABLE")
        print(f"  {status} ({pct:.0f}%): {claim}")

    # Compensation percentage range
    print("\n  Compensation percentage range across thresholds:")
    for mutation in ['G32A', 'R403C']:
        mut_data = sensitivity_df[sensitivity_df['mutation'] == mutation]
        min_pct = mut_data['Compensation_pct'].min()
        max_pct = mut_data['Compensation_pct'].max()
        median_pct = mut_data['Compensation_pct'].median()
        print(f"    {mutation}: {min_pct:.1f}% - {max_pct:.1f}% (median: {median_pct:.1f}%)")

    # Compensation as fraction of classifiable pathways
    print("\n  Compensation as % of classifiable (non-Complex) pathways:")
    for mutation in ['G32A', 'R403C']:
        mut_claims = claims_df[[f'comp_of_classifiable_{mutation}']].values.flatten()
        print(f"    {mutation}: {mut_claims.min():.1f}% - {mut_claims.max():.1f}% (median: {np.median(mut_claims):.1f}%)")

    print("\n" + "=" * 70)
    print(f"Results saved to: {OUTPUT_DIR}")
    print("=" * 70)


if __name__ == '__main__':
    main()
