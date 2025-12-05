#!/usr/bin/env python3
"""
Test Bump Chart Features
========================

Tests specific features requested in the review:
1. Custom highlighting of arbitrary pathways.
2. Exporting underlying data for information preservation (interactive readiness).
"""

import sys
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd

# Add module paths
sys.path.insert(0, str(Path(__file__).parent.parent / '01_Scripts'))

from Python.config import resolve_path
from Python.viz_bump_charts import (
    load_data,
    filter_by_scope,
    create_bump_chart,
    compute_global_line_weights,
    create_unified_legend
)

OUTPUT_DIR = resolve_path('03_Results/02_Analysis/Plots/Trajectory_Flow')
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

def test_custom_highlighting():
    print("\n=== Testing Custom Highlighting ===")
    df = load_data()
    
    # 1. Find pathways related to "mitophagy" or specific interesting ones
    # Let's search for some known interesting ones or just pick by string
    term = "Mitophagy"
    matches = df[df['Description'].str.contains(term, case=False, na=False)]
    
    if len(matches) == 0:
        print(f"  No pathways found matching '{term}'. Falling back to 'Oxidative'.")
        term = "Oxidative"
        matches = df[df['Description'].str.contains(term, case=False, na=False)]
    
    custom_ids = matches['pathway_id'].unique().tolist()
    print(f"  Found {len(custom_ids)} pathways matching '{term}':")
    for pid in custom_ids[:5]:
        desc = matches[matches['pathway_id'] == pid]['Description'].iloc[0]
        print(f"    - {pid}: {desc}")
        
    if len(custom_ids) > 10:
        print(f"  ...and {len(custom_ids)-10} more. Limiting to top 10 for clarity.")
        custom_ids = custom_ids[:10]

    # 2. Generate Chart with Custom Highlights
    print("  Generating chart with custom highlights...")
    
    fig, axes = plt.subplots(1, 3, figsize=(24, 8), 
                             gridspec_kw={'width_ratios': [5, 5, 1.5]})
    
    # Filter to 'all' scope to ensure our custom pathways are present
    df_scope = filter_by_scope(df, 'all')
    
    # Compute weights
    weights = compute_global_line_weights(df_scope)
    
    # Plot G32A
    create_bump_chart(df_scope, 'G32A', axes[0], y_type='nes', scope='all',
                      weights=weights, highlight_mode='labeled', 
                      custom_highlight_ids=custom_ids)
    
    # Plot R403C
    create_bump_chart(df_scope, 'R403C', axes[1], y_type='nes', scope='all',
                      weights=weights, highlight_mode='labeled',
                      custom_highlight_ids=custom_ids)
    
    # Legend
    create_unified_legend(axes[2], df_scope)
    
    fig.suptitle(f"Bump Chart Test: Custom Highlight ('{term}')", fontsize=16, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    
    out_file = OUTPUT_DIR / 'test_bump_custom_highlight.png'
    fig.savefig(out_file, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"  Saved: {out_file}")
    plt.close(fig)

    # 3. Export Data (Information Preservation Test)
    print("\n=== Testing Data Export (Information Preservation) ===")
    # Extract the data for the highlighted pathways
    df_export = df[df['pathway_id'].isin(custom_ids)].copy()
    
    # Select relevant columns for an "Explorer" tool
    export_cols = [
        'pathway_id', 'Description', 'database',
        'Pattern_G32A', 'Pattern_R403C',
        'NES_Early_G32A', 'NES_TrajDev_G32A', 'NES_Late_G32A',
        'NES_Early_R403C', 'NES_TrajDev_R403C', 'NES_Late_R403C',
        'p.adjust_G32A_vs_Ctrl_D35', 'p.adjust_G32A_vs_Ctrl_D65',
        'p.adjust_R403C_vs_Ctrl_D35', 'p.adjust_R403C_vs_Ctrl_D65'
    ]
    
    # Handle missing columns gracefully
    final_cols = [c for c in export_cols if c in df_export.columns]
    df_export = df_export[final_cols]
    
    csv_file = OUTPUT_DIR / 'test_bump_custom_highlight_data.csv'
    df_export.to_csv(csv_file, index=False)
    print(f"  Saved underlying data to: {csv_file}")
    print("  This demonstrates how visual data can be preserved for interactive exploration.")

if __name__ == '__main__':
    test_custom_highlighting()
