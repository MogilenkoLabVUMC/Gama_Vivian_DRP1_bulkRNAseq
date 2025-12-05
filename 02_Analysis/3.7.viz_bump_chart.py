#!/usr/bin/env python3
"""
Bump Chart Visualization Orchestrator
=====================================

Orchestrates the generation of pathway trajectory bump charts.
Uses the refactored modular logic in `Python.viz_bump_charts`.

Generated Outputs:
1. Uniform: Baseline visualization
2. Weighted: Line thickness based on pattern frequency
3. Highlighted: Top significant pathways labeled
4. Curved: Trajectory deviation visualization
5. Combined: Final paper version (Weighted + Curved + Highlighted)
"""

import sys
from pathlib import Path
import matplotlib.pyplot as plt

# Add module paths
sys.path.insert(0, str(Path(__file__).parent.parent / '01_Scripts'))

from Python.config import resolve_path
from Python.viz_bump_charts import (
    load_data,
    filter_by_scope,
    create_bump_chart,
    create_unified_legend,
    compute_weight_categories,
    BumpChartConfig
)

OUTPUT_ROOT = resolve_path('03_Results/02_Analysis/Plots/Trajectory_Flow')
OUTPUT_BUMP = OUTPUT_ROOT / 'bump'  # Non-key bump charts go to subfolder
OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
OUTPUT_BUMP.mkdir(parents=True, exist_ok=True)

# Key files that stay in root folder (for paper)
KEY_FILES = {
    ('focused', 'curved_nes'),
    ('focused', 'FINAL_paper_combined'),
    ('significant', 'curved_nes'),
}

def main():
    print("=" * 80)
    print("BUMP CHART VISUALIZATION (Refactored)")
    print("=" * 80)
    
    df = load_data()
    
    # Define scopes to run
    scopes = ['focused', 'significant'] 
    
    # 1. UNIFORM
    print("\n--- Generating Uniform Bump Charts ---")
    for scope in scopes:
        cfg = BumpChartConfig(scope=scope, mode='uniform', y_type='nes')
        generate_plot(df, cfg, suffix='uniform_nes')
        
        cfg_rank = BumpChartConfig(scope=scope, mode='uniform', y_type='rank')
        generate_plot(df, cfg_rank, suffix='uniform_rank')

    # 2. WEIGHTED
    print("\n--- Generating Weighted Bump Charts ---")
    for scope in scopes:
        cfg = BumpChartConfig(scope=scope, mode='weighted', y_type='nes')
        generate_plot(df, cfg, suffix='weighted_nes')
        
        cfg_rank = BumpChartConfig(scope=scope, mode='weighted', y_type='rank')
        generate_plot(df, cfg_rank, suffix='weighted_rank')

    # 3. HIGHLIGHTED
    print("\n--- Generating Highlighted Bump Charts ---")
    for scope in scopes:
        # Highlighted usually implies weighted for context
        cfg = BumpChartConfig(scope=scope, mode='weighted', show_highlights=True, y_type='nes')
        generate_plot(df, cfg, suffix='highlight_nes')

    # 4. CURVED
    print("\n--- Generating Curved Bump Charts ---")
    for scope in scopes:
        # Curved NES
        cfg = BumpChartConfig(scope=scope, mode='weighted', show_curves=True, y_type='nes')
        generate_plot(df, cfg, suffix='curved_nes')
        
        # Curved Rank (New Request)
        if scope == 'significant': # User specifically asked for this variation
            cfg_rank = BumpChartConfig(scope=scope, mode='weighted', show_curves=True, y_type='rank')
            generate_plot(df, cfg_rank, suffix='curved_rank')

    # 5. COMBINED (FINAL PAPER VERSION)
    print("\n--- Generating Combined (Paper) Bump Charts ---")
    # Weighted + Curved + Highlighted
    for scope in scopes:
        cfg = BumpChartConfig(
            scope=scope,
            mode='weighted',
            y_type='nes',
            show_highlights=True,
            show_curves=True,
            label_x_pos=1.5 # Tweak for paper tightness
        )
        generate_plot(df, cfg, suffix='FINAL_paper_combined')

    print("\nDone! Outputs in:", OUTPUT_DIR)

def generate_plot(df, config, suffix):
    """Generate a side-by-side plot for G32A and R403C."""
    print(f"  > {suffix} | {config.scope}")
    
    # Filter Data
    df_scope = filter_by_scope(df, config.scope)
    if len(df_scope) == 0:
        print("    Skipping: No data")
        return

    # Setup Figure
    # Wider if highlighted
    fig_width = 24 if config.show_highlights else 16
    width_ratios = [5, 5, 1.5] if config.show_highlights else [3.5, 3.5, 1.5]
    
    fig, axes = plt.subplots(1, 3, figsize=(fig_width, 8), 
                             gridspec_kw={'width_ratios': width_ratios})
    
    # Compute Weights globally for this scope
    # (Ensures consistent coloring/thickness across both subplots)
    weight_cats = None
    if config.mode == 'weighted':
        weight_cats = compute_weight_categories(df_scope)

    # Plot G32A
    create_bump_chart(df_scope, 'G32A', axes[0], config, weight_cats=weight_cats)
    axes[0].set_ylabel('NES' if config.y_type == 'nes' else 'Rank', fontsize=12)
    
    # Plot R403C
    create_bump_chart(df_scope, 'R403C', axes[1], config, weight_cats=weight_cats)
    
    # Share Y Limits (for NES)
    if config.y_type == 'nes':
        y_min = min(axes[0].get_ylim()[0], axes[1].get_ylim()[0])
        y_max = max(axes[0].get_ylim()[1], axes[1].get_ylim()[1])
        axes[0].set_ylim(y_min, y_max)
        axes[1].set_ylim(y_min, y_max)

    # Legend
    create_unified_legend(axes[2], df_scope, weight_cats=weight_cats)
    
    # Title
    # Use scope name from new config defaults if available, or fallback
    # Actually config object has scope_defaults but not the pretty name.
    # We can import SCOPE_CONFIGS or just format the string.
    title_suffix = suffix.replace('_', ' ').title()
    fig.suptitle(f"Bump Chart | {config.scope.capitalize()} | {title_suffix}", 
                 fontsize=16, fontweight='bold', y=0.98)
    
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    
    # Determine output directory: key files go to root, others to bump/
    filename = f"bump_{config.scope}_{suffix}"
    is_key_file = (config.scope, suffix) in KEY_FILES
    output_dir = OUTPUT_ROOT if is_key_file else OUTPUT_BUMP

    # Use try-catch for saving to avoid crashing on permission/path issues
    try:
        fig.savefig(output_dir / f"{filename}.png", dpi=150, bbox_inches='tight', facecolor='white')
        fig.savefig(output_dir / f"{filename}.pdf", dpi=300, bbox_inches='tight', facecolor='white')
        location = "(ROOT)" if is_key_file else "(bump/)"
        print(f"    Saved: {filename} {location}")
    except Exception as e:
        print(f"    Error saving {filename}: {e}")

    plt.close(fig)

if __name__ == '__main__':
    main()