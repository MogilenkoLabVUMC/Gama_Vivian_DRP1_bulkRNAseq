#!/usr/bin/env python3
"""
Trajectory Flow Visualization
==============================

Creates visualizations showing the FLOW of pathway dynamics from
early defects through maturation to late outcomes.

Complements 3.4.pattern_summary_normalized.py (pattern counts) with
trajectory-focused views that reveal:
1. How early defects flow into different outcomes
2. Active vs passive compensation mechanisms
3. Late_onset as a distinct phenomenon (no early input)

Figures:
- Fig A: Alluvial/Sankey diagram (binary left nodes)
- Fig B: Alluvial/Sankey diagram (graded severity)
- Fig C: Trajectory heatmap (per-pathway journey)
- Fig D: Violin plots (magnitude distributions by pattern/stage)
"""

import sys
from pathlib import Path

# Add module paths
sys.path.insert(0, str(Path(__file__).parent.parent / '01_Scripts'))

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

# Try to import plotly for interactive Sankey
try:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False
    print("Warning: plotly not available. Using matplotlib for Sankey diagrams.")

# Import project modules
from Python.config import resolve_path
from Python.semantic_categories import PATTERN_COLORS, MUTATION_COLORS
from Python.pattern_definitions import (
    PADJ_SIGNIFICANT, NES_EFFECT, NES_STRONG, IMPROVEMENT_RATIO, WORSENING_RATIO,
    MEANINGFUL_PATTERNS, ACTIVE_PATTERNS, PASSIVE_PATTERNS
)

# =============================================================================
# CONFIGURATION
# =============================================================================

OUTPUT_ROOT = resolve_path('03_Results/02_Analysis/Plots/Trajectory_Flow')
OUTPUT_DIR = OUTPUT_ROOT / 'alluvial'  # Alluvial diagrams go to alluvial subfolder
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Note: Thresholds (PADJ_SIGNIFICANT, NES_EFFECT, etc.) imported from pattern_definitions.py

# Color schemes
FLOW_COLORS = {
    # Early severity
    'Strong_defect': '#d62728',      # Red
    'Moderate_defect': '#ff7f0e',    # Orange
    'No_defect': '#2ca02c',          # Green
    'No_Early_Defect': '#e0e0e0',    # Light Gray (for Late_onset start)
    # Mechanisms
    'Active': '#1f77b4',             # Blue (darker)
    'Passive': '#aec7e8',            # Light blue
    'Flip_sign': '#9467bd',          # Purple (Sign reversal)
    'Late_onset': '#7f7f7f',         # Gray (Mechanism for Late_onset)
    # Late outcomes
    'Improved': '#98df8a',           # Light green
    'Resolved': '#2ca02c',           # Green
    'Worsened': '#d62728',           # Red
    'New_defect': '#d62728',         # Red (Late_onset outcome)
    # Special
    'Complex': '#999999',            # Gray
}

# =============================================================================
# DATA LOADING AND PREPROCESSING
# =============================================================================

def load_pathway_data():
    """Load pathway data from master GSEA table and compute derived columns."""
    print("Loading pathway data from master_gsea_table.csv...")

    master_file = resolve_path('03_Results/02_Analysis/master_gsea_table.csv')
    df = pd.read_csv(master_file)

    # Pivot to wide format (one row per pathway)
    id_cols = ['pathway_id', 'database', 'Description']
    pattern_cols = ['Pattern_G32A', 'Confidence_G32A', 'Pattern_R403C', 'Confidence_R403C']
    nes_cols = ['NES_Early_G32A', 'NES_Early_R403C',
                'NES_TrajDev_G32A', 'NES_TrajDev_R403C',
                'NES_Late_G32A', 'NES_Late_R403C']

    # Get p.adjust columns from original contrasts
    padj_cols = []
    for contrast in ['G32A_vs_Ctrl_D35', 'G32A_vs_Ctrl_D65', 'Maturation_G32A_specific',
                     'R403C_vs_Ctrl_D35', 'R403C_vs_Ctrl_D65', 'Maturation_R403C_specific']:
        padj_cols.append(f'p.adjust')

    keep_cols = id_cols + pattern_cols + nes_cols
    available_cols = [c for c in keep_cols if c in df.columns]

    df_wide = df[available_cols].drop_duplicates(subset=['pathway_id', 'database']).copy()

    print(f"  Loaded {len(df_wide)} unique pathways")

    # We need p.adjust columns - let's get them from the long format
    # Load separate file for p.adjust values
    padj_file = resolve_path('03_Results/02_Analysis/Python_exports/gsea_results_wide.csv')
    if padj_file.exists():
        df_padj = pd.read_csv(padj_file)
        # Merge p.adjust columns
        padj_cols_to_merge = [c for c in df_padj.columns if 'p.adjust' in c]
        merge_cols = ['pathway_id', 'database'] + padj_cols_to_merge
        available_merge_cols = [c for c in merge_cols if c in df_padj.columns]

        if len(available_merge_cols) > 2:  # More than just ID cols
            df_wide = df_wide.merge(
                df_padj[available_merge_cols].drop_duplicates(),
                on=['pathway_id', 'database'],
                how='left'
            )
            print(f"  Merged p.adjust columns: {len(padj_cols_to_merge)} columns")

    return df_wide


def compute_derived_columns(df, mutation='G32A'):
    """Compute derived columns for flow analysis."""
    print(f"\nComputing derived columns for {mutation}...")

    df = df.copy()

    # Column name mappings
    early_nes = f'NES_Early_{mutation}'
    trajdev_nes = f'NES_TrajDev_{mutation}'
    late_nes = f'NES_Late_{mutation}'
    pattern_col = f'Pattern_{mutation}'

    # Try different p.adjust column naming conventions
    early_padj_options = [
        f'p.adjust_Early_{mutation}',
        f'p.adjust_{mutation}_vs_Ctrl_D35',
        'p.adjust_G32A_vs_Ctrl_D35' if mutation == 'G32A' else 'p.adjust_R403C_vs_Ctrl_D35'
    ]
    trajdev_padj_options = [
        f'p.adjust_TrajDev_{mutation}',
        f'p.adjust_Maturation_{mutation}_specific'
    ]

    # Find available p.adjust columns
    early_padj = None
    for col in early_padj_options:
        if col in df.columns:
            early_padj = col
            break

    trajdev_padj = None
    for col in trajdev_padj_options:
        if col in df.columns:
            trajdev_padj = col
            break

    print(f"  Using Early NES: {early_nes}")
    print(f"  Using TrajDev NES: {trajdev_nes}")
    print(f"  Using Early p.adjust: {early_padj}")
    print(f"  Using TrajDev p.adjust: {trajdev_padj}")

    # --- 1. Has Early Defect ---
    if early_padj and early_padj in df.columns:
        df[f'Has_Early_Defect_{mutation}'] = (
            (df[early_padj] < PADJ_SIGNIFICANT) &
            (df[early_nes].abs() > NES_EFFECT)
        )
    else:
        # Fallback: use pattern classification
        df[f'Has_Early_Defect_{mutation}'] = df[pattern_col].isin([
            'Compensation', 'Progressive', 'Natural_improvement',
            'Natural_worsening', 'Transient'
        ])

    # --- 2. Early Severity ---
    def categorize_early_severity(row):
        nes = row[early_nes] if pd.notna(row[early_nes]) else 0
        has_defect = row[f'Has_Early_Defect_{mutation}']
        pattern = row[pattern_col]

        if pattern == 'Late_onset':
             return 'No_Early_Defect' # Special category for Late_onset start
        
        if not has_defect:
            return 'No_defect'
        elif abs(nes) > NES_STRONG:
            return 'Strong_defect'
        elif abs(nes) > NES_EFFECT:
            return 'Moderate_defect'
        else:
            return 'No_defect'

    df[f'Early_Severity_{mutation}'] = df.apply(categorize_early_severity, axis=1)

    # --- 3. Mechanism (Active vs Passive) ---
    # Need to handle Late_onset and Sign_reversal specifically
    def categorize_mechanism(row):
        pattern = row[pattern_col]
        
        if pattern == 'Late_onset':
            return 'Late_onset'
        elif pattern == 'Sign_reversal':
            return 'Flip_sign'
        
        # Standard Active/Passive logic
        is_active = False
        if trajdev_padj and trajdev_padj in df.columns:
            is_active = (
                (row[trajdev_padj] < PADJ_SIGNIFICANT) &
                (abs(row[trajdev_nes]) > NES_EFFECT)
            )
        else:
            is_active = pattern in ACTIVE_PATTERNS
            
        return 'Active' if is_active else 'Passive'

    df[f'Mechanism_{mutation}'] = df.apply(categorize_mechanism, axis=1)

    # --- 4. Late Outcome ---
    def categorize_late_outcome(row):
        pattern = row[pattern_col]
        early_nes_val = abs(row[early_nes]) if pd.notna(row[early_nes]) else 0
        late_nes_val = abs(row[late_nes]) if pd.notna(row[late_nes]) else 0

        if pattern == 'Late_onset':
            return 'New_defect'
        elif pattern == 'Transient':
            return 'Resolved'
        elif pattern in ['Natural_worsening', 'Progressive']:
            return 'Worsened'
        elif pattern in ['Compensation', 'Natural_improvement']:
            if late_nes_val < NES_EFFECT:
                return 'Resolved'
            else:
                return 'Improved'
        elif pattern == 'Sign_reversal':
            # Sign reversal: Direction flipped. Categorize by late outcome magnitude
            if late_nes_val < NES_EFFECT:
                return 'Resolved'
            elif late_nes_val < early_nes_val:
                return 'Improved'  # Magnitude reduced
            else:
                return 'Worsened'  # Magnitude increased
        else:
            return 'Other'

    df[f'Late_Outcome_{mutation}'] = df.apply(categorize_late_outcome, axis=1)

    # Print summary
    print(f"\n  Summary for {mutation}:")
    print(f"    Has Early Defect: {df[f'Has_Early_Defect_{mutation}'].sum()}")
    print(f"    Early Severity:")
    for val, cnt in df[f'Early_Severity_{mutation}'].value_counts().items():
        print(f"      {val}: {cnt}")
    print(f"    Mechanism:")
    for val, cnt in df[f'Mechanism_{mutation}'].value_counts().items():
        print(f"      {val}: {cnt}")
    print(f"    Late Outcome:")
    for val, cnt in df[f'Late_Outcome_{mutation}'].value_counts().items():
        print(f"      {val}: {cnt}")

    return df


# =============================================================================
# SANKEY/ALLUVIAL DIAGRAM
# =============================================================================

def build_sankey_data(df, mutation='G32A', use_graded=False):
    """
    Build node and link data for Sankey diagram.

    Parameters
    ----------
    df : pd.DataFrame
        Pathway data with derived columns
    mutation : str
        G32A or R403C
    use_graded : bool
        If True, use Strong/Moderate/No defect. If False, use binary Has/No defect.

    Returns
    -------
    dict
        Dictionary with 'nodes' and 'links' for Sankey diagram
    """
    pattern_col = f'Pattern_{mutation}'
    severity_col = f'Early_Severity_{mutation}'
    mechanism_col = f'Mechanism_{mutation}'
    outcome_col = f'Late_Outcome_{mutation}'
    has_defect_col = f'Has_Early_Defect_{mutation}'

    # Filter out Complex and Insufficient_data for clarity
    df_flow = df[~df[pattern_col].isin(['Complex', 'Insufficient_data'])].copy()
    
    # In previous version, Late_onset was filtered out or handled separately.
    # Now we include everyone. 
    # But for the "Early" node logic, we need to be careful.

    # --- Define nodes ---
    if use_graded:
        early_nodes = ['Strong_defect', 'Moderate_defect']
    else:
        early_nodes = ['Has_Early_Defect']
    
    # Add 'No_Early_Defect' node for Late_onset pathways to start from
    early_nodes.append('No_Early_Defect')

    mechanism_nodes = ['Active', 'Passive', 'Flip_sign', 'Late_onset']
    outcome_nodes = ['Improved', 'Resolved', 'Worsened', 'New_defect']

    all_nodes = early_nodes + mechanism_nodes + outcome_nodes
    node_indices = {node: i for i, node in enumerate(all_nodes)}

    # --- Build links ---
    links = []

    # Helper to get pathway text (truncated)
    def get_pathway_text(df_subset, limit=20):
        if len(df_subset) == 0:
            return ""
        
        # Use Description if available, else pathway_id
        if 'Description' in df_subset.columns:
            names = df_subset.apply(lambda x: x['Description'] if pd.notna(x['Description']) else x['pathway_id'], axis=1).tolist()
        else:
            names = df_subset['pathway_id'].tolist()
            
        total_count = len(names)
        names.sort()
        
        if total_count > limit:
            shown = names[:limit]
            return "<br>".join(shown) + f"<br>... and {total_count - limit} more"
        else:
            return "<br>".join(names)

    # Links from Early to Mechanism
    for early in early_nodes:
        # Select rows matching this early state
        if early == 'No_Early_Defect':
            # This captures Late_onset primarily
            df_early = df_flow[df_flow[severity_col] == 'No_Early_Defect']
        else:
            # Normal defects
            if use_graded:
                df_early = df_flow[df_flow[severity_col] == early]
            else:
                df_early = df_flow[df_flow[has_defect_col] == True]

        for mechanism in mechanism_nodes:
            df_subset = df_early[df_early[mechanism_col] == mechanism]
            count = len(df_subset)
            if count > 0:
                links.append({
                    'source': node_indices[early],
                    'target': node_indices[mechanism],
                    'value': count,
                    'label': f'{early} → {mechanism}',
                    'pathway_text': get_pathway_text(df_subset)
                })

    # Links from Mechanism to Outcome
    for mechanism in mechanism_nodes:
        df_mech = df_flow[df_flow[mechanism_col] == mechanism]
        for outcome in outcome_nodes:
            df_subset = df_mech[df_mech[outcome_col] == outcome]
            count = len(df_subset)
            if count > 0:
                links.append({
                    'source': node_indices[mechanism],
                    'target': node_indices[outcome],
                    'value': count,
                    'label': f'{mechanism} → {outcome}',
                    'pathway_text': get_pathway_text(df_subset)
                })

    # Node colors
    node_colors = []
    for node in all_nodes:
        if node in FLOW_COLORS:
            node_colors.append(FLOW_COLORS[node])
        elif node == 'Has_Early_Defect':
            node_colors.append(FLOW_COLORS['Strong_defect'])
        else:
            node_colors.append('#999999')

    return {
        'nodes': all_nodes,
        'node_colors': node_colors,
        'links': links,
        'late_onset_count': len(df_flow[df_flow[pattern_col] == 'Late_onset']),
        'node_indices': node_indices
    }


def create_sankey_plotly(sankey_data, mutation='G32A', use_graded=False, title_suffix=''):
    """Create interactive Sankey diagram using Plotly."""
    if not PLOTLY_AVAILABLE:
        print("  Plotly not available, skipping interactive Sankey")
        return None

    nodes = sankey_data['nodes']
    links = sankey_data['links']
    node_colors = sankey_data['node_colors']
    
    # We no longer need the separate annotation for Late_onset since it's in the graph
    # late_onset_count = sankey_data['late_onset_count'] 

    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=nodes,
            color=node_colors
        ),
        link=dict(
            source=[l['source'] for l in links],
            target=[l['target'] for l in links],
            value=[l['value'] for l in links],
            label=[l['label'] for l in links],
            customdata=[l.get('pathway_text', '') for l in links],
            hovertemplate='<b>%{label}</b><br>Count: %{value}<br><br>Pathways:<br>%{customdata}<extra></extra>',
            color=['rgba(100,100,100,0.3)' for _ in links]
        )
    )])

    version = 'Graded Severity' if use_graded else 'Binary'
    fig.update_layout(
        title=dict(
            text=f"{mutation} Trajectory Flow ({version}){title_suffix}",
            font=dict(size=16)
        ),
        font=dict(size=10),  # Reduced font size
        height=800, 
        width=1000
    )

    return fig


def create_sankey_matplotlib(sankey_data, mutation='G32A', use_graded=False):
    """Create static Sankey-like diagram using matplotlib."""

    nodes = sankey_data['nodes']
    links = sankey_data['links']
    late_onset_count = sankey_data['late_onset_count']

    # Create figure
    fig, ax = plt.subplots(figsize=(14, 10))

    # Define column positions
    col_positions = {
        'early': 0.1,
        'mechanism': 0.45,
        'outcome': 0.8
    }

    # Assign columns to nodes
    node_columns = {}
    if use_graded:
        node_columns['Strong_defect'] = 'early'
        node_columns['Moderate_defect'] = 'early'
    else:
        node_columns['Has_Early_Defect'] = 'early'
    node_columns['Active'] = 'mechanism'
    node_columns['Passive'] = 'mechanism'
    node_columns['Improved'] = 'outcome'
    node_columns['Resolved'] = 'outcome'
    node_columns['Worsened'] = 'outcome'
    node_columns['New_defect'] = 'outcome'

    # Calculate node sizes and positions
    # Group nodes by column
    cols = {'early': [], 'mechanism': [], 'outcome': []}
    for node in nodes:
        col = node_columns.get(node)
        if col:
            cols[col].append(node)

    # Calculate total flow per column
    node_totals = {node: 0 for node in nodes}
    for link in links:
        src_node = nodes[link['source']]
        tgt_node = nodes[link['target']]
        node_totals[src_node] += link['value']
        node_totals[tgt_node] += link['value']

    # Add Late_onset to New_defect
    node_totals['New_defect'] += late_onset_count

    # Normalize and position nodes
    node_positions = {}
    for col_name, col_nodes in cols.items():
        total = sum(node_totals[n] for n in col_nodes)
        if total == 0:
            total = 1

        y_pos = 0.9
        for node in col_nodes:
            height = max(0.05, node_totals[node] / total * 0.7)
            node_positions[node] = {
                'x': col_positions[col_name],
                'y': y_pos - height/2,
                'height': height,
                'width': 0.08
            }
            y_pos -= height + 0.05

    # Draw nodes
    for node, pos in node_positions.items():
        color = FLOW_COLORS.get(node, '#999999')
        if node == 'Has_Early_Defect':
            color = FLOW_COLORS['Strong_defect']

        rect = plt.Rectangle(
            (pos['x'], pos['y']), pos['width'], pos['height'],
            facecolor=color, edgecolor='black', linewidth=1, alpha=0.8
        )
        ax.add_patch(rect)

        # Node label
        label = node.replace('_', '\n')
        count = node_totals[node]
        ax.text(
            pos['x'] + pos['width']/2, pos['y'] + pos['height']/2,
            f"{label}\n({count})",
            ha='center', va='center', fontsize=9, fontweight='bold'
        )

    # Draw links as curved paths
    from matplotlib.patches import FancyArrowPatch
    from matplotlib.path import Path as MplPath
    import matplotlib.patches as patches

    for link in links:
        src_node = nodes[link['source']]
        tgt_node = nodes[link['target']]

        if src_node not in node_positions or tgt_node not in node_positions:
            continue

        src_pos = node_positions[src_node]
        tgt_pos = node_positions[tgt_node]

        # Start and end points
        x0 = src_pos['x'] + src_pos['width']
        y0 = src_pos['y'] + src_pos['height']/2
        x1 = tgt_pos['x']
        y1 = tgt_pos['y'] + tgt_pos['height']/2

        # Width proportional to flow
        max_flow = max(l['value'] for l in links)
        width = max(0.5, link['value'] / max_flow * 8)

        # Draw curved line
        verts = [
            (x0, y0),
            (x0 + 0.1, y0),
            (x1 - 0.1, y1),
            (x1, y1)
        ]
        codes = [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4]
        path = MplPath(verts, codes)
        patch = patches.PathPatch(path, facecolor='none', edgecolor='gray',
                                  alpha=0.4, linewidth=width)
        ax.add_patch(patch)

    # Draw Late_onset as separate annotation
    if late_onset_count > 0:
        ax.annotate(
            f'Late_onset\n({late_onset_count} pathways)\nNo early defect',
            xy=(0.7, 0.1), xytext=(0.3, 0.1),
            fontsize=10, ha='center', va='center',
            bbox=dict(boxstyle='round', facecolor=FLOW_COLORS['Late_onset'], alpha=0.7),
            arrowprops=dict(arrowstyle='->', color=FLOW_COLORS['Late_onset'], lw=2)
        )

    # Column labels
    ax.text(col_positions['early'] + 0.04, 0.98, 'Early Defect\nSeverity',
            ha='center', va='top', fontsize=12, fontweight='bold')
    ax.text(col_positions['mechanism'] + 0.04, 0.98, 'Mechanism\n(TrajDev)',
            ha='center', va='top', fontsize=12, fontweight='bold')
    ax.text(col_positions['outcome'] + 0.04, 0.98, 'Late\nOutcome',
            ha='center', va='top', fontsize=12, fontweight='bold')

    # Title
    version = 'Graded Severity' if use_graded else 'Binary'
    ax.set_title(f'{mutation} Trajectory Flow ({version})\n'
                 f'Showing flow from early defects → mechanism → late outcome',
                 fontsize=14, fontweight='bold', color=MUTATION_COLORS.get(mutation, 'black'))

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')

    plt.tight_layout()
    return fig


# =============================================================================
# TRAJECTORY HEATMAP
# =============================================================================

def create_trajectory_heatmap(df, mutation='G32A', patterns=None, filename_suffix=''):
    """
    Create trajectory heatmap showing per-pathway NES journey.

    Rows = pathways grouped by Pattern
    Columns = Early | TrajDev | Late
    Color = NES value
    """
    print(f"\nCreating trajectory heatmap for {mutation} (suffix='{filename_suffix}')...")

    pattern_col = f'Pattern_{mutation}'
    
    early_nes = f'NES_Early_{mutation}'
    trajdev_nes = f'NES_TrajDev_{mutation}'
    late_nes = f'NES_Late_{mutation}'

    # Filter to specific patterns
    if patterns is None:
        # Default to all classifiable patterns
        patterns_to_show = ['Compensation', 'Natural_improvement', 'Late_onset', 
                           'Natural_worsening', 'Progressive', 'Transient']
        title_extra = "grouped by pattern"
    else:
        patterns_to_show = patterns
        title_extra = "focused view"

    df_heat = df[df[pattern_col].isin(patterns_to_show)].copy()

    if len(df_heat) == 0:
        print(f"  No pathways with patterns: {patterns_to_show}")
        return None

    # Create matrix
    df_matrix = df_heat[[pattern_col, early_nes, trajdev_nes, late_nes]].copy()
    df_matrix = df_matrix.rename(columns={
        early_nes: 'Early',
        trajdev_nes: 'TrajDev',
        late_nes: 'Late'
    })

    # Sort by Pattern, then by Early NES within pattern
    # Define custom sort order for patterns based on input list order
    pattern_order = {k: v for v, k in enumerate(patterns_to_show)}
    df_matrix['sort_key'] = df_matrix[pattern_col].map(pattern_order)
    
    df_matrix = df_matrix.sort_values(['sort_key', 'Early'], ascending=[True, False])

    # Get pattern labels for row annotation
    pattern_labels = df_matrix[pattern_col].values

    # Create matrix for heatmap
    matrix = df_matrix[['Early', 'TrajDev', 'Late']].values

    # Dynamic figure height
    # Base height per row, plus some overhead
    num_rows = len(df_matrix)
    row_height = 0.005  # inches per row
    calc_height = max(4, num_rows * row_height + 2) # Min 4 inches, or scaled
    
    if calc_height > 20:
        calc_height = 20 # Cap max height
        
    fig, ax = plt.subplots(figsize=(8, calc_height))

    # Create diverging colormap
    cmap = sns.diverging_palette(240, 10, as_cmap=True)

    # Plot heatmap
    vmax = 3.0
    im = ax.imshow(matrix, aspect='auto', cmap=cmap, vmin=-vmax, vmax=vmax)

    # X-axis labels
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(['Early', 'TrajDev', 'Late'], fontsize=12, fontweight='bold')

    # Y-axis: show pattern group boundaries
    pattern_changes = [0]
    current_pattern = pattern_labels[0]
    for i, p in enumerate(pattern_labels):
        if p != current_pattern:
            pattern_changes.append(i)
            current_pattern = p
    pattern_changes.append(len(pattern_labels))

    # Draw horizontal lines at pattern boundaries
    # Use zorder to ensure visibility
    for i in pattern_changes[1:-1]:
        ax.axhline(y=i-0.5, color='white', linewidth=2, zorder=10)

    # Add pattern labels on the right
    for i in range(len(pattern_changes) - 1):
        start = pattern_changes[i]
        end = pattern_changes[i + 1]
        mid = (start + end) / 2
        pattern = pattern_labels[start]

        color = PATTERN_COLORS.get(pattern, 'black')
        
        # Add count to label
        count = end - start
        label = f"{pattern.replace('_', ' ')}\n(n={count})"
        
        ax.text(3.2, mid, label,
                ha='left', va='center', fontsize=10, fontweight='bold',
                color=color)

    # Colorbar
    cbar = plt.colorbar(im, ax=ax, shrink=0.5, aspect=20)
    cbar.set_label('NES', fontsize=11)

    # Title
    ax.set_title(f'{mutation} Trajectory Heatmap\n'
                 f'{len(df_heat)} pathways ({title_extra})',
                 fontsize=14, fontweight='bold',
                 color=MUTATION_COLORS.get(mutation, 'black'))

    ax.set_ylabel(f'Pathways (n={len(df_heat)})', fontsize=11)
    ax.set_yticks([])

    plt.tight_layout()
    
    # Save directly here since we have the suffix
    out_name = f'trajectory_heatmap_{mutation}{filename_suffix}'
    fig.savefig(OUTPUT_DIR / f'{out_name}.pdf', dpi=300, bbox_inches='tight')
    fig.savefig(OUTPUT_DIR / f'{out_name}.png', dpi=150, bbox_inches='tight')
    print(f"  Saved: {out_name}.pdf/png")
    
    plt.close(fig)
    return True


# =============================================================================
# VIOLIN PLOTS
# =============================================================================

def create_violin_plots(df, mutations=['G32A', 'R403C']):
    """
    Create violin plots showing NES magnitude distributions by pattern and stage.
    """
    print("\nCreating violin plots...")

    # Prepare long-format data for plotting
    plot_data = []

    for mutation in mutations:
        pattern_col = f'Pattern_{mutation}'
        early_nes = f'NES_Early_{mutation}'
        trajdev_nes = f'NES_TrajDev_{mutation}'
        late_nes = f'NES_Late_{mutation}'

        patterns_to_show = ['Compensation', 'Natural_improvement', 'Late_onset']
        df_sub = df[df[pattern_col].isin(patterns_to_show)].copy()

        for _, row in df_sub.iterrows():
            pattern = row[pattern_col]

            for stage, col in [('Early', early_nes), ('TrajDev', trajdev_nes), ('Late', late_nes)]:
                if pd.notna(row[col]):
                    plot_data.append({
                        'Mutation': mutation,
                        'Pattern': pattern,
                        'Stage': stage,
                        'abs_NES': abs(row[col]),
                        'NES': row[col]
                    })

    df_plot = pd.DataFrame(plot_data)

    if len(df_plot) == 0:
        print("  No data for violin plots")
        return None

    # Create figure with facets for each stage
    fig, axes = plt.subplots(1, 3, figsize=(15, 6), sharey=True)

    stages = ['Early', 'TrajDev', 'Late']
    patterns_order = ['Compensation', 'Natural_improvement', 'Late_onset']

    for i, stage in enumerate(stages):
        ax = axes[i]
        df_stage = df_plot[df_plot['Stage'] == stage]

        # Create violin plot
        sns.violinplot(
            data=df_stage,
            x='Pattern',
            y='abs_NES',
            hue='Mutation',
            split=True,
            ax=ax,
            palette=MUTATION_COLORS,
            order=patterns_order,
            inner='quartile'
        )

        # Add strip plot for individual points
        sns.stripplot(
            data=df_stage,
            x='Pattern',
            y='abs_NES',
            hue='Mutation',
            dodge=True,
            ax=ax,
            palette=MUTATION_COLORS,
            order=patterns_order,
            size=2,
            alpha=0.3,
            legend=False
        )

        ax.set_title(f'{stage} Stage', fontsize=12, fontweight='bold')
        ax.set_xlabel('')
        ax.set_ylabel('|NES|' if i == 0 else '')
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')

        # Add horizontal line at thresholds
        ax.axhline(y=NES_EFFECT, color='gray', linestyle='--', alpha=0.5, label='Effect threshold')
        ax.axhline(y=NES_STRONG, color='gray', linestyle=':', alpha=0.5, label='Strong threshold')

        if i == 2:
            ax.legend(loc='upper right', title='Mutation')
        else:
            ax.get_legend().remove() if ax.get_legend() else None

    fig.suptitle('NES Magnitude Distributions by Pattern and Stage\n'
                 'Compensation shows large Early + large TrajDev + small Late\n'
                 'Natural_improvement shows large Early + small TrajDev + small Late\n'
                 'Late_onset shows small Early + variable TrajDev + large Late',
                 fontsize=12, fontweight='bold')

    plt.tight_layout(rect=[0, 0, 1, 0.92])
    return fig


# =============================================================================
# MAIN EXECUTION
# =============================================================================

def main():
    """Generate all trajectory flow visualizations."""
    print("=" * 80)
    print("TRAJECTORY FLOW VISUALIZATION")
    print("=" * 80)
    print(f"Output directory: {OUTPUT_DIR}")

    # Load and preprocess data
    df = load_pathway_data()

    # Process both mutations
    for mutation in ['G32A', 'R403C']:
        df = compute_derived_columns(df, mutation)

    # --- Generate visualizations ---

    for mutation in ['G32A', 'R403C']:
        print(f"\n{'='*40}")
        print(f"Generating figures for {mutation}")
        print(f"{'='*40}")

        # 1. Sankey diagram (binary)
        print("\n1. Creating Sankey diagram (binary)...")
        sankey_data = build_sankey_data(df, mutation, use_graded=False)

        if PLOTLY_AVAILABLE:
            fig_sankey = create_sankey_plotly(sankey_data, mutation, use_graded=False)
            if fig_sankey:
                fig_sankey.write_html(OUTPUT_DIR / f'alluvial_binary_{mutation}.html')
                fig_sankey.write_image(OUTPUT_DIR / f'alluvial_binary_{mutation}.png', scale=2)
                print(f"  Saved: alluvial_binary_{mutation}.html/png")

        # fig_sankey_mpl = create_sankey_matplotlib(sankey_data, mutation, use_graded=False)
        # fig_sankey_mpl.savefig(OUTPUT_DIR / f'alluvial_binary_{mutation}.pdf', dpi=300, bbox_inches='tight')
        # plt.close(fig_sankey_mpl)
        # print(f"  Saved: alluvial_binary_{mutation}.pdf")

        # 2. Sankey diagram (graded severity)
        print("\n2. Creating Sankey diagram (graded severity)...")
        sankey_data_graded = build_sankey_data(df, mutation, use_graded=True)

        if PLOTLY_AVAILABLE:
            fig_sankey_gr = create_sankey_plotly(sankey_data_graded, mutation, use_graded=True)
            if fig_sankey_gr:
                fig_sankey_gr.write_html(OUTPUT_DIR / f'alluvial_graded_{mutation}.html')
                fig_sankey_gr.write_image(OUTPUT_DIR / f'alluvial_graded_{mutation}.png', scale=2)
                print(f"  Saved: alluvial_graded_{mutation}.html/png")

        # fig_sankey_gr_mpl = create_sankey_matplotlib(sankey_data_graded, mutation, use_graded=True)
        # fig_sankey_gr_mpl.savefig(OUTPUT_DIR / f'alluvial_graded_{mutation}.pdf', dpi=300, bbox_inches='tight')
        # plt.close(fig_sankey_gr_mpl)
        # print(f"  Saved: alluvial_graded_{mutation}.pdf")

        # 3. Trajectory heatmap (Full)
        print("\n3. Creating trajectory heatmap (Full)...")
        create_trajectory_heatmap(df, mutation)
        
        # 4. Trajectory heatmap (Zoomed - No Improved)
        print("\n4. Creating trajectory heatmap (Zoomed - No Improved)...")
        create_trajectory_heatmap(
            df, 
            mutation, 
            patterns=['Late_onset', 'Natural_worsening', 'Progressive', 'Transient'],
            filename_suffix='_zoom_no_improved'
        )

    # 5. Violin plots (both mutations together)
    print("\n5. Creating violin plots...")
    fig_violin = create_violin_plots(df, mutations=['G32A', 'R403C'])
    if fig_violin:
        fig_violin.savefig(OUTPUT_DIR / 'violin_magnitude_by_pattern.pdf',
                           dpi=300, bbox_inches='tight')
        fig_violin.savefig(OUTPUT_DIR / 'violin_magnitude_by_pattern.png',
                           dpi=150, bbox_inches='tight')
        plt.close(fig_violin)
        print("  Saved: violin_magnitude_by_pattern.pdf/png")

    # 6. Generate README
    generate_readme(df)

    print("\n" + "=" * 80)
    print("COMPLETE!")
    print("=" * 80)
    print(f"\nOutput directory: {OUTPUT_DIR}")
    print("\nGenerated files:")
    for f in sorted(OUTPUT_DIR.glob('*')):
        print(f"  - {f.name}")


def generate_readme(df):
    """Generate README documentation for the alluvial subfolder."""
    readme_content = """# Alluvial Diagrams (Sankey Flows)

## Purpose

These alluvial/Sankey diagrams show the **FLOW** of pathway dynamics from early defects
through maturation mechanisms to late outcomes.

## Files

### Interactive Alluvial (Plotly)
- `alluvial_binary_G32A.html` - Binary view: Has Early Defect vs No Defect (Interactive)
- `alluvial_binary_R403C.html` - Same for R403C mutation
- `alluvial_graded_G32A.html` - Graded view: Strong/Moderate/No defect
- `alluvial_graded_R403C.html` - Same for R403C mutation
- Corresponding `.png` files for static viewing

### Classical Alluvial (ggalluvial - R)
- `alluvial_ggalluvial_G32A.pdf/png` - Classical alluvial for G32A
- `alluvial_ggalluvial_R403C.pdf/png` - Classical alluvial for R403C
- `alluvial_combined.pdf/png` - Side-by-side comparison

## Structure
- **LEFT**: Early defect status/severity
- **MIDDLE**: Mechanism (Active TrajDev vs Passive)
- **RIGHT**: Late outcome (Improved, Resolved, Worsened)
- **Late_onset**: Shown as separate stream with no left-side input

## Scripts
- `02_Analysis/3.5.viz_trajectory_flow.py` - Interactive Plotly alluvial
- `02_Analysis/3.6.viz_alluvial_ggalluvial.R` - Classical ggalluvial

## See Also
- Parent folder README for bump charts and overall documentation
- `docs/PATTERN_CLASSIFICATION.md` for pattern definitions

---
Generated: """ + pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')

    with open(OUTPUT_DIR / 'README.md', 'w') as f:
        f.write(readme_content)
    print("\n  Saved: README.md")


if __name__ == '__main__':
    main()
