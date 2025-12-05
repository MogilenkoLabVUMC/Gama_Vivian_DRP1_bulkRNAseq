"""
Bump Chart Visualization Module
===============================

Core logic for generating bump charts (slope graphs) of pathway trajectories.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.path import Path
from adjustText import adjust_text

from Python.config import resolve_path
from Python.pattern_definitions import get_pattern_colors, NES_EFFECT, MEANINGFUL_PATTERNS
from Python.semantic_categories import (
    is_relevant_for_highlight,
    get_highlight_priority,
    assign_semantic_category,
    MAX_OTHER_CATEGORY_HIGHLIGHTS,
)
# Import unified color config for consistent mutation colors
from Python.color_config import MUTATION_COLORS

# =============================================================================
# CONFIGURATION & STYLING CLASSES
# =============================================================================

class BumpChartConfig:
    """
    Configuration container for Bump Chart visualization.
    Encapsulates all 'knobs' for easy tuning.
    """
    def __init__(self,
                 scope='focused',
                 y_type='nes',
                 mode='uniform',  # 'uniform', 'weighted'
                 show_highlights=False,
                 show_curves=False,
                 label_truncation=50,
                 label_x_pos=1.35, # Closer to chart
                 curve_strength=0.5):
        
        self.scope = scope
        self.y_type = y_type
        self.mode = mode
        self.show_highlights = show_highlights
        self.show_curves = show_curves
        self.label_truncation = label_truncation
        self.label_x_pos = label_x_pos
        self.curve_strength = curve_strength
        
        # Scope-specific defaults
        self.scope_defaults = {
            'focused': {'alpha': 0.7, 'lw': 1.5, 'expected_n': '~100'},
            'significant': {'alpha': 0.15, 'lw': 0.5, 'expected_n': '~2-4k'},
            'all': {'alpha': 0.05, 'lw': 0.3, 'expected_n': '~12k'},
        }

    def get_base_style(self):
        """Get base style dict for current scope."""
        return self.scope_defaults.get(self.scope, self.scope_defaults['focused'])

class BumpChartStyler:
    """
    Handles logic for line colors, weights, alphas, and zorders.
    """
    def __init__(self, config):
        self.config = config
        self.pattern_colors = get_pattern_colors()
        # Use unified color config (single source of truth)
        self.mutation_colors = MUTATION_COLORS.copy()

    def get_color(self, pattern):
        return self.pattern_colors.get(pattern, '#999999')

    def get_style(self, pattern, weight_category=None, is_highlighted=False):
        """
        Return visualization attributes (linewidth, alpha, zorder) for a line.
        
        weight_category: 'dominant', 'common', 'uncommon', 'rare', or None
        """
        base = self.config.get_base_style()
        base_alpha = base['alpha']
        base_lw = base['lw']
        
        # Default values
        lw = base_lw
        alpha = base_alpha
        zorder = 1
        
        # 1. Apply Weighting Logic (if enabled)
        # We scale relative to the base scope settings to handle different densities
        if self.config.mode == 'weighted' and weight_category:
            if weight_category == 'dominant':  # > 30%
                lw = base_lw * 0.8
                alpha = base_alpha 
                zorder = 1
            elif weight_category == 'common':  # > 10%
                lw = base_lw * 1.2
                alpha = min(1.0, base_alpha * 2.5) # Boost visibility
                zorder = 2
            elif weight_category == 'uncommon': # > 1%
                lw = base_lw * 1.8
                alpha = min(1.0, base_alpha * 4.0)
                zorder = 3
            elif weight_category == 'rare':    # < 1%
                lw = base_lw * 2.5
                alpha = min(1.0, base_alpha * 5.0)
                zorder = 4

        # 2. Apply Highlight Overrides
        if is_highlighted:
            lw = max(lw * 1.5, 2.0)
            alpha = 1.0
            zorder = 10  # Top layer
            
        return {'linewidth': lw, 'alpha': alpha, 'zorder': zorder}

# =============================================================================
# DATA PREPARATION
# =============================================================================

def load_data():
    """Load pathway data from master GSEA table."""
    master_file = resolve_path('03_Results/02_Analysis/master_gsea_table.csv')
    df = pd.read_csv(master_file)
    
    # Add TrajDev Significance for Curve plotting
    for mut in ['G32A', 'R403C']:
        traj_contrast = f'Maturation_{mut}_specific'
        col_name = f'Sig_TrajDev_{mut}'
        if traj_contrast in df['contrast'].values:
            sig_map = df[df['contrast'] == traj_contrast].set_index('pathway_id')['p.adjust'].apply(lambda x: x < 0.05).to_dict()
            df[col_name] = df['pathway_id'].map(sig_map).fillna(False)
        else:
            df[col_name] = False
    return df

def filter_by_scope(df, scope):
    """Filter dataframe based on scope configuration."""
    if scope == 'focused':
        return df[df['database'].isin(['MitoCarta', 'SynGO'])].copy()
    elif scope == 'significant':
        # Filter for meaningful patterns only
        patterns_to_show = [p for p in MEANINGFUL_PATTERNS if p != 'Complex']
        sig_mask = pd.Series(False, index=df.index)
        for mut in ['G32A', 'R403C']:
            pattern_col = f'Pattern_{mut}'
            if pattern_col in df.columns:
                sig_mask |= df[pattern_col].isin(patterns_to_show)
        return df[sig_mask].copy()
    return df.copy()

def compute_weight_categories(df, mutations=['G32A', 'R403C']):
    """
    Compute frequency-based weight categories globally across provided mutations.
    Returns dict: {pattern: 'category'}
    """
    pattern_counts = {}
    for mut in mutations:
        pattern_col = f'Pattern_{mut}'
        if pattern_col not in df.columns: continue
        
        # Count only complete rows
        nes_cols = [f'NES_Early_{mut}', f'NES_Late_{mut}']
        df_valid = df[df[nes_cols].notna().all(axis=1)]
        
        counts = df_valid[pattern_col].value_counts().to_dict()
        for p, c in counts.items():
            pattern_counts[p] = pattern_counts.get(p, 0) + c
            
    total = sum(pattern_counts.values())
    categories = {}
    if total == 0: return categories

    for pattern, count in pattern_counts.items():
        freq = count / total
        if freq > 0.3: categories[pattern] = 'dominant'
        elif freq > 0.1: categories[pattern] = 'common'
        elif freq > 0.01: categories[pattern] = 'uncommon'
        else: categories[pattern] = 'rare'
            
    return categories

def identify_highlight_pathways(df, mutation, max_per_pattern=5):
    """
    Identify pathways to highlight based on:
    1. Biological Relevance (Semantic Category)
    2. Priority
    3. Magnitude of Change
    """
    pattern_col = f'Pattern_{mutation}'
    early_col = f'NES_Early_{mutation}'
    late_col = f'NES_Late_{mutation}'
    
    highlight_ids = set()
    patterns_to_check = [p for p in MEANINGFUL_PATTERNS if p != 'Complex']

    for pattern in patterns_to_check:
        df_pat = df[df[pattern_col] == pattern].copy()
        if len(df_pat) == 0: continue
        
        # 1. Filter irrelevant
        df_pat = df_pat[df_pat['Description'].apply(is_relevant_for_highlight)]
        if len(df_pat) == 0: continue
        
        # 2. Prioritize
        df_pat['_priority'] = df_pat.apply(lambda r: get_highlight_priority(assign_semantic_category(r)), axis=1)
        
        # 3. Sort by Priority (asc) then Magnitude Change (desc)
        df_pat['abs_change'] = (df_pat[late_col] - df_pat[early_col]).abs()
        df_pat = df_pat.sort_values(['_priority', 'abs_change'], ascending=[True, False])
        
        # 4. Select top candidates
        # Take up to max_per_pattern from high priority
        # If not enough, fill with others up to limit (capped by MAX_OTHER...)
        
        # Simple selection: just take top N sorted by priority/magnitude
        # This implicitly favors priority categories
        selected = df_pat.head(max_per_pattern)['pathway_id'].tolist()
        highlight_ids.update(selected)
            
    return highlight_ids

# =============================================================================
# GEOMETRY HELPERS
# =============================================================================

def get_bezier_curve(p0, p1, p2, n=50):
    """Quadratic Bezier curve points."""
    t = np.linspace(0, 1, n)
    x0, y0 = p0; x1, y1 = p1; x2, y2 = p2
    x = (1 - t)**2 * x0 + 2 * (1 - t) * t * x1 + t**2 * x2
    y = (1 - t)**2 * y0 + 2 * (1 - t) * t * y1 + t**2 * y2
    return x, y

# =============================================================================
# PLOTTING CORE
# =============================================================================

def create_bump_chart(df, mutation, ax, config, weight_cats=None, custom_highlights=None):
    """
    Core plotting function using Config and Styler.
    """
    styler = BumpChartStyler(config)
    pattern_col = f'Pattern_{mutation}'
    nes_cols = [f'NES_Early_{mutation}', f'NES_Late_{mutation}']
    
    # Pre-filter data
    df_plot = df[df[nes_cols].notna().all(axis=1)].copy()
    if len(df_plot) == 0:
        ax.text(0.5, 0.5, "No Data", ha='center')
        return

    # Drop duplicates for this contrast
    df_plot = df_plot.drop_duplicates(subset=['pathway_id'])

    # Calculate Ranks if needed
    if config.y_type == 'rank':
        for col in nes_cols:
            rank_col = col.replace('NES_', 'Rank_')
            # ascending=False -> Max NES is Rank 1 (Top)
            df_plot[rank_col] = df_plot[col].rank(ascending=False, method='first')

    # Determine Highlights
    highlight_ids = set()
    if config.show_highlights:
        if custom_highlights:
            highlight_ids = set(custom_highlights).intersection(df_plot['pathway_id'].unique())
        else:
            highlight_ids = identify_highlight_pathways(df_plot, mutation)

    # Determine Drawing Order
    # Sort by zorder (weights) if available, then by frequency (Rare on top)
    # Default strategy: Frequency descending (Most frequent -> Bottom)
    # If weights exist: Rare (High zorder) -> Top
    counts = df_plot[pattern_col].value_counts()
    
    # Create list of (pattern, zorder) to sort
    unique_patterns = df_plot[pattern_col].unique()
    pattern_sort_keys = []
    for p in unique_patterns:
        cat = weight_cats.get(p) if weight_cats else None
        z = styler.get_style(p, cat)['zorder']
        # Secondary sort: frequency (more frequent = lower in same zorder group)
        freq = counts.get(p, 0)
        pattern_sort_keys.append((p, z, -freq)) # z asc, freq desc
        
    pattern_sort_keys.sort(key=lambda x: (x[1], x[2]))
    pattern_order = [x[0] for x in pattern_sort_keys]

    # --- Plotting Lines ---
    
    # Store label info
    labels_to_add = []
    
    for pattern in pattern_order:
        df_pat = df_plot[df_plot[pattern_col] == pattern]
        color = styler.get_color(pattern)
        w_cat = weight_cats.get(pattern) if weight_cats else None
        
        # Split into background and highlight
        df_hi = df_pat[df_pat['pathway_id'].isin(highlight_ids)]
        df_bg = df_pat[~df_pat['pathway_id'].isin(highlight_ids)]
        
        # Plot Background
        if len(df_bg) > 0:
            style = styler.get_style(pattern, w_cat, is_highlighted=False)
            _plot_lines(ax, df_bg, nes_cols, mutation, config, style, color, len(df_plot))
            
        # Plot Highlights
        if len(df_hi) > 0:
            style = styler.get_style(pattern, w_cat, is_highlighted=True)
            _plot_lines(ax, df_hi, nes_cols, mutation, config, style, color, len(df_plot))
            
            # Collect label data
            for _, row in df_hi.iterrows():
                y_end = row[nes_cols[1]] if config.y_type == 'nes' else row[nes_cols[1].replace('NES','Rank')]
                desc = str(row.get('Description', row.get('pathway_id')))
                trunc = config.label_truncation
                lbl = desc[:trunc] + '...' if len(desc) > trunc else desc
                labels_to_add.append((lbl, y_end, color))

    # --- Plotting Labels ---
    if labels_to_add and config.show_highlights:
        _add_labels(ax, labels_to_add, config)

    # --- Decoration ---
    _add_decorations(ax, df_plot, nes_cols, mutation, config, len(highlight_ids))

def _plot_lines(ax, df, cols, mutation, config, style, color, max_rank=100):
    """Helper to plot a batch of lines."""
    x_pos = [0, 1]
    
    for _, row in df.iterrows():
        # Get Y values
        if config.y_type == 'nes':
            y = [row[c] for c in cols]
        else:
            # Fallback for rank - assumes cols exist
            y = [row.get(c.replace('NES_', 'Rank_'), 0) for c in cols]

        # Curve Logic
        use_curve = False
        if config.show_curves:
             if row.get(f'Sig_TrajDev_{mutation}', False):
                use_curve = True
                
        if use_curve:
            traj_nes = row.get(f'NES_TrajDev_{mutation}', 0)
            if pd.isna(traj_nes): traj_nes = 0
            
            y_mid = (y[0] + y[1]) / 2
            
            if config.y_type == 'rank':
                # Rank Logic: 
                # 1. Scale: NES ~2 should be ~10-20% of chart height.
                #    Chart height is max_rank. Scale factor approx max_rank / 10.
                # 2. Direction: Positive NES (Upreg) = Visual Bump Up = Lower Rank Value.
                #    So we SUBTRACT from y_mid.
                scale = max_rank / 10.0
                offset = -1 * traj_nes * scale * config.curve_strength
            else:
                # NES Logic
                offset = traj_nes * config.curve_strength
                
            y_control = y_mid + offset
            bx, by = get_bezier_curve((0, y[0]), (0.5, y_control), (1, y[1]))
            ax.plot(bx, by, color=color, **style, solid_capstyle='round')
        else:
            ax.plot(x_pos, y, color=color, **style, solid_capstyle='round')

def _add_labels(ax, labels, config):
    """Add labels with smart placement."""
    texts = []
    y_origs = []
    
    x_label = config.label_x_pos
    
    # Create text objects
    for text, y, color in labels:
        t = ax.text(x_label, y, text, fontsize=7, color=color, 
                    ha='left', va='center', fontweight='bold', zorder=20)
        texts.append(t)
        y_origs.append(y)
        
    # Adjust positions
    # Use adjust_text to spread vertically
    # removing x argument allows text to find space naturally, prevents repulsion to the left
    adjust_text(
        texts,
        y=y_origs,
        ax=ax,
        only_move={'points': 'y', 'text': 'y'}, # primarily adjust Y
        force_text=(0.1, 0.5), # Allow slight X nudge if needed
        lim=100
    )
    
    # Draw connectors
    for t, y_start in zip(texts, y_origs):
        pos = t.get_position()
        y_text = pos[1]
        x_text = pos[0]
        
        # Safety: If text moved too far left (overlapping chart), push it back
        # This fixes the "text at wrong end" issue if repulsion happened
        if x_text < 1.05:
             x_text = max(x_text, x_label)
             t.set_position((x_text, y_text))

        # Draw line from Chart(x=1.0) to Label(x=x_text)
        ax.annotate(
            "", 
            xy=(1.0, y_start), 
            xytext=(x_text - 0.02, y_text),
            arrowprops=dict(arrowstyle="-", color='gray', alpha=0.4, lw=0.5, shrinkA=0, shrinkB=0),
            zorder=15
        )

def _add_decorations(ax, df, cols, mutation, config, n_highlights):
    """Add grid, axes, titles."""
    # Axes
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['Early\n(D35)', 'Late\n(D65)'], fontsize=10)
    
    if config.y_type == 'rank':
        ax.invert_yaxis()
        ax.set_ylabel('Rank')
    else:
        # Dynamic Limits
        vals = df[cols].values.flatten()
        vals = vals[~np.isnan(vals)]
        if len(vals) > 0:
            mx = max(4, np.max(vals) * 1.1)
            mn = min(-4, np.min(vals) * 1.1)
            ax.set_ylim(mn, mx)
        # Reference lines
        ax.axhline(0, c='black', lw=0.5, alpha=0.3)
        ax.axhline(NES_EFFECT, c='gray', ls='--', lw=0.5, alpha=0.3)
        ax.axhline(-NES_EFFECT, c='gray', ls='--', lw=0.5, alpha=0.3)
    
    # Title
    hl_text = f", {n_highlights} labeled" if n_highlights else ""
    ax.set_title(f"{mutation}\n(n={len(df)}{hl_text})", 
                 fontweight='bold', color={'G32A':'#0072B2', 'R403C':'#D55E00'}.get(mutation, 'black'))
    
    ax.grid(axis='y', alpha=0.2)
    
    # X Limits
    if config.show_highlights:
        ax.set_xlim(-0.1, config.label_x_pos + 1.0) # Space for labels
    else:
        ax.set_xlim(-0.1, 1.1)

# =============================================================================
# LEGEND HELPERS
# =============================================================================

def create_unified_legend(ax, df, weight_cats=None):
    """Create a unified legend."""
    from matplotlib.lines import Line2D
    
    # Gather Patterns present
    patterns = set()
    for m in ['G32A', 'R403C']:
        if f'Pattern_{m}' in df.columns:
            patterns.update(df[f'Pattern_{m}'].dropna().unique())
            
    if not patterns:
        ax.axis('off'); return

    # Sort by meaningful order (Frequency/Category)
    # If weights provided, sort by category (Dominant -> Rare)
    # Map category to sort index
    cat_order = {'dominant': 0, 'common': 1, 'uncommon': 2, 'rare': 3}
    
    def sort_key(p):
        cat = weight_cats.get(p, 'rare') if weight_cats else 'rare'
        return (cat_order.get(cat, 4), p)
        
    sorted_patterns = sorted(patterns, key=sort_key)
    
    handles = []
    labels = []
    # Dummy config for legend style
    styler = BumpChartStyler(BumpChartConfig(mode='weighted')) 
    
    for p in sorted_patterns:
        c = get_pattern_colors().get(p, 'gray')
        cat = weight_cats.get(p, 'common') if weight_cats else None
        style = styler.get_style(p, cat)
        
        # Create line handle
        h = Line2D([0], [0], color=c, lw=style['linewidth']*2) # Thicker for legend
        handles.append(h)
        
        # Label with count/freq?
        # For simplicity, just pattern name + category hint
        lbl = p
        if weight_cats and p in weight_cats:
            lbl += f" ({weight_cats[p]})"
        labels.append(lbl)
        
    ax.legend(handles, labels, loc='center left', fontsize=8, title='Patterns', framealpha=1.0)
    ax.axis('off')