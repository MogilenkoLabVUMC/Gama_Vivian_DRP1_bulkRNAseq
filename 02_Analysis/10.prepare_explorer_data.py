#!/usr/bin/env python3
"""
DRP1 Pathway Data Explorer Generator (v2)

This script creates a standalone, self-contained HTML file that allows
collaborators to explore GSEA pathway analysis and GSVA trajectory results.

Changes in v2:
- GSEA view: Vertical layout (dotplot on top, table below)
- GSVA view: Now supports ~12K pathways with search/filter
- Removed broken "All Modules Comparison" heatmap
- Consistent filter UI across both views

Usage:
    python3 02_Analysis/10.prepare_explorer_data.py

Output:
    03_Results/02_Analysis/Explorer/DRP1_Pathway_Explorer.html
"""

import os
import json
import pandas as pd
import numpy as np
from pathlib import Path
import urllib.request

# ============================================================================
# Configuration
# ============================================================================

PROJECT_ROOT = Path(__file__).parent.parent
DATA_DIR = PROJECT_ROOT / "03_Results" / "02_Analysis"
OUTPUT_DIR = DATA_DIR / "Explorer"

PATHWAY_FILE = DATA_DIR / "master_gsea_table.csv"
GSVA_FILE = DATA_DIR / "master_gsva_all_table.csv"  # Comprehensive GSVA
OUTPUT_HTML = OUTPUT_DIR / "DRP1_Pathway_Explorer.html"

# CDN URLs for libraries (will be embedded)
TABULATOR_JS_URL = "https://unpkg.com/tabulator-tables@5.5.0/dist/js/tabulator.min.js"
TABULATOR_CSS_URL = "https://unpkg.com/tabulator-tables@5.5.0/dist/css/tabulator.min.css"
PLOTLY_JS_URL = "https://cdn.plot.ly/plotly-basic-2.27.0.min.js"

# Colors (colorblind-safe)
COLORS = {
    "genotypes": {
        "Ctrl": "#808080",
        "G32A": "#0072B2",
        "R403C": "#D55E00"
    },
    "nes": {
        "negative": "#2166AC",
        "neutral": "#F7F7F7",
        "positive": "#B35806"
    },
    "patterns": {
        "Compensation": "#009E73",
        "Progressive": "#D55E00",
        "Natural_worsening": "#E69F00",
        "Natural_improvement": "#56B4E9",
        "Late_onset": "#CC79A7",
        "Transient": "#0072B2",
        "Persistent": "#999999",
        "Complex": "#F0E442"
    },
    "databases": {
        "MitoCarta": "#E41A1C",
        "SynGO": "#377EB8",
        "gobp": "#4DAF4A",
        "gocc": "#984EA3",
        "gomf": "#FF7F00",
        "hallmark": "#FFFF33",
        "kegg": "#A65628",
        "reactome": "#F781BF",
        "wiki": "#999999",
        "tf": "#66C2A5",
        "cgp": "#FC8D62",
        "canon": "#8DA0CB"
    }
}

# ============================================================================
# Data Loading and Preprocessing
# ============================================================================

def load_pathway_data():
    """Load and preprocess the master pathway table."""
    print("Loading pathway data...")
    df = pd.read_csv(PATHWAY_FILE)
    print(f"  Loaded {len(df)} rows")

    # Get unique pathways (one row per pathway instead of per contrast)
    pathways_wide = df.drop_duplicates(subset=['pathway_id']).copy()

    columns_to_keep = [
        'pathway_id', 'database', 'Description', 'setSize',
        'Pattern_G32A', 'Pattern_R403C', 'change_consistency',
        'NES_Early_G32A', 'NES_Early_R403C',
        'NES_TrajDev_G32A', 'NES_TrajDev_R403C',
        'NES_Late_G32A', 'NES_Late_R403C',
        'ever_significant', 'ever_significant_trajectory'
    ]

    pathways_wide = pathways_wide[columns_to_keep].copy()

    # Round numeric columns
    numeric_cols = pathways_wide.select_dtypes(include=[np.number]).columns
    pathways_wide[numeric_cols] = pathways_wide[numeric_cols].round(3)

    # Get p-values for each stage
    for contrast in ['G32A_vs_Ctrl_D35', 'R403C_vs_Ctrl_D35',
                     'Maturation_G32A_specific', 'Maturation_R403C_specific',
                     'G32A_vs_Ctrl_D65', 'R403C_vs_Ctrl_D65']:
        contrast_df = df[df['contrast'] == contrast][['pathway_id', 'p.adjust']].copy()
        stage_map = {
            'G32A_vs_Ctrl_D35': 'Early_G32A',
            'R403C_vs_Ctrl_D35': 'Early_R403C',
            'Maturation_G32A_specific': 'TrajDev_G32A',
            'Maturation_R403C_specific': 'TrajDev_R403C',
            'G32A_vs_Ctrl_D65': 'Late_G32A',
            'R403C_vs_Ctrl_D65': 'Late_R403C'
        }
        col_name = f'padj_{stage_map[contrast]}'
        contrast_df = contrast_df.rename(columns={'p.adjust': col_name})
        pathways_wide = pathways_wide.merge(contrast_df, on='pathway_id', how='left')

    padj_cols = [c for c in pathways_wide.columns if c.startswith('padj_')]
    pathways_wide[padj_cols] = pathways_wide[padj_cols].round(4)

    pathways_wide['ever_significant'] = pathways_wide['ever_significant'].map(
        {True: True, False: False, 'True': True, 'False': False})
    pathways_wide['ever_significant_trajectory'] = pathways_wide['ever_significant_trajectory'].map(
        {True: True, False: False, 'True': True, 'False': False})

    pathways_wide = pathways_wide.replace({np.nan: None})

    print(f"  Preprocessed to {len(pathways_wide)} unique pathways")
    return pathways_wide


def load_gsva_data():
    """Load the comprehensive GSVA trajectory data (~12K pathways)."""
    print("Loading comprehensive GSVA data...")
    df = pd.read_csv(GSVA_FILE)
    print(f"  Loaded {len(df)} rows")

    # Convert to wide format for the explorer (one row per pathway)
    # Keep the long format data for trajectory plotting
    gsva_long = df.copy()

    # Create wide format with GSVA scores by condition
    pivot_cols = ['pathway_id', 'database', 'pathway_name']

    # Pivot Mean_GSVA by Contrast_Label
    gsva_wide = gsva_long.pivot_table(
        index=['pathway_id', 'database', 'pathway_name'],
        columns='Contrast_Label',
        values='Mean_GSVA',
        aggfunc='first'
    ).reset_index()

    # Also get divergence and significance
    div_pivot = gsva_long.pivot_table(
        index='pathway_id',
        columns='Contrast_Label',
        values='Divergence_vs_Ctrl',
        aggfunc='first'
    ).reset_index()
    div_pivot.columns = ['pathway_id'] + [f'Div_{c}' for c in div_pivot.columns[1:]]

    sig_pivot = gsva_long.pivot_table(
        index='pathway_id',
        columns='Contrast_Label',
        values='significant',
        aggfunc='first'
    ).reset_index()
    sig_pivot.columns = ['pathway_id'] + [f'Sig_{c}' for c in sig_pivot.columns[1:]]

    pval_pivot = gsva_long.pivot_table(
        index='pathway_id',
        columns='Contrast_Label',
        values='p_adjusted',
        aggfunc='first'
    ).reset_index()
    pval_pivot.columns = ['pathway_id'] + [f'padj_{c}' for c in pval_pivot.columns[1:]]

    # Merge all pivots
    gsva_wide = gsva_wide.merge(div_pivot, on='pathway_id', how='left')
    gsva_wide = gsva_wide.merge(sig_pivot, on='pathway_id', how='left')
    gsva_wide = gsva_wide.merge(pval_pivot, on='pathway_id', how='left')

    # Add ever_significant flag
    sig_cols = [c for c in gsva_wide.columns if c.startswith('Sig_')]
    gsva_wide['ever_significant'] = gsva_wide[sig_cols].any(axis=1)

    # Round numeric columns
    numeric_cols = gsva_wide.select_dtypes(include=[np.number]).columns
    gsva_wide[numeric_cols] = gsva_wide[numeric_cols].round(3)

    gsva_wide = gsva_wide.replace({np.nan: None})

    print(f"  Preprocessed to {len(gsva_wide)} unique GSVA pathways")

    # Also return the long format for trajectory plotting
    gsva_long_clean = gsva_long[['pathway_id', 'database', 'pathway_name',
                                  'Genotype', 'Day', 'Contrast_Label',
                                  'Mean_GSVA', 'SD_GSVA', 'Divergence_vs_Ctrl',
                                  'p_value', 'p_adjusted', 'significant']].copy()
    numeric_cols = gsva_long_clean.select_dtypes(include=[np.number]).columns
    gsva_long_clean[numeric_cols] = gsva_long_clean[numeric_cols].round(3)
    gsva_long_clean = gsva_long_clean.replace({np.nan: None})

    return gsva_wide, gsva_long_clean


def get_metadata(pathway_df, gsva_wide_df):
    """Extract metadata for filters."""
    metadata = {
        "databases": sorted(pathway_df['database'].unique().tolist()),
        "patterns_g32a": sorted([p for p in pathway_df['Pattern_G32A'].unique() if p]),
        "patterns_r403c": sorted([p for p in pathway_df['Pattern_R403C'].unique() if p]),
        "consistency_types": sorted([c for c in pathway_df['change_consistency'].unique() if c]),
        "gsva_databases": sorted(gsva_wide_df['database'].unique().tolist()),
        "total_pathways": len(pathway_df),
        "significant_pathways": int(pathway_df['ever_significant'].sum()),
        "total_gsva_pathways": len(gsva_wide_df),
        "significant_gsva_pathways": int(gsva_wide_df['ever_significant'].sum()),
    }
    return metadata


# ============================================================================
# Library Fetching
# ============================================================================

def fetch_library(url, name):
    """Fetch a JS/CSS library from CDN."""
    print(f"  Fetching {name}...")
    try:
        with urllib.request.urlopen(url, timeout=30) as response:
            content = response.read().decode('utf-8')
            print(f"    Downloaded {len(content)} bytes")
            return content
    except Exception as e:
        print(f"    Warning: Could not fetch {name}: {e}")
        print(f"    Using CDN link instead (requires internet)")
        return None


def get_libraries():
    """Fetch all required libraries."""
    print("Fetching libraries...")
    libs = {
        'tabulator_js': fetch_library(TABULATOR_JS_URL, "Tabulator JS"),
        'tabulator_css': fetch_library(TABULATOR_CSS_URL, "Tabulator CSS"),
        'plotly_js': fetch_library(PLOTLY_JS_URL, "Plotly JS"),
    }
    return libs


# ============================================================================
# HTML Template
# ============================================================================

def generate_html(pathway_data, gsva_wide, gsva_long, metadata, libraries):
    """Generate the complete HTML file."""

    pathway_json = json.dumps(pathway_data.to_dict('records'), default=str)
    gsva_wide_json = json.dumps(gsva_wide.to_dict('records'), default=str)
    gsva_long_json = json.dumps(gsva_long.to_dict('records'), default=str)
    metadata_json = json.dumps(metadata)
    colors_json = json.dumps(COLORS)

    tabulator_js = f"<script>{libraries['tabulator_js']}</script>" if libraries['tabulator_js'] else f'<script src="{TABULATOR_JS_URL}"></script>'
    tabulator_css = f"<style>{libraries['tabulator_css']}</style>" if libraries['tabulator_css'] else f'<link href="{TABULATOR_CSS_URL}" rel="stylesheet">'
    plotly_js = f"<script>{libraries['plotly_js']}</script>" if libraries['plotly_js'] else f'<script src="{PLOTLY_JS_URL}"></script>'

    html = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>DRP1 Pathway Data Explorer</title>
    {tabulator_css}
    {plotly_js}
    {tabulator_js}
    <style>
        :root {{
            --primary-color: #2166AC;
            --secondary-color: #B35806;
            --bg-color: #f5f5f5;
            --card-bg: #ffffff;
            --text-color: #333333;
            --border-color: #ddd;
        }}

        * {{ box-sizing: border-box; margin: 0; padding: 0; }}

        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            background: var(--bg-color);
            color: var(--text-color);
            line-height: 1.6;
        }}

        .header {{
            background: linear-gradient(135deg, var(--primary-color), var(--secondary-color));
            color: white;
            padding: 1rem 2rem;
            display: flex;
            justify-content: space-between;
            align-items: center;
        }}

        .header h1 {{ font-size: 1.5rem; font-weight: 600; }}
        .header-stats {{ font-size: 0.9rem; opacity: 0.9; }}

        .nav-tabs {{
            display: flex;
            background: var(--card-bg);
            border-bottom: 2px solid var(--border-color);
        }}

        .nav-tab {{
            padding: 1rem 2rem;
            cursor: pointer;
            border: none;
            background: none;
            font-size: 1rem;
            border-bottom: 3px solid transparent;
            transition: all 0.2s;
        }}

        .nav-tab:hover {{ background: var(--bg-color); }}
        .nav-tab.active {{
            border-bottom-color: var(--primary-color);
            color: var(--primary-color);
            font-weight: 600;
        }}

        .filter-bar {{
            background: var(--card-bg);
            padding: 1rem 2rem;
            display: flex;
            flex-wrap: wrap;
            gap: 1rem;
            align-items: center;
            border-bottom: 1px solid var(--border-color);
        }}

        .filter-group {{
            display: flex;
            flex-direction: column;
            gap: 0.25rem;
        }}

        .filter-group label {{
            font-size: 0.75rem;
            color: #666;
            text-transform: uppercase;
        }}

        .filter-group select, .filter-group input {{
            padding: 0.5rem 1rem;
            border: 1px solid var(--border-color);
            border-radius: 4px;
            font-size: 0.9rem;
            min-width: 150px;
        }}

        .search-input {{ flex: 1; min-width: 250px; }}

        .checkbox-group {{
            display: flex;
            align-items: center;
            gap: 0.5rem;
        }}

        .checkbox-group input {{ min-width: auto; width: 18px; height: 18px; }}

        /* GSEA View - Vertical Layout */
        .gsea-content {{
            display: flex;
            flex-direction: column;
            padding: 1.5rem 2rem;
            gap: 2rem;
        }}

        .gsea-viz-panel {{
            background: var(--card-bg);
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            padding: 1.5rem;
            height: 600px;
        }}

        .gsea-detail-panel {{
            background: var(--card-bg);
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            padding: 1.5rem;
            display: none;
        }}

        .gsea-table-panel {{
            background: var(--card-bg);
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            padding: 1.5rem;
            height: calc(100vh - 800px);
            min-height: 450px;
        }}

        /* GSVA View - 70/30 split with table below */
        .gsva-content {{
            display: flex;
            flex-direction: column;
            padding: 1.5rem 2rem;
            gap: 2rem;
        }}

        .gsva-top-row {{
            display: flex;
            gap: 1.5rem;
            height: 500px;
        }}

        .gsva-viz-panel {{
            flex: 0 0 70%;
            background: var(--card-bg);
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            padding: 1.5rem;
        }}

        .gsva-stats-panel {{
            flex: 0 0 28%;
            background: var(--card-bg);
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            padding: 1.5rem;
            overflow: auto;
        }}

        .gsva-table-panel {{
            background: var(--card-bg);
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            padding: 1.5rem;
            height: calc(100vh - 720px);
            min-height: 450px;
        }}

        .panel-header {{
            display: flex;
            justify-content: space-between;
            align-items: center;
            margin-bottom: 0.5rem;
            padding-bottom: 0.5rem;
            border-bottom: 1px solid var(--border-color);
        }}

        .panel-title {{ font-size: 1rem; font-weight: 600; }}
        .panel-actions {{ display: flex; gap: 0.5rem; }}

        .btn {{
            padding: 0.4rem 0.8rem;
            border: 1px solid var(--border-color);
            border-radius: 4px;
            background: white;
            cursor: pointer;
            font-size: 0.8rem;
        }}

        .btn:hover {{ background: var(--bg-color); }}
        .btn-primary {{
            background: var(--primary-color);
            color: white;
            border-color: var(--primary-color);
        }}

        .plot-container {{ flex: 1; min-height: 350px; }}

        .view {{ display: none; }}
        .view.active {{ display: block; }}

        .status-bar {{
            background: var(--card-bg);
            padding: 0.5rem 2rem;
            font-size: 0.85rem;
            color: #666;
            border-top: 1px solid var(--border-color);
            display: flex;
            justify-content: space-between;
        }}

        .viz-controls {{
            display: flex;
            gap: 0.5rem;
            margin-bottom: 0.5rem;
            flex-wrap: wrap;
        }}

        .viz-controls select {{
            padding: 0.3rem 0.6rem;
            font-size: 0.8rem;
            border: 1px solid var(--border-color);
            border-radius: 4px;
        }}

        /* Pattern badges */
        .pattern-badge {{
            display: inline-block;
            padding: 0.2rem 0.5rem;
            border-radius: 3px;
            font-size: 0.75rem;
            font-weight: 500;
        }}

        .pattern-Compensation {{ background: #009E73; color: white; }}
        .pattern-Progressive {{ background: #D55E00; color: white; }}
        .pattern-Natural_worsening {{ background: #E69F00; color: white; }}
        .pattern-Natural_improvement {{ background: #56B4E9; color: black; }}
        .pattern-Late_onset {{ background: #CC79A7; color: white; }}
        .pattern-Transient {{ background: #0072B2; color: white; }}
        .pattern-Persistent {{ background: #999999; color: white; }}
        .pattern-Complex {{ background: #F0E442; color: black; }}

        .sig-indicator {{
            width: 10px;
            height: 10px;
            border-radius: 50%;
            display: inline-block;
        }}
        .sig-yes {{ background: #009E73; }}
        .sig-no {{ background: #ccc; }}

        .nes-bar {{ display: flex; align-items: center; gap: 0.5rem; }}
        .nes-bar-container {{
            width: 60px;
            height: 8px;
            background: #eee;
            border-radius: 4px;
            overflow: hidden;
            position: relative;
        }}
        .nes-bar-fill {{ height: 100%; position: absolute; }}
        .nes-positive {{ background: var(--secondary-color); right: 50%; }}
        .nes-negative {{ background: var(--primary-color); left: 50%; }}

        .tabulator {{ font-size: 0.85rem; }}
        .tabulator .tabulator-header {{ background: var(--bg-color); }}
        .tabulator-row.tabulator-selected {{ background: #e3f2fd !important; }}
        .tabulator-row:hover {{ background: #f5f5f5; }}

        /* Stats mini-table */
        .stats-table {{
            width: 100%;
            border-collapse: collapse;
            font-size: 0.85rem;
        }}
        .stats-table th, .stats-table td {{
            padding: 0.4rem;
            text-align: left;
            border-bottom: 1px solid var(--border-color);
        }}
        .stats-table th {{ background: var(--bg-color); font-weight: 600; }}
    </style>
</head>
<body>
    <header class="header">
        <h1>DRP1 Pathway Data Explorer</h1>
        <div class="header-stats">
            <span id="total-pathways"></span> pathways |
            <span id="significant-count"></span> significant
        </div>
    </header>

    <nav class="nav-tabs">
        <button class="nav-tab active" data-view="gsea">GSEA Pathways</button>
        <button class="nav-tab" data-view="gsva">GSVA Trajectories</button>
    </nav>

    <!-- GSEA View (Vertical Layout) -->
    <div id="gsea-view" class="view active">
        <div class="filter-bar">
            <div class="filter-group">
                <label>Database</label>
                <select id="filter-database"><option value="">All Databases</option></select>
            </div>
            <div class="filter-group">
                <label>G32A Pattern</label>
                <select id="filter-pattern-g32a"><option value="">All Patterns</option></select>
            </div>
            <div class="filter-group">
                <label>R403C Pattern</label>
                <select id="filter-pattern-r403c"><option value="">All Patterns</option></select>
            </div>
            <div class="filter-group checkbox-group">
                <label style="display: flex; align-items: center; gap: 0.5rem;">
                    <input type="checkbox" id="filter-significant"> Significant Only
                </label>
            </div>
            <div class="filter-group search-input">
                <label>Search Pathways</label>
                <input type="text" id="search-pathway" placeholder="e.g., calcium, ribosome, mitochondr...">
            </div>
        </div>

        <div class="gsea-content">
            <div class="gsea-viz-panel">
                <div class="panel-header">
                    <span class="panel-title">GSEA Dotplot</span>
                    <div class="panel-actions">
                        <div class="viz-controls">
                            <select id="dotplot-stage">
                                <option value="Early_G32A">Early G32A</option>
                                <option value="Early_R403C">Early R403C</option>
                                <option value="TrajDev_G32A">TrajDev G32A</option>
                                <option value="TrajDev_R403C">TrajDev R403C</option>
                                <option value="Late_G32A">Late G32A</option>
                                <option value="Late_R403C">Late R403C</option>
                            </select>
                            <select id="dotplot-color">
                                <option value="database">Color by Database</option>
                                <option value="pattern">Color by Pattern</option>
                            </select>
                            <select id="dotplot-n">
                                <option value="20">Top 20</option>
                                <option value="30" selected>Top 30</option>
                                <option value="50">Top 50</option>
                            </select>
                        </div>
                        <button class="btn" onclick="exportPlot('dotplot', 'svg')">SVG</button>
                        <button class="btn" onclick="exportPlot('dotplot', 'png')">PNG</button>
                    </div>
                </div>
                <div id="dotplot" class="plot-container"></div>
            </div>

            <div id="gsea-detail-panel" class="gsea-detail-panel">
                <div class="panel-header">
                    <span class="panel-title" id="detail-pathway-name">Selected Pathway</span>
                    <button class="btn" onclick="hideDetailPanel()">Close</button>
                </div>
                <div id="trajectory-plot" style="height: 220px;"></div>
            </div>

            <div class="gsea-table-panel">
                <div class="panel-header">
                    <span class="panel-title">Pathway Table</span>
                    <button class="btn btn-primary" onclick="exportTableCSV()">Export CSV</button>
                </div>
                <div id="pathway-table" style="height: calc(100% - 40px);"></div>
            </div>
        </div>
    </div>

    <!-- GSVA View (With Search/Filter) -->
    <div id="gsva-view" class="view">
        <div class="filter-bar">
            <div class="filter-group">
                <label>Database</label>
                <select id="gsva-filter-database"><option value="">All Databases</option></select>
            </div>
            <div class="filter-group checkbox-group">
                <label style="display: flex; align-items: center; gap: 0.5rem;">
                    <input type="checkbox" id="gsva-filter-significant"> Significant Only
                </label>
            </div>
            <div class="filter-group search-input">
                <label>Search Pathways</label>
                <input type="text" id="gsva-search-pathway" placeholder="e.g., ribosome, OXPHOS, translation...">
            </div>
        </div>

        <div class="gsva-content">
            <div class="gsva-top-row">
                <div class="gsva-viz-panel">
                    <div class="panel-header">
                        <span class="panel-title">GSVA Trajectory</span>
                        <div class="panel-actions">
                            <button class="btn" onclick="exportPlot('gsva-trajectory', 'svg')">SVG</button>
                            <button class="btn" onclick="exportPlot('gsva-trajectory', 'png')">PNG</button>
                        </div>
                    </div>
                    <div id="gsva-trajectory-plot" class="plot-container"></div>
                </div>

                <div class="gsva-stats-panel">
                    <div class="panel-header">
                        <span class="panel-title">Statistics</span>
                    </div>
                    <div id="gsva-stats-content"></div>
                </div>
            </div>

            <div class="gsva-table-panel">
                <div class="panel-header">
                    <span class="panel-title">Pathway Selector</span>
                    <button class="btn btn-primary" onclick="exportGSVACSV()">Export CSV</button>
                </div>
                <div id="gsva-table" style="height: calc(100% - 40px);"></div>
            </div>
        </div>
    </div>

    <footer class="status-bar">
        <span id="status-message">Ready</span>
        <span>DRP1 Pathway Analysis | Generated: <span id="gen-date"></span></span>
    </footer>

    <script>
        // ====================================================================
        // Data
        // ====================================================================
        const PATHWAY_DATA = {pathway_json};
        const GSVA_WIDE = {gsva_wide_json};
        const GSVA_LONG = {gsva_long_json};
        const METADATA = {metadata_json};
        const COLORS = {colors_json};

        // ====================================================================
        // State
        // ====================================================================
        let pathwayTable = null;
        let gsvaTable = null;
        let filteredPathways = [...PATHWAY_DATA];
        let filteredGSVA = [...GSVA_WIDE];
        let selectedGSVAPathway = GSVA_WIDE[0]?.pathway_id || null;
        let currentView = 'gsea';

        // ====================================================================
        // Initialization
        // ====================================================================
        document.addEventListener('DOMContentLoaded', () => {{
            document.getElementById('gen-date').textContent = new Date().toLocaleDateString();
            document.getElementById('total-pathways').textContent = METADATA.total_pathways;
            document.getElementById('significant-count').textContent = METADATA.significant_pathways;

            initFilters();
            initPathwayTable();
            initGSVATable();
            updateDotplot();
            if (selectedGSVAPathway) updateGSVATrajectory();

            initEventListeners();
        }});

        function initFilters() {{
            // GSEA filters
            const dbSelect = document.getElementById('filter-database');
            METADATA.databases.forEach(db => {{
                const opt = document.createElement('option');
                opt.value = db;
                opt.textContent = db;
                dbSelect.appendChild(opt);
            }});

            const g32aSelect = document.getElementById('filter-pattern-g32a');
            METADATA.patterns_g32a.forEach(p => {{
                const opt = document.createElement('option');
                opt.value = p;
                opt.textContent = p;
                g32aSelect.appendChild(opt);
            }});

            const r403cSelect = document.getElementById('filter-pattern-r403c');
            METADATA.patterns_r403c.forEach(p => {{
                const opt = document.createElement('option');
                opt.value = p;
                opt.textContent = p;
                r403cSelect.appendChild(opt);
            }});

            // GSVA filters
            const gsvaDbSelect = document.getElementById('gsva-filter-database');
            METADATA.gsva_databases.forEach(db => {{
                const opt = document.createElement('option');
                opt.value = db;
                opt.textContent = db;
                gsvaDbSelect.appendChild(opt);
            }});
        }}

        function initEventListeners() {{
            // View tabs
            document.querySelectorAll('.nav-tab').forEach(tab => {{
                tab.addEventListener('click', () => {{
                    document.querySelectorAll('.nav-tab').forEach(t => t.classList.remove('active'));
                    tab.classList.add('active');
                    currentView = tab.dataset.view;
                    document.querySelectorAll('.view').forEach(v => v.classList.remove('active'));
                    document.getElementById(currentView + '-view').classList.add('active');
                    setTimeout(() => {{
                        if (currentView === 'gsva') {{
                            Plotly.relayout('gsva-trajectory-plot', {{}});
                        }} else {{
                            Plotly.relayout('dotplot', {{}});
                        }}
                    }}, 100);
                }});
            }});

            // GSEA filter listeners
            document.getElementById('filter-database').addEventListener('change', applyFilters);
            document.getElementById('filter-pattern-g32a').addEventListener('change', applyFilters);
            document.getElementById('filter-pattern-r403c').addEventListener('change', applyFilters);
            document.getElementById('filter-significant').addEventListener('change', applyFilters);

            let searchTimeout;
            document.getElementById('search-pathway').addEventListener('input', () => {{
                clearTimeout(searchTimeout);
                searchTimeout = setTimeout(() => applyFilters(), 300);
            }});

            document.getElementById('dotplot-stage').addEventListener('change', updateDotplot);
            document.getElementById('dotplot-color').addEventListener('change', updateDotplot);
            document.getElementById('dotplot-n').addEventListener('change', updateDotplot);

            // GSVA filter listeners
            document.getElementById('gsva-filter-database').addEventListener('change', applyGSVAFilters);
            document.getElementById('gsva-filter-significant').addEventListener('change', applyGSVAFilters);

            let gsvaSearchTimeout;
            document.getElementById('gsva-search-pathway').addEventListener('input', () => {{
                clearTimeout(gsvaSearchTimeout);
                gsvaSearchTimeout = setTimeout(() => applyGSVAFilters(), 300);
            }});
        }}

        // ====================================================================
        // GSEA Filtering
        // ====================================================================
        function applyFilters() {{
            const database = document.getElementById('filter-database').value;
            const patternG32A = document.getElementById('filter-pattern-g32a').value;
            const patternR403C = document.getElementById('filter-pattern-r403c').value;
            const sigOnly = document.getElementById('filter-significant').checked;
            const search = document.getElementById('search-pathway').value.toLowerCase();

            filteredPathways = PATHWAY_DATA.filter(p => {{
                if (database && p.database !== database) return false;
                if (patternG32A && p.Pattern_G32A !== patternG32A) return false;
                if (patternR403C && p.Pattern_R403C !== patternR403C) return false;
                if (sigOnly && !p.ever_significant) return false;
                if (search && !p.Description.toLowerCase().includes(search)) return false;
                return true;
            }});

            pathwayTable.setData(filteredPathways);
            updateDotplot();
            updateStatus(`GSEA: ${{filteredPathways.length}} of ${{PATHWAY_DATA.length}} pathways`);
        }}

        // ====================================================================
        // GSVA Filtering
        // ====================================================================
        function applyGSVAFilters() {{
            const database = document.getElementById('gsva-filter-database').value;
            const sigOnly = document.getElementById('gsva-filter-significant').checked;
            const search = document.getElementById('gsva-search-pathway').value.toLowerCase();

            filteredGSVA = GSVA_WIDE.filter(p => {{
                if (database && p.database !== database) return false;
                if (sigOnly && !p.ever_significant) return false;
                if (search && !p.pathway_name.toLowerCase().includes(search)) return false;
                return true;
            }});

            gsvaTable.setData(filteredGSVA);
            updateStatus(`GSVA: ${{filteredGSVA.length}} of ${{GSVA_WIDE.length}} pathways`);

            // Select first pathway if current selection is filtered out
            if (filteredGSVA.length > 0 && !filteredGSVA.find(p => p.pathway_id === selectedGSVAPathway)) {{
                selectGSVAPathway(filteredGSVA[0].pathway_id);
            }}
        }}

        // ====================================================================
        // Pathway Table
        // ====================================================================
        function initPathwayTable() {{
            pathwayTable = new Tabulator("#pathway-table", {{
                data: PATHWAY_DATA,
                height: "100%",
                layout: "fitDataStretch",
                virtualDom: true,
                virtualDomBuffer: 300,
                selectable: 1,
                columns: [
                    {{title: "Pathway", field: "Description", width: 280, formatter: "textarea"}},
                    {{title: "DB", field: "database", width: 90,
                     formatter: cell => `<span style="color: ${{COLORS.databases[cell.getValue()] || '#999'}}; font-weight: 600;">${{cell.getValue()}}</span>`}},
                    {{title: "Size", field: "setSize", width: 60, hozAlign: "center"}},
                    {{title: "G32A", field: "Pattern_G32A", width: 100,
                     formatter: cell => cell.getValue() ? `<span class="pattern-badge pattern-${{cell.getValue()}}">${{cell.getValue().replace('_', ' ')}}</span>` : ''}},
                    {{title: "R403C", field: "Pattern_R403C", width: 100,
                     formatter: cell => cell.getValue() ? `<span class="pattern-badge pattern-${{cell.getValue()}}">${{cell.getValue().replace('_', ' ')}}</span>` : ''}},
                    {{title: "NES E.G32A", field: "NES_Early_G32A", width: 100, hozAlign: "right", formatter: nesFormatter}},
                    {{title: "NES L.G32A", field: "NES_Late_G32A", width: 100, hozAlign: "right", formatter: nesFormatter}},
                    {{title: "NES E.R403C", field: "NES_Early_R403C", width: 100, hozAlign: "right", formatter: nesFormatter}},
                    {{title: "NES L.R403C", field: "NES_Late_R403C", width: 100, hozAlign: "right", formatter: nesFormatter}},
                    {{title: "Sig", field: "ever_significant", width: 50, hozAlign: "center",
                     formatter: cell => `<span class="sig-indicator sig-${{cell.getValue() ? 'yes' : 'no'}}"></span>`}},
                ],
                rowClick: (e, row) => showPathwayDetail(row.getData())
            }});
        }}

        function nesFormatter(cell) {{
            const val = cell.getValue();
            if (val === null || val === undefined) return '';
            const color = val >= 0 ? COLORS.nes.positive : COLORS.nes.negative;
            const width = Math.min(Math.abs(val) / 3 * 50, 50);
            const dir = val >= 0 ? 'right' : 'left';
            return `<div class="nes-bar"><span style="width:40px;text-align:right;">${{val.toFixed(2)}}</span><div class="nes-bar-container"><div class="nes-bar-fill" style="background:${{color}};width:${{width}}%;${{dir}}:50%;"></div></div></div>`;
        }}

        // ====================================================================
        // GSVA Table
        // ====================================================================
        function initGSVATable() {{
            gsvaTable = new Tabulator("#gsva-table", {{
                data: GSVA_WIDE,
                height: "100%",
                layout: "fitDataStretch",
                virtualDom: true,
                virtualDomBuffer: 300,
                selectable: 1,
                columns: [
                    {{title: "Pathway", field: "pathway_name", width: 300, formatter: "textarea"}},
                    {{title: "DB", field: "database", width: 90,
                     formatter: cell => `<span style="color: ${{COLORS.databases[cell.getValue()] || '#999'}}; font-weight: 600;">${{cell.getValue()}}</span>`}},
                    {{title: "Ctrl D35", field: "Ctrl_D35", width: 80, hozAlign: "right",
                     formatter: cell => cell.getValue()?.toFixed(2) || ''}},
                    {{title: "Early G32A", field: "Early_G32A", width: 90, hozAlign: "right",
                     formatter: cell => cell.getValue()?.toFixed(2) || ''}},
                    {{title: "Late G32A", field: "Late_G32A", width: 90, hozAlign: "right",
                     formatter: cell => cell.getValue()?.toFixed(2) || ''}},
                    {{title: "Early R403C", field: "Early_R403C", width: 90, hozAlign: "right",
                     formatter: cell => cell.getValue()?.toFixed(2) || ''}},
                    {{title: "Late R403C", field: "Late_R403C", width: 90, hozAlign: "right",
                     formatter: cell => cell.getValue()?.toFixed(2) || ''}},
                    {{title: "Sig", field: "ever_significant", width: 50, hozAlign: "center",
                     formatter: cell => `<span class="sig-indicator sig-${{cell.getValue() ? 'yes' : 'no'}}"></span>`}},
                ],
                rowClick: (e, row) => selectGSVAPathway(row.getData().pathway_id)
            }});
        }}

        function selectGSVAPathway(pathwayId) {{
            selectedGSVAPathway = pathwayId;
            updateGSVATrajectory();
            updateGSVAStats();
        }}

        // ====================================================================
        // Dotplot
        // ====================================================================
        function updateDotplot() {{
            const stage = document.getElementById('dotplot-stage').value;
            const colorBy = document.getElementById('dotplot-color').value;
            const topN = parseInt(document.getElementById('dotplot-n').value);

            const nesCol = 'NES_' + stage;
            const padjCol = 'padj_' + stage;

            const sorted = [...filteredPathways]
                .filter(p => p[nesCol] !== null)
                .sort((a, b) => Math.abs(b[nesCol]) - Math.abs(a[nesCol]))
                .slice(0, topN);

            if (sorted.length === 0) {{
                Plotly.react('dotplot', [], {{title: 'No data', xaxis: {{title: 'NES'}}, yaxis: {{title: 'Pathway'}}}});
                return;
            }}

            const x = sorted.map(p => p[nesCol]);
            const y = sorted.map(p => truncate(p.Description, 45));
            const sizes = sorted.map(p => {{
                const padj = p[padjCol];
                if (!padj || padj === 0) return 20;
                return Math.min(Math.max(-Math.log10(padj) * 4, 5), 30);
            }});

            let colors;
            if (colorBy === 'database') {{
                colors = sorted.map(p => COLORS.databases[p.database] || '#999');
            }} else {{
                const patternCol = stage.includes('G32A') ? 'Pattern_G32A' : 'Pattern_R403C';
                colors = sorted.map(p => COLORS.patterns[p[patternCol]] || '#999');
            }}

            const trace = {{
                x, y,
                mode: 'markers',
                type: 'scatter',
                marker: {{ size: sizes, color: colors, opacity: 0.8, line: {{ width: 1, color: '#333' }} }},
                text: sorted.map(p => `${{p.Description}}<br>NES: ${{p[nesCol]?.toFixed(3)}}<br>p.adj: ${{p[padjCol]?.toFixed(4) || 'N/A'}}<br>DB: ${{p.database}}`),
                hoverinfo: 'text'
            }};

            const layout = {{
                title: `GSEA Dotplot - ${{stage.replace('_', ' ')}}`,
                xaxis: {{ title: 'NES', zeroline: true, zerolinewidth: 2, zerolinecolor: '#999' }},
                yaxis: {{ title: '', automargin: true }},
                margin: {{ l: 250, r: 20, t: 40, b: 40 }},
                hovermode: 'closest',
                showlegend: false
            }};

            Plotly.react('dotplot', [trace], layout, {{responsive: true}});
        }}

        // ====================================================================
        // Pathway Detail
        // ====================================================================
        function showPathwayDetail(pathway) {{
            document.getElementById('gsea-detail-panel').style.display = 'block';
            document.getElementById('detail-pathway-name').textContent = pathway.Description;

            const stages = ['Early_G32A', 'TrajDev_G32A', 'Late_G32A'];
            const stagesR = ['Early_R403C', 'TrajDev_R403C', 'Late_R403C'];

            const traceG32A = {{
                x: ['Early', 'TrajDev', 'Late'],
                y: stages.map(s => pathway['NES_' + s]),
                mode: 'lines+markers',
                name: 'G32A',
                line: {{ color: COLORS.genotypes.G32A, width: 2 }},
                marker: {{ size: 10 }}
            }};

            const traceR403C = {{
                x: ['Early', 'TrajDev', 'Late'],
                y: stagesR.map(s => pathway['NES_' + s]),
                mode: 'lines+markers',
                name: 'R403C',
                line: {{ color: COLORS.genotypes.R403C, width: 2 }},
                marker: {{ size: 10 }}
            }};

            const layout = {{
                title: truncate(pathway.Description, 60),
                xaxis: {{ title: 'Stage' }},
                yaxis: {{ title: 'NES', zeroline: true }},
                margin: {{ l: 60, r: 20, t: 40, b: 40 }},
                legend: {{ x: 0, y: 1.15, orientation: 'h' }},
                hovermode: 'x unified'
            }};

            Plotly.react('trajectory-plot', [traceG32A, traceR403C], layout, {{responsive: true}});
        }}

        function hideDetailPanel() {{
            document.getElementById('gsea-detail-panel').style.display = 'none';
        }}

        // ====================================================================
        // GSVA Trajectory
        // ====================================================================
        function updateGSVATrajectory() {{
            const pathwayData = GSVA_LONG.filter(d => d.pathway_id === selectedGSVAPathway);
            if (pathwayData.length === 0) return;

            const pathwayName = pathwayData[0]?.pathway_name || selectedGSVAPathway;
            const genotypes = ['Ctrl', 'G32A', 'R403C'];

            const traces = genotypes.map(geno => {{
                const genoData = pathwayData.filter(d => d.Genotype === geno).sort((a, b) => a.Day - b.Day);
                return {{
                    x: genoData.map(d => `D${{d.Day}}`),
                    y: genoData.map(d => d.Mean_GSVA),
                    error_y: {{ type: 'data', array: genoData.map(d => d.SD_GSVA), visible: true }},
                    mode: 'lines+markers',
                    name: geno,
                    line: {{ color: COLORS.genotypes[geno], width: 2 }},
                    marker: {{ size: 10 }}
                }};
            }});

            const layout = {{
                title: truncate(pathwayName, 50),
                xaxis: {{ title: 'Timepoint' }},
                yaxis: {{ title: 'Mean GSVA Score' }},
                margin: {{ l: 60, r: 20, t: 40, b: 40 }},
                legend: {{ x: 0, y: 1.15, orientation: 'h' }},
                hovermode: 'x unified'
            }};

            Plotly.react('gsva-trajectory-plot', traces, layout, {{responsive: true}});
        }}

        function updateGSVAStats() {{
            const pathwayData = GSVA_LONG.filter(d => d.pathway_id === selectedGSVAPathway);
            const statsContent = document.getElementById('gsva-stats-content');

            if (pathwayData.length === 0) {{
                statsContent.innerHTML = '<p>No data</p>';
                return;
            }}

            let html = '<table class="stats-table"><thead><tr><th>Condition</th><th>Mean</th><th>Div</th><th>p.adj</th></tr></thead><tbody>';

            const conditions = ['Ctrl_D35', 'Early_G32A', 'Early_R403C', 'Ctrl_D65', 'Late_G32A', 'Late_R403C'];
            conditions.forEach(cond => {{
                const row = pathwayData.find(d => d.Contrast_Label === cond);
                if (row) {{
                    const pStyle = row.p_adjusted && row.p_adjusted < 0.05 ? 'color: #009E73; font-weight: 600;' : '';
                    html += `<tr>
                        <td>${{cond.replace('_', ' ')}}</td>
                        <td>${{row.Mean_GSVA?.toFixed(3) || '-'}}</td>
                        <td style="color: ${{(row.Divergence_vs_Ctrl || 0) >= 0 ? COLORS.nes.positive : COLORS.nes.negative}}">${{row.Divergence_vs_Ctrl?.toFixed(3) || '-'}}</td>
                        <td style="${{pStyle}}">${{row.p_adjusted?.toFixed(4) || 'N/A'}}</td>
                    </tr>`;
                }}
            }});

            html += '</tbody></table>';
            statsContent.innerHTML = html;
        }}

        // ====================================================================
        // Export
        // ====================================================================
        function exportTableCSV() {{
            downloadCSV(pathwayTable.getData("active"), 'DRP1_pathways_filtered.csv');
        }}

        function exportGSVACSV() {{
            downloadCSV(gsvaTable.getData("active"), 'DRP1_gsva_filtered.csv');
        }}

        function downloadCSV(data, filename) {{
            if (data.length === 0) return;
            const headers = Object.keys(data[0]);
            const csv = [headers.join(','), ...data.map(row => headers.map(h => {{
                const val = row[h];
                if (val === null || val === undefined) return '';
                if (typeof val === 'string' && val.includes(',')) return `"${{val}}"`;
                return val;
            }}).join(','))].join('\\n');
            const blob = new Blob([csv], {{ type: 'text/csv' }});
            const link = document.createElement('a');
            link.href = URL.createObjectURL(blob);
            link.download = filename;
            link.click();
        }}

        function exportPlot(plotId, format) {{
            const id = plotId === 'dotplot' ? 'dotplot' : 'gsva-trajectory-plot';
            Plotly.downloadImage(id, {{ format, width: 1200, height: 800, filename: `DRP1_${{plotId}}` }});
        }}

        // ====================================================================
        // Utilities
        // ====================================================================
        function truncate(str, n) {{
            if (!str) return '';
            return str.length > n ? str.substr(0, n - 1) + '...' : str;
        }}

        function updateStatus(msg) {{
            document.getElementById('status-message').textContent = msg;
        }}
    </script>
</body>
</html>'''

    return html


# ============================================================================
# Main
# ============================================================================

def main():
    print("=" * 60)
    print("DRP1 Pathway Data Explorer Generator (v2)")
    print("=" * 60)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    pathway_data = load_pathway_data()
    gsva_wide, gsva_long = load_gsva_data()
    metadata = get_metadata(pathway_data, gsva_wide)

    libraries = get_libraries()

    print("Generating HTML...")
    html = generate_html(pathway_data, gsva_wide, gsva_long, metadata, libraries)

    print(f"Writing to {OUTPUT_HTML}...")
    with open(OUTPUT_HTML, 'w', encoding='utf-8') as f:
        f.write(html)

    file_size = OUTPUT_HTML.stat().st_size / (1024 * 1024)
    print("=" * 60)
    print(f"SUCCESS! Generated: {OUTPUT_HTML}")
    print(f"File size: {file_size:.2f} MB")
    print(f"GSEA pathways: {metadata['total_pathways']} ({metadata['significant_pathways']} significant)")
    print(f"GSVA pathways: {metadata['total_gsva_pathways']} ({metadata['significant_gsva_pathways']} significant)")
    print("=" * 60)


if __name__ == "__main__":
    main()
