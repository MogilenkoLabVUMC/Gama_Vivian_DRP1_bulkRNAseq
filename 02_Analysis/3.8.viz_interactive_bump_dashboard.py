#!/usr/bin/env python3
"""
Interactive Bump Chart Dashboard Generator
==========================================

Generates a single, self-contained HTML dashboard for exploring pathway trajectories.
Combines features of static bump charts (weighted, curved, rank/NES) with interactive filtering.

Usage:
    python3 02_Analysis/3.8.viz_interactive_bump_dashboard.py

Output:
    03_Results/02_Analysis/Plots/Trajectory_Flow/interactive_bump_dashboard.html
"""

import sys
import json
import numpy as np
import pandas as pd
from pathlib import Path

# Add module paths
sys.path.insert(0, str(Path(__file__).parent.parent / '01_Scripts'))

from Python.config import resolve_path
from Python.viz_bump_charts import (
    load_data, 
    filter_by_scope, 
    compute_weight_categories, 
    BumpChartConfig,
    MEANINGFUL_PATTERNS
)
from Python.pattern_definitions import get_pattern_colors, PATTERN_DEFINITIONS

OUTPUT_DIR = resolve_path('03_Results/02_Analysis/Plots/Trajectory_Flow')
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_FILE = OUTPUT_DIR / "interactive_bump_dashboard.html"

# CDN Links
PLOTLY_JS_URL = "https://cdn.plot.ly/plotly-2.27.0.min.js"

def prepare_data_with_pvals(df):
    """
    Pivot p-values and merge to create a wide-format dataframe with all stats.
    Adapted from 3.7.viz_interactive_bump.py.
    """
    # Create mapping for column renaming
    contrast_map = {
        'G32A_vs_Ctrl_D35': 'padj_Early_G32A',
        'G32A_vs_Ctrl_D65': 'padj_Late_G32A',
        'Maturation_G32A_specific': 'padj_TrajDev_G32A',
        'R403C_vs_Ctrl_D35': 'padj_Early_R403C',
        'R403C_vs_Ctrl_D65': 'padj_Late_R403C',
        'Maturation_R403C_specific': 'padj_TrajDev_R403C',
    }
    
    # Filter for relevant contrasts
    df_filtered = df[df['contrast'].isin(contrast_map.keys())]
    
    # Pivot p.adjust
    pvals = df_filtered.pivot(index='pathway_id', columns='contrast', values='p.adjust')
    pvals = pvals.rename(columns=contrast_map)
    
    # Base dataframe (unique pathways)
    df_base = df.drop_duplicates(subset=['pathway_id']).copy()
    
    # Join
    df_final = df_base.set_index('pathway_id').join(pvals).reset_index()
    
    return df_final

def prepare_dashboard_data():
    """
    Load, process, and structure data for the dashboard.
    Returns a dictionary of data and metadata.
    """
    print("Loading data...")
    df = load_data()
    
    # Enrich with p-values
    df = prepare_data_with_pvals(df)
    
    # 1. Filter Data
    # We want a superset of "Focused" and "Significant" to allow exploration
    df_focused = filter_by_scope(df, 'focused')
    df_sig = filter_by_scope(df, 'significant')
    
    # Combine and deduplicate
    df_combined = pd.concat([df_focused, df_sig]).drop_duplicates(subset=['pathway_id'])
    
    print(f"Data Loaded: {len(df)} total -> {len(df_combined)} in dashboard (Focused + Significant)")
    
    # 2. Compute Metadata & Attributes
    
    # Weight Categories (Global computation on the subset)
    weight_cats = compute_weight_categories(df_combined)
    
    # Pre-calculate Ranks for both mutations (NES based)
    # We want rank within this filtered dataset
    for mut in ['G32A', 'R403C']:
        # NES Ranks (Descending)
        for stage in ['Early', 'Late']:
            col = f'NES_{stage}_{mut}'
            rank_col = f'Rank_{stage}_{mut}'
            if col in df_combined.columns:
                df_combined[rank_col] = df_combined[col].rank(ascending=False, method='min')

    # 3. Structure Data for JSON Export
    # We need a list of pathway objects
    pathways = []
    
    # Columns to export
    base_cols = ['pathway_id', 'Description', 'database']
    
    for _, row in df_combined.iterrows():
        p = {col: row[col] for col in base_cols}
        
        # Per-mutation data
        for mut in ['G32A', 'R403C']:
            pattern = row.get(f'Pattern_{mut}')
            if pd.isna(pattern): pattern = None
            
            p[f'Pattern_{mut}'] = pattern
            p[f'NES_Early_{mut}'] = row.get(f'NES_Early_{mut}')
            p[f'NES_TrajDev_{mut}'] = row.get(f'NES_TrajDev_{mut}')
            p[f'NES_Late_{mut}'] = row.get(f'NES_Late_{mut}')
            
            p[f'padj_Early_{mut}'] = row.get(f'padj_Early_{mut}')
            p[f'padj_TrajDev_{mut}'] = row.get(f'padj_TrajDev_{mut}')
            p[f'padj_Late_{mut}'] = row.get(f'padj_Late_{mut}')
            
            p[f'Rank_Early_{mut}'] = row.get(f'Rank_Early_{mut}')
            p[f'Rank_Late_{mut}'] = row.get(f'Rank_Late_{mut}')
            
            p[f'Sig_TrajDev_{mut}'] = bool(row.get(f'Sig_TrajDev_{mut}', False))
            
        pathways.append(p)
        
    # Metadata
    metadata = {
        'weight_categories': weight_cats,
        'pattern_colors': get_pattern_colors(),
        'pattern_definitions': {k: v.get('interpretation', '') for k, v in PATTERN_DEFINITIONS.items()},
        'databases': sorted(df_combined['database'].dropna().unique().tolist()),
        'patterns': sorted([p for p in df_combined['Pattern_G32A'].dropna().unique() if p]), # G32A patterns as reference
    }
    
    return pathways, metadata

def generate_html(pathways, metadata):
    """Generate the HTML string."""
    
    pathways_json = json.dumps(pathways, default=str)
    metadata_json = json.dumps(metadata)
    
    # Use raw string and manual replacement to avoid f-string/JS conflict
    html_head = r"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Interactive Bump Chart Dashboard</title>
    <script src="__PLOTLY_JS_URL__"></script>
    <style>
        :root {
            --primary: #2c3e50;
            --accent: #3498db;
            --bg: #f8f9fa;
            --sidebar-width: 340px;
        }
        
        body { margin: 0; padding: 0; font-family: -apple-system, system-ui, sans-serif; background: var(--bg); display: flex; height: 100vh; overflow: hidden; }
        
        /* Sidebar */
        .sidebar {
            width: var(--sidebar-width);
            background: white;
            border-right: 1px solid #ddd;
            padding: 20px;
            overflow-y: auto;
            display: flex;
            flex-direction: column;
            gap: 20px;
            box-shadow: 2px 0 5px rgba(0,0,0,0.05);
            font-size: 0.9em;
            flex-shrink: 0;
        }
        
        .control-group { display: flex; flex-direction: column; gap: 8px; }
        .control-group label { font-weight: 600; font-size: 1em; color: var(--primary); display: flex; align-items: center; justify-content: space-between; }
        
        /* Custom Tooltip for Info Icons */
        .info-container {
            position: relative;
            display: inline-block;
        }
        
        .info-icon {
            cursor: help;
            color: #999;
            border: 1px solid #ccc;
            border-radius: 50%;
            width: 16px;
            height: 16px;
            font-size: 11px;
            display: inline-flex;
            align-items: center;
            justify-content: center;
            margin-left: 5px;
            background: #fff;
        }
        
        .info-icon:hover {
            color: var(--accent);
            border-color: var(--accent);
            font-weight: bold;
        }
        
        .tooltip-text {
            visibility: hidden;
            width: 220px;
            background-color: #333;
            color: #fff;
            text-align: left;
            border-radius: 6px;
            padding: 8px;
            position: absolute;
            z-index: 1000;
            bottom: 125%; /* Position above */
            right: 0; /* Align right */
            margin-bottom: 5px;
            opacity: 0;
            transition: none;
            font-size: 0.85em;
            font-weight: normal;
            line-height: 1.4;
            box-shadow: 0 2px 5px rgba(0,0,0,0.2);
            pointer-events: none;
        }
        
        .info-container:hover .tooltip-text {
            visibility: visible;
            opacity: 1;
        }
        
        .checkbox-list {
            display: flex;
            flex-direction: column;
            gap: 4px;
            max-height: 150px;
            overflow-y: auto;
            border: 1px solid #eee;
            padding: 5px;
            border-radius: 4px;
        }
        
        .checkbox-item { display: flex; align-items: center; gap: 6px; font-size: 0.9em; cursor: pointer; }
        .checkbox-item:hover { background: #f0f8ff; }
        
        input[type="text"] { padding: 8px; border: 1px solid #ddd; border-radius: 4px; width: 100%; box-sizing: border-box; }
        input[type="number"] { padding: 4px; border: 1px solid #ddd; border-radius: 4px; width: 60px; }
        
        .btn { padding: 8px 12px; background: var(--accent); color: white; border: none; border-radius: 4px; cursor: pointer; font-size: 0.9em; }
        .btn:hover { background: #2980b9; }
        
        /* Main Content */
        .main { flex: 1; display: flex; flex-direction: column; padding: 20px; gap: 10px; min-width: 0; }
        
        .header { display: flex; justify-content: space-between; align-items: center; padding-bottom: 10px; border-bottom: 1px solid #ddd; }
        .title { font-size: 1.2em; font-weight: bold; color: var(--primary); }
        .stats { font-size: 0.9em; color: #666; }
        
        .charts-container { flex: 1; display: flex; gap: 10px; min-height: 0; }
        .chart-box { flex: 1; background: white; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.05); display: flex; flex-direction: column; min-width: 0; }
        .chart-header { padding: 10px; font-weight: bold; text-align: center; border-bottom: 1px solid #eee; }
        .plotly-div { flex: 1; min-height: 0; width: 100%; }
        
        hr { border: 0; border-top: 1px solid #eee; margin: 10px 0; width: 100%; }
        
        .filter-row { display: flex; align-items: center; justify-content: space-between; font-size: 0.9em; margin-bottom: 4px; }
        
    </style>
</head>
"""

    html_body = r"""
<body>

<div class="sidebar">
    <div class="control-group">
        <label>Display Mode</label>
        <div class="checkbox-item"><input type="radio" name="mode" value="uniform" onchange="updateCharts()"> Uniform (All lines same width)</div>
        <div class="checkbox-item"><input type="radio" name="mode" value="weighted" checked onchange="updateCharts()"> Weighted (Thinner = Dominant)</div>
    </div>
    
    <div class="control-group">
        <label>
            Y-Axis Metric
            <div class="info-container">
                <span class="info-icon">?</span>
                <div class="tooltip-text"><b>NES (Normalized Enrichment Score):</b> Effect strength and direction. Positive = upregulated, negative = downregulated.<br><br><b>Rank:</b> Relative ordering based on NES (1 = highest NES). Shows position changes rather than magnitude.</div>
            </div>
        </label>
        <div class="checkbox-item"><input type="radio" name="ytype" value="nes" checked onchange="updateCharts()"> NES (Expression Strength)</div>
        <div class="checkbox-item"><input type="radio" name="ytype" value="rank" onchange="updateCharts()"> Rank (Relative Position)</div>
    </div>
    
    <div class="control-group">
        <label>
            Visual Style
            <div class="info-container">
                <span class="info-icon">?</span>
                <div class="tooltip-text"><b>Curved Lines:</b> Show trajectory deviation (TrajDev) at intermediate stage.<br><br><b>Works with both NES and Rank:</b><br>• NES: Curve magnitude = TrajDev NES value<br>• Rank: Curve scaled relative to chart height<br><br><b>Note:</b> Curves only appear for pathways with significant TrajDev (p.adjust &lt; 0.05).</div>
            </div>
        </label>
        <div class="checkbox-item"><input type="checkbox" id="chk-curved" onchange="updateCharts()"> Curved Lines (Trajectory Deviation)</div>
    </div>

    <div class="control-group">
        <label>
            Color By
            <div class="info-container">
                <span class="info-icon">?</span>
                <div class="tooltip-text">Choose how lines are colored.<br><br><b>Pattern:</b> Default trajectory classification colors.<br><b>NES options:</b> Blue-White-Orange diverging scale (Blue = downregulated, Orange = upregulated).</div>
            </div>
        </label>
        <select id="color-by" onchange="updateCharts()" style="padding: 6px; border: 1px solid #ddd; border-radius: 4px; width: 100%;">
            <option value="pattern" selected>Pattern (default)</option>
            <option value="nes_early">Early NES</option>
            <option value="nes_late">Late NES</option>
            <option value="nes_trajdev">TrajDev NES</option>
        </select>
    </div>

    <hr>
    
    <div class="control-group">
        <label>
            Advanced Filters
            <div class="info-container">
                <span class="info-icon">?</span>
                <div class="tooltip-text">Filter which lines are displayed based on significance and magnitude.</div>
            </div>
        </label>
        <div class="filter-row">
            <span>Min Abs(NES):
                <div class="info-container" style="display: inline;">
                    <span class="info-icon" style="font-size: 0.8em;">?</span>
                    <div class="tooltip-text">Pathway passes if <b>any</b> contrast (Early, TrajDev, or Late) has |NES| ≥ threshold.<br><br>Use Contrast-Specific Filters below for AND logic on individual contrasts.</div>
                </div>
            </span>
            <input type="number" id="min-nes" value="0" step="0.1" onchange="updateCharts()">
        </div>
        <div class="checkbox-item">
            <input type="checkbox" id="filter-sig-any" onchange="updateCharts()"> 
            Significant (Any Stage)
        </div>
        <div class="checkbox-item">
            <input type="checkbox" id="filter-sig-traj" onchange="updateCharts()">
            Significant TrajDev
        </div>
        <details style="margin-top: 10px;">
            <summary style="cursor: pointer; font-weight: 500; font-size: 0.9em; color: #555;">Contrast-Specific Filters
                <div class="info-container" style="display: inline;">
                    <span class="info-icon" style="font-size: 0.8em;">?</span>
                    <div class="tooltip-text"><b>AND logic:</b> All enabled filters must pass.<br><br>Example: Setting |Early NES| ≥ 1.5 AND |Late NES| ≥ 1.5 shows only pathways strong at both timepoints.</div>
                </div>
            </summary>
            <div style="padding: 8px 0; display: flex; flex-direction: column; gap: 6px;">
                <div class="filter-row">
                    <span>|Early NES| ≥</span>
                    <input type="number" id="filter-early-nes" value="" step="0.1" placeholder="0" onchange="updateCharts()">
                </div>
                <div class="filter-row">
                    <span>|Late NES| ≥</span>
                    <input type="number" id="filter-late-nes" value="" step="0.1" placeholder="0" onchange="updateCharts()">
                </div>
                <div class="filter-row">
                    <span>|TrajDev| ≥</span>
                    <input type="number" id="filter-trajdev-nes" value="" step="0.1" placeholder="0" onchange="updateCharts()">
                </div>
                <div class="checkbox-item">
                    <input type="checkbox" id="filter-sig-early" onchange="updateCharts()">
                    Sig. Early (p.adj &lt; 0.05)
                </div>
                <div class="checkbox-item">
                    <input type="checkbox" id="filter-sig-late" onchange="updateCharts()">
                    Sig. Late (p.adj &lt; 0.05)
                </div>
            </div>
        </details>
    </div>

    <hr>

    <div class="control-group">
        <label>Filter: Patterns</label>
        <div class="checkbox-list" id="pattern-list">
            <!-- Populated by JS -->
        </div>
        <button class="btn" style="margin-top: 5px; font-size: 0.8em; padding: 4px;" onclick="toggleAll('pattern-list')">Toggle All</button>
    </div>

    <div class="control-group">
        <label>Filter: Databases</label>
        <div class="checkbox-list" id="db-list">
             <!-- Populated by JS -->
        </div>
        <button class="btn" style="margin-top: 5px; font-size: 0.8em; padding: 4px;" onclick="toggleAll('db-list')">Toggle All</button>
    </div>
    
    <div class="control-group">
        <label>Highlight Pathway</label>
        <input type="text" id="search-box" placeholder="Search description..." oninput="updateHighlights()">
    </div>
    
    <div class="control-group" style="margin-top: auto;">
        <button class="btn" onclick="resetView()">Reset View</button>
    </div>
</div>

<div class="main">
    <div class="header">
        <div class="title">Pathway Trajectory Dashboard</div>
        <div class="stats" id="status-bar">Loading...</div>
    </div>
    
    <div class="charts-container">
        <div class="chart-box">
            <div class="chart-header" style="color: #0072B2;">G32A Mutation</div>
            <div id="chart-g32a" class="plotly-div"></div>
        </div>
        <div class="chart-box">
            <div class="chart-header" style="color: #D55E00;">R403C Mutation</div>
            <div id="chart-r403c" class="plotly-div"></div>
        </div>
    </div>
</div>
"""

    html_script = r"""
<script>
    // ========================================================================
    // DATA & CONFIG
    // ========================================================================
    const RAW_DATA = __PATHWAYS_JSON__;
    const METADATA = __METADATA_JSON__;
    
    let filteredData = [];
    let highlightSearch = "";
    
    // Global scale state
    let globalYMin = 0;
    let globalYMax = 0;
    
    // Config Defaults
    const PATTERN_COLORS = METADATA.pattern_colors;
    const WEIGHT_CATS = METADATA.weight_categories;
    const PATTERN_DEFS = METADATA.pattern_definitions;
    
    // Weight Style Mapping (LineWidth, Opacity)
    // Categories: dominant, common, uncommon, rare
    const WEIGHT_STYLES = {
        'dominant': { width: 1.0, opacity: 0.3 },
        'common':   { width: 2.0, opacity: 0.6 },
        'uncommon': { width: 3.0, opacity: 0.8 },
        'rare':     { width: 4.0, opacity: 1.0 }
    };
    
    const UNIFORM_STYLE = { width: 1.5, opacity: 0.5 };

    // Colorblind-safe diverging scale (matches R color_config.R)
    const NES_COLORS = {
        negative: '#2166AC',  // Blue
        neutral: '#F7F7F7',   // White
        positive: '#B35806'   // Orange
    };

    function nesToColor(nes, maxNes = 3.5) {
        if (nes === null || nes === undefined || isNaN(nes)) return '#999999';
        const t = Math.max(-1, Math.min(1, nes / maxNes)); // Clamp to [-1, 1]
        if (t <= 0) {
            // Negative: Blue to White
            return interpolateColor(NES_COLORS.negative, NES_COLORS.neutral, 1 + t);
        } else {
            // Positive: White to Orange
            return interpolateColor(NES_COLORS.neutral, NES_COLORS.positive, t);
        }
    }

    function interpolateColor(c1, c2, t) {
        const r1 = parseInt(c1.slice(1,3), 16);
        const g1 = parseInt(c1.slice(3,5), 16);
        const b1 = parseInt(c1.slice(5,7), 16);
        const r2 = parseInt(c2.slice(1,3), 16);
        const g2 = parseInt(c2.slice(3,5), 16);
        const b2 = parseInt(c2.slice(5,7), 16);
        const r = Math.round(r1 + (r2 - r1) * t);
        const g = Math.round(g1 + (g2 - g1) * t);
        const b = Math.round(b1 + (b2 - b1) * t);
        return `#${r.toString(16).padStart(2,'0')}${g.toString(16).padStart(2,'0')}${b.toString(16).padStart(2,'0')}`;
    }

    // ========================================================================
    // INITIALIZATION
    // ========================================================================
    function init() {
        populateCheckboxList('pattern-list', METADATA.patterns);
        populateCheckboxList('db-list', METADATA.databases);
        updateCharts();
    }
    
    function populateCheckboxList(id, items) {
        const container = document.getElementById(id);
        container.innerHTML = '';
        items.forEach(item => {
            const div = document.createElement('div');
            div.className = 'checkbox-item';
            div.innerHTML = `<input type="checkbox" value="${item}" checked onchange="updateCharts()"> ${item}`;
            if (id === 'pattern-list' && PATTERN_COLORS[item]) {
                div.innerHTML += `<span style="width:10px;height:10px;background:${PATTERN_COLORS[item]};border-radius:50%;display:inline-block;"></span>`;
            }
            container.appendChild(div);
        });
    }
    
    function toggleAll(id) {
        const inputs = document.querySelectorAll(`#${id} input`);
        const allChecked = Array.from(inputs).every(i => i.checked);
        inputs.forEach(i => i.checked = !allChecked);
        updateCharts();
    }
    
    function resetView() {
        document.querySelectorAll('input[type="checkbox"]').forEach(i => i.checked = true);
        document.getElementById('chk-curved').checked = false;
        document.querySelector('input[name="mode"][value="weighted"]').checked = true;
        document.querySelector('input[name="ytype"][value="nes"]').checked = true;
        document.getElementById('search-box').value = "";
        document.getElementById('min-nes').value = "0";
        document.getElementById('filter-sig-any').checked = false;
        document.getElementById('filter-sig-traj').checked = false;
        // Reset new controls
        document.getElementById('color-by').value = 'pattern';
        document.getElementById('filter-early-nes').value = '';
        document.getElementById('filter-late-nes').value = '';
        document.getElementById('filter-trajdev-nes').value = '';
        document.getElementById('filter-sig-early').checked = false;
        document.getElementById('filter-sig-late').checked = false;
        highlightSearch = "";
        updateCharts();
    }
    
    // ========================================================================
    // CORE LOGIC
    // ========================================================================
    
    function getSettings() {
        return {
            mode: document.querySelector('input[name="mode"]:checked').value,
            yType: document.querySelector('input[name="ytype"]:checked').value,
            curved: document.getElementById('chk-curved').checked,
            colorBy: document.getElementById('color-by').value,
            patterns: Array.from(document.querySelectorAll('#pattern-list input:checked')).map(i => i.value),
            dbs: Array.from(document.querySelectorAll('#db-list input:checked')).map(i => i.value),
            minNes: parseFloat(document.getElementById('min-nes').value) || 0,
            sigAny: document.getElementById('filter-sig-any').checked,
            sigTraj: document.getElementById('filter-sig-traj').checked,
            // Contrast-specific filters
            filterEarlyNes: parseFloat(document.getElementById('filter-early-nes').value) || 0,
            filterLateNes: parseFloat(document.getElementById('filter-late-nes').value) || 0,
            filterTrajdevNes: parseFloat(document.getElementById('filter-trajdev-nes').value) || 0,
            sigEarly: document.getElementById('filter-sig-early').checked,
            sigLate: document.getElementById('filter-sig-late').checked
        };
    }
    
    function updateHighlights() {
        highlightSearch = document.getElementById('search-box').value.toLowerCase();
        updateCharts(true); 
    }
    
    function updateCharts(isFastUpdate = false) {
        const settings = getSettings();
        
        // 1. Filter Data
        filteredData = RAW_DATA.filter(d => 
            settings.dbs.includes(d.database) &&
            (settings.patterns.includes(d.Pattern_G32A) || settings.patterns.includes(d.Pattern_R403C))
        );
        
        // 2. Calculate Global Scale
        if (settings.yType === 'nes') {
            let maxAbs = 0;
            filteredData.forEach(d => {
                [d.NES_Early_G32A, d.NES_Late_G32A, d.NES_TrajDev_G32A,
                 d.NES_Early_R403C, d.NES_Late_R403C, d.NES_TrajDev_R403C].forEach(v => {
                   if (v !== null && !isNaN(v)) maxAbs = Math.max(maxAbs, Math.abs(v));
                });
            });
            maxAbs = Math.max(maxAbs * 1.05, 1.0); 
            globalYMax = maxAbs;
            globalYMin = -maxAbs;
        } else {
            // Rank
            globalYMin = 1;
            globalYMax = filteredData.length || 100;
        }
        
        document.getElementById('status-bar').textContent = `Showing ${filteredData.length} pathways`;
        
        renderChart('chart-g32a', 'G32A', settings);
        renderChart('chart-r403c', 'R403C', settings);
    }
    
    function renderChart(divId, mutation, settings) {
        const traces = [];
        const fmt = (n) => n !== null && n !== undefined ? n.toFixed(2) : 'N/A';
        const fmtP = (n) => n !== null && n !== undefined ? n.toExponential(2) : 'N/A';
        const sig = (p) => (p !== null && p < 0.05) ? '*' : '';
        
        // 1. Filter for visibility (Mutation Specific)
        const patDataAll = filteredData.filter(d => {
            // Pattern check
            const pat = d[`Pattern_${mutation}`];
            if (!pat || !settings.patterns.includes(pat)) return false;
            
            // NES Threshold Check (Any stage > threshold)
            if (settings.minNes > 0) {
                const maxN = Math.max(
                    Math.abs(d[`NES_Early_${mutation}`] || 0), 
                    Math.abs(d[`NES_Late_${mutation}`] || 0),
                    Math.abs(d[`NES_TrajDev_${mutation}`] || 0)
                );
                if (maxN < settings.minNes) return false;
            }
            
            // Sig Check
            if (settings.sigAny) {
                const isSig = (d[`padj_Early_${mutation}`] < 0.05) || 
                              (d[`padj_Late_${mutation}`] < 0.05) || 
                              (d[`padj_TrajDev_${mutation}`] < 0.05);
                if (!isSig) return false;
            }
            
            if (settings.sigTraj && !d[`Sig_TrajDev_${mutation}`]) return false;

            // Contrast-specific NES thresholds
            if (settings.filterEarlyNes > 0) {
                if (Math.abs(d[`NES_Early_${mutation}`] || 0) < settings.filterEarlyNes) return false;
            }
            if (settings.filterLateNes > 0) {
                if (Math.abs(d[`NES_Late_${mutation}`] || 0) < settings.filterLateNes) return false;
            }
            if (settings.filterTrajdevNes > 0) {
                if (Math.abs(d[`NES_TrajDev_${mutation}`] || 0) < settings.filterTrajdevNes) return false;
            }

            // Contrast-specific significance filters
            if (settings.sigEarly && (d[`padj_Early_${mutation}`] === null || d[`padj_Early_${mutation}`] >= 0.05)) return false;
            if (settings.sigLate && (d[`padj_Late_${mutation}`] === null || d[`padj_Late_${mutation}`] >= 0.05)) return false;

            return true;
        });

        // Group by Pattern
        const patternCounts = {};
        patDataAll.forEach(d => {
            const p = d[`Pattern_${mutation}`];
            patternCounts[p] = (patternCounts[p] || 0) + 1;
        });
        
        // Sort Patterns: Primary = weight category, Secondary = frequency (higher freq first = bottom)
        const sortedPatterns = Object.keys(patternCounts).sort((a, b) => {
             const rank = p => {
                 const cat = WEIGHT_CATS[p] || 'rare';
                 if(cat === 'dominant') return 0;
                 if(cat === 'common') return 1;
                 if(cat === 'uncommon') return 2;
                 return 3;
             };
             // Primary sort by weight category
             const rankDiff = rank(a) - rank(b);
             if (rankDiff !== 0) return rankDiff;
             // Secondary sort by frequency (higher frequency = render first = bottom)
             return patternCounts[b] - patternCounts[a];
        });

        // DEBUG: Render order verification
        if (typeof console !== 'undefined') {
            console.group(`Render Order Debug - ${mutation}`);
            console.log('Pattern counts (per-mutation):', JSON.stringify(patternCounts, null, 2));
            console.log('Weight categories (global):', JSON.stringify(
                Object.fromEntries(sortedPatterns.map(p => [p, WEIGHT_CATS[p] || 'unknown']))
            , null, 2));
            console.log('Sorted patterns (first=bottom, last=top):', sortedPatterns);

            // Verify sort is correct
            const commonPatterns = sortedPatterns.filter(p => WEIGHT_CATS[p] === 'common');
            console.log('Common category patterns in order:', commonPatterns);
            console.log('Their counts:', commonPatterns.map(p => `${p}: ${patternCounts[p]}`));
            console.groupEnd();
        }

        // Render Traces - Strategy depends on colorBy setting
        let traceIndex = 0;  // DEBUG: track trace order
        if (settings.colorBy === 'pattern') {
            // Group by pattern - original approach
            sortedPatterns.forEach(pattern => {
                const rows = patDataAll.filter(d => d[`Pattern_${mutation}`] === pattern);
                console.log(`[${mutation}] Trace ${traceIndex++}: ${pattern} (${rows.length} paths, cat: ${WEIGHT_CATS[pattern]})`);

                // Style
                const color = PATTERN_COLORS[pattern] || '#999';
                let lw = UNIFORM_STYLE.width;
                let op = UNIFORM_STYLE.opacity;
                if (settings.mode === 'weighted') {
                    const cat = WEIGHT_CATS[pattern] || 'rare';
                    const style = WEIGHT_STYLES[cat];
                    lw = style.width;
                    op = style.opacity;
                }

                const xs = [];
                const ys = [];
                const texts = [];

                rows.forEach(row => {
                    const y1 = settings.yType === 'rank' ? row[`Rank_Early_${mutation}`] : row[`NES_Early_${mutation}`];
                    const y2 = settings.yType === 'rank' ? row[`Rank_Late_${mutation}`] : row[`NES_Late_${mutation}`];

                    if (y1 === null || y2 === null) return;

                    // Tooltip Wrapping
                    const desc = wrapText(row.Description, 35);
                    const patDesc = wrapText(PATTERN_DEFS[pattern] || "", 45);

                    const hover = `<b>${desc}</b><br>` +
                                  `DB: ${row.database} | Pat: ${pattern}<br><br>` +
                                  `<b>NES & Sig:</b><br>` +
                                  `E: ${fmt(row[`NES_Early_${mutation}`])}${sig(row[`padj_Early_${mutation}`])}<br>` +
                                  `T: ${fmt(row[`NES_TrajDev_${mutation}`])}${sig(row[`padj_TrajDev_${mutation}`])}<br>` +
                                  `L: ${fmt(row[`NES_Late_${mutation}`])}${sig(row[`padj_Late_${mutation}`])}<br>` +
                                  `<br><i>${patDesc}</i>`;

                    // Curve Logic (Enabled for both NES and Rank now)
                    if (settings.curved && row[`Sig_TrajDev_${mutation}`]) {
                        const traj = row[`NES_TrajDev_${mutation}`] || 0;
                        const yMid = (y1 + y2) / 2;
                        let offset = 0;

                        if (settings.yType === 'rank') {
                            // Rank deviation scaling
                            const chartHeight = globalYMax - globalYMin;
                            offset = -1 * traj * (chartHeight * 0.1);
                        } else {
                            offset = traj * 0.5;
                        }

                        const pts = getBezier(0, y1, 0.5, yMid + offset, 1, y2, 20);
                        xs.push(...pts.x, null);
                        ys.push(...pts.y, null);
                        texts.push(...Array(pts.x.length).fill(hover), null);
                    } else {
                        xs.push(0, 1, null);
                        ys.push(y1, y2, null);
                        texts.push(hover, hover, null);
                    }
                });

                if (xs.length > 0) {
                    traces.push({
                        x: xs, y: ys, mode: 'lines',
                        line: { color: color, width: lw }, opacity: op,
                        name: pattern, text: texts, hoverinfo: 'text',
                        type: 'scattergl'
                    });
                }
            });
        } else {
            // NES-based coloring - one trace per row for individual colors
            patDataAll.forEach(row => {
                const pattern = row[`Pattern_${mutation}`];
                const y1 = settings.yType === 'rank' ? row[`Rank_Early_${mutation}`] : row[`NES_Early_${mutation}`];
                const y2 = settings.yType === 'rank' ? row[`Rank_Late_${mutation}`] : row[`NES_Late_${mutation}`];

                if (y1 === null || y2 === null) return;

                // Determine NES value for coloring
                let nesForColor = 0;
                if (settings.colorBy === 'nes_early') {
                    nesForColor = row[`NES_Early_${mutation}`] || 0;
                } else if (settings.colorBy === 'nes_late') {
                    nesForColor = row[`NES_Late_${mutation}`] || 0;
                } else if (settings.colorBy === 'nes_trajdev') {
                    nesForColor = row[`NES_TrajDev_${mutation}`] || 0;
                }
                const color = nesToColor(nesForColor);

                // Style
                let lw = UNIFORM_STYLE.width;
                let op = UNIFORM_STYLE.opacity;
                if (settings.mode === 'weighted') {
                    const cat = WEIGHT_CATS[pattern] || 'rare';
                    const style = WEIGHT_STYLES[cat];
                    lw = style.width;
                    op = style.opacity;
                }

                // Tooltip
                const desc = wrapText(row.Description, 35);
                const patDesc = wrapText(PATTERN_DEFS[pattern] || "", 45);

                const hover = `<b>${desc}</b><br>` +
                              `DB: ${row.database} | Pat: ${pattern}<br><br>` +
                              `<b>NES & Sig:</b><br>` +
                              `E: ${fmt(row[`NES_Early_${mutation}`])}${sig(row[`padj_Early_${mutation}`])}<br>` +
                              `T: ${fmt(row[`NES_TrajDev_${mutation}`])}${sig(row[`padj_TrajDev_${mutation}`])}<br>` +
                              `L: ${fmt(row[`NES_Late_${mutation}`])}${sig(row[`padj_Late_${mutation}`])}<br>` +
                              `<br><i>${patDesc}</i>`;

                const xs = [];
                const ys = [];
                const texts = [];

                // Curve Logic
                if (settings.curved && row[`Sig_TrajDev_${mutation}`]) {
                    const traj = row[`NES_TrajDev_${mutation}`] || 0;
                    const yMid = (y1 + y2) / 2;
                    let offset = 0;

                    if (settings.yType === 'rank') {
                        const chartHeight = globalYMax - globalYMin;
                        offset = -1 * traj * (chartHeight * 0.1);
                    } else {
                        offset = traj * 0.5;
                    }

                    const pts = getBezier(0, y1, 0.5, yMid + offset, 1, y2, 20);
                    xs.push(...pts.x);
                    ys.push(...pts.y);
                    texts.push(...Array(pts.x.length).fill(hover));
                } else {
                    xs.push(0, 1);
                    ys.push(y1, y2);
                    texts.push(hover, hover);
                }

                traces.push({
                    x: xs, y: ys, mode: 'lines',
                    line: { color: color, width: lw }, opacity: op,
                    name: row.Description.slice(0, 30), text: texts, hoverinfo: 'text',
                    type: 'scattergl', showlegend: false
                });
            });
        }
        
        // Highlights Trace
        if (highlightSearch) {
            const hlXs = [];
            const hlYs = [];
            const hlTexts = [];
            
            // Search in visible data
            patDataAll.forEach(row => {
                 if (!row.Description.toLowerCase().includes(highlightSearch)) return;
                 const y1 = settings.yType === 'rank' ? row[`Rank_Early_${mutation}`] : row[`NES_Early_${mutation}`];
                 const y2 = settings.yType === 'rank' ? row[`Rank_Late_${mutation}`] : row[`NES_Late_${mutation}`];
                 if (y1 === null || y2 === null) return;

                 // Full rich tooltip for highlights (same format as regular lines)
                 const pattern = row[`Pattern_${mutation}`];
                 const desc = wrapText(row.Description, 35);
                 const patDesc = wrapText(PATTERN_DEFS[pattern] || "", 45);

                 const hover = `<b>${desc}</b><br>` +
                               `DB: ${row.database} | Pat: ${pattern}<br><br>` +
                               `<b>NES & Sig:</b><br>` +
                               `E: ${fmt(row[`NES_Early_${mutation}`])}${sig(row[`padj_Early_${mutation}`])} (p=${fmtP(row[`padj_Early_${mutation}`])})<br>` +
                               `T: ${fmt(row[`NES_TrajDev_${mutation}`])}${sig(row[`padj_TrajDev_${mutation}`])} (p=${fmtP(row[`padj_TrajDev_${mutation}`])})<br>` +
                               `L: ${fmt(row[`NES_Late_${mutation}`])}${sig(row[`padj_Late_${mutation}`])} (p=${fmtP(row[`padj_Late_${mutation}`])})<br>` +
                               `<br><i>${patDesc}</i>`;

                 hlXs.push(0, 1, null);
                 hlYs.push(y1, y2, null);
                 hlTexts.push(hover, hover, null);
            });
            
            if (hlXs.length > 0) {
                traces.push({
                    x: hlXs, y: hlYs, mode: 'lines+markers',
                    line: { color: 'black', width: 3 }, marker: { size: 6, color: 'black' },
                    text: hlTexts, hoverinfo: 'text', name: 'Highlights', type: 'scatter'
                });
            }
        }
        
        // Layout
        const layout = {
            margin: { t: 30, b: 30, l: 50, r: 20 },
            showlegend: false,
            xaxis: { tickvals: [0, 1], ticktext: ['Early', 'Late'], range: [-0.1, 1.1], zeroline: false },
            yaxis: {
                title: settings.yType === 'rank' ? 'Rank' : 'NES',
                range: [globalYMin, globalYMax],
                autorange: settings.yType === 'rank' ? 'reversed' : false
            },
            hovermode: 'closest',
            hoverlabel: {
                bgcolor: 'white',
                bordercolor: '#666',
                font: { size: 11, color: '#333', family: 'system-ui, sans-serif' },
                namelength: -1,
                align: 'left'
            }
        };
        
        Plotly.react(divId, traces, layout, { config: {displayModeBar: false} });
    }
    
    function getBezier(x0, y0, x1, y1, x2, y2, n) {
        const x = []; const y = [];
        for (let i = 0; i <= n; i++) {
            const t = i / n;
            const a = (1 - t) * (1 - t);
            const b = 2 * (1 - t) * t;
            const c = t * t;
            x.push(a * x0 + b * x1 + c * x2);
            y.push(a * y0 + b * y1 + c * y2);
        }
        return {x, y};
    }
    
    function wrapText(str, width = 35) {
        if (!str) return "";
        const words = str.split(' ');
        let currentLine = "";
        let result = "";
        words.forEach(word => {
            if ((currentLine + word).length > width) {
                result += currentLine + "<br>";
                currentLine = word + " ";
            } else {
                currentLine += word + " ";
            }
        });
        return result + currentLine;
    }

    init();

</script>
</body>
</html>"""

    # Assemble and Replace
    html = html_head + html_body + html_script
    
    html = html.replace('__PLOTLY_JS_URL__', PLOTLY_JS_URL)
    html = html.replace('__PATHWAYS_JSON__', pathways_json)
    html = html.replace('__METADATA_JSON__', metadata_json)
    
    return html


def main():
    """Generate the interactive bump chart dashboard."""
    print("=" * 80)
    print("INTERACTIVE BUMP CHART DASHBOARD GENERATOR")
    print("=" * 80)

    # 1. Load and Prepare Data
    pathways, metadata = prepare_dashboard_data()

    # 2. Generate HTML
    print("Generating dashboard HTML...")
    html = generate_html(pathways, metadata)

    # 3. Write to file
    print(f"Writing to {OUTPUT_FILE}...")
    OUTPUT_FILE.write_text(html, encoding='utf-8')

    print(f"\nDone! Dashboard saved to:\n  {OUTPUT_FILE}")
    print("\nOpen in browser to explore pathway trajectories interactively.")


if __name__ == '__main__':
    main()