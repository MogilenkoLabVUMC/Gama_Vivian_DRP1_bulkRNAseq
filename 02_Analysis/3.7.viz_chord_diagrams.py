#!/usr/bin/env python3
"""
GOChord-style Chord Diagram: Gene-Pathway Leading Edge Membership

Creates circular chord diagrams in the style of GOplot's GOChord function showing:
- Genes on the left with stacked rectangles for D35/D65 fold changes
- Pathways on the right as colored arcs with distinct colors
- Ribbon connections indicating gene membership in pathway's leading edge
- FDR > 0.25 shown in gray (non-significant) # threshold removed in the end, as no genes pass it

Generates separate figures for G32A and R403C mutations.

Usage:
    python3 02_Analysis/3.7.viz_chord_diagrams.py

Output:
    03_Results/02_Analysis/Plots/Chord_Diagrams/
"""

import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(project_root))
sys.path.insert(0, str(project_root / '01_Scripts'))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.path import Path as MplPath
from matplotlib.patches import PathPatch, Wedge, Rectangle, FancyBboxPatch
from matplotlib.colors import Normalize, LinearSegmentedColormap, to_rgba
from matplotlib.cm import ScalarMappable
import matplotlib.patheffects as path_effects
import warnings
warnings.filterwarnings('ignore')

# Import project color configuration
from Python.color_config import create_diverging_cmap, MUTATION_COLORS

# =============================================================================
# CONFIGURATION
# =============================================================================

# Updated pathways of interest - synaptic compartments
PATHWAYS_OF_INTEREST = {
    'syngo': [
        'SYNGO:presyn_ribosome',      # presynaptic ribosome
        'SYNGO:postsyn_ribosome',     # postsynaptic ribosome
        'GO:0099523',                 # presynaptic cytosol
        'GO:0099524',                 # postsynaptic cytosol
        'GO:0014069',                 # postsynaptic density
        'GO:0045211',                 # postsynaptic membrane
    ],
}

# Pathway display names (shorter for better readability)
PATHWAY_DISPLAY_NAMES = {
    'SYNGO:presyn_ribosome': 'Presynaptic\nRibosome',
    'SYNGO:postsyn_ribosome': 'Postsynaptic\nRibosome',
    'GO:0099523': 'Presynaptic\nCytosol',
    'GO:0099524': 'Postsynaptic\nCytosol',
    'GO:0014069': 'Postsynaptic\nDensity',
    'GO:0045211': 'Postsynaptic\nMembrane',
}

# Colors - colorblind-friendly palette (Wong 2011 Nature Methods + IBM Design)
# Optimized for deuteranopia, protanopia, and tritanopia
PATHWAY_COLORS = [
    '#882255',  # Wine/Burgundy - Presynaptic Ribosome
    '#332288',  # Indigo - Postsynaptic Ribosome
    '#117733',  # Forest Green - Presynaptic Cytosol
    '#88CCEE',  # Sky Blue - Postsynaptic Cytosol
    '#CC6677',  # Dusty Rose - Postsynaptic Density
    '#AA4499',  # Plum - Postsynaptic Membrane
    '#DDCC77',  # Sand/Gold
    '#999933',  # Olive
]

# Significance thresholds
PADJ_CUTOFF = 0.05
FDR_NONSIG = None  # Disabled - show all genes with fold change colors regardless of significance

# Output directory
OUTPUT_DIR = project_root / '03_Results/02_Analysis/Plots/Chord_Diagrams'

# =============================================================================
# DATA LOADING FUNCTIONS
# =============================================================================

def load_gsea_data_from_rds():
    """Load GSEA results from RDS files."""
    import subprocess

    # Build pathway list
    all_pathways = []
    for db, paths in PATHWAYS_OF_INTEREST.items():
        all_pathways.extend(paths)

    pathways_r = 'c("' + '","'.join(all_pathways) + '")'

    r_script = f'''
    library(jsonlite)

    syngo_results <- readRDS("03_Results/02_Analysis/checkpoints/syngo_gsea_results.rds")

    extract_pathway_data <- function(gsea_obj, pathway_ids) {{
        if (is.null(gsea_obj)) return(NULL)
        result <- gsea_obj@result[gsea_obj@result$ID %in% pathway_ids,
                                   c("ID", "Description", "NES", "pvalue", "p.adjust", "core_enrichment")]
        if (nrow(result) == 0) return(NULL)
        return(result)
    }}

    contrasts <- c("G32A_vs_Ctrl_D35", "G32A_vs_Ctrl_D65",
                   "R403C_vs_Ctrl_D35", "R403C_vs_Ctrl_D65")

    # Focus on synaptic compartment pathways
    syngo_pathways <- c("SYNGO:presyn_ribosome", "SYNGO:postsyn_ribosome",
                        "GO:0099523", "GO:0099524", "GO:0014069",
                        "GO:0045211")

    all_data <- list()

    for (contrast in contrasts) {{
        contrast_data <- list()

        if (contrast %in% names(syngo_results)) {{
            syngo_data <- extract_pathway_data(syngo_results[[contrast]], syngo_pathways)
            if (!is.null(syngo_data)) {{
                syngo_data$database <- "syngo"
                contrast_data$syngo <- syngo_data
            }}
        }}

        all_data[[contrast]] <- contrast_data
    }}

    combined <- do.call(rbind, lapply(names(all_data), function(contrast) {{
        contrast_df <- do.call(rbind, all_data[[contrast]])
        if (!is.null(contrast_df)) {{
            contrast_df$contrast <- contrast
            return(contrast_df)
        }}
        return(NULL)
    }}))

    write.csv(combined, "03_Results/02_Analysis/Plots/Chord_Diagrams/gsea_pathway_data.csv", row.names=FALSE)
    cat("Data exported successfully")
    '''

    result = subprocess.run(
        ['Rscript', '-e', r_script],
        cwd=str(project_root),
        capture_output=True,
        text=True
    )

    if result.returncode != 0:
        print(f"R script error: {result.stderr}")
        raise RuntimeError("Failed to extract GSEA data")

    gsea_df = pd.read_csv(OUTPUT_DIR / 'gsea_pathway_data.csv')
    return gsea_df


def load_de_results():
    """Load differential expression results."""
    de_dir = project_root / '03_Results/02_Analysis/DE_results'
    contrasts = ['G32A_vs_Ctrl_D35', 'G32A_vs_Ctrl_D65',
                 'R403C_vs_Ctrl_D35', 'R403C_vs_Ctrl_D65']

    de_data = {}
    for contrast in contrasts:
        file_path = de_dir / f'{contrast}_results.csv'
        if file_path.exists():
            df = pd.read_csv(file_path, index_col=0)
            de_data[contrast] = df

    return de_data


def prepare_chord_data(gsea_df, de_data, mutation='G32A', max_genes=70):
    """Prepare data for GOChord-style diagram."""
    d35_contrast = f'{mutation}_vs_Ctrl_D35'
    d65_contrast = f'{mutation}_vs_Ctrl_D65'

    mutation_df = gsea_df[gsea_df['contrast'].isin([d35_contrast, d65_contrast])].copy()

    # Get all pathways we want to show (not just significant ones)
    all_pathway_ids = []
    for db, paths in PATHWAYS_OF_INTEREST.items():
        all_pathway_ids.extend(paths)

    pathways = []
    all_genes = set()
    gene_pathway_map = {}

    for pathway_id in all_pathway_ids:
        pathway_rows = mutation_df[mutation_df['ID'] == pathway_id]
        if len(pathway_rows) == 0:
            continue

        nes_d35 = pathway_rows[pathway_rows['contrast'] == d35_contrast]['NES'].values
        nes_d65 = pathway_rows[pathway_rows['contrast'] == d65_contrast]['NES'].values
        padj_d35 = pathway_rows[pathway_rows['contrast'] == d35_contrast]['p.adjust'].values
        padj_d65 = pathway_rows[pathway_rows['contrast'] == d65_contrast]['p.adjust'].values

        core_genes_str = pathway_rows[pathway_rows['contrast'] == d35_contrast]['core_enrichment'].values
        if len(core_genes_str) == 0 or pd.isna(core_genes_str[0]):
            core_genes_str = pathway_rows[pathway_rows['contrast'] == d65_contrast]['core_enrichment'].values

        if len(core_genes_str) == 0 or pd.isna(core_genes_str[0]):
            continue

        core_genes = str(core_genes_str[0]).split('/')
        db = pathway_rows['database'].values[0]

        pathways.append({
            'id': pathway_id,
            'display_name': PATHWAY_DISPLAY_NAMES.get(pathway_id, pathway_id),
            'database': db,
            'nes_d35': nes_d35[0] if len(nes_d35) > 0 else np.nan,
            'nes_d65': nes_d65[0] if len(nes_d65) > 0 else np.nan,
            'padj_d35': padj_d35[0] if len(padj_d35) > 0 else 1.0,
            'padj_d65': padj_d65[0] if len(padj_d65) > 0 else 1.0,
            'core_genes': core_genes,
            'n_genes': len(core_genes)
        })

        all_genes.update(core_genes)
        for gene in core_genes:
            if gene not in gene_pathway_map:
                gene_pathway_map[gene] = []
            gene_pathway_map[gene].append(pathway_id)

    # Prioritize genes
    def gene_priority(gene):
        is_ribosomal = gene.startswith('RPL') or gene.startswith('RPS')
        n_pathways = len(gene_pathway_map.get(gene, []))
        return (0 if is_ribosomal else 1, -n_pathways, gene)

    sorted_genes = sorted(all_genes, key=gene_priority)

    if len(sorted_genes) > max_genes:
        print(f"   Filtering from {len(sorted_genes)} to {max_genes} genes")
        selected_genes = set(sorted_genes[:max_genes])
    else:
        selected_genes = set(sorted_genes)

    connections = []
    for gene in selected_genes:
        for pathway_id in gene_pathway_map.get(gene, []):
            connections.append((gene, pathway_id))

    genes = []
    de_d35 = de_data.get(d35_contrast)
    de_d65 = de_data.get(d65_contrast)

    for gene in selected_genes:
        logfc_d35 = de_d35.loc[gene, 'logFC'] if de_d35 is not None and gene in de_d35.index else np.nan
        logfc_d65 = de_d65.loc[gene, 'logFC'] if de_d65 is not None and gene in de_d65.index else np.nan
        padj_d35 = de_d35.loc[gene, 'adj.P.Val'] if de_d35 is not None and gene in de_d35.index else 1.0
        padj_d65 = de_d65.loc[gene, 'adj.P.Val'] if de_d65 is not None and gene in de_d65.index else 1.0

        genes.append({
            'name': gene,
            'logfc_d35': logfc_d35,
            'logfc_d65': logfc_d65,
            'padj_d35': padj_d35,
            'padj_d65': padj_d65,
            'n_pathways': len(gene_pathway_map.get(gene, [])),
            'is_ribosomal': gene.startswith('RPL') or gene.startswith('RPS')
        })

    genes = sorted(genes, key=lambda x: (0 if x['is_ribosomal'] else 1, -x['n_pathways'], x['name']))

    return {
        'pathways': pathways,
        'genes': genes,
        'connections': connections,
        'mutation': mutation,
        'total_genes': len(all_genes),
        'filtered_genes': len(selected_genes)
    }


# =============================================================================
# GOCHORD-STYLE DRAWING FUNCTIONS
# =============================================================================

def draw_ribbon(ax, gene_angle, gene_width, pathway_start, pathway_end, color, alpha=0.6):
    """Draw a ribbon connecting gene rectangle to pathway arc (GOChord style)."""
    # Gene rectangle position (on outer circle - same radius as pathways)
    gene_rad_start = np.radians(gene_angle - gene_width/2)
    gene_rad_end = np.radians(gene_angle + gene_width/2)
    inner_r = 0.85  # Same as pathway inner radius

    # Pathway arc position (on outer circle)
    pathway_rad_start = np.radians(pathway_start)
    pathway_rad_end = np.radians(pathway_end)
    outer_r = 0.85

    # Create ribbon path with bezier curves
    n_pts = 30

    # Gene side (bottom of ribbon)
    gene_pts_x = [inner_r * np.cos(gene_rad_start), inner_r * np.cos(gene_rad_end)]
    gene_pts_y = [inner_r * np.sin(gene_rad_start), inner_r * np.sin(gene_rad_end)]

    # Pathway side (top of ribbon)
    pathway_theta = np.linspace(pathway_rad_start, pathway_rad_end, n_pts)
    pathway_pts_x = outer_r * np.cos(pathway_theta)
    pathway_pts_y = outer_r * np.sin(pathway_theta)

    # Bezier control points
    ctrl_r = 0.3  # Control point radius (toward center)

    # Build path: gene_start -> pathway_start -> pathway_end -> gene_end -> close
    verts = []
    codes = []

    # Start at gene left edge
    verts.append((gene_pts_x[0], gene_pts_y[0]))
    codes.append(MplPath.MOVETO)

    # Bezier to pathway start
    ctrl1_x = ctrl_r * np.cos(gene_rad_start)
    ctrl1_y = ctrl_r * np.sin(gene_rad_start)
    ctrl2_x = ctrl_r * np.cos(pathway_rad_start)
    ctrl2_y = ctrl_r * np.sin(pathway_rad_start)

    verts.append((ctrl1_x, ctrl1_y))
    codes.append(MplPath.CURVE4)
    verts.append((ctrl2_x, ctrl2_y))
    codes.append(MplPath.CURVE4)
    verts.append((pathway_pts_x[0], pathway_pts_y[0]))
    codes.append(MplPath.CURVE4)

    # Arc along pathway
    for i in range(1, len(pathway_pts_x)):
        verts.append((pathway_pts_x[i], pathway_pts_y[i]))
        codes.append(MplPath.LINETO)

    # Bezier back to gene right edge
    ctrl3_x = ctrl_r * np.cos(pathway_rad_end)
    ctrl3_y = ctrl_r * np.sin(pathway_rad_end)
    ctrl4_x = ctrl_r * np.cos(gene_rad_end)
    ctrl4_y = ctrl_r * np.sin(gene_rad_end)

    verts.append((ctrl3_x, ctrl3_y))
    codes.append(MplPath.CURVE4)
    verts.append((ctrl4_x, ctrl4_y))
    codes.append(MplPath.CURVE4)
    verts.append((gene_pts_x[1], gene_pts_y[1]))
    codes.append(MplPath.CURVE4)

    # Close path
    verts.append((gene_pts_x[0], gene_pts_y[0]))
    codes.append(MplPath.CLOSEPOLY)

    path = MplPath(verts, codes)
    patch = PathPatch(path, facecolor=color, edgecolor=color, alpha=alpha, linewidth=0.3)
    ax.add_patch(patch)


def draw_gene_rectangles(ax, gene_angle, gene_width, logfc_d35, logfc_d65, padj_d35, padj_d65, cmap, norm):
    """Draw stacked rectangles for gene (inner=D35, outer=D65) like reference image."""
    # Position genes at same radius as pathways (0.85-1.0)
    outer_r = 1.0
    rect_height = 0.07

    # D35 rectangle (inner)
    r1_inner = outer_r - rect_height * 2
    r1_outer = outer_r - rect_height

    # D65 rectangle (outer)
    r2_inner = outer_r - rect_height
    r2_outer = outer_r

    # Get colors based on logFC (no significance filtering - show all fold changes)
    if pd.isna(logfc_d35):
        color_d35 = '#d3d3d3'  # Light gray for missing data
    else:
        color_d35 = cmap(norm(logfc_d35))

    if pd.isna(logfc_d65):
        color_d65 = '#d3d3d3'  # Light gray for missing data
    else:
        color_d65 = cmap(norm(logfc_d65))

    # Draw as wedges
    theta1 = gene_angle - gene_width/2
    theta2 = gene_angle + gene_width/2

    # D35 (inner)
    wedge1 = Wedge((0, 0), r1_outer, theta1, theta2, width=rect_height,
                   facecolor=color_d35, edgecolor='white', linewidth=0.3)
    ax.add_patch(wedge1)

    # D65 (outer)
    wedge2 = Wedge((0, 0), r2_outer, theta1, theta2, width=rect_height,
                   facecolor=color_d65, edgecolor='white', linewidth=0.3)
    ax.add_patch(wedge2)


def draw_pathway_arc(ax, start_angle, end_angle, color, label, inner_r=0.86, outer_r=1.0):
    """Draw pathway arc segment with label."""
    # Draw the arc - same radius as gene rectangles
    # Lower alpha (0.6) for pathway arcs to match ribbon transparency
    wedge = Wedge((0, 0), outer_r, start_angle, end_angle, width=outer_r-inner_r,
                  facecolor=color, edgecolor='white', linewidth=1.5, alpha=0.6)
    ax.add_patch(wedge)

    # Add label outside the arc
    mid_angle = (start_angle + end_angle) / 2
    label_r = outer_r + 0.08
    label_rad = np.radians(mid_angle)
    label_x = label_r * np.cos(label_rad)
    label_y = label_r * np.sin(label_rad)

    # Determine text rotation and alignment
    if -90 <= mid_angle <= 90:
        rotation = mid_angle
        ha = 'left'
    else:
        rotation = mid_angle + 180
        ha = 'right'

    ax.text(label_x, label_y, label, fontsize=9, ha=ha, va='center',
            rotation=rotation, rotation_mode='anchor', fontweight='bold',
            color='#333333')


def draw_gochord_diagram(data, output_path, figsize=(14, 12)):
    """Draw GOChord-style chord diagram."""
    pathways = data['pathways']
    genes = data['genes']
    connections = data['connections']
    mutation = data['mutation']

    fig, ax = plt.subplots(figsize=figsize, subplot_kw={'aspect': 'equal'})
    ax.set_xlim(-1.8, 1.8)
    ax.set_ylim(-1.7, 1.7)
    ax.axis('off')

    # Color setup for logFC
    cmap = create_diverging_cmap()
    norm = Normalize(vmin=-3, vmax=3)

    # Layout: genes on left (90-270 degrees), pathways on right (-75 to 75 degrees)
    n_genes = len(genes)
    n_pathways = len(pathways)

    gene_arc_total = 160  # degrees for genes
    pathway_arc_total = 140  # degrees for pathways

    gene_gap = 0.3
    pathway_gap = 3

    # Calculate gene positions
    total_gene_gap = gene_gap * (n_genes - 1)
    gene_width = (gene_arc_total - total_gene_gap) / n_genes if n_genes > 0 else 0

    gene_positions = {}
    current_angle = 100  # Start from top-left

    for gene_info in genes:
        gene = gene_info['name']
        mid_angle = current_angle + gene_width / 2
        gene_positions[gene] = {
            'start': current_angle,
            'end': current_angle + gene_width,
            'mid': mid_angle,
            'logfc_d35': gene_info['logfc_d35'],
            'logfc_d65': gene_info['logfc_d65'],
            'padj_d35': gene_info['padj_d35'],
            'padj_d65': gene_info['padj_d65'],
        }
        current_angle += gene_width + gene_gap

    # Calculate pathway positions
    total_pathway_gap = pathway_gap * (n_pathways - 1) if n_pathways > 1 else 0
    pathway_arc_width = (pathway_arc_total - total_pathway_gap) / n_pathways if n_pathways > 0 else 0

    pathway_positions = {}
    current_angle = -70  # Start from right side

    for i, pathway_info in enumerate(pathways):
        pathway_id = pathway_info['id']
        pathway_positions[pathway_id] = {
            'start': current_angle,
            'end': current_angle + pathway_arc_width,
            'mid': current_angle + pathway_arc_width / 2,
            'color': PATHWAY_COLORS[i % len(PATHWAY_COLORS)],
            'display_name': pathway_info['display_name'],
            'nes_d35': pathway_info['nes_d35'],
            'nes_d65': pathway_info['nes_d65'],
        }
        current_angle += pathway_arc_width + pathway_gap

    # Draw ribbons first (background)
    for gene, pathway_id in connections:
        if gene in gene_positions and pathway_id in pathway_positions:
            gene_pos = gene_positions[gene]
            pathway_pos = pathway_positions[pathway_id]

            draw_ribbon(ax, gene_pos['mid'], gene_width,
                       pathway_pos['start'], pathway_pos['end'],
                       pathway_pos['color'], alpha=0.5)

    # Draw gene rectangles
    for gene, pos in gene_positions.items():
        draw_gene_rectangles(ax, pos['mid'], gene_width,
                            pos['logfc_d35'], pos['logfc_d65'],
                            pos['padj_d35'], pos['padj_d65'],
                            cmap, norm)

        # Gene label - positioned outside the rectangles, radial orientation
        label_r = 1.08  # Outside the outer rectangle
        label_rad = np.radians(pos['mid'])
        label_x = label_r * np.cos(label_rad)
        label_y = label_r * np.sin(label_rad)

        # Radial text: rotate so text reads from center outward
        # For left side of circle (90-270 degrees), text should read right-to-left
        if 90 <= pos['mid'] <= 270:
            rotation = pos['mid'] - 180  # Flip so readable from outside
            ha = 'right'
        else:
            rotation = pos['mid']
            ha = 'left'

        fontsize = 7 if n_genes < 50 else 6
        ax.text(label_x, label_y, gene, fontsize=fontsize, ha=ha, va='center',
                rotation=rotation, rotation_mode='anchor', fontfamily='sans-serif')

    # Draw pathway arcs
    for pathway_id, pos in pathway_positions.items():
        draw_pathway_arc(ax, pos['start'], pos['end'], pos['color'], pos['display_name'])

    # Title
    mutation_color = MUTATION_COLORS.get(mutation, '#333333')
    ax.set_title(f'{mutation} Mutation: Gene-Pathway Leading Edge',
                fontsize=14, fontweight='bold', color=mutation_color, pad=20, y=1.02)

    # Add colorbar legend
    cbar_ax = fig.add_axes([0.02, 0.25, 0.02, 0.25])
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cbar_ax)
    cbar.set_label('log₂ Fold Change', fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    # Add ring legend box
    legend_ax = fig.add_axes([0.02, 0.55, 0.15, 0.15])
    legend_ax.axis('off')
    legend_ax.set_xlim(0, 1)
    legend_ax.set_ylim(0, 1)

    legend_ax.text(0.5, 0.95, 'Gene Rings', fontsize=9, fontweight='bold', ha='center', va='top')
    legend_ax.text(0.5, 0.70, 'Inner = D35', fontsize=8, ha='center', va='center')
    legend_ax.text(0.5, 0.50, 'Outer = D65', fontsize=8, ha='center', va='center')

    # Pathway legend
    pathway_legend_ax = fig.add_axes([0.82, 0.15, 0.16, 0.35])
    pathway_legend_ax.axis('off')
    pathway_legend_ax.set_xlim(0, 1)
    pathway_legend_ax.set_ylim(0, 1)

    pathway_legend_ax.text(0.5, 0.98, 'Pathways', fontsize=9, fontweight='bold', ha='center', va='top')

    y_pos = 0.88
    for i, pathway_info in enumerate(pathways):
        color = PATHWAY_COLORS[i % len(PATHWAY_COLORS)]
        pathway_legend_ax.add_patch(Rectangle((0.0, y_pos - 0.04), 0.08, 0.06,
                                               facecolor=color, edgecolor='white', linewidth=0.5))
        # Shorter name for legend
        short_name = pathway_info['display_name'].replace('\n', ' ')
        if len(short_name) > 25:
            short_name = short_name[:23] + '...'
        pathway_legend_ax.text(0.12, y_pos, short_name, fontsize=6.5, va='center')
        y_pos -= 0.11

    plt.tight_layout()
    plt.savefig(output_path.with_suffix('.pdf'), dpi=300, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.savefig(output_path.with_suffix('.png'), dpi=300, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.close()

    print(f"Saved: {output_path.with_suffix('.pdf')}")
    print(f"Saved: {output_path.with_suffix('.png')}")


def save_chord_data(data, output_path):
    """Save chord data to CSV."""
    gene_df = pd.DataFrame(data['genes'])
    gene_df['mutation'] = data['mutation']

    pathway_df = pd.DataFrame(data['pathways'])
    pathway_df['mutation'] = data['mutation']

    conn_df = pd.DataFrame(data['connections'], columns=['gene', 'pathway'])
    conn_df['mutation'] = data['mutation']

    gene_df.to_csv(output_path.parent / f"chord_genes_{data['mutation']}.csv", index=False)
    pathway_df.to_csv(output_path.parent / f"chord_pathways_{data['mutation']}.csv", index=False)
    conn_df.to_csv(output_path.parent / f"chord_connections_{data['mutation']}.csv", index=False)


def create_readme(output_dir):
    """Create README documentation."""
    readme_content = '''# GOChord-Style Chord Diagrams: Synaptic Ribosome Gene-Pathway Membership

## Overview

These chord diagrams visualize the relationship between synaptic ribosome pathways
and their core (leading-edge) genes from GSEA, styled after the GOplot R package's
GOChord function.

## Figure Interpretation

### Layout
- **Left side (Genes):** Core genes arranged as stacked colored rectangles
  - **Inner rectangle:** D35 fold change (log2FC)
  - **Outer rectangle:** D65 fold change (log2FC)
  - Color scale: Blue = downregulated, White = neutral, Orange = upregulated
  - **Gray rectangles:** FDR > 0.25 (non-significant)

- **Right side (Pathways):** Colored arcs representing enriched pathways
  - Each pathway has a distinct colorblind-friendly color
  - Labels positioned outside the arc

- **Ribbons:** Colored connections indicating that a gene is in the leading edge
  of that pathway (contributes to the enrichment signal)

### Color Scales
- **Gene fold change:** Diverging blue-white-orange scale (range: -3 to +3 log2FC)
- **Pathway colors:** Distinct colorblind-friendly palette per pathway

## Pathways Included

### SynGO (Synaptic Gene Ontology)
- Presynaptic Ribosome (SYNGO:presyn_ribosome)
- Postsynaptic Ribosome (SYNGO:postsyn_ribosome)
- Synapse (GO:0045202)

## Output Files

### Figures
- `chord_diagram_G32A.pdf/.png` - G32A mutation chord diagram
- `chord_diagram_R403C.pdf/.png` - R403C mutation chord diagram

### Data Files
- `chord_genes_*.csv` - Gene information (fold changes, significance)
- `chord_pathways_*.csv` - Pathway information (NES, significance)
- `chord_connections_*.csv` - Gene-pathway membership connections
- `gsea_pathway_data.csv` - Raw GSEA data extracted from RDS files

## Generating Script

**Script:** `02_Analysis/3.7.viz_chord_diagrams.py`

**Usage:**
```bash
python3 02_Analysis/3.7.viz_chord_diagrams.py
```

## Style Reference

Design inspired by GOplot R package GOChord function:
- Wencke Walter, et al. (2015) GOplot: an R package for visually combining
  expression data with functional analysis. Bioinformatics, 31(17), 2912–2914.

---
Generated by: `02_Analysis/3.7.viz_chord_diagrams.py`
'''

    with open(output_dir / 'README.md', 'w') as f:
        f.write(readme_content)


# =============================================================================
# MAIN
# =============================================================================

def main():
    print("=" * 60)
    print("GOChord-Style Chord Diagram Visualization")
    print("=" * 60)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("\n1. Loading GSEA data...")
    gsea_df = load_gsea_data_from_rds()
    print(f"   Loaded {len(gsea_df)} pathway-contrast combinations")

    print("\n2. Loading DE results...")
    de_data = load_de_results()
    print(f"   Loaded DE results for {len(de_data)} contrasts")

    for mutation in ['G32A', 'R403C']:
        print(f"\n3. Processing {mutation} mutation...")

        data = prepare_chord_data(gsea_df, de_data, mutation, max_genes=60)
        print(f"   Found {len(data['pathways'])} pathways")
        print(f"   Using {len(data['genes'])} genes")
        print(f"   Found {len(data['connections'])} connections")

        output_path = OUTPUT_DIR / f'chord_diagram_{mutation}'
        draw_gochord_diagram(data, output_path)
        save_chord_data(data, output_path)

    print("\n4. Creating README...")
    create_readme(OUTPUT_DIR)

    print("\n" + "=" * 60)
    print("Chord diagram generation complete!")
    print(f"Output: {OUTPUT_DIR}")
    print("=" * 60)


if __name__ == '__main__':
    main()
