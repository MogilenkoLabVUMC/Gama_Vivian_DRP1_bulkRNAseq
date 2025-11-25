#!/usr/bin/env python3
"""
UpSet Plot: Ribosome Pathway Gene Overlap
==========================================

Creates an UpSet plot showing gene overlap between ribosome-related pathways from:
- SynGO: Synaptic ribosome pathways (cellular component)
- MitoCarta: Mitochondrial ribosome pathways
- GO: Cytoplasmic ribosome biogenesis pathways

This complements Fig1_Ribosome_Paradox by showing the molecular basis
for why these pathways behave differently.
"""

import sys
from pathlib import Path

# Add module paths
sys.path.insert(0, str(Path(__file__).parent.parent / '01_Scripts'))

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from upsetplot import UpSet, from_contents
import warnings
warnings.filterwarnings('ignore')

from Python.config import resolve_path, ensure_output_dir

# Output directory
OUTPUT_DIR = ensure_output_dir()


def load_mitocarta_genes(pathway_names):
    """
    Load genes from MitoCarta GMX file for specified pathways.

    GMX format: columns are pathways, rows are:
    - Row 0: Short pathway name
    - Row 1: Full pathway name (hierarchy)
    - Rows 2+: Gene symbols
    """
    gmx_file = resolve_path('00_Data/MitoCarta_3.0/MitoPathways3.0.gmx')

    # Read as tab-separated, no header
    df = pd.read_csv(gmx_file, sep='\t', header=None)

    # First row contains pathway names
    pathway_header = df.iloc[0].values

    genes_dict = {}
    for pname in pathway_names:
        # Find column index for this pathway
        matches = [i for i, p in enumerate(pathway_header) if pname.lower() in str(p).lower()]
        if matches:
            col_idx = matches[0]
            # Genes start from row 2 (skip name and hierarchy rows)
            genes = df.iloc[2:, col_idx].dropna().tolist()
            genes = [g for g in genes if isinstance(g, str) and g.strip()]
            genes_dict[pname] = set(genes)
            print(f"  MitoCarta '{pname}': {len(genes)} genes")

    return genes_dict


def load_syngo_genes():
    """
    Load genes from SynGO for ribosome-related cellular component terms.

    SynGO uses:
    - syngo_ontologies.xlsx: Term definitions
    - syngo_annotations.xlsx: Gene-to-term mappings
    """
    syngo_dir = resolve_path('00_Data/SynGO_bulk_20231201')

    # Load ontology terms to find ribosome-related ones
    ontology_file = syngo_dir / 'syngo_ontologies.xlsx'
    annotations_file = syngo_dir / 'syngo_annotations.xlsx'

    ontologies = pd.read_excel(ontology_file, engine='openpyxl')
    annotations = pd.read_excel(annotations_file, engine='openpyxl')

    # Find ribosome-related CC terms
    ribosome_terms = ontologies[
        (ontologies['name'].str.contains('ribosome', case=False, na=False)) &
        (ontologies['domain'] == 'CC')  # Cellular Component
    ]

    print(f"  Found {len(ribosome_terms)} SynGO CC ribosome terms:")
    for _, row in ribosome_terms.iterrows():
        print(f"    - {row['id']}: {row['name']}")

    # Get gene symbols for these terms
    syngo_ribo_genes = set()
    for term_id in ribosome_terms['id']:
        term_genes = annotations[annotations['go_id'] == term_id]['hgnc_symbol'].dropna().unique()
        syngo_ribo_genes.update(term_genes)

    print(f"  SynGO synaptic ribosome genes: {len(syngo_ribo_genes)}")
    return syngo_ribo_genes


def load_go_ribosome_genes():
    """
    Load genes from GO ribosome biogenesis pathways.

    Uses curated list of core ribosome biogenesis genes from GO:0042254
    and related cytoplasmic ribosome terms.
    """
    # Core cytoplasmic ribosome biogenesis genes (from GO:0042254 and related terms)
    # These are genes involved in ribosome assembly in the cytoplasm
    cytoplasmic_ribo_genes = {
        # Large subunit proteins (RPL family)
        'RPL3', 'RPL4', 'RPL5', 'RPL6', 'RPL7', 'RPL7A', 'RPL8', 'RPL9', 'RPL10',
        'RPL10A', 'RPL11', 'RPL12', 'RPL13', 'RPL13A', 'RPL14', 'RPL15', 'RPL17',
        'RPL18', 'RPL18A', 'RPL19', 'RPL21', 'RPL22', 'RPL23', 'RPL23A', 'RPL24',
        'RPL26', 'RPL27', 'RPL27A', 'RPL28', 'RPL29', 'RPL30', 'RPL31', 'RPL32',
        'RPL34', 'RPL35', 'RPL35A', 'RPL36', 'RPL36A', 'RPL37', 'RPL37A', 'RPL38',
        'RPL39', 'RPL40', 'RPL41', 'RPLP0', 'RPLP1', 'RPLP2',
        # Small subunit proteins (RPS family)
        'RPS2', 'RPS3', 'RPS3A', 'RPS4X', 'RPS4Y1', 'RPS5', 'RPS6', 'RPS7', 'RPS8',
        'RPS9', 'RPS10', 'RPS11', 'RPS12', 'RPS13', 'RPS14', 'RPS15', 'RPS15A',
        'RPS16', 'RPS17', 'RPS18', 'RPS19', 'RPS20', 'RPS21', 'RPS23', 'RPS24',
        'RPS25', 'RPS26', 'RPS27', 'RPS27A', 'RPS28', 'RPS29', 'RPSA',
        # Ribosome biogenesis factors
        'BOP1', 'PES1', 'WDR12', 'NOP56', 'NOP58', 'FBL', 'DKC1', 'GAR1', 'NHP2',
        'NOP10', 'NCL', 'NPM1', 'GNL3', 'RRS1', 'RSL1D1', 'NLE1', 'RRP12', 'NOL11',
        'UTP14A', 'UTP14C', 'UTP15', 'UTP18', 'UTP20', 'WDR43', 'WDR3', 'DDX21',
        'DDX27', 'DDX47', 'DDX51', 'DDX52', 'DDX54', 'DDX56', 'DHX37', 'GTPBP4',
        'TSR1', 'TSR2', 'RIOK1', 'RIOK2', 'SBDS', 'EFL1', 'NMD3', 'LSG1', 'MDN1'
    }

    print(f"  GO cytoplasmic ribosome genes (curated): {len(cytoplasmic_ribo_genes)}")
    return cytoplasmic_ribo_genes


def create_upset_plot(gene_sets, output_file):
    """
    Create UpSet plot showing gene overlaps between ribosome pathway categories.
    """
    # Filter out empty sets
    gene_sets = {k: v for k, v in gene_sets.items() if len(v) > 0}

    if len(gene_sets) < 2:
        print("  Warning: Need at least 2 non-empty gene sets for UpSet plot")
        return

    # Create UpSet data structure
    upset_data = from_contents(gene_sets)

    # Create figure
    fig = plt.figure(figsize=(12, 8))

    upset = UpSet(
        upset_data,
        subset_size='count',
        show_counts=True,
        sort_by='cardinality',
        sort_categories_by='cardinality',
        facecolor='steelblue',
        element_size=40
    )

    upset.plot(fig=fig)

    # Add title
    fig.suptitle('Gene Overlap: Ribosome Pathways Across Databases\n'
                 '(Synaptic vs Mitochondrial vs Cytoplasmic)',
                 fontsize=14, fontweight='bold', y=1.02)

    # Calculate overlap statistics
    total_genes = sum(len(s) for s in gene_sets.values())
    unique_genes = len(set.union(*gene_sets.values()))
    overlap_pct = (total_genes - unique_genes) / total_genes * 100 if total_genes > 0 else 0

    # Add interpretation note
    fig.text(0.5, -0.02,
             f'Key insight: Only {overlap_pct:.1f}% overlap between ribosome compartments.\n'
             'This molecular distinctness explains why synaptic and mitochondrial ribosome\n'
             'pathways show opposite trajectory patterns in DRP1 mutations.',
             ha='center', fontsize=10, style='italic',
             bbox=dict(boxstyle='round', facecolor='#F5F5F5', edgecolor='gray', alpha=0.8))

    plt.tight_layout()

    # Save
    for ext, dpi in [('pdf', 300), ('png', 150)]:
        outfile = OUTPUT_DIR / f'{output_file}.{ext}'
        fig.savefig(outfile, dpi=dpi, bbox_inches='tight', facecolor='white')
        print(f"  Saved: {outfile}")

    plt.close(fig)


def main():
    """Generate UpSet plot for ribosome gene overlaps."""
    print("="*80)
    print("GENERATING RIBOSOME GENE OVERLAP UPSET PLOT")
    print("="*80)

    gene_sets = {}

    # 1. Load MitoCarta mitochondrial ribosome genes
    print("\n1. Loading MitoCarta mitochondrial ribosome genes...")
    mito_pathways = ['Mitochondrial_ribosome', 'Mitochondrial_ribosome_assembly',
                     'Translation_factors']
    mito_genes = load_mitocarta_genes(mito_pathways)

    # Combine all MitoCarta ribosome genes
    all_mito_genes = set()
    for genes in mito_genes.values():
        all_mito_genes.update(genes)
    gene_sets['Mitochondrial (MitoCarta)'] = all_mito_genes
    print(f"  Total MitoCarta ribosome genes: {len(all_mito_genes)}")

    # 2. Load SynGO synaptic ribosome genes
    print("\n2. Loading SynGO synaptic ribosome genes...")
    syngo_genes = load_syngo_genes()
    gene_sets['Synaptic (SynGO)'] = syngo_genes

    # 3. Load GO cytoplasmic ribosome genes
    print("\n3. Loading GO cytoplasmic ribosome genes...")
    go_genes = load_go_ribosome_genes()
    gene_sets['Cytoplasmic (GO)'] = go_genes

    # Print summary
    print("\n" + "-"*80)
    print("GENE SET SUMMARY:")
    for name, genes in gene_sets.items():
        print(f"  {name}: {len(genes)} genes")

    # Calculate and print pairwise overlaps
    print("\nPAIRWISE OVERLAPS:")
    names = list(gene_sets.keys())
    for i, n1 in enumerate(names):
        for n2 in names[i+1:]:
            overlap = gene_sets[n1] & gene_sets[n2]
            print(f"  {n1} & {n2}: {len(overlap)} genes")
            if overlap:
                print(f"    Genes: {', '.join(sorted(overlap)[:10])}{'...' if len(overlap) > 10 else ''}")

    # Triple overlap
    if len(names) >= 3:
        triple_overlap = gene_sets[names[0]] & gene_sets[names[1]] & gene_sets[names[2]]
        print(f"  All three: {len(triple_overlap)} genes")
        if triple_overlap:
            print(f"    Genes: {', '.join(sorted(triple_overlap))}")

    # Create plot
    print("\n" + "-"*80)
    print("Creating UpSet plot...")
    create_upset_plot(gene_sets, 'Fig1b_Ribosome_Gene_Overlap_UpSet')

    print("\n" + "="*80)
    print("COMPLETE!")
    print("="*80)


if __name__ == '__main__':
    main()

