"""
Semantic category definitions for DRP1 mutation analysis.

Defines biologically meaningful categories for organizing GSEA pathways,
with special focus on:
- Neuronal function and development
- Mitochondrial dynamics (DRP1's core function)
- Calcium signaling
- Translation machinery
"""

import pandas as pd

# Semantic category order (from neuronal to generic)
# Order reflects biological priority for DRP1 mutation analysis
# Note: Cytoskeletal removed as not relevant to DRP1 mutation narrative
SEMANTIC_CATEGORY_ORDER = [
    'Synapse',                       # SynGO synaptic pathways
    'Neuronal Development',          # Axon/synapse development
    'Mitochondrial Dynamics',        # Fission/fusion/trafficking (DRP1 core)
    'Electron Transport Chain',      # Complexes I-IV
    'ATP Synthase (Complex V)',      # OXPHOS Complex V
    'Mitochondrial Metabolism',      # TCA cycle, fatty acid oxidation
    'Mitochondrial Function',        # General mitochondrial
    'Mitochondrial Ribosome',        # Mito ribosome structure
    'Mitochondrial Translation',     # Mito protein synthesis
    'Ribosome Biogenesis',           # Ribosome assembly (not structure)
    'Cytoplasmic Ribosome',          # Cytosolic ribosome structure
    'Cytoplasmic Translation',       # Cytosolic protein synthesis
    'Calcium Signaling',             # Ca2+ homeostasis
    'Other'                          # Catch-all (renamed from RNA Processing)
]

# Colors for semantic category labels (for figure annotations)
# Note: Cytoskeletal removed as not relevant to DRP1 mutation narrative
SEMANTIC_COLORS = {
    'Synapse': '#8B4513',                     # Saddle brown
    'Neuronal Development': '#A0522D',        # Sienna
    'Mitochondrial Dynamics': '#006400',      # Dark green (DRP1 core!)
    'Electron Transport Chain': '#2E8B57',    # Sea green
    'ATP Synthase (Complex V)': '#DAA520',    # Goldenrod
    'Mitochondrial Metabolism': '#228B22',    # Forest green
    'Mitochondrial Function': '#32CD32',      # Lime green
    'Mitochondrial Ribosome': '#4169E1',      # Royal blue
    'Mitochondrial Translation': '#1E90FF',   # Dodger blue
    'Ribosome Biogenesis': '#6A5ACD',         # Slate blue
    'Cytoplasmic Ribosome': '#9932CC',        # Dark orchid
    'Cytoplasmic Translation': '#BA55D3',     # Medium orchid
    'Calcium Signaling': '#DC143C',           # Crimson
    'Other': '#696969',                       # Dim gray
}

# Pattern colors for trajectory classification
# IMPORTANT: Import from canonical source to maintain single source of truth
from .pattern_definitions import get_pattern_colors
PATTERN_COLORS = get_pattern_colors()  # Re-export for backward compatibility

# Mutation colors - import from unified color config (single source of truth)
# Note: For heatmap annotations, use HEATMAP_ANNOTATION_COLORS from color_config.py
from .color_config import MUTATION_COLORS  # Re-export for backward compatibility

# =============================================================================
# BIOLOGICAL RELEVANCE FILTERING FOR HIGHLIGHT SELECTION
# =============================================================================
# These constants control which pathways get highlighted/labeled in figures.
# Context: iPSC-derived cortical neurons with DRP1 mutations affecting
# mitochondrial dynamics, translation, and synaptic development.

# Keywords to EXCLUDE from highlighting (not relevant to cortical neurons)
# Case-insensitive matching against pathway Description
EXCLUDE_FROM_HIGHLIGHT_KEYWORDS = [
    # Developmental/tissue-specific (wrong cell type)
    'MATERNAL', 'ZYGOTIC', 'EMBRYONIC_AXIS',
    'CARDIOMYOCYTE', 'CARDIAC_MUSCLE', 'HEART_FIELD', 'CARDIAC_CHAMBER',
    'SKELETAL_MUSCLE', 'STRIATED_MUSCLE', 'MUSCLE_CONTRACTION', 'MYOFIBRIL',
    'OSTEOBLAST', 'OSTEOCLAST', 'BONE_REMODEL', 'CHONDROCYTE',
    'NEPHRON', 'KIDNEY', 'RENAL', 'GLOMERUL',
    'HEPATOCYTE', 'LIVER', 'HEPAT', 'BILE',
    'ADIPOCYTE', 'ADIPOGEN', 'FAT_CELL',
    'PANCREA', 'ISLET', 'INSULIN_SECRET',
    'SPERM', 'SPERMAT', 'TESTIS', 'OVARY', 'OOCYTE',
    # Cancer-specific (disease context mismatch)
    'CANCER', 'TUMOR', 'METASTA', 'ONCOGEN', 'CARCINOMA',
    'LEUKEMIA', 'LYMPHOMA', 'MELANOMA', 'GLIOBLASTOMA',
    # Immune-specific (not primary focus)
    'T_CELL_RECEPTOR', 'B_CELL_RECEPTOR', 'ANTIBODY', 'IMMUNOGLOBULIN',
    'ANTIGEN_PRESENT', 'MHC_CLASS',
    # Plant/non-mammalian
    'PLANT', 'CHLOROPLAST', 'PHOTOSYNTH',
]

# Priority semantic categories for DRP1 neuronal analysis
# Pathways in these categories get priority for highlighting
PRIORITY_CATEGORIES_FOR_HIGHLIGHT = [
    'Synapse',                       # Core: synaptic function
    'Neuronal Development',          # Core: neuronal maturation
    'Mitochondrial Dynamics',        # Core: DRP1's primary function
    'Mitochondrial Ribosome',        # Core: mito translation machinery
    'Mitochondrial Translation',     # Core: mito protein synthesis
    'Electron Transport Chain',      # Core: OXPHOS complexes I-IV
    'ATP Synthase (Complex V)',      # Core: ATP production
    'Cytoplasmic Ribosome',          # Important: cytosolic translation
    'Cytoplasmic Translation',       # Important: protein synthesis
    'Ribosome Biogenesis',           # Important: ribosome assembly
    'Calcium Signaling',             # Important: Ca2+ homeostasis
    'Mitochondrial Metabolism',      # Relevant: TCA, fatty acid oxidation
    'Mitochondrial Function',        # Relevant: general mito
]

# Maximum pathways from non-priority categories (fallback)
MAX_OTHER_CATEGORY_HIGHLIGHTS = 2


def is_relevant_for_highlight(description):
    """
    Check if a pathway is biologically relevant for highlighting.

    Parameters
    ----------
    description : str
        Pathway description/name

    Returns
    -------
    bool
        True if pathway should be considered for highlighting
    """
    if not description:
        return False

    desc_upper = str(description).upper()

    # Check against exclusion keywords
    for keyword in EXCLUDE_FROM_HIGHLIGHT_KEYWORDS:
        if keyword in desc_upper:
            return False

    return True


def get_highlight_priority(semantic_category):
    """
    Get priority score for a semantic category (lower = higher priority).

    Parameters
    ----------
    semantic_category : str
        The semantic category from assign_semantic_category()

    Returns
    -------
    int
        Priority score (0 = highest priority, 99 = lowest/Other)
    """
    if semantic_category in PRIORITY_CATEGORIES_FOR_HIGHLIGHT:
        return PRIORITY_CATEGORIES_FOR_HIGHLIGHT.index(semantic_category)
    return 99  # 'Other' or unknown categories


def assign_semantic_category(row, exclude_pathways=None):
    """
    Assign semantic category based on pathway description and database.

    Priority order matters - more specific categories checked first.
    Excluded pathways return None and should be filtered out.

    Parameters
    ----------
    row : pd.Series
        Row containing 'Description' and 'database' columns
    exclude_pathways : list, optional
        List of pathway names to exclude (return None)

    Returns
    -------
    str or None
        Semantic category name, or None if pathway should be excluded
    """
    desc_raw = str(row.get('Description', ''))
    desc = desc_raw.lower()
    db = row.get('database', '')

    # 0. EXCLUSION CHECK FIRST
    if exclude_pathways and desc_raw in exclude_pathways:
        return None  # Will be filtered out

    # Check for irrelevant pathway keywords (negative filter)
    # Cilium/sperm pathways are not relevant to cortical neurons
    # Note: 'axonem' matches both 'axoneme' and 'axonemal'
    irrelevant_keywords = ['cilium', 'sperm', 'axonem', 'flagell', 'motile cilia']
    if any(kw in desc for kw in irrelevant_keywords):
        return None

    # Cytoskeletal - EXCLUDED (not relevant to DRP1 mutation narrative)
    # These pathways are filtered out by returning None
    # Moved to top to ensure complete removal
    cytoskel_keywords = ['microtubule', 'actin', 'cytoskeleton', 'tubulin',
                        'microfilament', 'intermediate filament', 'spectrin']
    if any(kw in desc for kw in cytoskel_keywords):
        return None  # Exclude cytoskeletal pathways

    # Priority order matters!

    # 1. SynGO pathways are always Synapse (synaptic function)
    if db == 'SynGO':
        return 'Synapse'

    # 2. Mitochondrial Dynamics (check BEFORE general mito!)
    # This captures DRP1's core function: fission/fusion/trafficking
    mito_dynamics_keywords = [
        'fission', 'fusion', 'dynamin', 'drp1', 'dnm1l',
        'opa1', 'mfn1', 'mfn2', 'fis1', 'mff',
        'mitochondrial dynamics', 'mitochondrial morphology',
        'miro', 'rhot', 'trafficking'
    ]
    if any(kw in desc for kw in mito_dynamics_keywords):
        if 'mitochondri' in desc or db == 'MitoCarta':
            return 'Mitochondrial Dynamics'

    # 3. Neuronal Development
    # Captures axonogenesis, synaptogenesis, neurite development
    neuro_dev_keywords = [
        'axon', 'axonal', 'axonogenesis', 'axon guidance',
        'dendrit', 'dendritic morphogenesis',
        'synaptogenesis', 'synapse formation', 'synapse maturation',
        'neurogenesis', 'neuron differentiation', 'neuron development',
        'neurite', 'neurite outgrowth', 'growth cone',
        'neurodevelopment', 'neuronal development', 'neural development'
    ]
    if any(kw in desc for kw in neuro_dev_keywords):
        return 'Neuronal Development'

    # 4. (Removed - Cytoskeletal check moved to top)

    # 5. Electron Transport Chain (Complexes I-IV)
    etc_keywords = ['complex i', 'complex ii', 'complex iii', 'complex iv',
                   'nadh dehydrogenase', 'nadh-ubiquinone', 'succinate dehydrogenase',
                   'ubiquinol-cytochrome', 'cytochrome c oxidase', 'cytochrome bc1',
                   'respiratory chain complex', 'electron transport chain']
    if any(kw in desc for kw in etc_keywords):
        return 'Electron Transport Chain'

    # 6. ATP Synthase / Complex V
    if ('atp' in desc and ('synth' in desc or 'complex' in desc)) or 'complex v' in desc:
        return 'ATP Synthase (Complex V)'

    # 7. Ribosome Biogenesis (assembly, NOT structure)
    # Must check BEFORE ribosome structure categories
    biogenesis_keywords = ['ribosome biogenesis', 'ribosomal subunit assembly',
                          'ribosome assembly', 'preribosome', 'rrna processing',
                          'ribosomal rna', 'rrna metabol', 'gobp_ribosome_biogenesis']
    if any(kw in desc for kw in biogenesis_keywords):
        return 'Ribosome Biogenesis'

    # 8. Mitochondrial Ribosome (structure)
    if 'ribosome' in desc or 'ribosomal' in desc:
        if 'mitochondri' in desc or db == 'MitoCarta':
            return 'Mitochondrial Ribosome'

    # 9. Mitochondrial Translation
    if 'translation' in desc:
        if 'mitochondri' in desc or db == 'MitoCarta':
            return 'Mitochondrial Translation'

    # 10. Cytoplasmic Ribosome (non-mitochondrial, structure)
    if 'ribosome' in desc or 'ribosomal' in desc:
        if 'mitochondri' not in desc:
            return 'Cytoplasmic Ribosome'

    # 11. Cytoplasmic Translation (non-mitochondrial)
    if 'translation' in desc:
        if 'mitochondri' not in desc:
            return 'Cytoplasmic Translation'

    # 12. Calcium Signaling
    # Exclude "calcium-independent" pathways (they're actually NOT about calcium signaling)
    if 'calcium' in desc and 'independent' in desc:
        pass  # Skip calcium-independent pathways, they'll go to Other
    else:
        calcium_keywords = ['calcium', 'ca2+', 'calmodulin', 'calcineurin',
                           'ryanodine', 'serca', 'store-operated', 'cytosolic_ca2']
        if any(kw in desc for kw in calcium_keywords):
            return 'Calcium Signaling'

    # 13. Mitochondrial Metabolism (TCA, fatty acid, etc.)
    # Note: electron transport now has its own category
    mito_metab_keywords = ['tca', 'citrate', 'krebs', 'fatty acid',
                          'beta-oxidation', 'oxidative phosphorylation', 'oxphos']
    if any(kw in desc for kw in mito_metab_keywords):
        if 'mitochondri' in desc or db == 'MitoCarta':
            return 'Mitochondrial Metabolism'

    # 14. General Mitochondrial Function
    if 'mitochondri' in desc or db == 'MitoCarta':
        return 'Mitochondrial Function'

    # 15. Default: Other (catch-all)
    return 'Other'


def categorize_pathways(df, category_col='Semantic_Category', exclude_pathways=None):
    """
    Add semantic category column to pathway dataframe.

    Pathways that should be excluded (return None from assign_semantic_category)
    are removed from the output dataframe.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with 'Description' and 'database' columns
    category_col : str
        Name for the new category column
    exclude_pathways : list, optional
        List of pathway names to exclude

    Returns
    -------
    pd.DataFrame
        DataFrame with added category column, excluded pathways removed
    """
    df = df.copy()
    df[category_col] = df.apply(
        lambda row: assign_semantic_category(row, exclude_pathways), axis=1
    )
    # Filter out excluded pathways (those with None category)
    df = df[df[category_col].notna()].copy()
    return df


def deduplicate_pathways(df, nes_cols=None):
    """
    Remove duplicate pathways, keeping the one with highest signal.

    When the same biological pathway appears from multiple databases
    (e.g., GOBP_MITOCHONDRIAL_TRANSLATION and MitoCarta Translation),
    keep only the one with the highest max |NES| across trajectory stages.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with 'Description' column and NES columns
    nes_cols : list, optional
        List of NES column names. If None, looks for columns starting with 'NES_'

    Returns
    -------
    pd.DataFrame
        Deduplicated DataFrame
    """
    import re

    df = df.copy()

    # Find NES columns if not specified
    if nes_cols is None:
        nes_cols = [c for c in df.columns if c.startswith('NES_')]

    if not nes_cols:
        return df

    # Create normalized pathway name by removing database prefixes
    prefix_pattern = r'^(GOBP_|GOCC_|GOMF_|REACTOME_|KEGG_|WP_|HALLMARK_)'
    df['_pathway_base'] = df['Description'].str.replace(prefix_pattern, '', regex=True)
    df['_pathway_base'] = df['_pathway_base'].str.lower()

    # Calculate max |NES| for ranking
    df['_max_nes'] = df[nes_cols].abs().max(axis=1)

    # Sort by max NES descending, then deduplicate
    df = df.sort_values('_max_nes', ascending=False)
    df = df.drop_duplicates(subset=['_pathway_base'], keep='first')

    # Clean up temporary columns
    df = df.drop(columns=['_pathway_base', '_max_nes'])

    return df


def get_category_order_map():
    """Get mapping of category to sort order."""
    return {cat: i for i, cat in enumerate(SEMANTIC_CATEGORY_ORDER)}


def sort_by_category(df, nes_col='max_nes', category_col='Semantic_Category'):
    """
    Sort dataframe by semantic category, then by NES magnitude.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with category column
    nes_col : str
        Column name for NES-based secondary sorting
    category_col : str
        Column name for category

    Returns
    -------
    pd.DataFrame
        Sorted DataFrame
    """
    df = df.copy()
    order_map = get_category_order_map()
    df['_category_order'] = df[category_col].map(order_map).fillna(len(order_map))

    sort_cols = ['_category_order']
    ascending = [True]

    if nes_col in df.columns:
        sort_cols.append(nes_col)
        ascending.append(False)  # Higher |NES| first within category

    df = df.sort_values(sort_cols, ascending=ascending)
    df = df.drop(columns=['_category_order'])

    return df
