"""
DRP1 bulk RNA-seq project configuration.

Centralized configuration for paths, thresholds, and column definitions
used across all analysis scripts.
"""

from pathlib import Path

# Base paths (relative to project root)
CONFIG = {
    # Data directories
    'data_dir': Path('03_Results/02_Analysis/Python_exports'),
    'classified_data': Path('03_Results/02_Analysis/Plots/Cross_database_validation/pathways_classified.csv'),
    'output_dir': Path('03_Results/02_Analysis/Plots/Publication_Figures'),

    # Analysis thresholds
    'padj_cutoff': 0.05,           # Significance cutoff
    'nes_threshold': 1.5,          # Minimum |NES| for defect
    'low_threshold': 0.5,          # Threshold for "no effect"
    'min_data_points': 3,          # Minimum non-NA values across 6 columns
    'max_per_category': 10,        # Display limit per semantic category

    # Figure settings
    'dpi': 300,
    'vmax': 3.5,                   # Max |NES| for colorscale
    'nonsig_style': 'show_values', # 'hatching' or 'show_values' for non-sig cells

    # Databases to EXCLUDE from analyses
    # CGP: Cancer-focused pathways, not relevant for neuronal biology
    # Canon: Redundant with other curated pathway sets
    'excluded_databases': ['cgp', 'canon'],

    # Trajectory column definitions
    # These map the contrast names to the trajectory framework
    'nes_cols': [
        'NES_Early_G32A', 'NES_TrajDev_G32A', 'NES_Late_G32A',
        'NES_Early_R403C', 'NES_TrajDev_R403C', 'NES_Late_R403C'
    ],
    'padj_cols': [
        'p.adjust_Early_G32A', 'p.adjust_TrajDev_G32A', 'p.adjust_Late_G32A',
        'p.adjust_Early_R403C', 'p.adjust_TrajDev_R403C', 'p.adjust_Late_R403C'
    ],

    # Contrast mapping for trajectory framework
    # Maps original contrast names to trajectory stages
    'contrast_mapping': {
        'G32A_vs_Ctrl_D35': 'Early_G32A',
        'Maturation_G32A_specific': 'TrajDev_G32A',
        'G32A_vs_Ctrl_D65': 'Late_G32A',
        'R403C_vs_Ctrl_D35': 'Early_R403C',
        'Maturation_R403C_specific': 'TrajDev_R403C',
        'R403C_vs_Ctrl_D65': 'Late_R403C',
    }
}

# Pathways to EXCLUDE from analyses
# These are not relevant to iPSC-derived cortical neurons
EXCLUDE_PATHWAYS = [
    # Viral pathways (not relevant to iPSC neurons)
    'REACTOME_SARS_COV_1_MODULATES_HOST_TRANSLATION_MACHINERY',
    'REACTOME_EXPORT_OF_VIRAL_RIBONUCLEOPROTEINS_FROM_NUCLEUS',
    'REACTOME_SARS_COV_2_MODULATES_HOST_TRANSLATION_MACHINERY',

    # Cancer pathways
    'WP_RETINOBLASTOMA_GENE_IN_CANCER',

    # Cardiac pathways (not relevant to cortical neurons)
    'WP_CALCIUM_REGULATION_IN_CARDIAC_CELLS',
    'REACTOME_RESPONSE_TO_ELEVATED_PLATELET_CYTOSOLIC_CA2',

    # Sperm/cilium pathways (not cortical neuron relevant)
    'GOBP_SPERM_AXONEME_ASSEMBLY',
    'GOBP_CILIUM_MOVEMENT',
    'GOBP_AXONEME_ASSEMBLY',
    'GOBP_CILIUM_ASSEMBLY',
    'GOBP_CILIUM_ORGANIZATION',
    'GOBP_MOTILE_CILIUM_ASSEMBLY',
    'GOBP_SPERM_MOTILITY',
    'GOCC_AXONEMAL_DYNEIN_COMPLEX',
    'GOCC_AXONEMAL_DOUBLET_MICROTUBULE',
    'GOCC_CILIARY_PLASM',
    'GOCC_MOTILE_CILIUM',
    'GOCC_9_PLUS_2_MOTILE_CILIUM',

    # ECM pathways (not neuronal)
    'GOCC_COLLAGEN_CONTAINING_EXTRACELLULAR_MATRIX',

    # Immune pathways (not relevant to iPSC neurons)
    'GOBP_ALPHA_BETA_T_CELL_ACTIVATION',
    'GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION',
    'GOBP_T_CELL_ACTIVATION',
    'GOBP_IMMUNE_RESPONSE',
]


def get_project_root():
    """Get the project root directory."""
    # Walk up from this file to find the project root
    current = Path(__file__).resolve()
    for parent in current.parents:
        if (parent / 'CLAUDE.md').exists() or (parent / '02_Analysis').exists():
            return parent
    return Path.cwd()


def resolve_path(relative_path):
    """Resolve a path relative to project root."""
    root = get_project_root()
    if isinstance(relative_path, str):
        relative_path = Path(relative_path)
    return root / relative_path


def ensure_output_dir():
    """Ensure output directory exists."""
    output_dir = resolve_path(CONFIG['output_dir'])
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir
