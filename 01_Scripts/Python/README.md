# Python Modules for DRP1 Bulk RNA-seq Analysis

This directory contains core Python modules used across the analysis pipeline, providing centralized configuration, pattern classification, data loading, and visualization utilities.

## Module Overview

### Core Modules

| Module | Purpose | Status |
|--------|---------|--------|
| `pattern_definitions.py` | **Canonical pattern classification** (single source of truth) | Primary |
| `config.py` | Centralized paths, thresholds, and project configuration | Primary |
| `data_loader.py` | Data loading with database filtering | Primary |
| `semantic_categories.py` | Biological pathway categorization and colors | Primary |
| `color_config.py` | Unified color schemes (Python mirror of R colors) | Primary |
| `viz_bump_charts.py` | Bump chart visualization module | Primary |
| `patterns.py` | Legacy wrapper (redirects to pattern_definitions.py) | Deprecated |
| `__init__.py` | Package initialization and exports | Support |

---

## 1. pattern_definitions.py

**CANONICAL SOURCE OF TRUTH for pattern classification**

### Purpose
Defines the trajectory pattern classification system used to categorize how pathways respond to DRP1 mutations across neuronal maturation. This module is the **single authoritative source** for pattern logic, synchronized with the R implementation in `02_Analysis/1.7.create_master_gsva_table.R`.

### Trajectory Framework

```
Early:   Mutation vs Control at D35 (initial state)
TrajDev: Mutation-specific maturation trajectory deviation
         = (D65_mut - D35_mut) - (D65_ctrl - D35_ctrl)
Late:    Mutation vs Control at D65 (mature state)
```

**Key Insight:** TrajDev represents active transcriptional plasticity, not passive developmental changes.

### Pattern Taxonomy (8 Patterns)

| Pattern | Type | Description | Criteria |
|---------|------|-------------|----------|
| **Compensation** | Active | Adaptive response opposing defect | TrajDev sig + opposes Early + improved |
| **Sign_reversal** | Active | Trajectory reversal (sign flipped) | TrajDev sig + opposes Early + Late opposite sign |
| **Progressive** | Active | Active worsening | TrajDev sig + amplifies Early |
| **Natural_improvement** | Passive | Recovery without active TrajDev | TrajDev NS + improved |
| **Natural_worsening** | Passive | Deterioration without active TrajDev | TrajDev NS + worsened |
| **Late_onset** | Emerging | Maturation-dependent dysfunction | Early NS + Late sig |
| **Transient** | Resolving | Developmental delay (fully resolved) | Early strong + Late NS |
| **Complex** | Mixed | Multiphasic or inconsistent | Doesn't fit above |

**Complex Subtypes:**
- `Complex_opposing`: TrajDev opposes Early but insufficient magnitude
- `Complex_amplifying`: TrajDev amplifies Early but insufficient magnitude
- `Complex_multiphasic`: Inconsistent directionality across stages

### Super-Categories (Simplified for Main Figures)

```python
SUPER_CATEGORY_MAP = {
    'Compensation': 'Active_Compensation',
    'Sign_reversal': 'Active_Reversal',
    'Progressive': 'Active_Progression',
    'Natural_improvement': 'Passive',
    'Natural_worsening': 'Passive',
    'Late_onset': 'Late_onset',
    'Transient': 'Other',
    'Complex': 'Other'
}
```

### Thresholds

```python
# Significance
PADJ_SIGNIFICANT = 0.05      # p.adjust for "significant"
PADJ_TRENDING = 0.10         # p.adjust for "trending"

# Effect size (GSEA NES)
NES_EFFECT = 0.5             # Minimum |NES| for any effect
NES_STRONG = 1.0             # Strong defect threshold

# Improvement/worsening
IMPROVEMENT_RATIO = 0.7      # |Late|/|Early| < 0.7 = improved (≥30% reduction)
WORSENING_RATIO = 1.3        # |Late|/|Early| > 1.3 = worsened (≥30% increase)

# GSVA-specific (scaled for -1 to 1 range)
GSVA_EFFECT = 0.15           # Equivalent to NES 0.5
GSVA_STRONG = 0.30           # Equivalent to NES 1.0
```

### Key Functions

```python
# Single pathway classification
classify_pattern(early_nes, early_padj, trajdev_nes, trajdev_padj, late_nes, late_padj)
# Returns: (pattern, confidence)

# Batch classification
add_pattern_classification(df, mutations=['G32A', 'R403C'])
# Returns: df with Pattern_*, Confidence_*, Super_Category_* columns

# Add super-categories
add_super_category_columns(df, mutations=['G32A', 'R403C'])

# Pattern statistics
get_pattern_summary(df, mutation='G32A')
get_strict_counts(df, mutation='G32A')      # Strict definition (High + Medium confidence)
get_potential_counts(df, mutation='G32A')   # Potential (includes Low confidence)
```

### Usage by Downstream Scripts

**Master table generation:**
- `02_Analysis/1.5.create_master_pathway_table.py` - Applies pattern classification to all GSEA results

**Visualization scripts:**
- `02_Analysis/3.1.publication_figures.py` - Main publication figures
- `02_Analysis/3.2.publication_figures_dotplot.py` - Dotplot variants
- `02_Analysis/3.4.pattern_summary_normalized.py` - Pattern distribution plots
- `02_Analysis/3.5.viz_trajectory_flow.py` - Trajectory flow diagrams
- `02_Analysis/3.7.viz_bump_chart.py` - Bump charts (via viz_bump_charts module)
- `02_Analysis/3.8.viz_interactive_bump*.py` - Interactive bump chart viewers
- `02_Analysis/Supp4.sensitivity_analysis.py` - Threshold sensitivity analysis
- `02_Analysis/Supp6.app_bump_chart_explorer.py` - Streamlit bump chart explorer

### Synchronization with R Implementation

**CRITICAL:** Pattern logic must stay synchronized with R implementation:
- **R file:** `02_Analysis/1.7.create_master_gsva_table.R` (lines 348-486, 527-543)
- **Thresholds:** Lines 47-60 (Python) ↔ Lines 348-354 (R)
- **Classification logic:** Lines 276-411 (Python) ↔ Lines 379-486 (R)
- **Super-categories:** Lines 80-90 (Python) ↔ Lines 527-543 (R)

**Key difference:**
- **Python (GSEA):** Uses p-values for TrajDev significance (p.adjust < 0.05)
- **R (GSVA):** Uses magnitude threshold (|TrajDev| > GSVA_EFFECT)
- **Reason:** GSVA TrajDev is calculated (difference of differences), not statistically tested

**Verification checklist:**
```bash
# Regenerate master tables
python3 02_Analysis/1.5.create_master_pathway_table.py
Rscript 02_Analysis/1.7.create_master_gsva_table.R

# Compare pattern distributions (should be consistent)
# Verify thresholds match exactly in both files
```

---

## 2. config.py

**Centralized project configuration**

### Purpose
Single source for file paths, thresholds, and analysis parameters used across all Python scripts.

### Key Configuration

```python
CONFIG = {
    # Data paths
    'master_gsea_table': Path('03_Results/02_Analysis/master_gsea_table.csv'),
    'output_dir': Path('03_Results/02_Analysis/Plots/Publication_Figures'),

    # Analysis thresholds
    'padj_cutoff': 0.05,           # Significance cutoff
    'nes_threshold': 1.5,          # Minimum |NES| for defect
    'low_threshold': 0.5,          # Threshold for "no effect"

    # Figure settings
    'dpi': 300,
    'vmax': 3.5,                   # Max |NES| for colorscale
    'nonsig_style': 'show_values', # 'hatching' or 'show_values'

    # Database exclusions
    'excluded_databases': ['cgp', 'canon'],  # Not relevant for neurons

    # Contrast mapping to trajectory framework
    'contrast_mapping': {
        'G32A_vs_Ctrl_D35': 'Early_G32A',
        'Maturation_G32A_specific': 'TrajDev_G32A',
        'G32A_vs_Ctrl_D65': 'Late_G32A',
        'R403C_vs_Ctrl_D35': 'Early_R403C',
        'Maturation_R403C_specific': 'TrajDev_R403C',
        'R403C_vs_Ctrl_D65': 'Late_R403C',
    }
}
```

### Pathway Exclusions

`EXCLUDE_PATHWAYS` list filters out pathways irrelevant to iPSC-derived cortical neurons:
- Viral pathways (SARS-CoV-2, etc.)
- Cancer pathways
- Cardiac-specific pathways
- Sperm/cilium pathways
- ECM pathways
- Immune pathways

### Utility Functions

```python
get_project_root()           # Auto-detect project root
resolve_path(relative_path)  # Resolve paths relative to project root
ensure_output_dir()          # Create output directory if missing
```

---

## 3. data_loader.py

**Data loading utilities with filtering**

### Purpose
Provides standardized functions to load GSEA results and pathway classifications with automatic database filtering.

### Key Functions

```python
# Load master GSEA table with pattern classifications
load_classified_pathways(exclude_databases=None, min_data_points=None)
# Returns: DataFrame with Pattern_*, Confidence_*, Super_Category_* columns

# Filter pathways by relevance and quality
filter_pathways(df, min_effect=0.5, exclude_patterns=None)
```

### Features

- **Automatic database exclusion:** Removes CGP and canon by default
- **Pattern classification integration:** Loads pre-classified patterns from master table
- **p.adjust merging:** Integrates significance values from wide-format exports

---

## 4. semantic_categories.py

**Biological pathway categorization and visualization colors**

### Purpose
Defines biologically meaningful categories for organizing GSEA pathways with focus on neuronal function, mitochondrial dynamics, calcium signaling, and translation machinery.

### Semantic Category Hierarchy

```python
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
    'Ribosome Biogenesis',           # Ribosome assembly
    'Cytoplasmic Ribosome',          # Cytosolic ribosome structure
    'Cytoplasmic Translation',       # Cytosolic protein synthesis
    'Calcium Signaling',             # Ca2+ homeostasis
    'Other'                          # Catch-all
]
```

### Color Schemes

```python
SEMANTIC_COLORS = {...}              # Colors for semantic categories
PATTERN_COLORS = {...}               # Re-exported from pattern_definitions.py
MUTATION_COLORS = {...}              # Re-exported from color_config.py
```

### Biological Relevance Filtering

`EXCLUDE_FROM_HIGHLIGHT_KEYWORDS` list defines pathways to exclude from figure highlights:
- Wrong cell type (cardiac, skeletal muscle, kidney, liver)
- Developmental (maternal, zygotic)
- Immune pathways
- Viral/cancer pathways

### Key Functions

```python
assign_semantic_category(pathway_description)
# Returns: category name from SEMANTIC_CATEGORY_ORDER

is_relevant_for_highlight(pathway_description)
# Returns: True if pathway is relevant for cortical neurons

get_highlight_priority(pathway_description, semantic_category)
# Returns: Priority score for pathway highlighting
```

---

## 5. color_config.py

**Unified color configuration (Python mirror of R colors)**

### Purpose
Single source of truth for color schemes used across all Python visualizations. Mirrors the R color definitions in `01_Scripts/R_scripts/color_config.R`.

### Color System Overview

1. **DIVERGING_COLORS:** Blue-White-Orange for NES/logFC heatmaps
   - Use for continuous diverging metrics (NES, log2FC)
   - **Do NOT use for categorical annotations**

2. **HEATMAP_ANNOTATION_COLORS:** Distinct colors for heatmap annotations
   - Use for row/column annotations on diverging heatmaps
   - Avoids blue/orange to prevent confusion with gradient

3. **MUTATION_COLORS:** Standard mutation colors for general use

4. **PATTERN_COLORS:** Pattern-specific colors (imported from pattern_definitions.py)

### Diverging Color Palette

```python
DIVERGING_COLORS = {
    'negative': '#2166AC',  # Blue
    'neutral': '#F7F7F7',   # White
    'positive': '#B35806'   # Orange/Brown
}

# Create matplotlib colormap
create_diverging_cmap(name='BlueWhiteOrange')
```

### Heatmap Annotation Colors

```python
HEATMAP_ANNOTATION_COLORS = {
    'mutation': {
        'G32A': '#7B68EE',     # Medium Slate Blue (purple-ish)
        'R403C': '#DC143C',    # Crimson (red)
    },
    'timepoint': {
        'D35': '#2E8B57',      # Sea green
        'D65': '#8B4513',      # Saddle brown
    },
    'pattern': {...}           # Re-exported from pattern_definitions.py
}
```

### Mutation Colors (General Use)

```python
MUTATION_COLORS = {
    'G32A': '#1f77b4',    # Blue (GTPase domain)
    'R403C': '#ff7f0e',   # Orange (Stalk domain)
    'Ctrl': '#7f7f7f'     # Gray
}
```

**Design:** Colorblind-safe palette using deuteranopia/protanopia-safe colors.

---

## 6. viz_bump_charts.py

**Bump chart (slope graph) visualization module**

### Purpose
Core logic for generating bump charts of pathway trajectories (Early → TrajDev → Late). Refactored for modularity, configurability, and clean separation of concerns.

### Architecture

```
BumpChartConfig       → Configuration container (all tunable parameters)
BumpChartStyler       → Line colors, weights, alphas, zorders
BumpChartHighlighter  → Logic for selecting pathways to label
BumpChartRenderer     → Core plotting logic (lines, curves, labels)
BumpChartBuilder      → High-level interface (orchestrates above classes)
```

### Key Classes

**BumpChartConfig:**
```python
config = BumpChartConfig(
    scope='focused',           # 'focused', 'significant', 'all'
    y_type='nes',             # 'nes' or 'rank'
    mode='uniform',           # 'uniform' or 'weighted'
    show_highlights=True,
    show_curves=False,
    label_truncation=50,
    curve_strength=0.5
)
```

**BumpChartBuilder:**
```python
builder = BumpChartBuilder(config, df_subset, mutation='G32A')
fig, ax = builder.build(title="G32A Trajectories")
```

### Scope Modes

| Scope | Description | Expected N | Alpha | Linewidth |
|-------|-------------|-----------|-------|-----------|
| `focused` | Meaningful patterns only | ~100 | 0.7 | 1.5 |
| `significant` | All significant pathways | ~2-4k | 0.15 | 0.5 |
| `all` | All pathways (including NS) | ~12k | 0.05 | 0.3 |

### Highlighting System

Uses `semantic_categories.py` functions to select biologically relevant pathways:
- Filters by `is_relevant_for_highlight()` (excludes cardiac, immune, etc.)
- Ranks by `get_highlight_priority()` (prioritizes mitochondrial, synaptic)
- Limits `Other` category to `MAX_OTHER_CATEGORY_HIGHLIGHTS` (default: 3)

### Pattern Color Integration

```python
from Python.pattern_definitions import get_pattern_colors
pattern_colors = get_pattern_colors()
# Automatically applies pattern-specific colors from canonical source
```

### Usage by Visualization Scripts

**Bump chart generation:**
- `02_Analysis/3.7.viz_bump_chart.py` - Main bump chart script
- `02_Analysis/3.8.viz_interactive_bump_dashboard.py` - Dashboard viewer
- `02_Analysis/3.8.viz_interactive_bump.py` - Interactive viewer
- `02_Analysis/Supp6.app_bump_chart_explorer.py` - Streamlit explorer

---

## 7. patterns.py (Deprecated)

**Legacy wrapper for backward compatibility**

### Status
**DEPRECATED:** This module is kept for backward compatibility only. All new code should import directly from `pattern_definitions.py`.

### Purpose
Provides legacy function names that redirect to canonical implementations in `pattern_definitions.py`. Emits deprecation warnings when used.

### Migration Path

```python
# Old (deprecated)
from Python.patterns import classify_trajectory_pattern, is_compensation

# New (canonical)
from Python.pattern_definitions import classify_pattern, PATTERN_DEFINITIONS
```

---

## 8. __init__.py

**Package initialization and exports**

### Purpose
Defines the public API for the `Python` package, making key functions available via:

```python
from Python import load_classified_pathways, classify_pattern
```

### Exports

```python
__all__ = [
    'CONFIG',
    'SEMANTIC_CATEGORY_ORDER',
    'SEMANTIC_COLORS',
    'PATTERN_COLORS',
    'MUTATION_COLORS',
    'assign_semantic_category',
    'load_classified_pathways',
    'filter_pathways',
    'classify_trajectory_pattern',    # Legacy
    'is_compensation',                # Legacy
    'add_pattern_classification'
]
```

**Note:** Exports legacy function names from `patterns.py` for backward compatibility.

---

## Dependency Graph

```
config.py (base configuration)
    ↓
pattern_definitions.py (canonical pattern logic)
    ↓
data_loader.py → loads classified data
    ↓
semantic_categories.py → categorizes pathways
    ↓
color_config.py → visualization colors
    ↓
viz_bump_charts.py → bump chart rendering
    ↓
[Analysis scripts: 3.7, 3.8, Supp6, etc.]
```

---

## Usage Examples

### Load and Filter Classified Pathways

```python
from Python.data_loader import load_classified_pathways
from Python.config import CONFIG

# Load with default filtering (excludes CGP, canon)
df = load_classified_pathways()

# Load with custom exclusions
df = load_classified_pathways(
    exclude_databases=['cgp', 'canon', 'tf'],
    min_data_points=3
)

# Filter to meaningful patterns only
from Python.pattern_definitions import MEANINGFUL_PATTERNS
df_focused = df[df['Pattern_G32A'].isin(MEANINGFUL_PATTERNS)]
```

### Classify New Data

```python
from Python.pattern_definitions import classify_pattern, add_pattern_classification

# Single pathway
pattern, confidence = classify_pattern(
    early_nes=1.5,
    early_padj=0.01,
    trajdev_nes=-1.2,
    trajdev_padj=0.03,
    late_nes=0.4,
    late_padj=0.15
)
# Returns: ('Compensation', 'High')

# Batch classification
df = add_pattern_classification(df, mutations=['G32A', 'R403C'])
# Adds: Pattern_G32A, Confidence_G32A, Pattern_R403C, Confidence_R403C
```

### Generate Bump Charts

```python
from Python.viz_bump_charts import BumpChartConfig, BumpChartBuilder
from Python.data_loader import load_classified_pathways
from Python.pattern_definitions import MEANINGFUL_PATTERNS

# Load data
df = load_classified_pathways()
df_focused = df[df['Pattern_G32A'].isin(MEANINGFUL_PATTERNS)]

# Configure
config = BumpChartConfig(
    scope='focused',
    show_highlights=True,
    show_curves=False
)

# Build chart
builder = BumpChartBuilder(config, df_focused, mutation='G32A')
fig, ax = builder.build(title="G32A Compensation Pathways")
fig.savefig('g32a_bump_chart.pdf', bbox_inches='tight')
```

### Apply Semantic Categories

```python
from Python.semantic_categories import assign_semantic_category, is_relevant_for_highlight

# Categorize pathways
df['Semantic_Category'] = df['Description'].apply(assign_semantic_category)

# Filter to neuronally-relevant pathways only
df['Relevant'] = df['Description'].apply(is_relevant_for_highlight)
df_relevant = df[df['Relevant']]
```

---

## Testing Pattern Synchronization

To verify pattern classification consistency between Python and R:

```bash
# 1. Regenerate master tables
python3 02_Analysis/1.5.create_master_pathway_table.py
Rscript 02_Analysis/1.7.create_master_gsva_table.R

# 2. Compare pattern distributions
# Python:
python3 -c "
from Python.data_loader import load_classified_pathways
from Python.pattern_definitions import get_pattern_summary
df = load_classified_pathways()
print(get_pattern_summary(df, 'G32A'))
"

# R:
Rscript -e "
df <- read.csv('03_Results/02_Analysis/master_gsva_focused_table.csv')
table(df\$Pattern_G32A)
"

# 3. Verify thresholds match:
#    Python: pattern_definitions.py lines 47-60
#    R:      1.7.create_master_gsva_table.R lines 348-354
```

---

## File Change Log

| Date | File | Change |
|------|------|--------|
| 2025-11-26 | pattern_definitions.py | Initial canonical implementation |
| 2025-12-01 | pattern_definitions.py | Added Sign_reversal pattern, Complex subtypes |
| 2025-12-01 | viz_bump_charts.py | Refactored to modular class-based architecture |
| 2025-12-03 | semantic_categories.py | Updated biological relevance filters |
| 2025-12-03 | color_config.py | Added unified color configuration |

---

## Related Documentation

- **Pattern system details:** `/workspaces/Gama_Vivian_DRP1_bulkRNAseq/docs/PATTERN_CLASSIFICATION.md`
- **Pattern synchronization:** `/workspaces/Gama_Vivian_DRP1_bulkRNAseq/CLAUDE.md` (Pattern System Synchronization Checklist section)
- **R color config:** `/workspaces/Gama_Vivian_DRP1_bulkRNAseq/01_Scripts/R_scripts/color_config.R`
- **R pattern implementation:** `/workspaces/Gama_Vivian_DRP1_bulkRNAseq/02_Analysis/1.7.create_master_gsva_table.R`
- **Project overview:** `/workspaces/Gama_Vivian_DRP1_bulkRNAseq/CLAUDE.md`

---

## Technical Notes

### Pattern Classification: Python vs R

**Python (GSEA):**
- Uses actual p-values for significance testing
- TrajDev significance: `p.adjust < 0.05`
- Applied to GSEA results (NES-based)

**R (GSVA):**
- Uses magnitude thresholds for TrajDev
- TrajDev significance: `|TrajDev| > GSVA_EFFECT`
- Applied to GSVA enrichment scores

**Both approaches are valid:**
- GSVA TrajDev is calculated (difference of differences), not statistically tested
- Magnitude threshold provides equivalent filtering in GSVA context
- Pattern distributions should be consistent across both implementations

### Import Hierarchy

To avoid circular imports, follow this order:
1. `config.py` (no internal dependencies)
2. `pattern_definitions.py` (depends on config)
3. `data_loader.py` (depends on config, pattern_definitions)
4. `semantic_categories.py` (depends on pattern_definitions, color_config)
5. `color_config.py` (standalone, but imports from pattern_definitions)
6. `viz_bump_charts.py` (depends on all above)

### Color Scheme Philosophy

**Diverging vs Categorical:**
- **Diverging (Blue-White-Orange):** Use for continuous metrics (NES, log2FC)
- **Categorical (distinct colors):** Use for discrete groups (mutations, patterns)
- **Never mix:** Don't use blue/orange from diverging palette for categorical annotations

**Colorblind Accessibility:**
- All palettes tested for deuteranopia/protanopia
- Avoid red-green distinctions
- Use shape/pattern in addition to color when possible

---

**Last Updated:** 2025-12-04
**Maintainer:** Claude Code
**Questions:** Refer to `/workspaces/Gama_Vivian_DRP1_bulkRNAseq/CLAUDE.md`
