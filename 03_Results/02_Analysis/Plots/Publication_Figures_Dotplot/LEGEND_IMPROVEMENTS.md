# Size Legend Improvements

## Issues Fixed (2025-11-25)

### Problem 1: Text Overlap
**Before**: Text labels were overlapping each other in the legend
**After**:
- Increased horizontal spacing between dots (using `np.linspace(0.15, 0.85, 4)`)
- Separated labels into two lines with better vertical positioning:
  - Line 1 (bold): `padj = 0.001` at y=0.28
  - Line 2 (italic): `(highly sig)` at y=0.15
- Added color differentiation for description text (`color='#555'`)

### Problem 2: Incorrect Dot Sizes
**Before**: Dot sizes in legend were scaled (`dot_sizes * 2`) and didn't match actual figures
**After**:
- Using **actual** calculated sizes directly from `renderer._calculate_dot_sizes()`
- No artificial scaling applied
- Legend dots now perfectly match the dots in the figures

### Problem 3: Cramped Legend Area
**Before**: Legend area was too short (height_ratios: 0.15, 1.5, 2)
**After**:
- Fig1 (Ribosome Paradox): height_ratio increased from 0.15 to 0.25
- Fig2 (MitoCarta): height_ratio increased from 1.5 to 2.0
- Fig3b (SynGO): height_ratio increased from 1.5 to 2.0
- Fig4 (Semantic Overview): height_ratio increased from 2 to 3
- Overall figure heights increased by 1-2 inches

## Legend Layout

```
┌─────────────────────────────────────────────────────────┐
│  Dot Size Legend: Larger dots = more significant        │  ← Title (y=0.90)
│                                                          │
│    ●      ●●      ●●●     ●●                            │  ← Dots (y=0.50)
│  0.001   0.01    0.05    0.1                            │  ← padj values (y=0.28)
│ (highly) (very)  (thresh) (not)                         │  ← Labels (y=0.15)
│   sig     sig     old)    sig)                          │
│                                                          │
│  Black edge = significant | Gray edge = not significant │  ← Edge explanation (y=0.02)
└─────────────────────────────────────────────────────────┘
```

## Technical Details

### Dot Size Calculation
The legend uses the same calculation as the main figures:
```python
# From DotplotRenderer._calculate_dot_sizes()
neg_log_padj = -np.log10(padj_matrix)
neg_log_padj_clipped = np.clip(neg_log_padj, 0, 10)
normalized = neg_log_padj_clipped / 10
dot_sizes = min_dot_size + normalized * (max_dot_size - min_dot_size)
```

### Example Sizes (for min_dot_size=20, max_dot_size=300)
- padj = 0.001 → -log10 = 3.0 → size = 104
- padj = 0.01  → -log10 = 2.0 → size = 76
- padj = 0.05  → -log10 = 1.3 → size = 56
- padj = 0.1   → -log10 = 1.0 → size = 48

### Edge Colors
- Black outline + linewidth=2: padj < 0.05 (significant)
- Gray outline + linewidth=1: padj ≥ 0.05 (not significant)

## Verification
To verify dot sizes match between legend and figures:
1. Open any dotplot PDF
2. Compare dots in the legend with dots in the main figure
3. For a given padj value, sizes should be identical

## Files Updated
- `02_Analysis/3.2.publication_figures_dotplot.py`
  - `add_size_legend_visual()` function completely rewritten
  - All figure GridSpec layouts updated with taller legend areas
- All dotplot PDFs regenerated with improved legends

## Result
The size legend now clearly and accurately explains the dot size-to-padj relationship without any text overlap or size mismatches, making it easy for readers to interpret the significance of pathways in the dotplots.
