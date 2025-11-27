# Comparison of Pattern Classification Between `main` and `dev`

This report summarizes how trajectory patterns are defined and counted in the
`main` vs `dev` branches, explains the implementation differences that cause
the discrepancies (especially the loss of “Progressive” pathways in `dev`),
and evaluates which implementation has stronger biological grounding.

All numbers below are computed directly from:

- `03_Results/02_Analysis/Plots/Cross_database_validation/pathways_classified.csv`
- `03_Results/02_Analysis/master_gsea_table.csv`

using the threshold values encoded in:

- `01_Scripts/Python/pattern_definitions.py`
- `01_Scripts/Python/config.py`

Paths and line numbers in this report refer to the `dev` branch unless
otherwise noted.

---

## 1. Numeric Comparison of `main` vs `dev`

We focused on the 12,221 unique pathways represented in
`03_Results/02_Analysis/master_gsea_table.csv:1` (one pattern per pathway).
For each pathway and mutation (G32A, R403C) we applied:

1. The **`main` (magnitude-only) classifier** reimplemented from
   `main:01_Scripts/Python/patterns.py` (pure NES-based, ignoring p.adjust).
2. The **`dev` significance-based classifier** as stored in the `Pattern_*`
   columns of `master_gsea_table.csv`, computed via
   `pattern_definitions.add_pattern_classification` in
   `02_Analysis/9.create_master_pathway_table.py:60`–92.

Results by mutation:

### G32A (12,221 unique pathways)

**Main (magnitude-only)**  
Recomputed from the legacy logic in `main:01_Scripts/Python/patterns.py`:

- `Compensation`: 5,908  
- `Progressive`: 980  
- `Natural_worsening`: 1,970  
- `Natural_improvement`: 222  
- `Late_onset`: 19  
- `Transient`: 90  
- `Complex`: 30  
- `Insufficient_data`: 3,002  

Worsening patterns (Progressive + Natural_worsening): **2,950**

**Dev (significance-based)**  
Using the canonical classifier from `01_Scripts/Python/pattern_definitions.py`
with thresholds:

- `NES_EFFECT = 0.5`  
- `NES_STRONG = 1.0`  
- `IMPROVEMENT_RATIO = 0.7`  
- `WORSENING_RATIO = 1.3`  
- `PADJ_SIGNIFICANT = 0.05`  
- `PADJ_TRENDING = 0.10`

Pattern counts (all confidence levels):

- `Compensation`: 1,020  
- `Progressive`: **0**  
- `Natural_worsening`: 2  
- `Natural_improvement`: 667  
- `Late_onset`: 46  
- `Transient`: 1  
- `Complex`: 7,483  
- `Insufficient_data`: 3,002  

Worsening (Progressive + Natural_worsening): **2**

### R403C (12,221 unique pathways)

**Main (magnitude-only)**:

- `Compensation`: 6,196  
- `Progressive`: 848  
- `Natural_worsening`: 1,791  
- `Natural_improvement`: 299  
- `Late_onset`: 15  
- `Transient`: 53  
- `Complex`: 17  
- `Insufficient_data`: 3,002  

Worsening (Progressive + Natural_worsening): **2,639**

**Dev (significance-based)**:

- `Compensation`: 1,087  
- `Progressive`: **0**  
- `Natural_worsening`: 2  
- `Natural_improvement`: 952  
- `Late_onset`: 61  
- `Transient`: 11  
- `Complex`: 7,106  
- `Insufficient_data`: 3,002  

Worsening (Progressive + Natural_worsening): **2**

**Summary of numerical shifts:**

- All **Progressive** pathways disappear in `dev`
  (G32A: 980 → 0; R403C: 848 → 0).
- Worsening patterns (Progressive + Natural_worsening) drop from
  2,950 → 2 (G32A) and 2,639 → 2 (R403C).
- `Compensation` shrinks from ~6K pathways per mutation to ~1K.
- A large fraction of pathways become `Complex` in `dev`.

The dev counts correspond to the “Pattern Distribution Summary” you reported:
Progressive goes to 0% for both mutants, with most pathways falling into
Compensation, Passive patterns, or Complex.

---

## 2. Implementation Differences That Drive These Changes

There are two main changes in how patterns are defined and applied.

### 2.1 New canonical significance-based classifier (dev)

The `dev` branch introduces a **canonical, significance-based classifier** in:

- `01_Scripts/Python/pattern_definitions.py`

Key features:

- **Thresholds** (`pattern_definitions.py:223`–241):
  - `PADJ_SIGNIFICANT = 0.05`
  - `PADJ_TRENDING = 0.10`
  - `NES_EFFECT = 0.5`
  - `NES_STRONG = 1.0`
  - `IMPROVEMENT_RATIO = 0.7` (≥30% improvement)
  - `WORSENING_RATIO = 1.3` (≥30% worsening)

- **Classification function**:
  - `classify_pattern(early_nes, early_padj, trajdev_nes, trajdev_padj, late_nes, late_padj)`
    at `pattern_definitions.py:168`–294.

- **Early defect requirements** (`pattern_definitions.py:220`–226):
  - High confidence: `p.adjust_Early < 0.05` AND `|NES_Early| > 0.5`
  - Medium confidence (trending): `p.adjust_Early < 0.10` AND `|NES_Early| > 0.5`

- **Late outcome** (`pattern_definitions.py:231`–241):
  - Improved: `|Late|/|Early| < 0.7` OR `|Late| < 0.5`
  - Worsened: `|Late|/|Early| > 1.3`
  - Resolved: `|Late| < 0.5`

- **TrajDev significance and directionality** (`pattern_definitions.py:243`–247):
  - TrajDev significant: `p.adjust_TrajDev < 0.05` AND `|NES_TrajDev| > 0.5`
  - TrajDev opposes Early: `sign(TrajDev) != sign(Early)`
  - TrajDev amplifies Early: `sign(TrajDev) == sign(Early)`

- **Progressive definition** (`pattern_definitions.py:271`–275):
  - Early significant defect (High or Medium),
  - TrajDev significant and amplifying,
  - Late worsened (ratio > 1.3).

- **Compensation definition** (`pattern_definitions.py:271`–272):
  - Early significant defect,
  - TrajDev significant and opposing,
  - Late improved/resolved.

- **Passive patterns** (`pattern_definitions.py:282`–287):
  - When TrajDev is **not** significant, trajectories with improvement or
    worsening are assigned to `Natural_improvement` or `Natural_worsening`.

The canonical pattern logic is heavily documented in
`docs/PATTERN_CLASSIFICATION.md`, which explicitly ties patterns to
statistically supported Early and TrajDev components.

`01_Scripts/Python/patterns.py` now acts as a **wrapper** over this canonical
logic (`patterns.py:55`–141):

- Default mode: `use_significance=True` → calls `classify_pattern(...)`.
- Legacy mode: `use_significance=False` → calls `_classify_trajectory_pattern_legacy(...)`
  which reproduces the main-branch magnitude-only behavior.

The master table generator (`02_Analysis/9.create_master_pathway_table.py`)
now **recomputes patterns using the canonical classifier**:

- `load_pattern_classifications()` at `9.create_master_pathway_table.py:60`–93
  loads NES and p.adjust values from `pathways_classified.csv`.
- It then calls `add_pattern_classification(df_patterns, mutations=['G32A', 'R403C'])`,
  which comes from `pattern_definitions.py` and produces:
  - `Pattern_G32A`, `Confidence_G32A`
  - `Pattern_R403C`, `Confidence_R403C`

These recomputed patterns are what you see in `master_gsea_table.csv`.

### 2.2 Magnitude-only classifier in `main`

In the `main` branch, pattern classification is done entirely in:

- `01_Scripts/Python/patterns.py` (old version; see `git show main:01_Scripts/Python/patterns.py`)

Key characteristics of the **legacy, magnitude-only classifier**:

- Uses only **NES values**, no p.adjust, via `classify_trajectory_pattern(row, ...)`.
- Thresholds:
  - `low_threshold = CONFIG['low_threshold']` (0.5)
  - `defect_threshold = 1.0`
  - From `01_Scripts/Python/config.py` (same in main and dev):
    - `low_threshold: 0.5`
    - `nes_threshold: 1.5` (used elsewhere, not in the legacy pattern function).

- Pattern logic (now preserved in dev as `_classify_trajectory_pattern_legacy`
  at `patterns.py:144`–190):
  - `Late_onset`: `|Early| < 0.5` and `|Late| > 1.0`
  - `Transient`: `|Early| > 1.0` and `|Late| < 0.5`
  - `Compensation`: `|Early| > 0.5`, `Late < |Early|`, and
    TrajDev opposes Early with `|TrajDev| > 0.5`
  - `Natural_improvement`: same as above but TrajDev magnitude ≤ 0.5
  - `Progressive`: `|Early| > 0.5`, `Late > |Early|`, and
    TrajDev amplifies Early with `|TrajDev| > 0.5`
  - `Natural_worsening`: same as above but TrajDev magnitude ≤ 0.5
  - `Persistent`: stable defect with little change and small TrajDev
  - `Complex`: everything else

Crucially, **pattern labels in main are assigned without any control for
statistical significance** (p.adjust is not consulted at all for pattern
classification). Even mildly noisy trajectories with small p.adjust values
or borderline NES values can be labelled as Progressive or Compensation.

In `main`, the master table did **not** recompute patterns; it simply copied
what was in `pathways_classified.csv`, which had been generated using this
magnitude-only logic.

---

## 3. Why do Progressive pathways disappear in `dev`?

To understand the disappearance of Progressive patterns in `dev`, we examined
the set of pathways that main calls Progressive and asked how many of them
would meet dev’s stricter criteria.

Using `pathways_classified.csv`, for each mutant:

### G32A Progressive pathways (main)

- Total Progressive (main): **980** pathways
- Among these 980:
  - Early significant defects by dev criteria (`p.adjust_Early < 0.05`,
    `|NES_Early| > 0.5`): **2**
  - TrajDev significant by dev criteria (`p.adjust_TrajDev < 0.05`,
    `|NES_TrajDev| > 0.5`): **10**
  - **Both** Early and TrajDev significant: **0**
  - Late worsening ratio > 1.3 (`|Late|/|Early| > 1.3`): **619**
  - Early+TrajDev significant AND ratio > 1.3: **0**

### R403C Progressive pathways (main)

- Total Progressive (main): **848** pathways
- Among these 848:
  - Early significant defects (dev criteria): **6**
  - TrajDev significant (dev criteria): **3**
  - Both Early and TrajDev significant: **0**
  - Late worsening ratio > 1.3: **462**
  - Early+TrajDev significant AND ratio > 1.3: **0**

So, for both mutants:

- **None** of the main-branch Progressive pathways satisfy dev’s
  **Progressive** definition:
  - You never see a pathway with:
    - a statistically robust Early defect,
    - a statistically robust, amplifying TrajDev,
    - and a clearly worsened Late state (≥30% increase),
    all at once.

That is exactly why `dev` has **0 Progressive pathways** for both G32A and
R403C. The transition from main → dev does not remove biological signals;
it applies more stringent, biologically motivated criteria that the data
simply do not satisfy for Progressive patterns.

The same analysis for main **Compensation** pathways shows that many of them
do satisfy dev’s stricter criteria, which is why Compensation remains
substantial in dev:

### G32A Compensation pathways (main)

- Total Compensation (main): **5,908**
- Among these:
  - Early significant defects (dev criteria): **1,579**
  - TrajDev significant (dev criteria): **1,221**
  - Both Early and TrajDev significant: **1,057**
  - Late improvement (ratio < 0.7 or `|Late| < 0.5`): **3,320**
  - Early+TrajDev significant AND improvement: **963**

In dev, G32A has 1,020 Compensation pathways, very close to the 963
“fully satisfying” main Compensation trajectories above, as expected.

### R403C Compensation pathways (main)

- Total Compensation (main): **6,196**
- Among these:
  - Early significant defects: **2,166**
  - TrajDev significant: **1,423**
  - Both Early and TrajDev significant: **1,266**
  - Late improvement: **3,389**
  - Early+TrajDev significant AND improvement: **1,044**

In dev, R403C has 1,087 Compensation pathways, again very close to the
“fully satisfying” subset of main Compensation calls.

So:

- The dev classifier **retains** many of main’s Compensation pathways
  that have strong statistical and effect-size support.
- It **removes** essentially all Progressive pathways because none of them
  meet the stronger criteria.

---

## 4. Changes in the Pattern Landscape

Putting this together, for the 9,219 pathways with complete trajectory data
(12,221 total minus 3,002 Insufficient_data), the landscape changes as follows.

### Main (magnitude-only)

- Many pathways are labeled as **Progressive** or **Natural_worsening**
  (2,950 G32A; 2,639 R403C), i.e. worsening is very common.
- **Compensation** is also abundant (~6,000 per mutant).
- Very few pathways are labeled `Complex` (30 G32A; 17 R403C).
- FDR/p.adjust is not used to gate these labels; any pathway with the right
  shape and NES magnitudes is assigned to a named pattern.

### Dev (significance-based)

- Worsening patterns nearly vanish:
  - `Progressive = 0`, `Natural_worsening = 2` for each mutant.
- **Compensation** is reduced but still substantial (~1,000 per mutant),
  representing pathways with clear evidence of early defect, significant
  trajectory deviation, and late improvement.
- **Natural_improvement** increases, capturing passive recovery when
  TrajDev is not significant.
- `Complex` explodes to ~7,000 pathways, absorbing:
  - trajectories with ambiguous or inconsistent dynamics,
  - marginal/noisy changes that fail strict Early/TrajDev/Late criteria,
  - patterns that do not fit cleanly into the 7 canonical categories.

These shifts are exactly what your summary table shows: dev is much more
stringent and conservative about calling active worsening.

---

## 5. Biological Interpretation: Which Implementation is Better Grounded?

Given the biological framing in `docs/PATTERN_CLASSIFICATION.md` and
`pattern_definitions.py`, the **dev** implementation is more biologically and
statistically sound than the **main** implementation.

### 5.1 TrajDev as active plasticity

The methods document defines TrajDev as:

> TrajDev = (D65_mut - D35_mut) - (D65_ctrl - D35_ctrl)

and emphasizes:

- Significant TrajDev implies **active transcriptional plasticity**.
- Active patterns (Compensation, Progressive) should reflect **active
  deviation** from normal maturation, not just passive drift.

The dev classifier respects this:

- Active patterns (Compensation, Progressive) require **significant TrajDev**
  (p.adjust < 0.05, |NES| > 0.5).

The main classifier does **not**: it labels many pathways as Progressive or
Compensation solely based on the sign and magnitude of NES_TrajDev, without
any statistical support.

### 5.2 Early defects must be statistically supported

Biologically, calling a pathway “defective early” implies that its
perturbation is reliably detected at D35. The dev classifier enforces this:

- Early defect requires significance or at least trend-level evidence:
  - High: `p.adjust_Early < 0.05`
  - Medium: `p.adjust_Early < 0.10`

Main’s magnitude-only approach treats any |NES_Early| > 0.5 as a meaningful
defect, regardless of FDR, which is unsafe when screening ~10k pathways.

### 5.3 Late outcome thresholds respected

The dev classifier uses relative change thresholds
(`IMPROVEMENT_RATIO`, `WORSENING_RATIO`) plus absolute NES cutoffs. This:

- Avoids over-interpreting tiny Late–Early differences as “worsening”
  or “improvement”.
- Aligns better with the idea that “progression” should be a **robust**
  worsening with effect size, not just a fractional increase.

Main simply compares Late vs Early, so any Late > Early trajectory with
TrajDev pointing in the same direction becomes Progressive, even for
borderline cases.

### 5.4 Honest handling of ambiguity and missing data

Dev explicitly marks:

- `Insufficient_data` when any required NES/p.adjust is missing
  (`pattern_definitions.py:212`–215, 154–159).
- `Complex` when the trajectory does not cleanly match other patterns
  (`pattern_definitions.py:144`–151, 294).

Main puts almost every pathway into a named category if NES values are
non-missing, even when data are noisy or inconsistent. This can create an
illusion of precision.

### 5.5 Interpreting “no Progressive pathways”

Under the dev criteria, “no Progressive pathways” means:

- In this dataset, there are no pathways with:
  - robust, significant Early defects,
  - robust, significant TrajDev deviations that amplify those defects,
  - and clearly worsened Late states (≥30% increase),
  considered together.

This does **not** mean:

- There are no pathways with any suggestion of worsening.
- There is no trajectory-level evidence of vulnerability.

It simply means:

- The evidence for strong, **actively progressive** transcriptional trajectories
  does not pass stringent effect-size and FDR filters across all stages.

From a biological standpoint, this is plausible and preferable to claiming
thousands of strongly progressive pathways based purely on NES shape without
statistical backing.

---

## 6. Conclusions and Recommendations

1. **Implementation differences**  
   - `main` uses a **magnitude-only** classifier (NES thresholds, no p.adjust),
     implemented in `patterns.py` and applied to generate the original pattern
     labels.
   - `dev` uses a **significance-based** classifier, implemented in
     `pattern_definitions.py` and re-applied to trajectory data in
     `9.create_master_pathway_table.py` to produce `Pattern_*` and
     `Confidence_*` columns.

2. **Effect on pattern counts**  
   - Progressive and Natural_worsening patterns are numerous in main but
     nearly vanish in dev because they lack the required combination of
     significant Early and TrajDev components plus strong Late worsening.
   - Compensation shrinks but remains substantial, corresponding to the subset
     of main Compensation calls that are statistically well-supported.
   - Many pathways become Complex under dev, reflecting trajectories that do
     not meet strict canonical criteria.

3. **Biological grounding**  
   - For mechanistic claims (e.g., “this pathway is actively progressive”),
     the dev implementation is **much better grounded**:
     - It enforces significance for both Early and TrajDev.
     - It uses explicit thresholds for meaningful late changes.
     - It encodes the biological role of TrajDev as active plasticity.
   - The main implementation is useful for **exploratory** pattern discovery
     but over-calls progression and compensation when interpreted as
     “hard” biological categories.

4. **Practical recommendation**  
   - For manuscript figures, conclusions, and biological interpretation,
     rely on **dev** pattern calls (or at least treat main calls as
     hypothesis-generating that require dev-style validation).
   - The absence of Progressive pathways in dev should be reported as a
     **negative but informative finding**: within this dataset, robust,
     actively progressive transcriptional trajectories are rare or below the
     detection threshold given the adopted statistical and effect size
     criteria.

If needed, we can extend this report with:

- Per-database pattern distributions under dev vs main.
- A list of “near-progressive” pathways (strong NES changes but failing one
  dev criterion) for targeted biological review.

