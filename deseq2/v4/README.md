# TICTOC V4 — factorial interaction model (reproducible successor to iDEP v3)

`run_factorial_v4.R` reproduces the multifactorial analysis from the original
`TICTOC_3_factor_model/markdown reports/TICTOC_markdown_v3.html` (iDEP 0.95) in plain, versioned
DESeq2 — and carries its two most useful results into the manuscript. **This supersedes the v3 HTML**
for the interaction analysis (no hard-coded paths, pinned versions in `../../ENVIRONMENT.txt`).

## Model
`~ Treatment * Genotype * Tissue`, reference levels **Ground / WT / Root** (same as v3). Input: the
filtered count matrix (59,918 × 48). DEGs at P.adj < 0.05, |log₂FC| ≥ 1.

## Key results (executed 2026-07-04)

**Interaction-term DEG counts** (`interaction_DEG_counts.csv`):

| Interaction term | up | down | total |
|---|--:|--:|--:|
| **Treatment(Flight) × Tissue(Shoot)** | 1670 | 1704 | **3374** |
| Treatment(Flight) × Genotype(A68) | 2 | 60 | 62 |
| Treatment(Flight) × Genotype(D130) | 0 | 20 | 20 |
| Genotype(D130) × Tissue(Shoot) | 1 | 4 | 5 |
| Genotype(A68) × Tissue(Shoot) | 1 | 0 | 1 |
| 3-way terms | — | — | ≤1 |

→ **The spaceflight response is overwhelmingly tissue-specific** (Treatment×Tissue = 3,374 DEGs, dwarfing
all else). The **AVP-OX × Flight interactions are small and strongly down-biased** — the engineered lines
*attenuate* the wild-type flight programme (converges with GO R5 / PhysioSpace R6). Genotype×Tissue and
3-way interactions are negligible.

**PC–factor correlation** (`pc_factor_correlation.csv`, Kruskal–Wallis):

| PC | % var | tracks |
|---|--:|---|
| PC1 | 92.5% | **Tissue** (p = 2.9×10⁻⁹) |
| PC2 | 1.7% | **Treatment** (p = 6.3×10⁻⁶) |
| PC3 | 1.0% | Treatment (p = 3.4×10⁻³) |

→ Tissue is the dominant transcriptomic axis; spaceflight is the second. Genotype is not significantly
associated with any top PC (a subtle, gene-specific effect — consistent with the small genotype-interaction DEG counts).

**Depth QC** (`v4_summary.txt`): total read depth differs by **Tissue** (Kruskal p ≈ 0) but **not by
Treatment** (p = 0.41) or Genotype (p = 0.98). DESeq2 size-factor normalisation handles the tissue depth
difference, and critically the **spaceflight (Treatment) comparisons are not depth-confounded** — supporting
the Flight-vs-Ground claims.

## What V4 improves over v3
- Reproducible & versioned (no hard-coded `setwd`; `../../ENVIRONMENT.txt`).
- Same factorial model, now feeding the integrated thesis alongside the crosswalk-based GO, PhysioSpace,
  and morphometrics (all absent from v3).
- Quantifies the depth confound explicitly and shows Treatment is clean.
