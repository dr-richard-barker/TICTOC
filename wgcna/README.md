# WGCNA — root co-expression modules ↔ Flight / genotype / root traits

Rebuilds the co-expression side of the iDEP v3 report reproducibly and adds the integrative step
(Fig 4 / R4): correlating module eigengenes with treatment, genotype, **and the RSML root
morphometric traits** — linking transcriptome to root architecture.

## Run
```bash
cd wgcna
Rscript run_wgcna_root.R      # -> results/
```
Deps: `WGCNA`, `DESeq2`, `clusterProfiler`, `org.At.tair.db`. Root samples only (24; tissue dominates
PC1, so co-expression is built within the tissue where both the flight response and the RSML traits live).

## Method
VST → top 5,000 variable root genes → signed network, soft power 12, `minModuleSize` 30 →
6 modules. Module eigengenes correlated (Pearson) with Flight, AVP-OX/genotype indicators, and
group-level day-6 RSML traits (total length, primary length, lateral count).

## Result (executed 2026-07-05) — `module_trait_correlation.csv`

| Module (genes) | Flight | AVP-OX | total len | primary len | n_lateral |
|---|--:|--:|--:|--:|--:|
| **turquoise (1963)** | **+0.89** | +0.06 | **+0.86** | +0.50 | **+0.83** |
| **blue (1767)** | **−0.83** | −0.08 | **−0.83** | −0.54 | **−0.81** |
| **brown (1046)** | **+0.66** | **−0.43** | +0.41 | −0.19 | +0.33 |
| yellow (68) | +0.32 | −0.05 | +0.35 | +0.24 | +0.35 |
| grey/unassigned (118) | −0.09 | +0.28 (D130 +0.53) | | | |

(bold = P < 0.05.)

**Interpretation:**
- **turquoise / blue** are a flight-responsive, growth-coupled axis: turquoise is up in flight and
  positively correlated with a larger, more-branched root system; blue is its suppressed mirror. GO of
  turquoise: signaling / signal transduction, isoprenoid–terpenoid metabolism, transcription.
- **brown** is flight-induced but **lower in the AVP-OX lines** (AVP-OX r = −0.43) — the module-level
  transcriptomic correlate of the AVP-OX dampening seen in GO (R5) and PhysioSpace (R6).

## Caveats
- 24 samples is modest for WGCNA — treat modules as exploratory; the Flight axis is the robust signal.
- Morphometric correlations are **group-level** (6 genotype × treatment groups drive them) and inherit the
  FL/GC image-calibration caveat (see `../morphometrics/README.md`): they show that flight-responsive
  modules co-vary with the *traced* root phenotype, not an individual-plant match.

## Module names (from GO BP/MF/CC — `results/module_names.csv`)
turquoise = **Signalling & isoprenoid metabolism** · blue = **Translation & ribosome biogenesis** · brown = **Defence & ubiquitin signalling** · yellow = **Metal transport & phenylpropanoid metabolism** · green = **Photosynthesis (light reactions)** · grey = Unassigned. Per-module GO tables: `results/GO_*_module_*.csv`; integrated model + limitations: `../integration/README.md`.
