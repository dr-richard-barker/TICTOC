# DESeq2 — full contrast set

Fits one DESeq2 model over all 48 samples and extracts the canonical contrasts, so functional analysis
(`../go_analysis/`) can run across **all** genotypes × tissues (not just the two legacy iDEP tables).
This is roadmap §4.1–4.2, executed.

## Run
```bash
cd deseq2
Rscript run_deseq2_contrasts.R          # -> contrasts/Diff_genes_heatmap_*.csv + DEG_counts_summary.csv
```
Deps: `DESeq2` (+ `ashr` for LFC shrinkage; falls back to MLE LFC if absent).

## Model
- Design `~ 0 + group`, where `group = Genotype_Tissue_Treatment` (12 groups × 4 reps).
- Input: `../TICTOC_run1_filteredCounts_v3.csv` (59,918 genes × 48). Replicate columns share a header
  label; the script uniquifies them.
- DEG call: `P.adj < 0.05` and `|log₂FC| ≥ 1` (ashr-shrunk LFC).

## Contrasts (10)
- **Spaceflight** (Flight − Ground) within each Genotype × Tissue — 6.
- **AVP-OX interaction** `(line FL−GC) − (WT FL−GC)` per tissue — 4 (A68, D130 × Root, Shoot).

## Executed result (2026-07-04) — `contrasts/DEG_counts_summary.csv`
| Contrast | up | down |
|---|--:|--:|
| WT_Root_FlightVsGround | 2503 | 578 |
| A68_Root_FlightVsGround | 2806 | 3610 |
| D130_Root_FlightVsGround | 420 | 299 |
| WT_Shoot_FlightVsGround | 335 | 259 |
| A68_Shoot_FlightVsGround | 79 | 419 |
| D130_Shoot_FlightVsGround | 327 | 890 |
| A68vsWT_Root_FlightInteraction | 2 | 55 |
| D130vsWT_Root_FlightInteraction | 0 | 20 |
| A68vsWT_Shoot_FlightInteraction | 2 | 3 |
| D130vsWT_Shoot_FlightInteraction | 40 | 155 |

Root ≫ shoot response; interaction contrasts down-biased (AVP-OX attenuates the WT flight programme).

## Notes
- DEG tables use the same columns as `../TICTOC_3_factor_model/Diff_genes_heatmap_*.csv`, so they feed
  `../go_analysis/run_go_clusterprofiler.R --deg-glob "../deseq2/contrasts/Diff_genes_heatmap_*.csv"`.
- This model supersedes the two legacy iDEP tables for the full-coverage analysis; keep both until the
  manuscript figures are locked.
