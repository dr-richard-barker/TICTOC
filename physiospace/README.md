# PhysioSpace stress-pattern decoding

Projects the TICTOC cotton Flight-vs-Ground contrasts onto **Arabidopsis stress reference spaces**
to quantify which known stress programs the spaceflight response resembles (PhysioScores). This is
roadmap step **4.5** and mirrors the validated OSD-767 tomato pipeline.

## Run
```bash
cd physiospace
Rscript run_physiospace.R          # loads AT stress spaces from the PlantPhysioSpace package
```
Install (once):
```r
install.packages(c("BiocManager","remotes"))
BiocManager::install(c("DESeq2","org.At.tair.db"))
remotes::install_github("JRC-COMBINE/PhysioSpaceMethods")   # method (calculatePhysioMap, PhysioHeatmap)
remotes::install_github("JRC-COMBINE/PlantPhysioSpace")     # data: AT/rice/soy/wheat stress spaces
```

## What it does
1. VST-transforms the count matrix and computes, per **genotype × tissue**, the mean
   **Flight − Ground** fold change (input matrix, genes × contrasts).
2. Re-indexes cotton → Arabidopsis: `Gohir → AT locus` (via
   [`../crosswalk/gohir_to_arabidopsis.tsv`](../crosswalk/gohir_to_arabidopsis.tsv)) →
   `AT Entrez` (`org.At.tair.db`), collapsing many-Gohir → one-Entrez by mean.
3. Runs `PhysioSpaceMethods::calculatePhysioMap(..., GenesRatio = 0.05, TTEST = FALSE,
   ImputationMethod = "PCA")` against each reference space — the **same parameters as OSD-767**.
4. Writes `PhysioScores_<space>.csv` + a heatmap per reference space to `results/`.

## Reference spaces (from the `PlantPhysioSpace` data package)
Loaded automatically via `data(<name>, package = "PlantPhysioSpace")` — no manual files needed.
Default = the three used by OSD-767 (all Arabidopsis, Entrez-indexed):

- `AT_Stress_Space` (22 stress axes, ATH1 microarray compendium, 85 GEO datasets)
- `AT_Stress_Space_Meta` (meta-grouped stress axes)
- `AT_Stress_Space_RNASeq` (RNA-seq-derived axes)

Override with `--spaces "A,B,C"`. Also available in the package: `*_Detailed` variants, and
non-Arabidopsis spaces `OS_Stress_Space` (rice), `GM_Stress_Space` (soybean), `TA_Stress_Space` (wheat).
To use local copies instead, pass `--ref-dir /folder` of `<name>.rds` matrices.

## Executed run (2026-07-04, R 4.6)
Ran end-to-end. 59,918 genes → 6 Flight−Ground contrasts (WT/A68/D130 × root/shoot) → 14,790 Arabidopsis
Entrez genes, ~13,367 shared with `AT_Stress_Space`. Outputs in `results/` (default signed-p) and
`results_static/` (`--static TRUE`).

> **Use `--static TRUE` for this dataset.** In default signed-p mode every **root** contrast saturates to
> `±Inf` (p-value underflow — the root Flight response is very strong), so roots are unusable there.
> `STATICResponse=TRUE` returns the bounded statistic (~−0.3…+0.8) and is fully interpretable.

**Preliminary observations (descriptive; from `results_static/AT_Stress_Space`):**
- Spaceflight **roots of all genotypes activate a coordinated abiotic/defence stress signature** —
  Osmotic, Drought, FarRed, Genotoxic, Wounding, Drought.Light most strongly.
- **Shoots are muted/mixed** (WT shoot: Salt, Hormone, UV, Biotic positive; A68/D130 shoots mostly negative).
- The **AVP-OX roots show somewhat *lower* stress-program activation than WT** (e.g. Osmotic WT 0.78 >
  A68 0.69 > D130 0.61; Wounding WT 0.72 > A68 0.51 > D130 0.44) — a possible transcriptomic correlate of
  engineered stress tolerance. ⚠ modest, descriptive; needs a formal contrast + the frozen-model input.

## Caveats
- The `±Inf` saturation above: prefer `--static TRUE`, or cap infinities before plotting signed-p.
- The cotton→Arabidopsis bridge is best-BLASTP-hit orthology (approximate); raise `--min-pid` for a
  stricter mapping. OSD-767 used Ensembl Compara 1:1 orthologs — expect somewhat noisier cotton scores.
- Fold changes here are VST group-mean differences (blind VST). If you freeze a DESeq2 model
  (roadmap §4.1) with shrunken LFCs, swap those in as the input matrix for consistency with the DEG tables.

## Options
`--spaces` (default `AT_Stress_Space,AT_Stress_Space_Meta,AT_Stress_Space_RNASeq`),
`--counts` (default `../TICTOC_run1_filteredCounts_v3.csv`), `--crosswalk`, `--out` (default `results`),
`--min-pid` (default 30), `--genes-ratio` (default 0.05), `--ref-dir` (optional local `.rds` override).

## References
- Lenz M. *et al.* (2013) *PhysioSpace: relating gene expression experiments from heterogeneous sources
  using shared physiological processes.* PLoS ONE 8(10):e77627. doi:10.1371/journal.pone.0077627
- Hadizadeh Esfahani A. *et al.* (2021) *Plant PhysioSpace: transferring stress-response knowledge across
  species.* Plant Physiology 187(3):1795. (PMC8566250)
- Packages: `PhysioSpaceMethods` and `PlantPhysioSpace` (JRC-COMBINE), GNU GPLv3.
