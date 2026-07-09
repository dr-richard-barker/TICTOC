# Reproducing the TICTOC analysis

Every result, figure and table in the manuscript regenerates from the committed scripts and the
count matrix / RSML tracings in this repo. This is the ordered pipeline; `run_all.sh` runs it end to end.

## Prerequisites
- **R ≥ 4.5** (tested on 4.6.0). Packages: `DESeq2`, `clusterProfiler`, `org.At.tair.db`, `enrichplot`,
  `WGCNA`, `ashr`, `ggplot2`, `pheatmap`, and (from JRC-COMBINE GitHub) `PhysioSpaceMethods`,
  `PlantPhysioSpace`. Exact versions in [`ENVIRONMENT.txt`](ENVIRONMENT.txt).
- **Python 3** with `numpy`, `pandas`, `statsmodels`, `matplotlib` (stdlib `csv/xml` for the rest).

> ⚠ **No-compiler machines (no Rtools/`make`):** several Bioconductor packages ship source-only and
> fail via `BiocManager::install` (segfaults on the annotation-DB build). Install those from tarball with
> `R CMD INSTALL` instead — they are pure R/data and need no compilation:
> `GO.db`, `HDO.db`, `DOSE`, `enrichplot` (download via `download.packages(..., type="source")`, then
> `R CMD INSTALL --no-multiarch pkg.tar.gz`). CRAN/Bioc **binaries** install normally.

## Pipeline (order matters)
| # | Stage | Command | Key output |
|---|---|---|---|
| 1 | Ortholog crosswalk | `python crosswalk/build_gohir_to_arabidopsis.py` | `crosswalk/gohir_to_arabidopsis.tsv` |
| 2 | Morphometrics | `python morphometrics/extract_rsml_traits.py` · `summarise_traits.py` · `morphometric_stats.py` | `morphometrics/rsml_traits*.csv`, stats |
| 3 | DESeq2 contrasts | `Rscript deseq2/run_deseq2_contrasts.R` | `deseq2/contrasts/` (10 contrasts + DEG summary) |
| 4 | DESeq2 factorial (V4) | `Rscript deseq2/run_factorial_v4.R` | `deseq2/v4/` (interactions, PC–factor) |
| 5 | GO / KEGG | `Rscript go_analysis/run_go_clusterprofiler.R --deg-glob "deseq2/contrasts/Diff_genes_heatmap_*.csv" --out go_analysis/results_full` | `go_analysis/results_full/` |
| 6 | PhysioSpace | `Rscript physiospace/run_physiospace.R --static TRUE --out physiospace/results_static` | `physiospace/results_static/` |
| 7 | WGCNA | `Rscript wgcna/run_wgcna_root.R` · `Rscript wgcna/annotate_modules_go.R` | `wgcna/results/` (+ module names) |
| 8 | Integration | `Rscript integration/integrate_root_expression.R` · `Rscript integration/make_named_heatmap.R` | `integration/results/`, Fig 4 |
| 9 | Figures | `Rscript deseq2/make_pca.R` · `python manuscript/figures/make_figures.py` | `manuscript/figures/` (Fig 1,2,3,5,6,7) |
| 10 | Word draft | `cd manuscript && npm install docx && node build_docx.js` | `manuscript/TICTOC_manuscript_draft.docx` |

## Notes
- **PhysioSpace:** use `--static TRUE`; the default signed-p score saturates root contrasts to ±Inf.
- **WGCNA:** if `blockwiseModules` errors with "unused arguments (weights.x …)", another package masked
  `cor` — the script sets `cor <- WGCNA::cor` before it; keep that line.
- **Integration is group-level (n=6).** For individual-level (n≈24) analysis, supply a library→plant
  manifest and run `integration/pair_rnaseq_to_images.py` (see `integration/PAIRING_WHATS_NEEDED.md`).
- Raw reads and images are on NASA OSDR/GeneLab (released on publication); this repo carries the curated
  count matrix, design table, and curated RSML set needed to reproduce everything downstream.
