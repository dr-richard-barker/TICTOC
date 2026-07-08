# Supplementary materials index — TICTOC

All supplementary items exist in the repository and are archived with the Zenodo release. Paths are
repo-relative.

## Supplementary data
| ID | Item | Path |
|---|---|---|
| S1 | Full DEG tables, all 10 contrasts (per-gene LFC + padj) | `deseq2/contrasts/Diff_genes_heatmap_*.csv` |
| S2 | DEG-count summary | `deseq2/contrasts/DEG_counts_summary.csv` |
| S3 | Factorial interaction model (interaction DEG counts, PC–factor) | `deseq2/v4/` |
| S4 | GO/KEGG enrichment tables, all contrasts | `go_analysis/results_full/*.csv` |
| S5 | Per-module GO (BP/MF/CC) + module names | `wgcna/results/GO_*_module_*.csv`, `module_names.csv` |
| S6 | WGCNA module assignments + module–trait correlations | `wgcna/results/module_*.csv` |
| S7 | PhysioScore matrices (3 reference spaces) | `physiospace/results_static/*.csv` |
| S8 | RSML per-plant morphometric traits + group means | `morphometrics/rsml_traits*.csv` |
| S9 | Morphometric mixed-model coefficients (total + primary) | `morphometrics/*_stats_*` |
| S10 | Integration: GS/MM, hub-growth genes, module summary | `integration/results/*.csv` |
| S11 | *Gohir*→Arabidopsis ortholog crosswalk | `crosswalk/gohir_to_arabidopsis.tsv` |
| S12 | RSML QC inventory | `Data/RSML_QC_summary.md` |
| S13 | RNA-seq QC (MultiQC) | *[raw/trimmed MultiQC reports — add from run-1 pipeline]* |

## Supplementary figures (candidates)
- Named module–trait heatmap (also main Fig 4): `manuscript/Fig4_module_trait_named.pdf`
- Per-contrast GO dotplots: `go_analysis/results_full/GO_*.pdf`
- PhysioScore heatmaps per reference space: `physiospace/results_static/*.pdf`
- Primary-root sensitivity model summary: `morphometrics/primary_stats_summary.txt`

## Methods reproducibility
- Environment / package versions: `ENVIRONMENT.txt`
- Every result is regenerable from the committed scripts (`deseq2/`, `go_analysis/`, `physiospace/`,
  `wgcna/`, `morphometrics/`, `integration/`); data dictionary: `DATA_DICTIONARY.md`.
