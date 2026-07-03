# GO / KEGG enrichment (`clusterProfiler`)

Enriches the TICTOC cotton DEGs for GO terms and KEGG pathways by translating `Gohir` genes to
Arabidopsis (via [`../crosswalk/gohir_to_arabidopsis.tsv`](../crosswalk/gohir_to_arabidopsis.tsv))
and running `clusterProfiler` against the richly-annotated Arabidopsis databases. This is roadmap
step **4.4** and the on-ramp to PhysioSpace (4.5).

## Run
```bash
cd go_analysis
Rscript run_go_clusterprofiler.R          # uses sensible defaults (see below)
```
First-time install:
```r
install.packages("BiocManager")
BiocManager::install(c("clusterProfiler","org.At.tair.db","enrichplot"))
install.packages("ggplot2")
```

## What it does
1. Loads the `Gohir → AT` crosswalk (optionally filtered by `--min-pid`).
2. Builds the **background universe** = AT loci reachable from the *expressed* cotton genes
   (from the count matrix) — not the whole genome. This is essential for correct enrichment.
3. For every `Diff_genes_heatmap_*.csv` DEG table, splits genes into **up / down / all**, maps to
   unique AT loci, and runs `enrichGO` (BP, MF, CC) + `enrichKEGG` (organism `ath`).
4. Writes per-contrast, per-ontology CSVs and dotplot PDFs to `results/`.

## Inputs (defaults, overridable via flags)
| Flag | Default | Notes |
|---|---|---|
| `--deg-glob` | `../TICTOC_3_factor_model/Diff_genes_heatmap_*.csv` | DEG tables to process |
| `--crosswalk` | `../crosswalk/gohir_to_arabidopsis.tsv` | Gohir→AT map |
| `--counts` | `../TICTOC_run1_filteredCounts_v3.csv` | defines the universe |
| `--out` | `results` | output dir |
| `--use-regulation` | `TRUE` | trust the DEG table's `Regulation` (Up/Down) column; set `FALSE` to re-threshold |
| `--padj` / `--lfc` | `0.05` / `1` | thresholds used only when `--use-regulation FALSE` |
| `--min-pid` | `30` | drop crosswalk hits below this % identity |
| `--kegg` | `TRUE` | KEGG needs internet; set `FALSE` for offline GO-only |

## Caveats (also in the script header)
- **Best-BLASTP-hit orthology is approximate.** Many Gohir → one AT (A/D homoeologs + paralogs), so
  gene counts collapse on mapping; enrichment is over *unique* AT loci. Raise `--min-pid` for a
  stricter mapping.
- **Keep the universe.** Do not remove the `universe=` argument — enriching a mapped subset against
  the whole Arabidopsis genome inflates significance.
- Not yet executed against a live R install in this repo — run once and sanity-check the mapping-rate
  messages it prints before trusting the tables.

## Alternative (GO directly in cotton space)
If you prefer to avoid the ortholog step for GO, CottonGen ships `Gh_TM1_UTX_v2.1_genes2Go.xlsx.gz`
(Gohir→GO). The Arabidopsis route is still required for **PhysioSpace**, which is Arabidopsis-based.
