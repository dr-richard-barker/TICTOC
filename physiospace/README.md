# PhysioSpace stress-pattern decoding

Projects the TICTOC cotton Flight-vs-Ground contrasts onto **Arabidopsis stress reference spaces**
to quantify which known stress programs the spaceflight response resembles (PhysioScores). This is
roadmap step **4.5** and mirrors the validated OSD-767 tomato pipeline.

## Run
```bash
cd physiospace
Rscript run_physiospace.R --ref-dir /path/to/reference_spaces
```
Install (once):
```r
install.packages(c("BiocManager","remotes","pheatmap"))
BiocManager::install(c("DESeq2","org.At.tair.db"))
remotes::install_github("JRC-COMBINE/PhysioSpaceMethods")
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

## You must supply the reference spaces
Not committed here (they are the shared Arabidopsis stress compendia from Hadizadeh Esfahani et al.
2021, the same files used for OSD-767):

- `AT_Stress_Space`, `AT_Stress_Space_Meta`, `AT_Stress_Space_RNASeq`

Provide them as `.rds` matrices (rownames = Arabidopsis **Entrez** IDs, columns = stress axes) in a
folder passed via `--ref-dir`. Reuse the OSD-767 copies rather than rebuilding.

## Caveats
- Not yet executed against a live R install — validate the printed gene-mapping counts and the
  "shared Entrez genes with input" line (expect a few thousand) before interpreting scores.
- The cotton→Arabidopsis bridge is best-BLASTP-hit orthology (approximate); raise `--min-pid` for a
  stricter mapping. OSD-767 used Ensembl Compara 1:1 orthologs — expect somewhat noisier cotton scores.
- Fold changes here are VST group-mean differences (blind VST). If you freeze a DESeq2 model
  (roadmap §4.1) with shrunken LFCs, swap those in as the input matrix for consistency with the DEG tables.

## Options
`--counts` (default `../TICTOC_run1_filteredCounts_v3.csv`), `--crosswalk`, `--out` (default `results`),
`--min-pid` (default 30), `--genes-ratio` (default 0.05).
