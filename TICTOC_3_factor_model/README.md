# TICTOC_3_factor_model/ — RNA-seq analysis (3-factor DESeq2 / iDEP)

Differential-expression analysis of the TICTOC cotton RNA-seq under the 3-factor design
**Treatment (Flight/Ground) × Genotype (WT/A68/D130) × Tissue (Root/Shoot)**, adapted from the
iDEP 0.95 workflow. Gene IDs are Phytozome `Gohir.*` (see `../DATA_DICTIONARY.md`).

## Analysis source (canonical)
| File | Role |
|---|---|
| `TICTOC_markdown_evolved_v2.Rmd` | **canonical pipeline** — titled "TICTOC R Markdown v3", adapted from iDEP 0.95 (dated 2022-10-21). Renders to the v3 report below. Paths de-hardcoded 2026-07-03. |

Required external inputs (not committed — see the note inside the `.Rmd` and `../DATA_DICTIONARY.md`):
`Downloaded_Converted_Data.csv`, `iDEP_core_functions.R`, `Cotton__ghirsutum_eg_gene_GeneInfo.csv`,
`Cotton__ghirsutum_eg_gene.db`.

## Rendered reports (version history)
| File | Version | Status |
|---|---|---|
| `markdown reports/TICTOC_markdown_v3.html` | **v3** (current) | ✅ canonical rendered report (from the `.Rmd` above) |
| `markdown reports/TICTOC_iDEP_R_Markdown_v2.html` | v2 (iDEP) | superseded — earlier iDEP export, kept for provenance |

> A byte-identical duplicate of the v3 report previously sat at the repo root; removed 2026-07-03.
> The README's report link points here, to `markdown reports/`.

## Derived output tables
| File | Content |
|---|---|
| `Diff_genes_heatmap_*.csv` | per-contrast DEG tables with expression. Contrast is in the filename: `A68-WT (Root_GC)` = A68 vs WT in Root/Ground; `I_Treatment_Flight.Genotype_A68` = the Treatment[Flight]×Genotype[A68] interaction term. Columns per `../DATA_DICTIONARY.md` §4. |
| `Enriched.csv` | over-representation enrichment (Direction, adj.Pval, nGenes, Pathways, Genes) |
| `D130-WT(GAGE).csv` | GAGE gene-set analysis for the D130-vs-WT contrast |
| `AllGeneListsGMT.gmt` | gene sets (GMT) for enrichment; IDs are upper-case `GOHIR.*` |
| `enrichmentPlotDEG2.eps` | enrichment figure (vector) |

## Notes
- These outputs cover a subset of the full contrast grid. When the models are frozen (roadmap §4.1),
  regenerate a complete, documented set of contrasts from the single canonical `.Rmd`.
- For GO/pathway work in Arabidopsis space, map `Gohir.*` via `../crosswalk/gohir_to_arabidopsis.tsv`
  (the older raimondii-based `../Ara_vs_Cotton_biomart_export` does not join).
