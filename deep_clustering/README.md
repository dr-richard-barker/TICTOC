# Deep-learning expression-program analysis

Autoencoder-based gene co-expression **programs** from the bulk RNA-seq, annotated by GO and Arabidopsis
cell-type markers, with a per-program spaceflight stress response compared **WT vs AVP-OX** (manuscript
R8, Fig 8 & Fig 9).

> ⚠ **This is BULK RNA-seq (48 libraries), not single-cell data.** We do **not** cluster cells. An
> autoencoder learns a latent gene embedding and *k*-means clusters **genes** into co-expression programs
> (a neural analogue of WGCNA). Cell-type labels are **marker-overlap associations** from a curated
> Arabidopsis panel — an interpretive aid, never a claim of single-cell identity. A true cell-resolved
> analysis would require single-nucleus/cell RNA-seq or spatial data.

## Pipeline (run from repo root)
```bash
python deep_clustering/01_autoencoder_programs.py      # AE -> 12 gene programs (torch)      -> results/
python deep_clustering/02_celltype_markers.py          # curated AT cell-type marker overlap -> results/
Rscript deep_clustering/03_program_go_export.R         # GO per program + export space/entrez-> results/
python deep_clustering/04_program_stress.py            # per-program stress response, WT vs AVP-OX
python deep_clustering/05_figures.py                   # Fig 8 (atlas) + Fig 9 (WT vs AVP-OX)
```
Deps: `torch`, `scikit-learn`, `umap-learn`, `numpy`, `pandas`, `scipy` (Python); `clusterProfiler`,
`org.At.tair.db`, `PlantPhysioSpace` (R). Seeds fixed for reproducibility.

## Method
1. **Autoencoder** (48→32→12→32→48) on the top-8,000 variable genes (log-CPM, per-gene z), latent codes
   *k*-means-clustered into 12 programs; UMAP for the map.
2. **GO** (BP) per program; **cell-type** marker overlap (curated canonical Arabidopsis root/leaf markers).
3. **Stress response**: correlation of each program's Flight−Ground fold change (dominant tissue,
   Gohir→AT→Entrez) with each Arabidopsis `AT_Stress_Space` axis over shared genes (≥15), per genotype.

## Result
Programs partition by tissue and map to interpretable functions (P1 lignin/endodermis, P5 cell cycle,
P2/P11 defence, P10 photosynthesis/mesophyll, P0 guard-cell/osmotic). The **AVP-OX lines attenuate the
stress/defence programs** (P2, P11, P4, shoot osmotic/hormone; Δ < 0) while **sustaining growth/metabolic
programs** (cell-cycle P5, nucleotide P9; Δ > 0) — localising the whole-transcriptome "less stress, more
growth" thesis to specific functional programs. Key tables in `results/` (`program_GO.csv`,
`program_celltypes.csv`, `program_stress_scores.csv`, `program_stress_WTvsAVPOX.csv`).

## Caveats
- Bulk data → programs, not cells (above). Cell-type panel is small/curated; treat as suggestive.
- Per-program stress score is a PhysioSpace-style correlation over gene sets, not the exact
  `calculatePhysioMap` output; small programs (few shared genes) are less reliable.
