# TICTOC — figures & tables plan (with draft legends)

Planned display items for the npj Microgravity manuscript (matches roadmap §7 and the Results stubs in
`manuscript_tictoc.md`). Each entry lists **what it shows**, the **source** (script/data that produces
it), and a **draft self-contained legend**. Legends describe intended content — finalise numbers only
from executed outputs. npj spec: TIFF, ≥300 dpi (roadmap §7).

## Main figures

**Fig 1 — Experiment overview & design.** *(schematic + summary)*
Source: hand-drawn schematic + `TICTOC_target_v5.csv`, `Data/RSML_QC_summary.md`.
Legend (draft): *TIC-TOC experimental design. (a) Wild-type and two AVP1-overexpressing cotton lines
(A68, D130) were grown aboard the ISS (SpaceX CRS-22) alongside ground controls. (b) Root systems were
imaged on days 3–6 and traced (RSML); root and shoot tissues were profiled by RNA-seq in a full
Treatment × Genotype × Tissue design (n = 4/group). (c) Sample and plant counts per group.*

**Fig 2 — Root architecture over days 3–6.** *(RSML morphometrics)* ✅ **generated: `figures/Fig2_root_architecture.pdf`**
Source: `../morphometrics/rsml_traits.csv`. Total length (top) vs primary length (bottom) × genotype; Flight/Ground ± SE.
Legend (draft): *Root-system architecture of AVP-OX vs wild-type cotton under spaceflight. Trait
trajectories (e.g. total root length, mean diameter) across days 3–6 for Flight vs Ground, by genotype;
points = plants, lines = model fits. Statistics: [mixed model]; unbalanced n per QC.*

**Fig 3 — Transcriptome overview.** *(PCA + DEG counts)*
Source: `../deseq2/` (executed model). ✅ **PCA panel generated: `Fig3_PCA.pdf`** (PC1=96%=tissue); DEG counts in `../deseq2/contrasts/DEG_counts_summary.csv`.
✅ **Fig 3b generated: `figures/Fig3b_DEG_counts.pdf`** (diverging up/down bars per genotype×tissue).
Legend (draft): *Spaceflight transcriptome. (a) PCA (PC1=tissue, PC2=treatment). (b) DEG counts per genotype×tissue, up (orange) vs down (blue).*

**Fig 4 — Root co-expression modules ↔ traits.** *(WGCNA)* ✅ **generated: `Fig4_module_trait_named.pdf`**
Source: `../wgcna/` + `../integration/`; modules named by GO (`../wgcna/results/module_names.csv`).
Legend (draft): *Co-expression modules and their trait associations. Module eigengene–trait correlations
for Flight, genotype, and root morphometric traits; heatmap of representative module expression.*

**Fig 5 — Functional enrichment.** *(GO/KEGG per module)*
Source: `go_analysis/run_go_clusterprofiler.R` (dotplots).
Legend (draft): *GO/KEGG enrichment of spaceflight-responsive gene sets, ortholog-mapped to Arabidopsis
(CottonGen BLASTP; enriched against the expressed-gene universe). Dot size = gene count, colour = q-value.*

**Fig 6 — PhysioSpace stress decoding.** *(PhysioScore heatmap)*
Source: `physiospace/run_physiospace.R` (`PhysioHeatmap`).
Legend (draft): *PhysioSpace stress-pattern decoding. PhysioScores (signed −log₁₀ p) for each
genotype × tissue Flight−Ground contrast across Arabidopsis stress axes (AT_Stress_Space); positive =
activation of that stress program. Compares WT vs AVP-OX programs.*

**Fig 7 — Integrative model.** *(synthesis schematic)*
Source: hand-drawn from R2–R6 results.
Legend (draft): *Integrative model linking root-trait modules, expression modules, and stress programs,
summarising how AVP-OX engineering modulates the cotton spaceflight response.*

## Main tables
- **T1 — Sample/design summary** (genotype × treatment × tissue × reps; plant/day counts). Source: `TICTOC_target_v5.csv`, RSML QC.
- **T2 — DEG counts per contrast** (up/down, thresholds). Source: frozen DESeq2 model.
- **T3 — Top enriched GO/KEGG terms** per contrast. ✅ **generated: `Table3_top_GO_terms.csv`** (from `../go_analysis/results_full/`).
- **T4 — PhysioScores** per genotype × tissue × stress axis. Source: `physiospace/results/`.

## Supplementary (roadmap §7)
Full DEG tables · GO/KEGG tables · PhysioScore matrices · RSML raw traits · MultiQC/QC · crosswalk TSV.

## Production notes
- Keep a single figure-build script per panel where possible; commit source + rendered PNG, put heavy
  TIFFs in the Zenodo release (see `.gitignore`).
- Cross-check every figure/table is cited in the text before submission (an npj/reviewer flag).
