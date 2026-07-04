# TICTOC — Project Roadmap, FAIR & Publication Tracker

**Targeting Improved Cotton Through Orbital Cultivation (TIC-TOC)**
Gilroy Lab, University of Wisconsin–Madison · CASIS grant UA-2018-276 · flown on SpaceX CRS-22 to the ISS.

> This is the **living working document** for finishing the analysis, making the repo FAIR, and
> preparing the Zenodo deposit + npj Microgravity manuscript. It sits alongside the public GitBook
> landing page (`README.md`), which is left untouched. Update the checkboxes as items land.
>
> Last reviewed: 2026-07-03 · Maintainer: R. Barker

---

## 0. Progress log

**2026-07-03 — FAIR pass + analysis-pipeline scaffold** (this session):
- ✅ Added this roadmap, [`DATA_DICTIONARY.md`](DATA_DICTIONARY.md), `.gitignore`, `CITATION.cff`;
  README now links here.
- ✅ **Built the corrected `Gohir → Arabidopsis` ortholog crosswalk** ([`crosswalk/`](crosswalk/)) —
  the GO/PhysioSpace blocker. CottonGen BLASTP best-hit, **92.3%** count-matrix coverage, reproducible builder.
  (The old `Ara_vs_Cotton_biomart_export` was *G. raimondii* and didn't join.)
- ✅ Scaffolded the analysis pipeline: [`go_analysis/`](go_analysis/) (clusterProfiler GO/KEGG) and
  [`physiospace/`](physiospace/) (PhysioSpace decoding, wired to the `PlantPhysioSpace` data package).
- ✅ Cleanup: de-duped leaf zips + duplicate v3 report (~40 MB), pruned 32 redundant RSML scratch files
  (kept a unique timelapse + `A68_2A/2B` images), de-hardcoded the analysis `.Rmd`, documented the
  DEG reports and RSML QC ([`Data/RSML_QC_summary.md`](Data/RSML_QC_summary.md)), pinned genome versions.
- ✅ Drafted the [`manuscript/`](manuscript/) skeleton (npj Microgravity): Results stubs Fig 1–7,
  concrete Methods, and a figure/table plan with draft legends.
- ✅ **Extracted RSML morphometric traits + ran the mixed model** ([`morphometrics/`](morphometrics/)) —
  real per-plant-day table (198 rows) + group means, and a plant-random-effect model on log root length.
  Traced-data result: Flight roots elongate faster (day×condition p≈1e-31) and **AVP-OX lines amplify the
  flight response vs WT** (A68 ×2.34 p=0.032; D130 ×2.04 p=0.037) — the paper's thesis. Faster-growth and
  AVP-OX effects are robust to *uniform* miscalibration; ⚠ still verify image scale + lateral-tracing artifact.

**Still needs the maintainer:** a live R run of the GO + PhysioSpace scripts; OSDR accession; Zenodo DOI;
author ORCIDs / co-authors; `sessionInfo`; static figure exports for the dead Slides links; and the
scientific call on canonical DEG contrasts (§4.1) + whether `A68_2A/2B` belong in the analysis.

---

## 1. What this experiment is

Cotton (*Gossypium hirsutum*, `Gohir` genome) was grown on the ISS and on the ground to test whether
the **AVP1-overexpressing (AVP-OX)** lines — engineered for salt/drought tolerance and larger root
systems (Park et al. 2005) — alter their root architecture and transcriptome in microgravity.

**Full factorial design captured in the data:**

| Factor | Levels |
|---|---|
| **Treatment** | Flight (FL, ISS) vs Ground Control (GC) |
| **Genotype** | WT · A68 (AVP-OX) · D130 (AVP-OX) |
| **Tissue** | Root · Shoot/Leaf |
| **Replicates** | 4 per group → **48 RNA-seq libraries** |
| **Time series** | imaging days 3, 4, 5, 6 (root photography + RSML tracings) |

**Three measurement modalities:**
1. **Transcriptomics** — bulk RNA-seq, ~60 k `Gohir` genes, analysed with a 3-factor model (DESeq2 / iDEP).
2. **Root morphometrics** — SmartRoot tracings exported as **RSML** ("root markup language") from time-series ISS/ground photographs.
3. **Time-series root photography** — grayscale, inverted, 300 dpi TIFFs, days 3–6, incl. clarity-enhanced ISS downlink frames.

**Raw data:** deposited on NASA **OSDR/GeneLab** (embargoed; released on manuscript publication).

---

## 2. Repository map (Findability)

```
TICTOC/
├── README.md ....................... GitBook public landing page (abstract, ASGSR deck, video)
├── PROJECT_ROADMAP.md .............. ← THIS FILE: goals, FAIR & publication tracker
├── DATA_DICTIONARY.md .............. column/file schemas for every dataset
├── CITATION.cff · .gitignore · LICENSE (CC0-1.0)
│
├── RNA-seq — counts & design
│   ├── TICTOC_run1_filteredCounts_v3.csv .... filtered count matrix (~59.9k genes × 48)
│   ├── TICTOC_cotton_run2_csv/ ............... UNfiltered counts (+ README: name-mismatch note)
│   ├── TICTOC_target_v5.csv ................. sample/design table (Treatment×Genotype×Tissue)
│   └── Ara_vs_Cotton_biomart_export.txt.zip . ⚠ legacy raimondii map — superseded by crosswalk/
│
├── crosswalk/ ...................... ✅ Gohir→Arabidopsis ortholog map (TSV + builder + README)
├── morphometrics/ ................. ✅ RSML→traits extraction (rsml_traits.csv + summary + scripts)
├── go_analysis/ ................... ✅ clusterProfiler GO/KEGG scaffold (run_go_clusterprofiler.R)
├── physiospace/ ................... ✅ PhysioSpace decoding scaffold (run_physiospace.R)
├── manuscript/ .................... ✅ npj draft skeleton + figure/table plan with legends
│
├── RNA-seq — analysis
│   └── TICTOC_3_factor_model/ ...... DESeq2/iDEP pipeline + DEG/enrichment outputs (see its README;
│                                     canonical source = TICTOC_markdown_evolved_v2.Rmd → v3 report)
│
├── Root morphometrics (RSML)
│   └── Data/
│       ├── Final_RSML_format/ ...... ★ curated RSML set ("outlier removed") — USE THIS
│       ├── RSML_QC_summary.md ...... ✅ QC inventory (198 tracings, 53 plants, 44 full series)
│       └── Cotton {FL,GC} day{3..6} .../ .... raw per-day TIFFs + tracings (scratch pruned;
│                                     retains unique timelapse .avi + A68_2A/2B images — see QC doc)
│
├── Leaves/ ......................... leaf/shoot DESeq2 (interaction + linear), WGCNA, OSD cross-refs
│                                     (duplicate zips in TITCO_leaves_2024/ removed)
├── Figures/ ........................ root-phenotype figures (readme + external links)
└── *.pdf, drbs-*.md, data-1.md ..... ASGSR/CASIS decks & narrative stubs
```

Analysis pipeline order: `crosswalk/` → DEGs (`TICTOC_3_factor_model/`) → `go_analysis/` → `physiospace/`.

---

## 3. Status at a glance

| Workstream | State | Notes |
|---|---|---|
| RNA-seq counts (run 1) | ✅ done | filtered + unfiltered matrices present |
| 3-factor DESeq2 model | ✅ done (needs consolidation) | multiple report versions (v2/v3/iDEP) coexist |
| Root RSML tracings | ✅ done | `Final_RSML_format/` curated; raw days 3–6 retained |
| Root morphometric statistics | 🟡 partial | correlation plots exist; formal FL-vs-GC × genotype stats to finalise |
| Time-series imagery | ✅ captured | days 3–6, FL + GC |
| Leaf/shoot DEG + WGCNA | 🟡 partial | outputs present, not yet written up |
| Tissue-specific DEG clustering | ⬜ to do | §4 |
| GO/pathway enrichment (ortholog-mapped) | 🟡 partial | GAGE/GMT outputs exist; redo cleanly per tissue |
| PhysioSpace stress decoding | ⬜ to do | §4 — mirror the OSD-767 tomato pipeline |
| FAIR repo cleanup | ⬜ to do | §5 |
| Zenodo deposit | ⬜ to do | §6 |
| npj Microgravity manuscript | ⬜ to do | §7 |

Legend: ✅ done · 🟡 partial · ⬜ not started

---

## 4. Analysis roadmap (the science still to do)

Target endpoint: an integrated **root-architecture + transcriptome** story, decoded into stress
programs — the same arc validated on the OSD-767 tomato paper, now applied to AVP-OX cotton.

- [ ] **4.1 Freeze the DEG models.** Pick ONE canonical DESeq2 pipeline from the v2/v3/iDEP variants;
      archive the rest under `_archive/`. Export a documented design matrix and contrast list:
      - Main effects: Flight vs Ground; Genotype (A68, D130 vs WT); Root vs Shoot.
      - Interactions: **Treatment × Genotype** (does AVP-OX change the spaceflight response?) and
        **Treatment × Tissue** (root- vs shoot-specific flight response).
- [ ] **4.2 Organ-specific DEG sets.** Derive root-only and shoot-only DEG lists per contrast, with
      shrunken LFCs and padj, as tidy CSVs (one row per gene, columns per contrast).
- [ ] **4.3 Tissue-specific clustering.** Cluster organ-specific DEGs (WGCNA is already started for
      leaves; extend to root and to a joint set). Annotate modules by eigengene–trait correlation
      (Flight, genotype, and **root morphometric traits**). ✅ RSML traits now extracted to
      [`morphometrics/rsml_traits.csv`](morphometrics/rsml_traits.csv) (per plant-day) ready for the
      module–trait correlation; ⚠ verify FL/GC length calibration first (see `morphometrics/README.md`). *(2026-07-03)*
- [ ] **4.4 GO / pathway analysis per cluster.** Map `Gohir` → Arabidopsis and run GO/KEGG enrichment on
      each module/DEG set. Standardise on one tool (clusterProfiler recommended for reproducibility over
      the iDEP GAGE outputs). ✅ **Blocker resolved:** built [`crosswalk/gohir_to_arabidopsis.tsv`](crosswalk/gohir_to_arabidopsis.tsv)
      from the CottonGen *G. hirsutum* v2.1 → Arabidopsis BLASTP best-hit (92.3% count-matrix coverage;
      reproducible via `crosswalk/build_gohir_to_arabidopsis.py`). Use this instead of the *raimondii*
      `Ara_vs_Cotton_biomart_export`. ✅ **Script scaffolded:**
      [`go_analysis/run_go_clusterprofiler.R`](go_analysis/run_go_clusterprofiler.R) runs enrichGO
      (BP/MF/CC) + enrichKEGG over the DEG tables with a proper expressed-gene universe. Needs a first
      run against a live R install to validate. *(2026-07-03)*
- [~] **4.5 PhysioSpace stress-pattern decoding.** Project cotton contrasts onto the plant PhysioSpace
      stress compendium (Hadizadeh Esfahani et al.) to produce PhysioScores per genotype/tissue.
      ✅ **Script scaffolded:** [`physiospace/run_physiospace.R`](physiospace/run_physiospace.R) builds
      VST Flight−Ground contrasts, re-indexes Gohir→AT→Entrez via the crosswalk, and runs
      `calculatePhysioMap` (GenesRatio=0.05, TTEST=FALSE, ImputationMethod="PCA") — same params as OSD-767.
      Reference spaces load directly from the **`PlantPhysioSpace`** data package
      (`AT_Stress_Space` / `_Meta` / `_RNASeq`); no manual `.rds` needed. **Needs:** a live R run to
      validate (install `PhysioSpaceMethods` + `PlantPhysioSpace` from JRC-COMBINE GitHub). *(2026-07-03)*
- [ ] **4.6 Root morphometrics ↔ transcriptome integration.** Formal stats on RSML traits
      (mixed model: trait ~ Treatment × Genotype, day as time covariate) and correlate trait modules
      with expression modules — the headline integrative figure.

---

## 5. FAIR & repository hygiene

**Findable / Accessible**
- [ ] Add this roadmap link + a one-paragraph "repo contents" section to the public `README.md`.
- [ ] Reserve a **Zenodo DOI** and add the badge to `README.md` once minted (§6).
- [ ] Record the OSDR/GeneLab accession number(s) and cross-link both ways.

**Interoperable**
- [x] Add a **data dictionary** ([`DATA_DICTIONARY.md`](DATA_DICTIONARY.md)): column-level docs for the
      design/count/DEG/enrichment/WGCNA tables, RSML trait schema, imagery naming, ortholog map, and the
      sample-ID grammar (`{Genotype}_{plate}{well}_{FL|GC}_{day}`). *(2026-07-02)*
- [ ] State genome build/annotation version for `Gohir` IDs and the BioMart release used for orthologs.
- [ ] Keep RSML as-is (open XML standard) — good; just document which tracer/version produced them.

**Reusable**
- [ ] Add an **environment spec** (`renv.lock` or `sessionInfo.txt` + R version) for the pipeline.
- [x] De-hard-code paths: removed the absolute `setwd('/Users/drbhomeoffice/...')` from
      `TICTOC_markdown_evolved_v2.Rmd` (R Markdown knits from the file's own dir) and documented the
      external inputs it needs (`Downloaded_Converted_Data.csv`, `iDEP_core_functions.R`,
      `Cotton__ghirsutum_eg_gene_GeneInfo.csv`, `Cotton__ghirsutum_eg_gene.db`). Still to do: commit or
      link those missing input files. *(2026-07-02)*
- [ ] Clarify licensing: `LICENSE` is CC0. Confirm CC0 (data) + a code license (e.g. MIT) split, or state CC0 covers all.

**Cleanup (the "lots of folders/code/graphs" pass)**
- [~] Consolidate DEG models → keep one, archive the rest (§4.1). ✅ Documented in
      [`TICTOC_3_factor_model/README.md`](TICTOC_3_factor_model/README.md): canonical source is
      `TICTOC_markdown_evolved_v2.Rmd` → v3 report; iDEP v2 kept for provenance. Removed a byte-identical
      duplicate of the v3 report (~17 MB) from the repo root. Full contrast regen still pending (§4.1). *(2026-07-03)*
- [~] Resolve the mislabeled `TICTOC_cotton_run2_csv/` (holds `run1_allCounts`) — added an explanatory
      `README.md` in the folder; **rename decision pending maintainer** (was there ever a real run 2?).
- [x] De-duplicate `Leaves/` vs `Leaves/TITCO_leaves_2024/` — removed the two byte-identical DEG-value
      zips from the `TITCO_leaves_2024/` data-dump subfolder (~22 MB); canonical copies kept in
      top-level `Leaves/`. *(2026-07-02)*
- [x] Prune RSML scratch dirs under `Data/` and settle on `Final_RSML_format/` as canonical. ✅ QC in
      [`Data/RSML_QC_summary.md`](Data/RSML_QC_summary.md) (198 tracings, 53 plants, 44 full day-3→6
      series; unbalanced design). ✅ Removed 32 blob-SHA-verified redundant files (~35 MB) from the
      `test/`, `Convert late/`, `delete me_templates/` etc. sub-folders. ⚠ **Kept** 10 unique-content
      files (~107 MB) that were mislabeled as scratch — a timelapse `.avi` and `A68_2A/2B` day-5 images;
      see the QC doc for the maintainer decision on those. *(2026-07-03)*
- [~] Move the large rendered `*.html` reports + slide-deck PDFs to the Zenodo release or Git LFS;
      **added a `.gitignore`** (2026-07-02) for OS/R/Python cruft, with commented patterns ready to
      enable once the existing rendered reports are relocated. Relocation still pending.
- [ ] Replace dead/edit-only Google Slides links in `README.md`/`data-1.md`/`Figures/readme.md` with
      exported static figures committed to `Figures/`.

---

## 6. Zenodo deposit checklist

- [ ] Reserve DOI (Zenodo → "reserve" before publishing) and add badge to `README.md`.
- [ ] Assemble the release bundle: canonical count matrices, design table, curated RSML set,
      final DEG/GO/PhysioSpace tables, figure source files, and the frozen analysis code.
- [x] Write `CITATION.cff` at repo root *(2026-07-02)* — authors/keywords/license/repo populated;
      **TODO in-file:** add ORCIDs, co-authors, reserved DOI, and version/date on release.
- [ ] Fill Zenodo metadata: authors + ORCIDs, funding (**CASIS UA-2018-276**), keywords
      (spaceflight, cotton, AVP1, root architecture, RSML, transcriptomics), related identifiers
      (OSDR accession, manuscript DOI once known).
- [ ] Tag a GitHub release (`v1.0.0`) to trigger the Zenodo archive.
- [ ] Confirm data license (CC0/CC-BY) and embargo lift aligns with manuscript publication + OSDR release.

---

## 7. npj Microgravity manuscript checklist

**Narrative spine:** AVP-OX cotton in spaceflight — does the salt/drought-tolerance engineering change
the root-architecture and transcriptomic response to microgravity, and what stress programs does it
engage? Integrates morphometrics (RSML) + transcriptome (DEG→cluster→GO→PhysioSpace).

- [~] **Skeleton drafted** ([`manuscript/manuscript_tictoc.md`](manuscript/manuscript_tictoc.md)):
      title, abstract template, intro, Results stubs (Fig 1–7), concrete Methods, data-availability.
      Figure/table plan with draft legends in [`manuscript/figures_and_tables.md`](manuscript/figures_and_tables.md). *(2026-07-03)*
- [ ] **Outline & authorship** — author list, affiliations, ORCIDs, contributions, corresponding author.
- [ ] **Figures** (draft list — refine as analysis lands):
  - [ ] F1. Experiment overview + design schematic (genotype × treatment × tissue × time).
  - [ ] F2. Root architecture over days 3–6: RSML traits, Flight vs Ground × genotype.
  - [ ] F3. Transcriptome overview: PCA, DEG counts per contrast, root vs shoot.
  - [ ] F4. Tissue-specific DEG clusters / WGCNA modules ↔ traits.
  - [ ] F5. GO/pathway enrichment per module (ortholog-mapped).
  - [ ] F6. PhysioSpace stress-pattern decoding across genotype × tissue.
  - [ ] F7. Integrative model: root-trait modules ↔ expression modules ↔ stress programs.
- [ ] **Tables:** T1 sample/design summary; T2 DEG counts per contrast; T3 top enriched terms; T4 PhysioScores.
- [ ] **Figure legends** — write one per figure/table (self-contained, methods-cited).
- [ ] **Supplementary:** full DEG tables, GO/KEGG tables, PhysioScore matrices, RSML raw traits, QC/MultiQC.
- [ ] **Methods:** growth/hardware, RNA-seq pipeline + genome/annotation versions, DESeq2 model &
      contrasts, ortholog mapping, WGCNA params, GO tool, PhysioSpace reference & version, RSML/SmartRoot
      params, statistical models for morphometrics.
- [ ] **Data & code availability:** OSDR accession + Zenodo DOI + this GitHub repo.
- [ ] **Cover letter** + suggested reviewers; confirm npj figure spec (TIFF, resolution) and reference style.

---

## 8. Data availability statement (draft)

> Raw sequencing reads, processed count matrices, root images, and RSML tracings are archived at the
> NASA Open Science Data Repository (OSDR/GeneLab), accession **OSD-XXX** (released on publication).
> Analysis code and derived tables are archived at Zenodo (**DOI: 10.5281/zenodo.XXXXXXX**) and
> developed at https://github.com/dr-richard-barker/TICTOC.

---

### Related projects (cross-reference)
The DEG → tissue-clustering → GO → **PhysioSpace** arc mirrors the OSD-767 tomato pipeline; reuse those
`physioworker`/enrichment scripts as templates here rather than rebuilding from scratch.
