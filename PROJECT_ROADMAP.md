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
- ✅ **Sensitivity check (primary root only):** AVP-OX × flight interaction **holds on primary-root length**
  (A68 ×1.84 P=0.04; D130 ×1.84 P=0.017) — robust, since the primary root is equally imaged in FL/GC. The
  big total-length flight gain is **lateral-driven** (primary elongation not accelerated) → verify before claiming.

**2026-07-04 — analyses executed (R 4.6 installed):** ✅ Ran **GO/KEGG** and **PhysioSpace** end-to-end
(after working around a no-`make` toolchain via `R CMD INSTALL`). Both converge: **AVP-OX dampens the
spaceflight jasmonate/defense + osmotic/wounding stress response vs WT.** Real outputs committed
(`go_analysis/results/`, `physiospace/results_static/`); manuscript R5 + R6 filled. Also ran the
**morphometric mixed model** (AVP-OX amplifies the primary-root flight response). Remaining analysis:
freeze the full DESeq2 contrast set (§4.1) to extend GO/PhysioSpace to all genotypes × tissues.

**2026-07-05 — module naming + transcriptome↔architecture integration:** ✅ Named every WGCNA module by
GO (turquoise=signalling & isoprenoid metab.; blue=translation & ribosome; brown=defence & ubiquitin;
yellow=metal/phenylpropanoid; green=photosynthesis; grey=unassigned). Built [`integration/`](integration/)
(WGCNA GS/MM framework): turquoise↑/blue↓ track root growth, brown = AVP-OX-attenuated defence; candidate
mechanism = isoprenoid→growth-hormone precursor metabolism. ⚠ group-level (n=6), Flight/growth collinear →
correlational, not causal. Named Fig 4, updated R4/R7/Methods/Discussion.

**Still needs the maintainer:** a **library→plant manifest** to pair RNA-seq to imaged plants (unlocks individual-level integration, n≈24 — see [`integration/PAIRING_WHATS_NEEDED.md`](integration/PAIRING_WHATS_NEEDED.md)); verify FL/GC root-image calibration; OSDR accession; Zenodo DOI;
author ORCIDs / co-authors; the **RNA-seq read-processing pipeline** (trimmer/aligner/quantifier + versions — not in repo); static figure exports for the dead Slides links; and the
scientific call on canonical DEG contrasts (§4.1) + whether `A68_2A/2B` belong in the analysis.

---

## 0b. Public website (GitHub Pages)

A self-contained landing site lives in [`docs/`](docs/) (`index.html` + `assets/` figures, `.nojekyll` to
bypass Jekyll so the GitBook README can't break the build). **To publish:** GitHub → *Settings → Pages →
Source: Deploy from a branch → `main` / `docs`*. Site URL: `https://dr-richard-barker.github.io/TICTOC/`.
Regenerate figures then re-copy PNGs into `docs/assets/` if results change.

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

- [x] **4.1 DEG model executed (2026-07-04).** Ran a single DESeq2 model (`~0+group`, 12 groups) over the
      count matrix → [`deseq2/`](deseq2/): 10 contrasts (6 spaceflight + 4 AVP-OX interaction) with ashr
      shrinkage + DEG-count summary. Root≫shoot; interactions down-biased (AVP-OX attenuates WT programme).
      Supersedes the legacy iDEP tables for full coverage.
- [~] **4.1b (legacy) Pick ONE canonical DESeq2 pipeline from the v2/v3/iDEP variants;
      archive the rest under `_archive/`. Export a documented design matrix and contrast list:
      - Main effects: Flight vs Ground; Genotype (A68, D130 vs WT); Root vs Shoot.
      - Interactions: **Treatment × Genotype** (does AVP-OX change the spaceflight response?) and
        **Treatment × Tissue** (root- vs shoot-specific flight response).
- [x] **4.2 Organ-specific DEG sets.** ✅ Root- and shoot-specific Flight-vs-Ground DEG tables per genotype
      written to `deseq2/contrasts/` (shrunken LFC + padj, tidy CSV per contrast). *(2026-07-04)*
- [x] **4.3 Tissue-specific clustering — EXECUTED.** Root WGCNA ([`wgcna/`](wgcna/)): 6 modules; turquoise
      (flight r=+0.89, root-growth-coupled) / blue mirror; **brown flight-induced but down in AVP-OX** — module-level
      correlate of the dampening. Fed manuscript R4/Fig4. *(2026-07-05)* | _legacy note:_ (WGCNA was already started for
      leaves; extend to root and to a joint set). Annotate modules by eigengene–trait correlation
      (Flight, genotype, and **root morphometric traits**). ✅ RSML traits now extracted to
      [`morphometrics/rsml_traits.csv`](morphometrics/rsml_traits.csv) (per plant-day) ready for the
      module–trait correlation; ⚠ verify FL/GC length calibration first (see `morphometrics/README.md`). *(2026-07-03)*
- [x] **4.4 GO / pathway analysis — EXECUTED.** Map `Gohir` → Arabidopsis and run GO/KEGG enrichment on
      each module/DEG set. Standardise on one tool (clusterProfiler recommended for reproducibility over
      the iDEP GAGE outputs). ✅ **Blocker resolved:** built [`crosswalk/gohir_to_arabidopsis.tsv`](crosswalk/gohir_to_arabidopsis.tsv)
      from the CottonGen *G. hirsutum* v2.1 → Arabidopsis BLASTP best-hit (92.3% count-matrix coverage;
      reproducible via `crosswalk/build_gohir_to_arabidopsis.py`). Use this instead of the *raimondii*
      `Ara_vs_Cotton_biomart_export`. ✅ **EXECUTED (2026-07-04, R 4.6):**
      [`go_analysis/run_go_clusterprofiler.R`](go_analysis/run_go_clusterprofiler.R) → `go_analysis/results/`.
      Flight×A68 interaction: spaceflight **suppresses jasmonate/wounding/defense** in AVP-OX (response to
      wounding p.adj=3.6e-11; α-linolenic KEGG) + up-regulates *negative* defense regulators —
      **converges with PhysioSpace**. Fed into manuscript R5. Full contrast set awaits §4.1. *(2026-07-04)*
- [x] **4.5 PhysioSpace stress-pattern decoding — EXECUTED.** Project cotton contrasts onto the plant PhysioSpace
      stress compendium (Hadizadeh Esfahani et al.) to produce PhysioScores per genotype/tissue.
      ✅ **Script scaffolded:** [`physiospace/run_physiospace.R`](physiospace/run_physiospace.R) builds
      VST Flight−Ground contrasts, re-indexes Gohir→AT→Entrez via the crosswalk, and runs
      `calculatePhysioMap` (GenesRatio=0.05, TTEST=FALSE, ImputationMethod="PCA") — same params as OSD-767.
      Reference spaces load directly from the **`PlantPhysioSpace`** data package. ✅ **EXECUTED
      (2026-07-04, R 4.6):** real PhysioScores in `physiospace/results_static/`. Roots of all genotypes
      activate a coordinated osmotic/drought/genotoxic/wounding stress signature; AVP-OX roots show
      *attenuated* activation vs WT. Use `--static TRUE` (signed-p saturates roots to Inf). Fed into
      manuscript R6. Remaining: formal genotype contrast + frozen-model input. *(2026-07-04)*
- [x] **4.6 Root morphometrics ↔ transcriptome integration — EXECUTED.** Mixed model on RSML traits
      ([`morphometrics/`](morphometrics/)) + root WGCNA module–trait correlation ([`wgcna/`](wgcna/)):
      flight-responsive modules (turquoise +/blue −) co-vary with root architecture; brown module ties the
      AVP-OX attenuation to expression. Synthesised in manuscript R7. ⚠ trait link is group-level +
      inherits the FL/GC calibration caveat. *(2026-07-05)*

---

## 5. FAIR & repository hygiene

**Findable / Accessible**
- [x] Roadmap link + "repo contents" section added to public `README.md`. *(done)*
- [ ] Reserve a **Zenodo DOI** and add the badge to `README.md` once minted (§6).
- [ ] Record the OSDR/GeneLab accession number(s) and cross-link both ways.

**Interoperable**
- [x] Add a **data dictionary** ([`DATA_DICTIONARY.md`](DATA_DICTIONARY.md)): column-level docs for the
      design/count/DEG/enrichment/WGCNA tables, RSML trait schema, imagery naming, ortholog map, and the
      sample-ID grammar (`{Genotype}_{plate}{well}_{FL|GC}_{day}`). *(2026-07-02)*
- [x] Genome build/annotation + ortholog source stated in `DATA_DICTIONARY.md` §10 (UTX-TM1 v2.1; CottonGen BLASTP). *(done)*
- [x] RSML kept as-is; tracer documented (SmartRoot; `DATA_DICTIONARY.md` §7). *(done)*

**Reusable**
- [x] **Environment spec** [`ENVIRONMENT.txt`](ENVIRONMENT.txt) — R 4.6.0 + pinned versions (DESeq2 1.52, clusterProfiler 4.20, PhysioSpaceMethods 0.99.77, PlantPhysioSpace 0.9.14, ashr). *(2026-07-05)*
- [x] De-hard-code paths: removed the absolute `setwd('/Users/drbhomeoffice/...')` from
      `TICTOC_markdown_evolved_v2.Rmd` (R Markdown knits from the file's own dir) and documented the
      external inputs it needs (`Downloaded_Converted_Data.csv`, `iDEP_core_functions.R`,
      `Cotton__ghirsutum_eg_gene_GeneInfo.csv`, `Cotton__ghirsutum_eg_gene.db`). Still to do: commit or
      link those missing input files. *(2026-07-02)*
- [x] Licensing split set up: data CC0 (`LICENSE`) + code MIT (`LICENSE-CODE`); maintainer to confirm. *(2026-07-08)*

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
- [~] `Figures/readme.md` now points at the committed reproducible manuscript figures (Slides link kept, flagged legacy). README/data-1 Slides links remain (need your exports).

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
  - [x] F1. Design schematic ✅ `manuscript/figures/Fig1_design.pdf`
  - [x] F2. Root architecture ✅ `Fig2_root_architecture.pdf`
  - [x] F3. PCA + DEG counts ✅ `Fig3_PCA.pdf` + `Fig3b_DEG_counts.pdf`
  - [x] F4. Named WGCNA module–trait heatmap ✅ `Fig4_module_trait_named.pdf`
  - [x] F5. GO dotplots ✅ `go_analysis/results_full/GO_*.pdf`
  - [x] F6. PhysioScore heatmap ✅ `physiospace/results_static/`
  - [x] F7. Integrative model ✅ `Fig7_integrative_model.pdf`
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
