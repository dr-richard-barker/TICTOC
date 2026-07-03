# TICTOC — Project Roadmap, FAIR & Publication Tracker

**Targeting Improved Cotton Through Orbital Cultivation (TIC-TOC)**
Gilroy Lab, University of Wisconsin–Madison · CASIS grant UA-2018-276 · flown on SpaceX CRS-22 to the ISS.

> This is the **living working document** for finishing the analysis, making the repo FAIR, and
> preparing the Zenodo deposit + npj Microgravity manuscript. It sits alongside the public GitBook
> landing page (`README.md`), which is left untouched. Update the checkboxes as items land.
>
> Last reviewed: 2026-07-02 · Maintainer: R. Barker

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
├── LICENSE ......................... CC0-1.0
│
├── RNA-seq — counts & design
│   ├── TICTOC_run1_filteredCounts_v3.csv .... filtered count matrix (~59.9k genes × 48)
│   ├── TICTOC_cotton_run2_csv/
│   │   └── TICTOC_run1_allCounts.csv ........ UNfiltered counts  ⚠ folder says "run2", file says "run1"
│   ├── TICTOC_target_v5.csv ................. sample/design table (Treatment×Genotype×Tissue)
│   └── Ara_vs_Cotton_biomart_export.txt.zip . cotton→Arabidopsis ortholog map (for GO/PhysioSpace)
│
├── RNA-seq — analysis (the "many models" to consolidate)
│   └── TICTOC_3_factor_model/
│       ├── TICTOC_markdown_evolved_v2.Rmd ... iDEP-derived DESeq2 pipeline (⚠ hard-coded local paths)
│       ├── AllGeneListsGMT.gmt .............. gene sets for enrichment
│       ├── *_GAGE/Enriched/Diff_genes_*.csv . DEG + enrichment outputs
│       └── markdown reports/*.html .......... rendered v2/v3/iDEP reports (large)
│
├── Root morphometrics (RSML)
│   └── Data/
│       ├── Final_RSML_format/ .............. ★ curated RSML set ("outlier removed") — USE THIS
│       └── Cotton {FL,GC} day{3..6} .../ .... raw per-day TIFFs + multiple tracings (.rsml, .rsml01–04)
│                                              ⚠ contains scratch dirs: "test", "Late", "delete me_templates"
│
├── Leaves/ ......................... leaf/shoot DESeq2 (interaction + linear), WGCNA modules, OSD cross-refs
│   └── TITCO_leaves_2024/ .......... ⚠ duplicates several Leaves/*.zip byte-for-byte
│
├── Figures/ ........................ root-phenotype figures (currently just a readme + external links)
└── *.pdf, drbs-*.md, data-1.md ..... ASGSR/CASIS decks & narrative stubs
```

A cleaned, canonical map is a FAIR deliverable — see §5 hygiene tasks.

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
      (Flight, genotype, and **root morphometric traits** — length, surface area, volume from RSML).
- [ ] **4.4 GO / pathway analysis per cluster.** Map `Gohir` → Arabidopsis and run GO/KEGG enrichment on
      each module/DEG set. Standardise on one tool (clusterProfiler recommended for reproducibility over
      the iDEP GAGE outputs). ✅ **Blocker resolved:** built [`crosswalk/gohir_to_arabidopsis.tsv`](crosswalk/gohir_to_arabidopsis.tsv)
      from the CottonGen *G. hirsutum* v2.1 → Arabidopsis BLASTP best-hit (92.3% count-matrix coverage;
      reproducible via `crosswalk/build_gohir_to_arabidopsis.py`). Use this instead of the *raimondii*
      `Ara_vs_Cotton_biomart_export`. *(2026-07-03)*
- [ ] **4.5 PhysioSpace stress-pattern decoding.** Project cotton contrasts onto the plant PhysioSpace
      stress compendium (Hadizadeh Esfahani et al.) to produce PhysioScores per genotype/tissue.
      Reuse the OSD-767 `physioworker` scripts as the template; verify the ortholog bridge for cotton.
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
- [ ] Consolidate DEG models → keep one, archive the rest (§4.1).
- [~] Resolve the mislabeled `TICTOC_cotton_run2_csv/` (holds `run1_allCounts`) — added an explanatory
      `README.md` in the folder; **rename decision pending maintainer** (was there ever a real run 2?).
- [x] De-duplicate `Leaves/` vs `Leaves/TITCO_leaves_2024/` — removed the two byte-identical DEG-value
      zips from the `TITCO_leaves_2024/` data-dump subfolder (~22 MB); canonical copies kept in
      top-level `Leaves/`. *(2026-07-02)*
- [~] Prune RSML scratch dirs under `Data/` (`test`, `Late`, `Convert late`, `delete me_templates`) and
      settle on `Final_RSML_format/` as canonical. ✅ QC documented in
      [`Data/RSML_QC_summary.md`](Data/RSML_QC_summary.md) — 198 curated tracings, 53 plants, 44 with the
      full day-3→6 series; design is unbalanced (A68 smallest). **Scratch-dir deletion still pending
      maintainer sign-off.** *(2026-07-03)*
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
