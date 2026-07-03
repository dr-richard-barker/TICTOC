# TICTOC — Data Dictionary

Column-level and file-level documentation for the TIC-TOC cotton spaceflight datasets, to support
reuse (the **I** and **R** of FAIR). See [`PROJECT_ROADMAP.md`](PROJECT_ROADMAP.md) for status and the
analysis plan, and `README.md` for the study overview.

- **Organism / genome:** *Gossypium hirsutum* (upland cotton), allotetraploid (A + D subgenomes).
- **Expression gene IDs:** `Gohir.<subgenome><chr>G<number>` (e.g. `Gohir.A01G000301`; `1Z…`/`D…`
  denote subgenome/scaffold). ~59.9k genes after filtering.
- **Raw data of record:** NASA OSDR/GeneLab, accession **OSD-XXX** (fill in; embargoed until publication).

---

## 1. Sample-ID grammar

Root/plant samples are named:

```
{Genotype}_{plate}{well}_{Condition}_{Day}
      WT  _  a3    _   FL       _  3
```

| Token | Values | Meaning |
|---|---|---|
| Genotype | `WT`, `A68`, `D130` | WT = wild type; **A68** and **D130** = AVP1-overexpressing (AVP-OX) transgenic lines |
| plate+well | e.g. `a3`, `b12`, `c11` | plate letter + well number (plant identifier) |
| Condition | `FL`, `GC` | **FL** = Flight (ISS); **GC** = Ground Control |
| Day | `3`, `4`, `5`, `6` | imaging day (root time series) |

RNA-seq sample/column labels instead use the compact form
`{Genotype}{Tissue}{Condition}` (e.g. `A68RootFlight`, `WTShootGround`), with 4 replicates per group.
Flight image labels embed the ISS downlink frame ID, e.g. `iss065e092623` (ISS Expedition 65).

---

## 2. RNA-seq — design / sample tables

**`TICTOC_target_v5.csv`** — full design (48 samples, root + shoot). Transposed layout: first row =
sample names (repeated per replicate), then one row per factor.

| Row label | Values | Notes |
|---|---|---|
| `Name` | `A68RootFlight`, … | group label, repeated ×4 replicates |
| `Treatment` | `Flight`, `Ground` | spaceflight condition |
| `Genotype` | `WT`, `A68`, `D130` | genotype |
| `Tissue` | `Root`, `Shoot` | organ |

**`Leaves/TICTOC_target_v5_shoot.csv`** — same schema, shoot/leaf subset only (24 samples).

---

## 3. RNA-seq — count matrices

Rows = genes, columns = samples (4 replicates per group; replicate columns share a header name or use
`…Flight`, `…Flight1`, `…Flight2`, `…Flight3`).

| File | First column | Content |
|---|---|---|
| `TICTOC_run1_filteredCounts_v3.csv` | `Gene ID` | **filtered** raw counts (~59.9k genes × 48) — analysis input |
| `TICTOC_cotton_run2_csv/TICTOC_run1_allCounts.csv` | `ID` | **unfiltered** raw counts (all genes) ⚠ folder mislabeled "run2"; file is run1 |
| `Leaves/TICTOC_run1_allCountsMulti_titles_v4_shoots.csv` | `ID` | shoot-only counts, multi-title header |

Values are integer raw read counts (RSEM/featureCounts-style), not normalized.

---

## 4. RNA-seq — differential expression & heatmap tables

**`TICTOC_3_factor_model/Diff_genes_heatmap_*.csv`** — per-contrast DEG tables with expression.
Contrast is encoded in the filename, e.g. `A68-WT (Root_GC)` = A68 vs WT in Root, Ground;
`I_Treatment_Flight.Genotype_A68` = interaction term (Treatment[Flight] × Genotype[A68]).

| Column | Meaning |
|---|---|
| `Regulation` | `Up` / `Down` call |
| `Ensembl ID` | gene ID (`Gohir…`) |
| `log2 Fold Change` | contrast LFC |
| `Adj.Pval` | BH-adjusted p-value |
| `Symbol` | gene symbol (if mapped) |
| `Chr`, `Type` | chromosome / gene biotype |
| *sample columns* | per-sample expression (e.g. `A68RootFlight`, `A68RootFlight1`, …) |

---

## 5. RNA-seq — pathway / enrichment outputs

**`TICTOC_3_factor_model/Enriched.csv`** (ORA-style):

| Column | Meaning |
|---|---|
| `Direction` | Up / Down set tested |
| `adj.Pval` | adjusted enrichment p-value |
| `nGenes` | genes in the term ∩ set |
| `Pathways` | pathway / GO term name |
| `Genes` | member gene IDs |

**`TICTOC_3_factor_model/*(GAGE).csv`** (e.g. `D130-WT(GAGE).csv`) — GAGE gene-set analysis:
`Direction`, `Pathways`, `statistic`, `nGenes`, `adj.Pval`, `Genes`.

**`TICTOC_3_factor_model/AllGeneListsGMT.gmt`** — gene sets in GMT format: each line is
`<set name>\t<count/description>\t<gene>\t<gene>…` (note: gene IDs here are upper-cased `GOHIR.…`).

---

## 6. Leaf/shoot co-expression (WGCNA)

**`Leaves/Titoc_cotton_leaves_WGCNA_modules.csv`** — one row per gene:

| Column | Meaning |
|---|---|
| `symbol__gene_id` | symbol + gene id composite key |
| `gene_id` | `Gohir…` gene id |
| `module_color`, `module` | WGCNA module assignment (color label + index) |
| *sample columns* | expression per shoot sample/replicate (`A68ShootFlight`, `A68ShootFlight1`, …) |

---

## 7. Root morphometrics — RSML

**`Data/Final_RSML_format/*.rsml`** — curated, outlier-removed SmartRoot tracings (canonical set).
Raw multi-tracing versions (`.rsml`, `.rsml01`–`.rsml04`) live under the per-day `Data/Cotton …/`
folders and represent repeated/alternate traces of the same image — prefer `Final_RSML_format/`.

RSML is an open XML standard (root system markup language). Schema (SmartRoot v1):

- **`<metadata>`** — `unit` (inch), `resolution` (px per unit; ~37.89 → 300 dpi class), `software`
  (`smartroot`), and `<property-definitions>` declaring trait labels & units.
- **Declared traits:** `diameter` (cm), `length` (cm), `angle` (degree), `insertion` (cm),
  `lauz`, `lbuz` (cm; length above/below unbranched zone), `node-orientation` (radian).
- **`<scene><plant><root>`** — each root has `po:accession` (Plant Ontology term, `PO:0009005` = root),
  `<properties>` (`length`, `orientation`, `rulerAtOrigin`), `<geometry><polyline>` of `x,y` points
  (image pixel coords), and per-node `<functions>` (`diameter`, `orientation` samples along the polyline).
- **`<image><label>`** — links the trace to its source photo (e.g. `WT_a3_FL_3_iss065e092623_clarity_enhanced_300dpi`).

> Derived phenotypes for analysis (total root length, surface area, volume, count) are computed from
> these polylines + diameter functions — document the extraction script when added.

---

## 8. Time-series photography

**`Data/Cotton {FL|GC} day{3..6} …_inverted/*.tif`** — grayscale, inverted, 300 dpi root scans/photos.
Flight frames are clarity-enhanced ISS downlink images (`iss065e…`); ground frames use camera IDs
(e.g. `DSC_0450`). One image per plant per day; filenames follow the §1 sample grammar plus the source
frame ID. These are the inputs the RSML traces were drawn on.

---

## 9. Cotton ↔ Arabidopsis ortholog map

**`Ara_vs_Cotton_biomart_export.txt.zip`** — BioMart export (Ensembl Plants, ~2022-10) bridging
Arabidopsis to cotton, used for GO annotation transfer and PhysioSpace mapping. Tab-delimited columns:

`Gene stable ID` (AT…) · `Transcript stable ID` · `Gene description` · `Gene name` · `Gene type` ·
`Gene Synonym` · `Gene stable ID` (cotton) · `Transcript stable ID` (cotton) · `Gene type` ·
`Ensembl Canonical` · `Gene description` · `Transcript type` · `STRING ID`.

> ⚠ **ID-space caveat (interoperability):** the cotton IDs in this export are `B456_…` / `Gorai…`
> (an older *G. raimondii* / JGI annotation), **not** the `Gohir.…` IDs used in the count matrices and
> DEG tables — so this export does **not** join to the expression data. ✅ **Use the corrected crosswalk
> instead:** [`crosswalk/gohir_to_arabidopsis.tsv`](crosswalk/gohir_to_arabidopsis.tsv) maps `Gohir.` →
> Arabidopsis directly (CottonGen *G. hirsutum* v2.1 BLASTP best-hit, 92.3% coverage; see `crosswalk/README.md`).
> This older raimondii export is retained only for provenance.

---

## 10. Versions to pin (fill in during cleanup)

- [ ] *G. hirsutum* genome assembly + annotation release (for `Gohir` IDs).
- [ ] Ensembl Plants / BioMart release used for the ortholog export.
- [ ] SmartRoot version; DESeq2 / iDEP versions; R version (`sessionInfo`).
- [ ] OSDR/GeneLab accession number(s).
