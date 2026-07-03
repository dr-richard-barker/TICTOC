# Root morphometrics (RSML) — QC inventory

Summary of the curated root-tracing set in [`Final_RSML_format/`](Final_RSML_format/) — the
outlier-removed canonical RSML files used for morphometric analysis. Raw multi-tracing versions
(`.rsml`, `.rsml01`–`.rsml04`) under the per-day `Cotton {FL|GC} day{3..6} …/` folders are working
files and are superseded by this set.

Generated from filenames following the grammar `{Genotype}_{plate}{well}_{FL|GC}_{Day}.rsml`
(see `../DATA_DICTIONARY.md` §1). All **198** files parsed cleanly (0 malformed).

## Totals
- **198** curated RSML tracings
- **53** unique plants (genotype × plate-well × condition)
- **44 / 53** plants have the **complete day 3→6 time series**
- Time points: day 3 (49), day 4 (50), day 5 (49), day 6 (50) — evenly sampled

## Tracings per genotype × condition × day

| Genotype | Condition | d3 | d4 | d5 | d6 | Total |
|---|---|---:|---:|---:|---:|---:|
| WT   | Ground (GC) | 11 | 9 | 9 | 11 | 40 |
| WT   | Flight (FL) | 13 | 13 | 12 | 11 | 49 |
| A68  | Ground (GC) | 4 | 5 | 5 | 5 | 19 |
| A68  | Flight (FL) | 5 | 5 | 5 | 5 | 20 |
| D130 | Ground (GC) | 10 | 12 | 12 | 12 | 46 |
| D130 | Flight (FL) | 6 | 6 | 6 | 6 | 24 |

## Unique plants per genotype × condition

| Genotype | Ground (GC) | Flight (FL) | Total |
|---|---:|---:|---:|
| WT   | 12 | 13 | 25 |
| A68  | 5 | 5 | 10 |
| D130 | 12 | 6 | 18 |
| **All** | **29** | **24** | **53** |

## Day-series completeness (per plant)

| Days tracked | # plants |
|---:|---:|
| 4 (full 3–6) | 44 |
| 3 | 6 |
| 2 | 1 |
| 1 | 2 |

## Scratch-folder prune (2026-07-03)
Removed **32 truly-redundant files (~35 MB)** from the working sub-folders (`test/`, `Convert late/`,
`Converted late/`, `delete me_templates/`) — verified safe by git blob-SHA match (identical content
kept elsewhere) or the explicit `delete me_templates/` set.

**Deliberately kept — unique content that was mislabeled as scratch (10 files, ~107 MB):**
- `Cotton FL day 3 …/a68/A68_Cotton_SVT_timelapse.avi` — a root-growth **timelapse video** (no copy elsewhere).
- `Cotton FL day5 …/Late/A68_2A_FL_5…` and `A68_2B_FL_5…` — day-5 images of A68 plants labelled `2A`/`2B`
  that are **not in the curated set** (note: `A68_2B`'s TIFF is byte-identical to `A68_b2_FL_5`, so `2B`
  is likely a relabel of `b2`; `A68_2A` appears genuinely extra).
- Unique `Convert/Convert late/Converted late` source images for `D130_c11_GC_4`, `D130_a11_GC_6`,
  `WT_c10_GC_5`.

> **Maintainer decision:** confirm whether `A68_2A`/`A68_2B` are extra plants (and whether they belong in
> the analysis), and whether the timelapse + these images are also on OSDR — if so they can be dropped
> from git; if not, move them to a clearly-named location (they are primary data, not scratch).

## Notes for analysis
- The design is **unbalanced** across genotype × condition (A68 has the fewest plants: 5 FL / 5 GC;
  D130 Flight has 6 vs 12 Ground). Use models that tolerate unbalanced groups (mixed models with
  plant as a random effect) rather than assuming equal n.
- 9 plants have an incomplete series (1–3 days). Decide up front whether trait-vs-day analyses use
  only the 44 complete-series plants or all plants with a day covariate — and state it in Methods.
- This summary is regenerable from `Final_RSML_format/` filenames; update it if the curated set changes.
