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

## Notes for analysis
- The design is **unbalanced** across genotype × condition (A68 has the fewest plants: 5 FL / 5 GC;
  D130 Flight has 6 vs 12 Ground). Use models that tolerate unbalanced groups (mixed models with
  plant as a random effect) rather than assuming equal n.
- 9 plants have an incomplete series (1–3 days). Decide up front whether trait-vs-day analyses use
  only the 44 complete-series plants or all plants with a day covariate — and state it in Methods.
- This summary is regenerable from `Final_RSML_format/` filenames; update it if the curated set changes.
