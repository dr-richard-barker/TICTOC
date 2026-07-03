# Root morphometrics — RSML trait extraction

Turns the curated SmartRoot tracings (`../Data/Final_RSML_format/*.rsml`) into analysis-ready
morphometric tables for **Results R2 / Figure 2** (roadmap §4.2–4.3, morphometric side).

## Run
```bash
cd morphometrics
python extract_rsml_traits.py     # -> rsml_traits.csv        (one row per plant-day, 198 rows)
python summarise_traits.py        # -> rsml_traits_summary.csv (genotype x condition x day means)
```
No third-party dependencies (Python stdlib only).

## Outputs
- **`rsml_traits.csv`** — per plant-day: `n_roots`, `n_primary`, `n_lateral`,
  `total_length_native`, `primary_length_native`, `total_length_px`, `mean_diameter_px`,
  plus parsed `genotype/condition/day` and the file's `resolution`/`unit`.
- **`rsml_traits_summary.csv`** — group means ± SD per genotype × condition × day.

## Column notes (see also `../DATA_DICTIONARY.md` §7)
- `*_native` = SmartRoot's own `<properties><length>` (RSML metadata unit); `*_px` = geometric
  polyline length in image pixels. The two are linearly proportional here (ratio ≈ 14.9 px/native
  across all files), so **relative** comparisons are sound.
- Root **order**: top-level `<root>` = primary (order 1); nested `<root>` = lateral (order ≥ 2).

## Descriptive pattern (NOT yet a validated result)
At day 6, traced **total root length is much higher in Flight than Ground** in all genotypes, and the
Flight/Ground ratio is larger in the AVP-OX lines (WT ≈ 4.0×; A68 ≈ 5.9×; D130 ≈ 6.4×). Ground plants
stay largely single-primary (few laterals) while Flight plants develop many laterals by day 6.

> ⚠ **Read before interpreting.** This magnitude gap must be checked for artifacts before any
> biological claim:
> 1. **Calibration** — Flight (ISS `iss065e…`, clarity-enhanced) and Ground (`DSC_…`) images come from
>    different cameras. All files record the same `resolution` (37.888332 px/inch), so if the *real*
>    spatial scales differ, absolute FL-vs-GC lengths are confounded. Verify a scale bar / known
>    dimension per image set first.
> 2. **Tracing depth** — enhanced flight images may permit tracing finer laterals than ground images,
>    inflating length/lateral counts. Check tracer consistency.
> 3. **Statistics** — the design is unbalanced (`../Data/RSML_QC_summary.md`); use a mixed model
>    (plant random effect; treatment × genotype × day) — the summary here is descriptive only.

Once calibration + stats are settled, `rsml_traits.csv` feeds Fig 2 and the module–trait correlations
in Fig 4 (roadmap §4.3).
