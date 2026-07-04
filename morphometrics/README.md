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

## Statistical model (executed — `morphometric_stats.py`)
Linear mixed model on `log(total_length_native)`, random intercept per plant (repeated days 3–6),
Ground/WT as references. Outputs: `morphometric_stats_summary.txt`, `morphometric_stats_coefficients.csv`.
Key fixed effects:

| Effect | exp(coef) | p | Reading |
|---|---|---|---|
| `day × condition[FL]` | ×1.65 / day | ~1.6e-31 | **Flight roots elongate faster** over days 3–6 |
| `condition[FL] × genotype[A68]` | ×2.34 | 0.032 | A68 flight response **larger than WT** |
| `condition[FL] × genotype[D130]` | ×2.04 | 0.037 | D130 flight response **larger than WT** |

> **Robustness vs the calibration caveat.** A *constant* Flight-image scale error loads onto the
> `condition` main effect (a fixed log offset), **not** the `day × condition` slope or the
> `condition × genotype` interaction. Stats quantify the traced-data pattern only.

## Sensitivity check — primary root only (`morphometric_stats.py --trait primary_length_native`)
The **primary root** is thick and equally visible in flight and ground images, so it is **immune to the
lateral-tracing-depth artifact**. Re-running the model on `primary_length_native` (outputs
`primary_stats_*`) separates a robust result from a fragile one:

| Effect | total length | **primary length** | Reading |
|---|---|---|---|
| `day × condition[FL]` | ×1.65/day (P≈1e-31) | **×0.95/day (P=0.009)** | The huge flight length gain is **lateral-driven**; primary elongation is *not* accelerated (slightly slower) |
| `condition[FL] × A68` | ×2.34 (P=0.032) | **×1.84 (P=0.040)** | AVP-OX amplification of the flight response **holds on primary root** |
| `condition[FL] × D130` | ×2.04 (P=0.037) | **×1.84 (P=0.017)** | …and is **robust to the tracing artifact** |

**Conclusions:**
1. ✅ **Robust:** the two AVP-OX lines mount a larger primary-root spaceflight response than WT — seen in
   the primary root, which both image sets resolve equally, so *not* a tracing-depth artifact.
2. ⚠ **Verify:** the large total-length flight increase is lateral-root proliferation. This is either a
   real spaceflight phenotype **or** inflated by finer laterals being traceable in enhanced flight
   images. Verify image scale and re-count laterals under a matched detection threshold before claiming it.
