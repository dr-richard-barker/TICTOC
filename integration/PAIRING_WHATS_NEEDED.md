# Pairing RNA-seq libraries to imaged plants — what's needed

*A request for collaborators / the Gilroy lab. Answering this would substantially strengthen the
transcriptome ↔ root-architecture integration in the manuscript.*

## The situation
Two datasets exist per group but are only linked at the **group** level:

| Dataset | Identifier it carries | Granularity |
|---|---|---|
| Root imaging → RSML tracings | `genotype_plateWell_condition` (e.g. `A68_a12_FL`) | **per plant** ✅ |
| RNA-seq (count matrix + design) | `Genotype × Tissue × Treatment` only, ×4 | **group only** ❌ |

The four replicate RNA-seq libraries per group are **anonymous** — the count-matrix columns were
renamed to group labels, so we cannot tell which plant each library came from. There are **53 imaged
root plants but only 24 root RNA-seq libraries** (imaging covered more plants than were sequenced).

Because of this, the current expression↔growth integration is **group-level (n = 6)**, where spaceflight
and root growth are collinear — so we can say a pathway is *consistent with* driving growth but can't
separate that from the flight response.

## Update (2026-07): RNA-seq sample metadata received — what it does and doesn't give

The TIC-TOC "from biotech" sample sheet ([`../metadata/TICTOC_rnaseq_sample_metadata.csv`](../metadata/TICTOC_rnaseq_sample_metadata.csv))
**confirms the sample structure and gives every library a definitive ID** (e.g. `A68R-F1_S5_L004` =
A68, Root, Flight, replicate 1), and its order **exactly matches the count-matrix columns** — so each of
the 48 count columns is now labelled ([`../deseq2/sample_sheet.csv`](../deseq2/sample_sheet.csv)).
**However it contains no well/position column** — the replicate number (1–4) is an ordinal, not a
documented location. So it does not, by itself, link a library to an imaged plant.

**The gap is now a single, tight lookup:** for each of the 48 libraries (or at least the 24 root ones),
**which well/plant (`a12`, `b2`, …) — or which row (a/b/c/d) — was it extracted from?** That is the RNA
extraction / planting-map log. With it, `pair_by_location.py` (row level) or `pair_rnaseq_to_images.py`
(well level) joins immediately. If replicates were simply taken in row order, telling us "rep 1=row a … rep
4=row d" is enough.

## The two questions we need answered

1. **Were the sequenced plants the same individuals that were imaged?** (i.e. each plant imaged days 3–6,
   then harvested for RNA.)
2. **Was each RNA-seq library one plant, or a pool of several plants?**

- **Same individuals, one plant per library →** pairing is *recoverable* (see below). 🎯
- **Different plants, or pooled libraries →** per-plant pairing is impossible in principle; please just
  send the pooling scheme so we can document it. Only group-level integration remains valid.

## ⭐ Easier alternative — pair by LOCATION (row block), not exact plant

We don't strictly need the exact plant. The experiment has spatial **locations** (well = row-letter a–d +
number). At **row** level, most genotype × condition cells hold **2–4 imaged plants**, which we average as
**technical replicates** (`location_traits.py` → `results/location_mean_traits.csv`; 23 location cells,
14 with ≥2 plants — ~3.8× finer than the 6 group means). Because RNA-seq has **4 replicates** per
genotype × tissue × treatment and imaging has **4 rows (a–d)**, they plausibly correspond (a common block
design).

**So the minimal ask is just:** *do the 4 RNA-seq replicates correspond to the 4 rows a, b, c, d — and in
what order?* That single answer (or a `library_id → row` CSV) lets us join RNA-seq to row-mean traits at
~23-cell resolution via `pair_by_location.py`. (We can preview it now under the assumption "reps = rows
a,b,c,d in column order" with `python integration/pair_by_location.py --rows-in-order`, but that assumption
must be confirmed before the result is trusted.) This is much easier to recover than exact-plant identity.

## What to send for exact-plant (1:1) pairing (if available)

A single **library → plant manifest**: one row per RNA-seq library mapping it to the plate+well of the
plant it came from. Any of these sources likely contains it:

- **The OSDR/GeneLab ISA-Tab submission** — the sample table (`s_*.txt`) and assay table (`a_*.txt`)
  usually record the *source plant / sample name*. Since the OSDR release is being prepared now, this is
  the most natural place to restore the link. ⭐
- **Original FASTQ filenames** or the **sequencing-core submission form** (sample names almost always
  encoded the plant/plate/well/tube).
- A **harvest or RNA-extraction log** linking each tube to a plant.

### Manifest format (CSV)
```
library_id,genotype,tissue,treatment,plate_well
lib01,A68,Root,Flight,a12        # plate_well matches the RSML naming: a12, b2, c11, d8, ...
lib02,A68,Root,Flight,b2
...                               # ideally all 48 libraries (24 root + 24 shoot)
```
`plate_well` must use the same lower-case plate+well token as the RSML filenames (see the `plant` column
in `../morphometrics/rsml_traits.csv`).

## What happens next (already built)
Drop that manifest in and run:
```bash
python integration/pair_rnaseq_to_images.py --manifest library_to_plant.csv
```
It joins each library to its plant's RSML traits and writes `results/paired_rnaseq_rsml.csv`. That table
enables an **individual-level** integrated model (n up to 24 root plants) that **breaks the Flight/growth
collinearity** — letting us test, within a treatment group, whether plants with higher expression of the
signalling/isoprenoid module actually have larger roots, and fit a causal mediation model
(Treatment → expression → root trait).

## The payoff
- Group-level (now): n = 6, correlational, Flight/growth confounded.
- Paired (with manifest): n ≈ 24, within-group correlations, causal/mediation modelling, multi-block
  methods (sPLS/DIABLO) — turning the isoprenoid → growth-hormone hypothesis into a testable claim.
