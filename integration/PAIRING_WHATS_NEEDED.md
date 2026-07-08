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

## The two questions we need answered

1. **Were the sequenced plants the same individuals that were imaged?** (i.e. each plant imaged days 3–6,
   then harvested for RNA.)
2. **Was each RNA-seq library one plant, or a pool of several plants?**

- **Same individuals, one plant per library →** pairing is *recoverable* (see below). 🎯
- **Different plants, or pooled libraries →** per-plant pairing is impossible in principle; please just
  send the pooling scheme so we can document it. Only group-level integration remains valid.

## What to send (if recoverable)

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
