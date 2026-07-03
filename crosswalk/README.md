# Gohir → Arabidopsis ortholog crosswalk

Bridges the TICTOC cotton gene IDs (`Gohir.*`, *Gossypium hirsutum* Phytozome UTX-TM1 v2.1) to
*Arabidopsis thaliana* (`AT*`) loci, so GO annotation transfer and **PhysioSpace** stress-pattern
decoding (both Arabidopsis-based) can run against the cotton DEGs.

## Why this was needed
The repo's older `../Ara_vs_Cotton_biomart_export.txt.zip` maps Arabidopsis to `B456_`/`Gorai` IDs —
that's *G. raimondii* (the D-genome progenitor), because Ensembl Plants never carried *G. hirsutum*.
Those IDs **do not join** to the `Gohir.*` count matrices / DEG tables. This crosswalk fixes that.

## Files
| File | What |
|---|---|
| `gohir_to_arabidopsis.tsv` | the crosswalk — one best-hit Arabidopsis gene per Gohir gene |
| `build_gohir_to_arabidopsis.py` | reproducible builder (downloads source, rebuilds the TSV) |

## Source (authoritative, public, no login)
CottonGen → *G. hirsutum* UTX-TM1 v2.1 → `homology/blastp_G.hirsutum_UTX_v2.1_vs_arabidopsis.xlsx.gz`
(BLASTP of every Gohir protein against TAIR Arabidopsis proteins; 4 sheets, ~97.6k hit rows).

## Method
Transcript hits are collapsed to gene level (drop `.<isoform>` and the `_UTX-TM1…` suffix); the single
**best hit per Gohir gene** is kept by min E-value → max bit-score → max % identity.

## Coverage
- **67,082** Gohir genes mapped → **15,650** unique Arabidopsis genes.
- **55,320 / 59,918 (92.3%)** of the filtered count-matrix genes have an Arabidopsis ortholog.
- Many Gohir → one AT is expected (A/D homoeologs + paralog expansion in the allotetraploid).

## Columns (`gohir_to_arabidopsis.tsv`, tab-separated)
| Column | Meaning |
|---|---|
| `gohir_gene` | Gohir gene ID, canonical `Gohir.` casing (matches the count matrix) |
| `gohir_gene_upper` | UPPER-CASE key — use to join the DEG tables, which use `GOHIR.*` |
| `arabidopsis_gene` | best-hit Arabidopsis locus (e.g. `AT5G17700`) |
| `pid` | % identity of the BLASTP hit |
| `evalue`, `bit_score`, `align_length` | hit quality metrics (filter here for stricter ortholog sets) |
| `gohir_transcript_hit`, `arabidopsis_transcript_hit` | the exact isoforms behind the hit (provenance) |
| `arabidopsis_description` | TAIR defline / coordinates for the AT hit |

## Usage notes
- **Join casing:** count matrix uses `Gohir.*` → join on `gohir_gene`; DEG heatmap tables use
  `GOHIR.*` → join on `gohir_gene_upper` (or uppercase both sides).
- **A BLASTP best hit ≠ a curated ortholog.** For a stricter set, filter e.g. `pid ≥ 40` and
  `evalue ≤ 1e-10`, or restrict to reciprocal-best-hit pairs.
- **Rebuild:** `python build_gohir_to_arabidopsis.py` (needs `pandas` + `openpyxl`; auto-downloads the
  source xlsx). See `../DATA_DICTIONARY.md` §9 and `../PROJECT_ROADMAP.md` §4.4.

## Also available on CottonGen (for the GO step, if you prefer direct cotton GO over AT transfer)
`functional/Gh_TM1_UTX_v2.1_genes2Go.xlsx.gz` — Gohir → GO term assignments (InterPro/KAAS based).
