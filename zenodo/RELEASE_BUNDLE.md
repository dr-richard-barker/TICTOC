# Zenodo release bundle — what to deposit

Checklist for the `v1.0.0` archive that mints the citable DOI. Tagging a GitHub release
(`v1.0.0`) auto-archives the whole repo to Zenodo once the GitHub↔Zenodo hook is enabled; this file
records what the bundle should contain and the metadata to set.

## Include (already in repo — archived automatically on release)
- **Analysis code:** `crosswalk/ deseq2/ go_analysis/ physiospace/ wgcna/ morphometrics/ integration/`
- **Derived results:** each module's `results/` (DEG tables, GO/KEGG, PhysioScores, module data, traits)
- **Docs:** `README.md PROJECT_ROADMAP.md DATA_DICTIONARY.md ENVIRONMENT.txt CITATION.cff`
- **Manuscript + figures:** `manuscript/` (draft, figures F1–F7, tables T1–T4, supplementary index)
- **Curated inputs:** `TICTOC_run1_filteredCounts_v3.csv`, `TICTOC_target_v5.csv`,
  `crosswalk/gohir_to_arabidopsis.tsv`, `Data/Final_RSML_format/` (curated RSML)

## Keep OUT of the code archive (large / lives elsewhere)
- Raw FASTQ + raw images → **NASA OSDR/GeneLab** (cite the OSD accession, do not duplicate on Zenodo)
- Bulk per-day TIFFs under `Data/Cotton …/`, rendered `*.html` reports, slide-deck PDFs (regenerable /
  heavy; see `.gitignore`)

## Zenodo metadata to set (from `CITATION.cff`, then complete)
- **Title:** Targeting Improved Cotton Through Orbital Cultivation (TICTOC) — data & analysis
- **Authors + ORCIDs:** *[to complete]*
- **License:** data CC0-1.0; code MIT (`LICENSE`, `LICENSE-CODE`)
- **Funding:** CASIS UA-2018-276
- **Related identifiers:** OSDR accession (isDerivedFrom/isSupplementTo), manuscript DOI (once known),
  GitHub repo URL
- **Version:** v1.0.0 · **Keywords:** spaceflight, microgravity, cotton, AVP1, RSML, RNA-seq, PhysioSpace

## Steps
1. Reserve a DOI on Zenodo (before publishing) → add the badge to `README.md`.
2. Complete author/ORCID + license fields in `CITATION.cff`.
3. Enable the GitHub↔Zenodo webhook for `dr-richard-barker/TICTOC`.
4. Tag & publish `v1.0.0` → Zenodo archives → DOI minted.
5. Add the DOI to the manuscript Data-availability statement.
