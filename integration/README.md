# Integrating the root transcriptome with root architecture

Connects the RNA-seq / WGCNA modules to the measured root-growth changes, to nominate the genes,
pathways and metabolic processes **consistent with** driving the spaceflight root phenotype.

## ⚠ Statistical framing (read first)
RNA-seq (4 replicate libraries per Genotype × Tissue × Treatment group) is **not individually paired**
to the imaged plants (RSML is per-plant). The only valid join is at the **group level: 6 root
genotype × treatment groups**. Two consequences:
1. **n = 6.** Expression–trait correlations are computed on group means (not 24 pseudo-replicated
   samples). With 6 points, |r| ≈ 0.9 is easy to reach — results are **hypothesis-generating, not
   confirmatory**, and per-gene p-values are not trusted.
2. **Flight and traced root-growth are collinear across the 6 groups** (flight groups have the large
   roots). So a module can be flagged as *consistent with* root growth, but its growth-association
   **cannot be statistically separated from its Flight-response**. (The traced growth trait also carries
   the FL/GC image-calibration caveat — see `../morphometrics/README.md`.)

A causal integrated model (multi-block sPLS/DIABLO, or a mediation model Treatment → expression → trait)
would require **individually-paired expression + imaging**, or more groups/timepoints.

> **Can we pair them?** The RNA-seq libraries currently carry no plant ID, so pairing needs one extra
> file (a library→plant manifest). See [`PAIRING_WHATS_NEEDED.md`](PAIRING_WHATS_NEEDED.md) for the exact
> ask to send collaborators, and [`pair_rnaseq_to_images.py`](pair_rnaseq_to_images.py) — a ready-to-run
> join that produces an individual-level (n≈24) table the moment a manifest arrives.
> **Easier path:** pair by LOCATION (row block) — `location_traits.py` averages imaged plants per row as technical reps (23 cells), and `pair_by_location.py` joins libraries to rows once the replicate↔row order is confirmed (see PAIRING_WHATS_NEEDED.md).

## Model (WGCNA gene-significance / module-membership framework)
`integrate_root_expression.R`:
1. Module eigengene ↔ trait correlation at group level (n = 6).
2. **Module Membership (MM)** = cor(gene, module eigengene) over 24 samples (within-modality → hubness).
3. **Gene Significance (GS)** = group-level cor(gene, total root length) → growth-correlated genes.
4. **Hub-and-growth genes** = |MM| > 0.7 in a growth-associated module **and** |GS| > 0.7.

## Result — named modules vs Flight / AVP-OX / root architecture (group-level)

| Module — name | Flight | AVP-OX | Total len | Lateral n | hub-growth genes |
|---|--:|--:|--:|--:|--:|
| **turquoise — Signalling & isoprenoid metabolism** | +0.99 | +0.06 | +0.95 | +0.92 | 1678 |
| **blue — Translation & ribosome biogenesis** | −0.96 | −0.09 | −0.96 | −0.95 | 1444 |
| yellow — Metal transport & phenylpropanoid metab. | +0.70 | −0.10 | +0.75 | +0.76 | 21 |
| **brown — Defence & ubiquitin signalling** | +0.75 | **−0.49** | +0.46 | +0.38 | 0 |
| green — Photosynthesis (light reactions) | −0.70 | +0.41 | −0.50 | −0.43 | 0 |
| grey — Unassigned | −0.16 | +0.48 | −0.11 | −0.10 | 0 |

## Candidate mechanism (hypothesis)
The spaceflight root-growth phenotype co-occurs with **induction of a signalling / isoprenoid-metabolism
module (turquoise)** and **suppression of translation/ribosome biogenesis (blue)**. Isoprenoid/terpenoid
metabolism supplies precursors for growth-regulating hormones (gibberellins, brassinosteroids,
strigolactones, cytokinin/ABA moieties) — a **testable candidate route** from spaceflight signalling to
altered root growth. The **defence & ubiquitin module (brown)** is the one specifically **attenuated in
AVP-OX** (r = −0.49), tying the engineered-genotype effect to reduced defence signalling rather than to
the growth axis. Top growth-correlated hub genes are listed in `results/hub_growth_genes.csv`.

## Outputs
`results/module_trait_grouplevel.csv`, `gene_GS_MM.csv`, `hub_growth_genes.csv`,
`integration_summary.csv`; figure `../manuscript/Fig4_module_trait_named.pdf`;
module names `../wgcna/results/module_names.csv`.
