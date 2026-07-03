<!--
  TICTOC — npj Microgravity manuscript SKELETON (draft, 2026-07-03).
  STATUS: structure + Methods are drafted; Results are STUBS with [PLACEHOLDER] markers.
  Do NOT treat any number here as real — fill Results only from executed analysis outputs.
  Figure/table numbering matches figures_and_tables.md and roadmap §7.
-->

# Engineering for stress does not exempt it: root-architecture and transcriptomic response of AVP-OX cotton to spaceflight

**Authors:** [AUTHOR LIST — incl. R. Barker; ORCIDs]
**Affiliations:** [Gilroy Life Science Lab, University of Wisconsin–Madison; …]
**Corresponding author:** [name, email]
**Funding:** CASIS UA-2018-276.

---

## Abstract
*(structured, ~150–200 words — write last, from final results)*
Background: roots drive cotton yield and stress resilience; AVP1-overexpressing (AVP-OX) lines show
enhanced salt/drought tolerance and larger root systems, but how this engineering behaves in
microgravity is unknown. We flew wild-type and two AVP-OX cotton lines (A68, D130) on the ISS with
matched ground controls, and profiled **root architecture** (time-series imaging → RSML tracing) and the
**root + shoot transcriptome** (RNA-seq, full Treatment × Genotype × Tissue design). [KEY RESULT 1 —
morphometrics]. [KEY RESULT 2 — DEG/tissue pattern]. [KEY RESULT 3 — GO/PhysioSpace stress programs].
[INTERPRETATION — does AVP-OX buffer or amplify the spaceflight response]. Implications for engineering
crops for space agriculture.

---

## Introduction
- Roots, gravity, and cotton productivity; auxin-directed root system architecture. [refs]
- AVP1/AVP-OX technology: vacuolar H⁺-PPase overexpression → salt/drought tolerance, larger root
  systems, +fiber yield under stress (Park et al. 2005). [refs]
- Spaceflight remodels plant roots and transcriptomes; microgravity disrupts gravitropic auxin
  transport and cytoskeletal/cell-wall processes. [refs]
- Gap: whether stress-tolerance engineering changes the *spaceflight* response has not been tested in a
  crop, at the intersection of architecture and transcriptome.
- This study: WT vs AVP-OX cotton, Flight vs Ground, root vs shoot, days 3–6 → morphometrics +
  transcriptome, decoded into stress programs via ortholog-mapped GO and PhysioSpace.

---

## Results

### R1. Experimental design and dataset overview  → **Fig 1**, **Table 1**
[STUB] Describe the flown material (WT, A68, D130), Flight/Ground, root/shoot, day-3–6 imaging, and the
RNA-seq sample counts. Point to Table 1 (design) and the RSML QC (n = 53 plants, 44 full series).

### R2. Root architecture over days 3–6  → **Fig 2**, **Table S?**
[STUB — from RSML morphometrics] Flight vs Ground trajectories of root traits (total length, diameter,
growth rate…) per genotype. State whether AVP-OX lines diverge from WT under spaceflight. [PLACEHOLDER
stats: mixed model, treatment × genotype × day; unbalanced design — see QC].

### R3. Transcriptome landscape: light of the design  → **Fig 3**, **Table 2**
[STUB] PCA; DEG counts per contrast (Flight vs Ground; genotype effects; root vs shoot); up/down
asymmetry. [PLACEHOLDER numbers from the frozen DESeq2 model, roadmap §4.1].

### R4. Tissue-specific DEG clusters and trait linkage  → **Fig 4**
[STUB] Organ-specific clustering / WGCNA modules; correlate module eigengenes with Flight, genotype, and
**root morphometric traits** (the RSML link). [PLACEHOLDER].

### R5. Functional enrichment of responsive modules  → **Fig 5**, **Table 3**
[STUB — from `go_analysis/`] GO/KEGG enrichment per module/DEG set, ortholog-mapped to Arabidopsis
(CottonGen BLASTP crosswalk; enriched against the expressed-gene universe). [PLACEHOLDER top terms].

### R6. PhysioSpace stress-pattern decoding  → **Fig 6**, **Table 4**
[STUB — from `physiospace/`] PhysioScores per genotype × tissue against the Arabidopsis stress spaces
(`AT_Stress_Space` / `_Meta` / `_RNASeq`). Which known stress programs the spaceflight response
resembles, and whether AVP-OX shifts them. [PLACEHOLDER].

### R7. Integrative model  → **Fig 7**
[STUB] Synthesis linking root-trait modules ↔ expression modules ↔ stress programs; a model for how (or
whether) AVP-OX engineering buffers the microgravity response.

---

## Discussion
- Interpret the architecture ↔ transcriptome ↔ stress-program integration. [PLACEHOLDER]
- AVP-OX in space: buffering vs amplification; auxin/gravitropism and cell-wall themes. [PLACEHOLDER]
- Limitations: best-hit orthology (approximate); unbalanced morphometric design; single mission.
- Implications for engineering crops for spaceflight agriculture.

---

## Methods
*(Concrete — fill bracketed hardware/version specifics.)*

**Plant material & spaceflight.** Cotton (*Gossypium hirsutum*): wild type and two AVP1-overexpressing
lines (A68, D130). Grown aboard the ISS (SpaceX CRS-22, Expedition 65) with matched ground controls in
[hardware]. Root systems imaged on days 3–6.

**Root imaging & morphometrics.** Time-series photographs (grayscale, 300 dpi; flight frames
clarity-enhanced) were traced in **SmartRoot** and exported as **RSML**. Traits (root length, diameter,
angle, insertion; see `../DATA_DICTIONARY.md` §7) were extracted per plant/day. Curated set: 198
tracings, 53 plants (44 with the full day-3–6 series). Because the design is unbalanced across
genotype × condition, traits were modelled with [linear mixed models, plant as random effect;
treatment × genotype × day fixed effects].

**RNA-seq & differential expression.** Root and shoot RNA-seq, WT/A68/D130 × Flight/Ground, 4
replicates. Reads processed [pipeline], quantified against *G. hirsutum* **UTX-TM1 v2.1**
(Phytozome Ghirsutum_527_v2.1; `Gohir` IDs). Filtered counts (~59.9k genes × 48) analysed with
**DESeq2** using a Treatment × Genotype × Tissue model; contrasts and shrinkage per roadmap §4.1.

**Ortholog mapping.** `Gohir` genes mapped to *Arabidopsis thaliana* by CottonGen BLASTP best-hit
(`../crosswalk/`; 92.3% of expressed genes mapped).

**Functional enrichment.** GO (BP/MF/CC) and KEGG (`ath`) enrichment in Arabidopsis space with
**clusterProfiler**, against the universe of AT loci reachable from expressed cotton genes
(`../go_analysis/`).

**Co-expression.** WGCNA modules [params]; module–trait correlation to Flight, genotype, and RSML traits.

**PhysioSpace.** Per-(genotype × tissue) Flight−Ground VST fold changes were re-indexed to Arabidopsis
Entrez and projected onto the Arabidopsis stress reference spaces with **PhysioSpaceMethods**
(`calculatePhysioMap`, GenesRatio = 0.05, TTEST = FALSE, ImputationMethod = "PCA") using the
**PlantPhysioSpace** data package (`../physiospace/`). Refs: Lenz et al. 2013; Hadizadeh Esfahani et al. 2021.

**Statistics.** [Thresholds: padj < 0.05, |LFC| ≥ 1 for DEGs; enrichment BH q < 0.05; morphometric model
details]. Software versions in `sessionInfo` (`../DATA_DICTIONARY.md` §10).

---

## Data & code availability
Raw reads, count matrices, root images, and RSML tracings: NASA OSDR/GeneLab **OSD-XXX** (released on
publication). Analysis code and derived tables: Zenodo **DOI 10.5281/zenodo.XXXXXXX** and
https://github.com/dr-richard-barker/TICTOC.

## Author contributions / Acknowledgements / Competing interests
[TODO]

## References
[TODO — numbered Vancouver, npj style]
