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
*(PRELIMINARY — from executed RSML analysis; total-length magnitude pending image-calibration check.)*
Linear mixed models of root length (plant random effect; days 3–6) revealed a genotype-dependent
spaceflight response. Both **AVP-OX lines mounted a larger primary-root spaceflight response than wild
type** (condition × genotype on primary-root length: A68 ×1.84, *P* = 0.040; D130 ×1.84, *P* = 0.017) —
an effect measured on the primary root, which is equally resolved in flight and ground images and is
therefore robust to image artifacts. Total traced root length increased far more under spaceflight
(day × condition ≈ 1.65×/day), driven by **lateral-root proliferation**: flight plants developed many
laterals by day 6 while ground plants remained largely single-primary (Fig 2). Notably, primary-root
elongation itself was *not* accelerated in flight (day × condition ×0.95, *P* = 0.009), localising the
spaceflight length gain to lateral roots. [The magnitude of the lateral effect is reported pending
verification of Flight/Ground image scale and lateral-detection thresholds; see `../morphometrics/README.md`.]

### R3. Transcriptome landscape of the spaceflight response  → **Fig 3**, **Table 2**
*(From the executed DESeq2 model over all 12 groups; DEGs at P.adj < 0.05, |log₂FC| ≥ 1, ashr-shrunk.)*
The spaceflight (Flight vs Ground) response was **markedly stronger in root than shoot** across all
genotypes. In roots, wild type showed 2,503 up / 578 down DEGs and the AVP-OX line A68 the largest
response of all (2,806 up / 3,610 down), whereas D130 roots responded more modestly (420 / 299). Shoot
responses were smaller and skewed toward down-regulation (e.g. D130 shoot 327 up / 890 down; A68 shoot
79 / 419). The genotype × treatment **interaction contrasts — the AVP-OX flight response relative to wild
type — were consistently down-biased** (A68-vs-WT root 55 down / 2 up; D130-vs-WT root 20 down / 0 up),
i.e. the engineered lines *attenuate* rather than amplify the wild-type spaceflight transcriptional
programme (Table 2), consistent with the dampened defence/stress signatures in R5–R6.

Principal-component analysis placed **tissue as the dominant transcriptomic axis** (PC1, 92.5% of
variance, associated with tissue at P = 2.9×10⁻⁹) and **spaceflight as the second** (PC2, associated with
treatment at P = 6.3×10⁻⁶; Fig 3a). A full factorial model (Treatment × Genotype × Tissue) confirmed that
the spaceflight response is overwhelmingly **tissue-specific**: the Treatment × Tissue interaction
comprised 3,374 DEGs — far larger than the Treatment × Genotype (AVP-OX) interactions (A68 62, D130 20;
both down-biased) or the negligible genotype × tissue and three-way terms. Total read depth differed by
tissue (P ≈ 0) but **not by treatment** (P = 0.41), so the Flight-vs-Ground comparisons are not
depth-confounded. [Fig 3: PCA + DEG-count barplot; per-contrast tables in `../deseq2/contrasts/`,
factorial model in `../deseq2/v4/`.]

### R4. Tissue-specific DEG clusters and trait linkage  → **Fig 4**
[STUB] Organ-specific clustering / WGCNA modules; correlate module eigengenes with Flight, genotype, and
**root morphometric traits** (the RSML link). [PLACEHOLDER].

### R5. Functional enrichment of the spaceflight response  → **Fig 5**, **Table 3**
*(From the executed clusterProfiler run over all 10 DESeq2 contrasts; Arabidopsis-mapped, expressed-gene universe.)*
In **wild-type roots**, spaceflight up-regulated a canonical low-oxygen/mechanical-stress programme —
*cellular response to hypoxia / decreased oxygen* (P.adj = 7.3×10⁻²⁰) and *response to wounding* — while
down-regulating *translation* (P.adj = 1.5×10⁻²²), *ribosome biogenesis*, and *rRNA processing*,
indicating suppression of growth-associated protein synthesis. The AVP-OX line **A68 shifted this response
toward defence/immunity** (up-regulated *defence response to bacterium*, *immune response*), and D130 showed
the same themes more weakly. Directly contrasting each AVP-OX line's flight response with wild type
(interaction), the genes **specifically down-regulated in AVP-OX were enriched for jasmonic-acid
biosynthesis/metabolism and response to wounding** (A68-vs-WT root P.adj = 3.8×10⁻⁴; KEGG *α-linolenic-acid
metabolism*, *plant–pathogen interaction*), with up-regulation of *negative regulators of defence*. Thus,
against a shared hypoxia/wounding spaceflight signature, the engineered lines **attenuate the jasmonate/
defence arm** — **converging with the independent PhysioSpace result (R6)**.

### R6. PhysioSpace stress-pattern decoding  → **Fig 6**, **Table 4**
*(PRELIMINARY — from executed PhysioSpace run; VST group-mean input, `STATICResponse` scoring.)*
Projecting the Flight−Ground contrasts onto the Arabidopsis stress reference spaces (Plant PhysioSpace;
13,367 shared orthologs) showed that **spaceflight roots of all three genotypes engage a coordinated
abiotic-stress / defence programme** — most strongly osmotic, drought, far-red/shade, genotoxic, and
wounding axes — whereas the **shoot response was weaker and more heterogeneous**. Notably, the two
**AVP-OX lines showed *attenuated* root stress-programme activation relative to wild type** (e.g. osmotic
axis WT > A68 > D130; wounding axis WT 0.72 vs D130 0.44), a transcriptomic signature consistent with the
engineered stress tolerance operating under spaceflight. [Descriptive; formal genotype contrast on
PhysioScores + the frozen-model input pending. Root scores saturate under signed-p scoring — the bounded
`STATICResponse` statistic is reported; see `../physiospace/README.md`.]

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
(Phytozome Ghirsutum_527_v2.1; `Gohir` IDs). Filtered counts (~59.9k genes × 48) analysed with **DESeq2** using a Treatment × Genotype × Tissue model (references Ground/WT/Root; ashr LFC shrinkage). Both per-group Flight−Ground contrasts (`deseq2/contrasts/`) and the factorial interaction model (`deseq2/v4/`) were fit. Read depth differed by tissue but not treatment, so size-factor normalisation is not confounded with the spaceflight comparison. Package versions in `ENVIRONMENT.txt`.

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
