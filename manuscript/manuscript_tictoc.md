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
**root + shoot transcriptome** (RNA-seq, full Treatment × Genotype × Tissue design). Spaceflight roots
proliferated lateral roots, and both AVP-OX lines mounted a larger primary-root spaceflight response than
wild type (P ≤ 0.04). The transcriptomic response was strongly **tissue-specific** (Treatment × Tissue the
dominant interaction; root ≫ shoot), and in wild-type roots comprised a canonical hypoxia/wounding
induction with translational suppression. Across three independent methods — differential-expression
interaction tests, ortholog-mapped GO/KEGG, and PhysioSpace stress-pattern decoding — the AVP-OX
engineering **attenuated** rather than amplified the spaceflight defence/stress programme (suppressed
jasmonate/wounding signalling; lower osmotic/wounding PhysioScores). Co-expression analysis tied a
signalling/isoprenoid-metabolism module to the root-growth phenotype and a defence module to the AVP-OX
effect. Thus AVP-OX cotton appears to experience microgravity as **less stressful while sustaining stronger
root growth** — a favourable trait for space agriculture, pending image-calibration and replication checks.
These results have implications for engineering crops for spaceflight environments.

---

## Introduction
- Roots, gravity, and cotton productivity; auxin-directed root system architecture. [refs]
- AVP1/AVP-OX technology: vacuolar H⁺-PPase overexpression → salt/drought tolerance, larger root
  systems, +20% fibre yield under field stress in cotton (Pasapula et al. 2011¹; Arabidopsis origin Park et al. 2005²).
- Spaceflight remodels plant roots and transcriptomes; microgravity disrupts gravitropic auxin
  transport and cytoskeletal/cell-wall processes. [refs]
- Gap: whether stress-tolerance engineering changes the *spaceflight* response has not been tested in a
  crop, at the intersection of architecture and transcriptome.
- This study: WT vs AVP-OX cotton, Flight vs Ground, root vs shoot, days 3–6 → morphometrics +
  transcriptome, decoded into stress programs via ortholog-mapped GO and PhysioSpace.

---

## Results

### R1. Experimental design and dataset overview  → **Fig 1**, **Table 1**
Wild-type cotton and two independent AVP1-overexpressing lines (A68, D130) were grown aboard the ISS
(SpaceX CRS-22) alongside matched ground controls, in a full factorial **Treatment (Flight/Ground) ×
Genotype (WT/A68/D130) × Tissue (root/shoot)** design with four replicates per group (**48 RNA-seq
libraries**; Table 1). Root systems were imaged daily on days 3–6 and traced (SmartRoot → RSML), yielding
198 curated tracings from 53 plants, 44 with the complete four-day series (Fig 1; `../Data/RSML_QC_summary.md`).
Transcriptomes were quantified against the *G. hirsutum* UTX-TM1 v2.1 assembly (59,918 expressed genes).
Together these provide paired **root-architecture and transcriptomic** readouts of the spaceflight response
across genotypes and tissues.

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

### R4. Co-expression modules link the spaceflight transcriptome to root architecture  → **Fig 4**
*(From the executed root WGCNA; module eigengenes vs treatment/genotype/RSML traits. Trait correlations are group-level — see caveat.)*
Weighted co-expression network analysis of root samples resolved six modules, named here by their GO
enrichment (BP/MF/CC): **turquoise — *signalling & isoprenoid metabolism*** (1,963 genes),
**blue — *translation & ribosome biogenesis*** (1,767), **brown — *defence & ubiquitin signalling***
(1,046), **yellow — *metal transport & phenylpropanoid metabolism*** (68), **green — *photosynthesis
(light reactions)*** (38) and an unassigned grey module. Two formed a **flight-responsive, growth-coupled
axis**: turquoise was strongly induced by spaceflight and positively correlated with a larger, more-branched
root system (module-eigengene r with Flight/total-length/lateral-count = +0.89/+0.86/+0.83), and blue was
its suppressed mirror (Flight r = −0.83). Critically, the **defence & ubiquitin (brown) module** was
flight-induced yet **down-regulated in the AVP-OX lines** (r = −0.43), providing a module-level
transcriptomic correlate of the attenuated defence/stress response seen independently by enrichment (R5)
and PhysioSpace (R6). Thus a single, interpretable co-expression structure ties the spaceflight response,
the engineered-genotype effect, and the root-architecture phenotype together (Fig 4).

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
Across four independent analyses a single, consistent model emerges. Spaceflight elicits a strongly
**tissue-specific** transcriptional response (Treatment × Tissue ≫ all other interactions; R3), which in
roots is a canonical **hypoxia/wounding induction with translational suppression** (R5) that maps onto a
coordinated **osmotic/drought/wounding stress programme** (R6) and a growth-coupled co-expression module
(turquoise; R4). Against this shared wild-type response, the **AVP-OX engineering acts as an attenuator,
not an amplifier, of stress signalling**: the genotype × treatment interaction is small and down-biased
(R3), specifically suppressing the jasmonate/defence arm (R5), lowering PhysioSpace stress-programme
activation (R6), and down-regulating the flight-induced *brown* module (R4). Yet at the organ level the
AVP-OX lines mount a *larger* primary-root spaceflight response (R2). The coherent interpretation: **AVP-OX
cotton experiences microgravity as less of a stress — mounting a calmer defence/stress transcriptional
programme while sustaining stronger root growth** — a favourable phenotype for space agriculture, pending
the calibration and full-replication checks noted throughout. [Fig 7: schematic integrating
root-trait ↔ module ↔ stress-programme layers.]

To connect expression to the growth phenotype explicitly, we correlated module eigengenes and genes with
root-architecture traits (group level) and ranked genes by module membership and gene significance for
total root length. The growth phenotype co-occurred with induction of the **signalling & isoprenoid-metabolism
(turquoise)** module and suppression of the **translation & ribosome (blue)** module (1,678 and 1,444
hub-and-growth genes). Because isoprenoid/terpenoid metabolism supplies precursors for growth hormones
(gibberellins, brassinosteroids, strigolactones), this nominates **hormone-precursor metabolism as a
testable candidate route** from the spaceflight signal to altered root growth. *Statistical limitation:*
RNA-seq libraries are not individually paired to imaged plants, so the integration is group-level (n = 6),
and across those groups spaceflight and traced root growth are **collinear** — it identifies pathways
*consistent with* driving growth but cannot separate growth-driving from flight-responsive effects. A causal
model (paired multi-omics/DIABLO or Treatment→expression→trait mediation) is proposed for follow-up
(see `../integration/README.md`).

---

## Discussion
- Interpret the architecture ↔ transcriptome ↔ stress-program integration. [PLACEHOLDER]
- AVP-OX in space: buffering vs amplification; auxin/gravitropism and cell-wall themes. [PLACEHOLDER]
- **Named co-expression modules** (turquoise = signalling & isoprenoid metabolism; blue = translation &
  ribosome; brown = defence & ubiquitin; yellow = metal/phenylpropanoid; green = photosynthesis) give an
  interpretable middle layer linking the spaceflight signal to root growth; isoprenoid → growth-hormone
  precursor metabolism is the leading mechanistic hypothesis.
- Limitations: best-hit orthology (approximate); unbalanced morphometric design; single mission;
  **RNA-seq and imaging are not individually paired, so transcriptome↔architecture integration is
  group-level (n=6) and Flight/growth are collinear — correlational, not causal.** A paired multi-omics
  design (or mediation model) is needed to test the isoprenoid-hormone route directly.
- Implications for engineering crops for spaceflight agriculture.

---

## Methods
*(Concrete — fill bracketed hardware/version specifics.)*

**Plant material & spaceflight.** Cotton (*Gossypium hirsutum*): wild type and two independent
AVP1-overexpressing lines (A68, D130; AVP1 confers drought/salt tolerance and increased fibre yield¹).
Plants were grown aboard the ISS (SpaceX CRS-22, launched 3 June 2021; Expedition 65) in the **Vegetable
Production System (Veggie)** facility, in custom-designed **Target Veggie Chambers** (≈10 × 3 in, clear
plastic with a translucent gel substrate in place of soil so roots were visible for imaging). Fifteen
seed samples were launched and twelve grown; crew captured high-resolution images over the ≈6-day
experiment (days 3–6 analysed here). Matched ground controls were grown in identical Veggie hardware at
NASA Kennedy Space Center within the ISS Environmental Simulator (matched temperature, CO₂ and O₂), on a
48-h offset. Sponsored by the ISS National Laboratory (CASIS UA-2018-276) with support from the Target
Corporation.

**Root imaging & morphometrics.** Time-series photographs (grayscale, 300 dpi; flight frames
clarity-enhanced) were traced in **SmartRoot** and exported as **RSML**. Traits (root length, diameter,
angle, insertion; see `../DATA_DICTIONARY.md` §7) were extracted per plant/day with a stdlib parser
(`../morphometrics/extract_rsml_traits.py`). Curated set: 198 tracings, 53 plants (44 with the full
day-3–6 series). Because the design is unbalanced across genotype × condition, log-transformed total and
primary root length were modelled with **linear mixed models** (`statsmodels`; random intercept per
plant; fixed effects condition × genotype + day + day:condition, Ground/WT references). A primary-root
sensitivity analysis (immune to lateral-tracing depth) tested robustness of the genotype × condition
effects. *[Growth medium, planting density, harvest day, and image spatial calibration to add from the
flight protocol.]*

**RNA-seq & differential expression.** Root and shoot RNA-seq, WT/A68/D130 × Flight/Ground, 4
replicates. Reads *[quality-trimming, aligner and quantifier to add from the run-1 processing pipeline —
the exact tools/versions are needed from the sequencing/bioinformatics collaborators or the OSDR
processed-data record]*, quantified against the *G. hirsutum* accession TM-1 **UTX-TM1 v2.1** reference
genome³ (CottonGen/Phytozome Ghirsutum_527_v2.1; `Gohir` IDs). Filtered counts (~59.9k genes × 48) analysed with **DESeq2** using a Treatment × Genotype × Tissue model (references Ground/WT/Root; ashr LFC shrinkage). Both per-group Flight−Ground contrasts (`deseq2/contrasts/`) and the factorial interaction model (`deseq2/v4/`) were fit. Read depth differed by tissue but not treatment, so size-factor normalisation is not confounded with the spaceflight comparison. Package versions in `ENVIRONMENT.txt`.

**Ortholog mapping.** `Gohir` genes mapped to *Arabidopsis thaliana* by CottonGen BLASTP best-hit
(`../crosswalk/`; 92.3% of expressed genes mapped).

**Functional enrichment.** GO (BP/MF/CC) and KEGG (`ath`) enrichment in Arabidopsis space with
**clusterProfiler**, against the universe of AT loci reachable from expressed cotton genes
(`../go_analysis/`).

**Co-expression & integration.** Signed WGCNA on the top 5,000 variable root genes (soft power 12,
minModuleSize 30) gave six modules, named by GO (BP/MF/CC) enrichment. Module eigengenes were correlated
with Flight, genotype and group-level day-6 RSML root traits; within growth-associated modules, genes were
ranked by module membership and gene significance for total root length (`integration/`). Because RNA-seq
libraries are not individually paired to imaged plants, expression–trait correlations are computed on the
six genotype × treatment group means (hypothesis-generating; Flight and traced growth are collinear at this level).

**PhysioSpace.** Per-(genotype × tissue) Flight−Ground VST fold changes were re-indexed to Arabidopsis
Entrez and projected onto the Arabidopsis stress reference spaces with **PhysioSpaceMethods**
(`calculatePhysioMap`, GenesRatio = 0.05, TTEST = FALSE, ImputationMethod = "PCA") using the
**PlantPhysioSpace** data package (`../physiospace/`). Refs: Lenz et al. 2013; Hadizadeh Esfahani et al. 2021.

**Statistics.** DEGs: BH-adjusted *P* < 0.05 and |log₂FC| ≥ 1 (ashr-shrunk). GO/KEGG: BH *q* < 0.05
against the expressed-gene universe. PC–factor and module–trait associations: Kruskal–Wallis / Pearson.
Morphometrics: linear mixed models (above). Exact package versions in `../ENVIRONMENT.txt`.

---

## Figure & table legends

**Figure 1. Experimental design and analysis workflow.** Wild-type and two AVP1-overexpressing cotton
lines (A68, D130) were grown aboard the ISS (SpaceX CRS-22) with matched ground controls in a full
factorial Treatment × Genotype × Tissue design (4 replicates/group; 48 RNA-seq libraries). Root systems
were imaged on days 3–6 and traced (SmartRoot → RSML); root and shoot were profiled by RNA-seq. The
analysis pipeline (DESeq2 → *Gohir*→Arabidopsis crosswalk → GO/KEGG → PhysioSpace → WGCNA → integration)
is shown.

**Figure 2. Root architecture over days 3–6.** Group mean ± SE of total root length (top) and primary
root length (bottom) for Flight (orange) vs Ground (blue), per genotype (native RSML units; day-6 n per
group in Table 1). Total length diverges strongly under spaceflight whereas primary-root length differs
little, localising the flight effect to lateral roots. *Absolute magnitude pending Flight/Ground image-scale calibration.*

**Figure 3. Transcriptomic landscape.** (a) PCA of VST-transformed counts — PC1 (92.5% variance)
separates tissue, PC2 separates treatment. (b) Differentially expressed genes per Flight-vs-Ground
contrast (P.adj < 0.05, |log₂FC| ≥ 1), up in Flight (orange) vs down (blue), by genotype × tissue.

**Figure 4. Root co-expression modules linked to Flight, genotype and root architecture.** Correlations
between WGCNA module eigengenes (rows, named by GO enrichment) and Flight, AVP-OX and day-6 RSML root
traits (columns), at group level (n = 6). Turquoise (signalling & isoprenoid) tracks Flight and root
growth; brown (defence & ubiquitin) is down-regulated in AVP-OX.

**Figure 5. Functional enrichment.** GO/KEGG dot-/bar-plots of spaceflight-responsive gene sets,
ortholog-mapped to Arabidopsis, per contrast (`../go_analysis/results_full/`).

**Figure 6. PhysioSpace stress-pattern decoding.** PhysioScores (bounded `STATICResponse`) for each
genotype × tissue Flight−Ground contrast across Arabidopsis stress axes; AVP-OX roots show attenuated
osmotic/wounding activation vs WT.

**Figure 7. Integrative model.** Spaceflight elicits a tissue-specific hypoxia/wounding response with
translational suppression; AVP-OX attenuates the defence/stress arm while sustaining stronger root growth.

**Table 1.** Sample/design summary. **Table 2.** DEG counts per contrast (up/down). **Table 3.** Top
enriched GO/KEGG terms per contrast. **Table 4.** PhysioScores per genotype × tissue × stress axis.
See `Supplementary_index.md` for all supplementary items.

---

## Data & code availability
Raw reads, count matrices, root images, and RSML tracings: NASA OSDR/GeneLab **OSD-XXX** (released on
publication). Analysis code and derived tables: Zenodo **DOI 10.5281/zenodo.XXXXXXX** (`v1.0.0` release)
and https://github.com/dr-richard-barker/TICTOC. Data under CC0-1.0; analysis code under MIT.

## Author contributions / Acknowledgements / Competing interests
*[Author contributions per CRediT — to complete.]* Supported by CASIS UA-2018-276; conducted in the
Gilroy Life Science Lab, University of Wisconsin–Madison. The authors declare no competing interests.
*[Confirm.]*

## References
*(Method/anchor citations — verified where noted; confirm exact page numbers in your reference manager,
and add spaceflight-plant-biology background refs for the Introduction/Discussion.)*

1. Pasapula V, Shen G, Kuppu S, *et al.* Expression of an Arabidopsis vacuolar H⁺-pyrophosphatase gene
   (*AVP1*) in cotton improves drought- and salt tolerance and increases fibre yield in the field.
   *Plant Biotechnol J.* 2011;9(1):88–99. doi:10.1111/j.1467-7652.2010.00535.x  *(cotton AVP1 — verified)*
2. Park S, Li J, Pittman JK, *et al.* Up-regulation of a H⁺-pyrophosphatase (*AVP1*) as a strategy to
   engineer drought-resistant crop plants. *Proc Natl Acad Sci USA.* 2005;102(52):18830–18835.
3. Hu Y, Chen J, Fang L, *et al.* *Gossypium barbadense* and *Gossypium hirsutum* genomes provide
   insights into the origin and evolution of allotetraploid cotton. *Nat Genet.* 2019;51(4):739–748.
   *(TM-1 UTX reference genome — verified)*
4. Yu J, Jung S, Cheng CH, *et al.* CottonGen: a genomics, genetics and breeding database for cotton
   research. *Nucleic Acids Res.* 2014;42(D1):D1229–D1236.
5. Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with
   DESeq2. *Genome Biol.* 2014;15(12):550.
6. Stephens M. False discovery rates: a new deal. *Biostatistics.* 2017;18(2):275–294. *(ashr)*
7. Wu T, Hu E, Xu S, *et al.* clusterProfiler 4.0: a universal enrichment tool for interpreting omics
   data. *Innovation (Camb).* 2021;2(3):100141.
8. Langfelder P, Horvath S. WGCNA: an R package for weighted correlation network analysis.
   *BMC Bioinformatics.* 2008;9:559.
9. Lenz M, Müller FJ, Zenke M, Schuppert A. PhysioSpace: relating gene expression experiments from
   heterogeneous sources using shared physiological processes. *PLoS ONE.* 2013;8(10):e77627.
10. Hadizadeh Esfahani A, *et al.* Plant PhysioSpace: a robust tool to compare stress response across
    plant species. *Plant Physiol.* 2021;187(3):1795–1809.
11. Lobet G, Pagès L, Draye X. A novel image-analysis toolbox enabling quantitative analysis of root
    system architecture. *Plant Physiol.* 2011;157(1):29–39. *(SmartRoot)*
12. Lobet G, Pound MP, Diener J, *et al.* Root System Markup Language: toward a unified root architecture
    description language. *Plant Physiol.* 2015;167(3):617–627. *(RSML)*
13. Massa GD, Dufour NF, Crabb VA, *et al.* VEG-01: Veggie hardware validation testing on the
    International Space Station. *Open Agric.* 2017;2(1):33–41. *(Veggie facility)*
