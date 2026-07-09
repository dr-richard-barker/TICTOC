#!/usr/bin/env python3
"""
Stage 5 — the two manuscript montages from the autoencoder program analysis:
  Fig 8  Expression-program atlas (gene UMAP + program x group-mean heatmap + annotations)
  Fig 9  Cellular-program stress response, WT vs transgenic (program x stress + AVP-OX-WT delta)
Programs are autoencoder gene co-expression clusters from BULK RNA-seq (not single cells).
Deps: numpy, pandas, matplotlib.
"""
import os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.abspath(os.path.join(HERE, ".."))
R = os.path.join(HERE, "results"); OUTDIR = os.path.join(REPO, "manuscript", "figures")

um = pd.read_csv(os.path.join(R, "gene_umap.csv"))
gm = pd.read_csv(os.path.join(R, "program_group_means.csv"))
go = pd.read_csv(os.path.join(R, "program_GO.csv")).set_index("program")
ct = pd.read_csv(os.path.join(R, "program_celltypes.csv")).set_index("program")
progs = sorted(gm.program)
groups = [f"{g}_{t}_{c}" for g in ["WT", "A68", "D130"] for t in ["Root", "Shoot"] for c in ["Flight", "Ground"]]

def dom_tissue(p):
    r = gm[gm.program == p].iloc[0]
    root = np.mean([r[f"{g}_Root_Flight"] + r[f"{g}_Root_Ground"] for g in ["WT", "A68", "D130"]])
    shoot = np.mean([r[f"{g}_Shoot_Flight"] + r[f"{g}_Shoot_Ground"] for g in ["WT", "A68", "D130"]])
    return "Root" if root >= shoot else "Shoot"

def short(s, n=26):
    s = str(s).replace(" process", "").replace("biosynthetic", "biosynth.").replace("metabolic", "metab.")
    return s[:n]

def plabel(p):
    fn = short(go.loc[p, "top1"]) if p in go.index and str(go.loc[p, "top1"]) != "nan" else "—"
    cell = ct.loc[p, "top_celltype"] if p in ct.index and isinstance(ct.loc[p, "top_celltype"], str) and ct.loc[p, "top_celltype"] else ""
    tag = f" · {cell}" if cell else ""
    return f"P{p} {dom_tissue(p)[0]} · {fn}{tag}"

cmap12 = plt.get_cmap("tab20")

# ================= Fig 8 =================
fig = plt.figure(figsize=(14, 6.4))
gs = fig.add_gridspec(1, 2, width_ratios=[1.35, 1], wspace=0.55)
# A (left): program x group heatmap — function labels on the far left margin
axB = fig.add_subplot(gs[0])
M = gm.set_index("program")[groups].loc[progs].to_numpy()
vmax = np.abs(M).max()
im = axB.imshow(M, aspect="auto", cmap="RdBu_r", norm=TwoSlopeNorm(0, -vmax, vmax))
axB.set_xticks(range(len(groups))); axB.set_xticklabels([g.replace("_", " ") for g in groups], rotation=55, ha="right", fontsize=7.5)
axB.set_yticks(range(len(progs))); axB.set_yticklabels([plabel(p) for p in progs], fontsize=8)
for x in (4, 8): axB.axvline(x - 0.5, color="#333", lw=0.8)
axB.axvline(11.5, color="#333", lw=1.4)
axB.set_title("a  Program mean expression across genotype × tissue × treatment (z)", fontsize=11, loc="left")
fig.colorbar(im, ax=axB, shrink=0.5, pad=0.02, label="mean z")
# B (right): UMAP with program numbers annotated at cluster centroids (no legend box)
axA = fig.add_subplot(gs[1])
for p in progs:
    s = um[um.program == p]
    axA.scatter(s.umap1, s.umap2, s=4, color=cmap12(p % 20), alpha=0.7, linewidths=0)
    cx, cy = s.umap1.median(), s.umap2.median()
    axA.text(cx, cy, f"P{p}", fontsize=9, fontweight="bold", ha="center", va="center",
             bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.75))
axA.set_title("b  Gene program map (autoencoder + UMAP)", fontsize=11, loc="left")
axA.set_xlabel("UMAP-1"); axA.set_ylabel("UMAP-2"); axA.set_xticks([]); axA.set_yticks([])
fig.suptitle("Fig 8. Autoencoder expression programs — tissue-specific co-expression clusters (bulk RNA-seq; programs, not single cells)", fontsize=12, y=1.02)
fig.savefig(os.path.join(OUTDIR, "Fig8_program_atlas.pdf"), bbox_inches="tight")
fig.savefig(os.path.join(OUTDIR, "Fig8_program_atlas.png"), bbox_inches="tight", dpi=140)
plt.close(fig); print("wrote Fig8_program_atlas")

# ================= Fig 9 =================
ss = pd.read_csv(os.path.join(R, "program_stress_scores.csv"))
axes = sorted(ss.axis.unique(), key=lambda a: ss[ss.axis == a].score.abs().mean(), reverse=True)
# keep programs that have any scored axis
scored = [p for p in progs if not ss[(ss.program == p) & ss.score.notna()].empty]
def mat(geno):
    piv = ss[ss.genotype == geno].pivot_table(index="program", columns="axis", values="score")
    return piv.reindex(index=scored, columns=axes)
WT = mat("WT"); AV = (mat("A68") + mat("D130")) / 2; DELTA = AV - WT
fig, (a1, a2) = plt.subplots(1, 2, figsize=(14, 6.2))
vmax = np.nanmax(np.abs(WT.to_numpy()))
im1 = a1.imshow(WT.to_numpy(), aspect="auto", cmap="RdBu_r", norm=TwoSlopeNorm(0, -vmax, vmax))
a1.set_title("a  WT — program flight response vs stress axes", fontsize=11, loc="left")
im2 = a2.imshow(DELTA.to_numpy(), aspect="auto", cmap="PuOr_r", norm=TwoSlopeNorm(0, -np.nanmax(np.abs(DELTA.to_numpy())), np.nanmax(np.abs(DELTA.to_numpy()))))
a2.set_title("b  AVP-OX − WT (Δ) — attenuated stress (orange) vs sustained (purple)", fontsize=11, loc="left")
for ax, im, lab in ((a1, im1, "PhysioSpace-style score"), (a2, im2, "AVP-OX − WT Δ")):
    ax.set_xticks(range(len(axes))); ax.set_xticklabels(axes, rotation=55, ha="right", fontsize=8)
    ax.set_yticks(range(len(scored))); ax.set_yticklabels([plabel(p) for p in scored], fontsize=8)
    fig.colorbar(im, ax=ax, shrink=0.5, pad=0.02, label=lab)
fig.suptitle("Fig 9. Program-level stress response to spaceflight, WT vs AVP-OX (bulk RNA-seq gene programs)", fontsize=12, y=1.01)
fig.tight_layout(rect=[0, 0, 1, 0.98])
fig.savefig(os.path.join(OUTDIR, "Fig9_program_stress_WTvsAVPOX.pdf"), bbox_inches="tight")
fig.savefig(os.path.join(OUTDIR, "Fig9_program_stress_WTvsAVPOX.png"), bbox_inches="tight", dpi=140)
plt.close(fig); print("wrote Fig9_program_stress_WTvsAVPOX")
