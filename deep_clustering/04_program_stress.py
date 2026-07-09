#!/usr/bin/env python3
"""
Stage 4 — per-program stress response, PhysioSpace-style, compared WT vs transgenic.

For each expression program and genotype we compute the Flight-Ground log2 fold change of the program's
genes (in the program's dominant tissue), map Gohir->AT->Entrez, and score it against each Arabidopsis
stress axis (AT_Stress_Space) as the Pearson correlation between the program's fold-changes and the
axis's gene weights over shared genes (>= 15). This is a PhysioSpace-style, gene-set-restricted stress
score: positive = the program's flight response resembles that stress signature. We then contrast WT vs
the AVP-OX lines per program.

BULK data caveat: programs are gene co-expression clusters, not cells. Outputs feed Fig 9.
Deps: numpy, pandas, scipy.
"""
import os, re, numpy as np, pandas as pd
from scipy.stats import pearsonr

HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.abspath(os.path.join(HERE, ".."))
R = os.path.join(HERE, "results")
MIN_SHARED = 15

# ---- expression -> per genotype/tissue Flight-Ground log2FC ----
df = pd.read_csv(os.path.join(REPO, "TICTOC_run1_filteredCounts_v3.csv"), index_col=0)
X = df.to_numpy(float); cpm = np.log2(X / X.sum(0, keepdims=True) * 1e6 + 1)
cols = df.columns.to_numpy()
def parse(g): return (re.sub(r'(Root|Shoot).*', '', g), re.sub(r'.*?(Root|Shoot).*', r'\1', g), re.sub(r'.*(Flight|Ground)$', r'\1', g))
meta = np.array([parse(c) for c in cols])
def fc(geno, tis):
    fl = (meta[:, 0] == geno) & (meta[:, 1] == tis) & (meta[:, 2] == "Flight")
    gr = (meta[:, 0] == geno) & (meta[:, 1] == tis) & (meta[:, 2] == "Ground")
    return cpm[:, fl].mean(1) - cpm[:, gr].mean(1)
genes = df.index.to_numpy()

# ---- maps: Gohir -> AT -> Entrez ----
cw = pd.read_csv(os.path.join(REPO, "crosswalk", "gohir_to_arabidopsis.tsv"), sep="\t")
g2at = dict(zip(cw["gohir_gene_upper"].str.upper(), cw["arabidopsis_gene"]))
a2e = pd.read_csv(os.path.join(R, "at_to_entrez.csv")).dropna()
at2e = dict(zip(a2e["arabidopsis_gene"], a2e["entrez"].astype(int).astype(str)))
def gohir_to_entrez(gs):
    out = {}
    for g in gs:
        at = g2at.get(g.upper()); e = at2e.get(at) if at else None
        if e: out[g] = e
    return out

# ---- stress space + programs ----
sp = pd.read_csv(os.path.join(R, "AT_Stress_Space.csv")); sp["entrez"] = sp["entrez"].astype(str)
sp = sp.set_index("entrez")
sp = sp[~sp.index.duplicated(keep="first")]      # collapse duplicate Entrez so sp.at[e,a] is scalar
axes = list(sp.columns)
pa = pd.read_csv(os.path.join(R, "program_assignments.csv"))
pgm = pd.read_csv(os.path.join(R, "program_group_means.csv"))
# dominant tissue per program
def dom_tissue(p):
    r = pgm[pgm.program == p].iloc[0]
    root = np.mean([r[f"{g}_Root_Flight"] + r[f"{g}_Root_Ground"] for g in ["WT", "A68", "D130"]])
    shoot = np.mean([r[f"{g}_Shoot_Flight"] + r[f"{g}_Shoot_Ground"] for g in ["WT", "A68", "D130"]])
    return "Root" if root >= shoot else "Shoot"

genoset = ["WT", "A68", "D130"]
rows = []
for p in sorted(pa.program.unique()):
    tis = dom_tissue(p)
    pg = pa.gene[pa.program == p].to_numpy()
    g2e = gohir_to_entrez(pg)                      # program gene -> entrez
    gidx = {g: i for i, g in enumerate(genes)}
    for geno in genoset:
        f = fc(geno, tis)
        # entrez -> mean FC (collapse many gohir per entrez)
        acc = {}
        for g, e in g2e.items():
            if g in gidx: acc.setdefault(e, []).append(f[gidx[g]])
        fce = {e: float(np.mean(v)) for e, v in acc.items()}
        common = [e for e in fce if e in sp.index]
        for a in axes:
            sh = [e for e in common if not np.isnan(sp.at[e, a])]
            if len(sh) >= MIN_SHARED:
                r_, _ = pearsonr([fce[e] for e in sh], [sp.at[e, a] for e in sh])
            else:
                r_ = np.nan
            rows.append({"program": p, "tissue": tis, "genotype": geno, "axis": a,
                         "score": round(r_, 4) if not np.isnan(r_) else np.nan, "n_shared": len(sh)})
res = pd.DataFrame(rows)
res.to_csv(os.path.join(R, "program_stress_scores.csv"), index=False)

# WT-vs-transgenic delta per program x axis (mean AVP-OX - WT)
piv = res.pivot_table(index=["program", "axis"], columns="genotype", values="score")
piv["AVPOX_mean"] = piv[["A68", "D130"]].mean(1); piv["delta_vs_WT"] = piv["AVPOX_mean"] - piv["WT"]
piv.reset_index().to_csv(os.path.join(R, "program_stress_WTvsAVPOX.csv"), index=False)
print("wrote program_stress_scores.csv and program_stress_WTvsAVPOX.csv")
print("programs x axes scored:", res.dropna().program.nunique(), "programs")
