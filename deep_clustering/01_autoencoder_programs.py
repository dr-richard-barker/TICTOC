#!/usr/bin/env python3
"""
Stage 1 — autoencoder-based gene expression PROGRAMS (bulk RNA-seq).

IMPORTANT: this is BULK RNA-seq (48 libraries), not single-cell data. We do NOT cluster cells.
We train an autoencoder on gene expression profiles (each gene = its 48-condition vector), learn a
latent embedding, and cluster GENES into co-expression PROGRAMS (a neural analogue of WGCNA). Programs
are later associated with cell types by marker enrichment (stage 2) — an in-silico, marker-based label,
never a claim of literal single cells.

Outputs (deep_clustering/results/):
  program_assignments.csv  gene, program
  program_group_means.csv  program x 12 (genotype x tissue x treatment) mean z-expression
  gene_umap.csv            gene, umap1, umap2, program   (for the figure)
  ae_summary.txt

Deps: torch, scikit-learn, umap-learn, numpy, pandas.
"""
import os, numpy as np, pandas as pd, torch, torch.nn as nn
from sklearn.cluster import KMeans
import umap

SEED = 0
np.random.seed(SEED); torch.manual_seed(SEED)
HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.abspath(os.path.join(HERE, ".."))
OUT = os.path.join(HERE, "results"); os.makedirs(OUT, exist_ok=True)
K = 12                      # number of expression programs
N_VAR = 8000                # top-variable genes used
LATENT = 12

# ---- load + normalise ----
df = pd.read_csv(os.path.join(REPO, "TICTOC_run1_filteredCounts_v3.csv"), index_col=0)
genes0 = df.index.to_numpy(); X = df.to_numpy(dtype=float)          # genes x 48
cols = df.columns.to_numpy()
cpm = X / X.sum(0, keepdims=True) * 1e6
logx = np.log2(cpm + 1)
# top-variable genes
v = logx.var(1); keep = np.argsort(v)[::-1][:N_VAR]
genes = genes0[keep]; Z = logx[keep]
Z = (Z - Z.mean(1, keepdims=True)) / (Z.std(1, keepdims=True) + 1e-8)   # per-gene z-score (pattern, not level)
Xt = torch.tensor(Z, dtype=torch.float32)

# ---- autoencoder ----
class AE(nn.Module):
    def __init__(s, d, lat):
        super().__init__()
        s.enc = nn.Sequential(nn.Linear(d, 32), nn.ReLU(), nn.Linear(32, lat))
        s.dec = nn.Sequential(nn.Linear(lat, 32), nn.ReLU(), nn.Linear(32, d))
    def forward(s, x): z = s.enc(x); return s.dec(z), z

dev = "cuda" if torch.cuda.is_available() else "cpu"
ae = AE(Xt.shape[1], LATENT).to(dev); Xd = Xt.to(dev)
opt = torch.optim.Adam(ae.parameters(), lr=1e-3, weight_decay=1e-5)
lossf = nn.MSELoss()
ae.train()
for ep in range(400):
    opt.zero_grad(); rec, _ = ae(Xd); loss = lossf(rec, Xd); loss.backward(); opt.step()
    if (ep + 1) % 100 == 0: print(f"  epoch {ep+1}  recon MSE {loss.item():.4f}")
ae.eval()
with torch.no_grad(): _, lat = ae(Xd)
lat = lat.cpu().numpy()

# ---- cluster genes into programs ----
km = KMeans(n_clusters=K, random_state=SEED, n_init=10).fit(lat)
prog = km.labels_
pd.DataFrame({"gene": genes, "program": prog}).to_csv(os.path.join(OUT, "program_assignments.csv"), index=False)

# ---- program x group-mean z-expression (12 groups) ----
grp = np.array([c[:-1] if c[-1].isdigit() else c for c in cols])  # already group labels (dup names) -> use parsed
def parse(g):
    import re
    return re.sub(r'(Root|Shoot).*', '', g), re.sub(r'.*?(Root|Shoot).*', r'\1', g), re.sub(r'.*(Flight|Ground)$', r'\1', g)
glab = np.array(["%s_%s_%s" % parse(c) for c in cols])
groups = ["%s_%s_%s" % (ge, ti, tr) for ge in ["WT", "A68", "D130"] for ti in ["Root", "Shoot"] for tr in ["Flight", "Ground"]]
rows = []
for p in range(K):
    prof = Z[prog == p].mean(0)                       # mean z across the top-var genes of this program
    gm = {g: prof[glab == g].mean() for g in groups}
    rows.append({"program": p, "n_genes": int((prog == p).sum()), **gm})
pd.DataFrame(rows).to_csv(os.path.join(OUT, "program_group_means.csv"), index=False)

# ---- UMAP for the figure ----
um = umap.UMAP(n_neighbors=15, min_dist=0.3, random_state=SEED).fit_transform(lat)
pd.DataFrame({"gene": genes, "umap1": um[:, 0], "umap2": um[:, 1], "program": prog}).to_csv(
    os.path.join(OUT, "gene_umap.csv"), index=False)

with open(os.path.join(OUT, "ae_summary.txt"), "w") as f:
    f.write(f"Autoencoder gene-program clustering (BULK RNA-seq — programs, not cells)\n")
    f.write(f"genes used: {N_VAR} top-variable of {len(genes0)}; latent dim {LATENT}; K={K} programs; device {dev}\n")
    f.write("program sizes: " + ", ".join(f"{p}:{int((prog==p).sum())}" for p in range(K)) + "\n")
print("done. programs:", np.bincount(prog).tolist())
