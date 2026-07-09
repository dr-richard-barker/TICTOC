#!/usr/bin/env python3
"""
Regenerate the matplotlib manuscript figures reproducibly:
  Fig 1 (design schematic), Fig 2 (root architecture), Fig 3b (DEG counts),
  Fig 5 (GO dotplot), Fig 6 (PhysioScore heatmap), Fig 7 (integrative model).
(Fig 3a PCA -> deseq2/make_pca.R ; Fig 4 module-trait heatmap -> integration/make_named_heatmap.R.)

Run from the repo root:  python manuscript/figures/make_figures.py
Deps: matplotlib, numpy (stdlib csv). Reads committed CSVs; writes PDF+PNG into manuscript/figures/
(and Fig 4 lives in manuscript/). Colours: Okabe-Ito-based, CVD-safe (Flight #E69F00, Ground #0072B2).
"""
import os, csv, math
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from matplotlib.colors import TwoSlopeNorm

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
OUT = HERE
FL, GC, TQ, BR = "#E69F00", "#0072B2", "#1abc9c", "#8a6d3b"
UP, DN = "#c0392b", "#0072B2"


def save(fig, name):
    fig.savefig(os.path.join(OUT, name + ".pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(OUT, name + ".png"), bbox_inches="tight", dpi=140)
    plt.close(fig); print("wrote", name)


def box(ax, x, y, w, h, fc="#fff", ec="#ccc", lw=1.2):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.2,rounding_size=1.5",
                 fc=fc, ec=ec, lw=lw, mutation_aspect=0.6))


def arrow(ax, x1, y1, x2, y2, c="#999", lw=1.6):
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1), arrowprops=dict(arrowstyle="-|>", color=c, lw=lw))


# ---------------- Fig 1: design ----------------
def fig1():
    fig, ax = plt.subplots(figsize=(11, 7)); ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off"); ax.invert_yaxis()
    ax.text(2, 5, "Fig 1. TIC-TOC experimental design and analysis workflow", fontsize=15, fontweight="bold")
    ax.text(2, 9.5, "AVP-OX cotton (Gossypium hirsutum) aboard the ISS (SpaceX CRS-22) vs matched ground controls", fontsize=9.5, color="#555")
    ax.text(2, 17, "A. Full factorial design — 4 replicates/group, 48 RNA-seq libraries", fontsize=11, fontweight="bold")
    for g, x in zip(["WT", "A68 (AVP-OX)", "D130 (AVP-OX)"], [30, 52, 74]): ax.text(x, 22, g, ha="center", fontsize=10, fontweight="bold")
    for ty, tl in zip([27, 34], ["Root", "Shoot"]):
        ax.text(18, ty + 2.4, tl, ha="right", fontsize=10)
        for x in [30, 52, 74]:
            box(ax, x - 8, ty, 7, 4.2, fc=FL); ax.text(x - 4.5, ty + 2.4, "FL×4", ha="center", va="center", fontsize=8, color="#fff")
            box(ax, x + 1, ty, 7, 4.2, fc=GC); ax.text(x + 4.5, ty + 2.4, "GC×4", ha="center", va="center", fontsize=8, color="#fff")
    box(ax, 88, 27, 3, 3, fc=FL); ax.text(92, 29.3, "Flight", fontsize=8.5, va="center")
    box(ax, 88, 32, 3, 3, fc=GC); ax.text(92, 34.3, "Ground", fontsize=8.5, va="center")
    ax.text(2, 45, "B. Two paired readouts", fontsize=11, fontweight="bold")
    box(ax, 2, 48, 44, 9); ax.text(4, 51.5, "Root architecture", fontsize=10, fontweight="bold"); ax.text(4, 54.5, "days 3–6 imaging → SmartRoot → RSML · 198 tracings, 53 plants", fontsize=8.5, color="#555")
    box(ax, 50, 48, 48, 9); ax.text(52, 51.5, "Transcriptome", fontsize=10, fontweight="bold"); ax.text(52, 54.5, "root + shoot RNA-seq · UTX-TM1 v2.1 · 59,918 genes · DESeq2", fontsize=8.5, color="#555")
    ax.text(2, 64, "C. Analysis workflow", fontsize=11, fontweight="bold")
    for s, x in zip(["DESeq2", "crosswalk\nGohir→AT", "GO / KEGG", "PhysioSpace", "WGCNA", "Integration\nexpr ↔ growth"], [2, 18, 34, 50, 66, 82]):
        box(ax, x, 67, 14, 7, fc="#f4f4f2"); ax.text(x + 7, 71, s, ha="center", va="center", fontsize=8)
    for x in [2, 18, 34, 50, 66]: arrow(ax, x + 14, 70.5, x + 16, 70.5)
    ax.text(2, 82, "D. Root imaging time series", fontsize=11, fontweight="bold")
    ax.plot([14, 58], [88, 88], color="#888", lw=1.5)
    for d, x in zip([3, 4, 5, 6], [16, 30, 44, 58]): ax.plot(x, 88, "o", color="#555", ms=7); ax.text(x, 92, f"Day {d}", ha="center", fontsize=8.5, color="#555")
    save(fig, "Fig1_design")


# ---------------- Fig 2: root architecture ----------------
def fig2():
    rows = list(csv.DictReader(open(os.path.join(REPO, "morphometrics", "rsml_traits_summary.csv"))))
    genos = ["WT", "A68", "D130"]; days = [3, 4, 5, 6]
    def series(g, cond, col):
        m, e = [], []
        for d in days:
            r = next((x for x in rows if x["genotype"] == g and x["condition"] == cond and int(x["day"]) == d), None)
            m.append(float(r[col + "_mean"]) if r and r[col + "_mean"] else np.nan)
            sd = float(r[col + "_sd"]) if r and r[col + "_sd"] else 0; n = int(r["n_plants"]) if r else 1
            e.append(sd / math.sqrt(max(n, 1)))
        return np.array(m), np.array(e)
    fig, axs = plt.subplots(2, 3, figsize=(11, 6.4), sharex=True)
    for row, (col, lbl) in enumerate([("total_length_native", "Total root length"), ("primary_length_native", "Primary root length")]):
        for j, g in enumerate(genos):
            ax = axs[row][j]
            for cond, c, nm in [("FL", FL, "Flight"), ("GC", GC, "Ground")]:
                m, e = series(g, cond, col)
                ax.errorbar(days, m, yerr=e, marker="o", ms=5, lw=2, color=c, capsize=3, label=nm)
            if row == 0: ax.set_title(g, fontsize=11, fontweight="bold")
            if j == 0: ax.set_ylabel(lbl + "\n(RSML units)", fontsize=9)
            if row == 1: ax.set_xlabel("Day", fontsize=9)
            ax.set_xticks(days); ax.grid(alpha=0.25)
    axs[0][2].legend(fontsize=8, frameon=True)
    fig.suptitle("Fig 2. Root architecture days 3–6 (mean ± SE) — total length is lateral-driven; primary length is the robust AVP-OX signal", fontsize=10.5)
    fig.tight_layout(rect=[0, 0, 1, 0.96]); save(fig, "Fig2_root_architecture")


# ---------------- Fig 3b: DEG counts ----------------
def fig3b():
    rows = [r for r in csv.DictReader(open(os.path.join(REPO, "deseq2", "contrasts", "DEG_counts_summary.csv"))) if r["type"] == "spaceflight"]
    labs = [r["contrast"].replace("_FlightVsGround", "").replace("_", " ") for r in rows]
    up = [int(r["up"]) for r in rows]; dn = [-int(r["down"]) for r in rows]
    y = np.arange(len(rows))
    fig, ax = plt.subplots(figsize=(7.6, 4.4))
    ax.barh(y, up, color=UP, label="Up in Flight"); ax.barh(y, dn, color=DN, label="Down in Flight")
    for yi, (u, d) in enumerate(zip(up, dn)):
        ax.text(u + 40, yi, str(u), va="center", fontsize=8, color=UP)
        ax.text(d - 40, yi, str(-d), va="center", ha="right", fontsize=8, color=DN)
    ax.set_yticks(y); ax.set_yticklabels(labs, fontsize=9); ax.invert_yaxis()
    ax.axvline(0, color="#333", lw=1); ax.set_xlabel("DEGs (P.adj < 0.05, |log₂FC| ≥ 1)", fontsize=9)
    ax.legend(fontsize=8, loc="lower right"); ax.grid(axis="x", alpha=0.25)
    ax.set_title("Fig 3b. Spaceflight DEG counts per genotype × tissue (root ≫ shoot)", fontsize=11)
    fig.tight_layout(); save(fig, "Fig3b_DEG_counts")


# ---------------- Fig 5: GO dotplot ----------------
def fig5():
    base = os.path.join(REPO, "go_analysis", "results_full")
    contrasts = [("GO_BP_up_WT_Root_FlightVsGround.csv", "WT root ↑"),
                 ("GO_BP_down_WT_Root_FlightVsGround.csv", "WT root ↓"),
                 ("GO_BP_up_A68_Root_FlightVsGround.csv", "A68 root ↑"),
                 ("GO_BP_down_A68vsWT_Root_FlightInteraction.csv", "A68×Flight ↓")]
    order, data = [], {}
    for f, lab in contrasts:
        for r in list(csv.DictReader(open(os.path.join(base, f))))[:4]:
            t = r["Description"]
            if t not in order: order.append(t)
            data[(t, lab)] = (float(r["p.adjust"]), int(r["Count"]))
    clabs = [l for _, l in contrasts]; order = order[:16]
    fig, ax = plt.subplots(figsize=(8.0, 0.44 * len(order) + 1.8))
    for yi, t in enumerate(order):
        for xi, l in enumerate(clabs):
            if (t, l) in data:
                p, c = data[(t, l)]
                ax.scatter(xi, yi, s=20 + c * 3.2, c=[-math.log10(p)], cmap="viridis_r", vmin=2, vmax=20, edgecolors="#444", linewidths=0.4, zorder=3)
    ax.set_xticks(range(len(clabs))); ax.set_xticklabels(clabs, rotation=25, ha="right", fontsize=9)
    ax.set_yticks(range(len(order))); ax.set_yticklabels([t[:44] for t in order], fontsize=8)
    ax.set_xlim(-0.6, len(clabs) - 0.4); ax.set_ylim(-0.6, len(order) - 0.4); ax.invert_yaxis(); ax.grid(axis="y", ls=":", lw=0.4, alpha=0.5)
    sm = plt.cm.ScalarMappable(cmap="viridis_r", norm=plt.Normalize(2, 20)); sm.set_array([])
    fig.colorbar(sm, ax=ax, shrink=0.4, pad=0.02).set_label("−log₁₀ P.adj", fontsize=9)
    ax.scatter([], [], s=20 + 5 * 3.2, c="#888", label="5"); ax.scatter([], [], s=20 + 40 * 3.2, c="#888", label="40")
    ax.legend(loc="lower left", fontsize=7, frameon=True, title="genes", labelspacing=1.1, borderpad=0.8)
    ax.set_title("Fig 5. GO biological-process enrichment\nof spaceflight-responsive genes", fontsize=11.5)
    fig.tight_layout(); save(fig, "Fig5_GO_enrichment")


# ---------------- Fig 6: PhysioScore heatmap ----------------
def fig6():
    rows = list(csv.reader(open(os.path.join(REPO, "physiospace", "results_static", "PhysioScores_AT_Stress_Space.csv"))))
    hdr = [h.replace("_FlightVsGround", "").replace("_", " ") for h in rows[0][1:]]
    axes_lbl = [r[0] for r in rows[1:]]; M = np.array([[float(x) for x in r[1:]] for r in rows[1:]])
    order = ["WT Root", "A68 Root", "D130 Root", "WT Shoot", "A68 Shoot", "D130 Shoot"]
    idx = [hdr.index(o) for o in order]; M = M[:, idx]; hdr = order
    fig, ax = plt.subplots(figsize=(6.6, 7.2)); vmax = np.abs(M).max()
    im = ax.imshow(M, aspect="auto", cmap="RdBu_r", norm=TwoSlopeNorm(0, -vmax, vmax))
    ax.set_xticks(range(len(hdr))); ax.set_xticklabels(hdr, rotation=40, ha="right", fontsize=9)
    ax.set_yticks(range(len(axes_lbl))); ax.set_yticklabels(axes_lbl, fontsize=8.5); ax.axvline(2.5, color="#333", lw=1.2)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            ax.text(j, i, f"{M[i,j]:.2f}", ha="center", va="center", fontsize=6.2, color="white" if abs(M[i, j]) > vmax * 0.6 else "#222")
    ax.set_title("Fig 6. PhysioSpace stress-pattern decoding\n(Flight−Ground; bounded statistic; roots ≫ shoots)", fontsize=11)
    fig.colorbar(im, ax=ax, shrink=0.5, pad=0.02).set_label("PhysioScore", fontsize=9)
    fig.tight_layout(); save(fig, "Fig6_physioscores")


# ---------------- Fig 7: integrative model ----------------
def fig7():
    fig, ax = plt.subplots(figsize=(11, 7.4)); ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off"); ax.invert_yaxis()
    ax.text(2, 5, "Fig 7. Integrative model — AVP-OX cotton: less stress, more root", fontsize=15, fontweight="bold")
    box(ax, 2, 14, 20, 8, fc="#fbe4c8", ec=FL); ax.text(12, 17, "Spaceflight (ISS)", ha="center", fontsize=10, fontweight="bold"); ax.text(12, 20, "root ≫ shoot", ha="center", fontsize=8.5, color="#555")
    ax.text(30, 12, "WT root transcriptome", fontsize=10, fontweight="bold")
    box(ax, 30, 15, 26, 12); ax.text(31.5, 19, "▲ hypoxia / wounding", fontsize=9, color=UP); ax.text(31.5, 23, "▼ translation / ribosome", fontsize=9, color=DN)
    ax.text(60, 12, "Co-expression modules", fontsize=10, fontweight="bold")
    box(ax, 60, 15, 24, 3.6, fc="#d6f0ee", ec=TQ); ax.text(61, 17.4, "▲ turquoise · signalling/isoprenoid", fontsize=8)
    box(ax, 60, 19.2, 24, 3.6, fc="#dbe8fb", ec=GC); ax.text(61, 21.6, "▼ blue · translation", fontsize=8)
    box(ax, 60, 23.4, 24, 3.6, fc="#e8ddcf", ec=BR); ax.text(61, 25.8, "brown · defence/ubiquitin", fontsize=8)
    ax.text(87, 12, "Phenotype", fontsize=10, fontweight="bold")
    box(ax, 87, 15, 12, 5.4); ax.text(93, 17.8, "osmotic/\nwounding", ha="center", va="center", fontsize=8)
    box(ax, 87, 21.4, 12, 5.4, fc="#fbe4c8", ec=FL); ax.text(93, 24.2, "lateral-root\nproliferation", ha="center", va="center", fontsize=8)
    arrow(ax, 22, 18, 29, 20); arrow(ax, 56, 20, 59, 20); arrow(ax, 84, 17.5, 86, 17.5); arrow(ax, 84, 25, 86, 24)
    ax.annotate("", xy=(90, 23), xytext=(84, 16.8), arrowprops=dict(arrowstyle="-|>", color="#159a82", lw=2, connectionstyle="arc3,rad=-0.35"))
    ax.text(72, 31, "turquoise module is growth-coupled", fontsize=8, color="#159a82", ha="center")
    box(ax, 2, 36, 96, 26, fc="#f4f8f4", ec="#bcd")
    ax.text(4, 41, "AVP-OX engineering acts as an attenuator, not an amplifier", fontsize=13, fontweight="bold")
    box(ax, 4, 44, 29, 15); ax.text(5.5, 47.5, "Transcriptome", fontsize=10, fontweight="bold"); ax.text(5.5, 51, "small T×G interaction", fontsize=8.5, color="#555"); ax.text(5.5, 54, "▼ jasmonate/defence (GO)", fontsize=8.5, color=DN); ax.text(5.5, 57, "▼ brown module (r=−0.43)", fontsize=8.5, color=DN)
    box(ax, 36, 44, 29, 15); ax.text(37.5, 47.5, "Stress state", fontsize=10, fontweight="bold"); ax.text(37.5, 51, "▼ osmotic/drought/wounding", fontsize=8.5, color=DN); ax.text(37.5, 54, "PhysioScores vs WT", fontsize=8.5, color="#555"); ax.text(37.5, 57, "(WT > A68 > D130)", fontsize=8.5, color="#555")
    box(ax, 68, 44, 30, 15, fc="#fbe4c8", ec=FL); ax.text(69.5, 47.5, "Root growth", fontsize=10, fontweight="bold"); ax.text(69.5, 51, "▲ larger primary-root", fontsize=8.5, color=UP); ax.text(69.5, 54, "flight response than WT", fontsize=8.5, color="#555"); ax.text(69.5, 57, "(P ≤ 0.04, robust)", fontsize=8.5, color="#555")
    arrow(ax, 33, 51.5, 36, 51.5, "#159a82"); arrow(ax, 65, 51.5, 68, 51.5, "#159a82")
    ax.text(2, 68, "Caveats: RNA-seq unpaired to imaged plants → integration group-level (n=6), Flight/growth collinear (correlational);", fontsize=8.5, color="#555")
    ax.text(2, 71, "root-growth magnitude pending FL/GC image-scale calibration. Candidate mechanism: isoprenoid → growth-hormone precursors.", fontsize=8.5, color="#555")
    ax.text(2, 76, "Implication: stress-tolerance engineering (AVP1-OX) yields a calmer-yet-more-vigorous root phenotype for space agriculture.", fontsize=11, fontweight="bold")
    save(fig, "Fig7_integrative_model")


if __name__ == "__main__":
    fig1(); fig2(); fig3b(); fig5(); fig6(); fig7()
    print("All figures regenerated in", OUT)
