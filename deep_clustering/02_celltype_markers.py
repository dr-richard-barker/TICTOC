#!/usr/bin/env python3
"""
Stage 2 — associate expression programs with cell types via a CURATED canonical Arabidopsis
marker panel (root + leaf). BULK data: this is a marker-overlap association, NOT single-cell
identity. Reports, per program, which cell-type markers it carries (mapped AT -> Gohir via crosswalk).
Deps: pandas.
"""
import os, pandas as pd
from collections import defaultdict
HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.abspath(os.path.join(HERE, ".."))
R = os.path.join(HERE, "results")

# curated canonical markers: AT locus -> (cell type, symbol)
MARKERS = {
 "AT1G12560": ("Root hair", "EXPA7"), "AT5G49270": ("Root hair", "COBL9"),
 "AT1G79840": ("Epidermis (non-hair)", "GL2"), "AT5G14750": ("Epidermis (non-hair)", "WER"),
 "AT3G54220": ("Endodermis", "SCR"), "AT5G57620": ("Endodermis", "MYB36"), "AT2G36100": ("Endodermis/Casparian", "CASP1"),
 "AT4G37650": ("Stele/vascular", "SHR"), "AT1G79430": ("Phloem", "APL"), "AT1G22710": ("Phloem", "SUC2"),
 "AT1G71930": ("Xylem", "VND7"), "AT5G12870": ("Xylem", "MYB46"),
 "AT3G11260": ("Quiescent centre", "WOX5"),
 "AT1G29930": ("Mesophyll", "CAB1"), "AT5G38430": ("Mesophyll", "RBCS1B"), "AT3G27690": ("Mesophyll", "LHCB2.3"),
 "AT1G08810": ("Guard cell", "MYB60"), "AT5G46240": ("Guard cell", "KAT1"), "AT3G24140": ("Guard cell", "FAMA"),
 "AT5G41315": ("Trichome", "GL3"),
}
cw = pd.read_csv(os.path.join(REPO, "crosswalk", "gohir_to_arabidopsis.tsv"), sep="\t")
at2gohir = defaultdict(set)
for _, r in cw.iterrows(): at2gohir[r["arabidopsis_gene"]].add(str(r["gohir_gene_upper"]).upper())
marker_gohir = {}  # gohir(upper) -> (celltype, symbol)
for at, (ct, sym) in MARKERS.items():
    for gh in at2gohir.get(at, []): marker_gohir[gh] = (ct, sym)

pa = pd.read_csv(os.path.join(R, "program_assignments.csv"))
rows = []
for p in sorted(pa.program.unique()):
    hits = defaultdict(list)
    for g in pa.gene[pa.program == p]:
        mk = marker_gohir.get(g.upper())
        if mk: hits[mk[0]].append(mk[1])
    if hits:
        summ = "; ".join(f"{ct} ({','.join(sorted(set(s)))})" for ct, s in sorted(hits.items(), key=lambda x: -len(x[1])))
    else:
        summ = "(no panel markers)"
    rows.append({"program": p, "n_marker_hits": sum(len(v) for v in hits.values()),
                 "top_celltype": max(hits, key=lambda k: len(hits[k])) if hits else "", "cell_type_markers": summ})
    print(f"  P{p:>2}: {summ}")
pd.DataFrame(rows).to_csv(os.path.join(R, "program_celltypes.csv"), index=False)
print("wrote program_celltypes.csv (curated marker panel; bulk association, not single-cell)")
