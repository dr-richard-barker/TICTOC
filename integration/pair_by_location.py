#!/usr/bin/env python3
"""
Location-level pairing: join RNA-seq libraries to ROW-aggregated root traits (technical reps per row).

Complements pair_rnaseq_to_images.py (which does exact-plant, 1:1). This one needs only a coarse
library->row map — either an explicit manifest, or the confirmation that the 4 replicates per
genotype x tissue x treatment correspond to rows a,b,c,d (in a given order).

MANIFEST (CSV), one row per RNA-seq library:
    library_id,genotype,tissue,treatment,row
    lib01,A68,Root,Flight,a       # row = a/b/c/d (well row-letter used in the RSML plant IDs)
    lib02,A68,Root,Flight,b
    ...
  If you instead only know "the 4 reps are rows a,b,c,d in order", generate this manifest with
  --rows-in-order (assigns a,b,c,d to replicate order within each genotype x tissue x treatment group,
  reading the count-matrix column order) — clearly an ASSUMPTION to confirm before trusting results.

USAGE:
  python pair_by_location.py --manifest lib_to_row.csv
  python pair_by_location.py --rows-in-order            # build assumed manifest from the count matrix
Output: results/paired_by_location.csv (each library + its row's mean root traits, root libraries only).
Deps: pandas.
"""
import os, re, argparse, sys, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.abspath(os.path.join(HERE, ".."))
COND = {"Flight": "FL", "Ground": "GC"}

def rows_in_order():
    df = pd.read_csv(os.path.join(REPO, "TICTOC_run1_filteredCounts_v3.csv"), nrows=0)
    cols = list(df.columns[1:])
    def parse(g):
        g = re.sub(r'\.\d+$', '', g)                 # strip pandas dup-column suffix (.1/.2/.3)
        return (re.sub(r'(Root|Shoot).*', '', g), re.sub(r'.*?(Root|Shoot).*', r'\1', g), re.sub(r'.*(Flight|Ground)$', r'\1', g))
    rows, seen = [], {}
    for i, c in enumerate(cols):
        ge, ti, tr = parse(c); k = (ge, ti, tr); idx = seen.get(k, 0); seen[k] = idx + 1
        rows.append({"library_id": f"{c}_col{i+1}", "genotype": ge, "tissue": ti, "treatment": tr,
                     "row": "abcd"[idx] if idx < 4 else f"r{idx}"})
    return pd.DataFrame(rows)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest"); ap.add_argument("--rows-in-order", action="store_true")
    ap.add_argument("--loc", default=os.path.join(HERE, "results", "location_mean_traits.csv"))
    ap.add_argument("--out", default=os.path.join(HERE, "results", "paired_by_location.csv"))
    a = ap.parse_args()
    if not os.path.exists(a.loc):
        sys.exit("Run location_traits.py first to produce location_mean_traits.csv")
    if a.manifest:
        man = pd.read_csv(a.manifest)
    elif a.rows_in_order:
        man = rows_in_order()
        print("ASSUMED manifest (reps -> rows a,b,c,d in column order) — CONFIRM before trusting results.")
    else:
        sys.exit("Provide --manifest lib_to_row.csv or --rows-in-order (see PAIRING_WHATS_NEEDED.md).")
    loc = pd.read_csv(a.loc)
    man = man[man.tissue.str.capitalize() == "Root"].copy()          # RSML is root-only
    man["condition"] = man.treatment.str.capitalize().map(COND)
    m = man.merge(loc, on=["genotype", "condition", "row"], how="left", indicator=True)
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    m.to_csv(a.out, index=False)
    ok = (m._merge == "both").sum()
    print(f"Root libraries: {len(m)} | joined to a row's mean traits: {ok} | unmatched: {len(m)-ok}")
    if ok < len(m):
        print("Unmatched (row not imaged / spelling):")
        print(m[m._merge != "both"][["library_id", "genotype", "condition", "row"]].to_string(index=False))
    print(f"Wrote {a.out}. Fit an individual/location-level model on this (n up to ~{loc.shape[0]}) "
          "to sharpen the transcriptome<->architecture link beyond the n=6 group means.")

if __name__ == "__main__":
    main()
