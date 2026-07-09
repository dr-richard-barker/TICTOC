#!/usr/bin/env python3
"""
Location-level root-trait aggregation (the imaging half of location-based pairing).

Rationale: RNA-seq libraries can't be matched 1:1 to imaged plants, but the experiment has spatial
LOCATIONS (well = row-letter + number). At ROW level (a/b/c/d) most genotype x condition cells hold
2-4 imaged plants, which we average as technical replicates. If the 4 RNA-seq replicates per
genotype x tissue x treatment correspond to the 4 rows (a-d) - a common block design - the resulting
per-(genotype x row x condition) trait means can be joined to the libraries at ROW resolution (up to
~23 cells vs 6 group means), sharpening the transcriptome<->architecture integration.

This script produces the location-mean trait table from the committed RSML traits. The join itself
needs a library->row (or replicate-order) confirmation from the collaborators (see PAIRING_WHATS_NEEDED.md).

Usage: python integration/location_traits.py [--day 6]
Output: integration/results/location_mean_traits.csv   (genotype, row, condition, n_plants, trait means/SD)
Deps: pandas.
"""
import os, re, argparse, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.abspath(os.path.join(HERE, ".."))
TRAITS = ["total_length_native", "primary_length_native", "n_lateral", "n_roots", "mean_diameter_px"]

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--day", default="6")
    ap.add_argument("--out", default=os.path.join(HERE, "results", "location_mean_traits.csv"))
    a = ap.parse_args()
    df = pd.read_csv(os.path.join(REPO, "morphometrics", "rsml_traits.csv"))
    df = df[df.genotype.astype(bool) & (df.day.astype(str) == str(a.day))].copy()
    df["row"] = df["plant"].map(lambda w: re.match(r"([a-zA-Z]+)", str(w)).group(1))
    for t in TRAITS: df[t] = pd.to_numeric(df[t], errors="coerce")
    g = df.groupby(["genotype", "row", "condition"])
    out = g[TRAITS].mean().round(3)
    out.columns = [f"{c}_mean" for c in out.columns]
    out["n_plants"] = g.size()
    for t in TRAITS: out[f"{t}_sd"] = g[t].std().round(3)
    out = out.reset_index()
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    out.to_csv(a.out, index=False)
    reps = (out.n_plants > 1).sum()
    print(f"Wrote {len(out)} (genotype x row x condition) location cells -> {a.out}")
    print(f"  cells with >1 plant (true technical reps): {reps}/{len(out)}")
    print(f"  vs 6 genotype x condition group means — ~{len(out)/6:.1f}x finer for the imaging side.")
    print("\nExample (A68):")
    print(out[out.genotype == "A68"][["row", "condition", "n_plants", "total_length_native_mean", "primary_length_native_mean"]].to_string(index=False))

if __name__ == "__main__":
    main()
