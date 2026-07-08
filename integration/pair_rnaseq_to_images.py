#!/usr/bin/env python3
"""
Pair RNA-seq libraries to imaged plants (RSML root traits), IF a library->plant manifest exists.

Background: the committed count matrix labels the 4 replicate libraries per group only by
Genotype x Tissue x Treatment (no plant ID), while the RSML tracings are keyed per plant
(genotype_plateWell_condition, e.g. A68_a12_FL). Pairing therefore needs ONE extra file — a manifest
mapping each RNA-seq library to the plant (plate+well) it came from. See PAIRING_WHATS_NEEDED.md.

This script is ready to run the moment that manifest is supplied. Until then it validates inputs and
tells you exactly what is missing.

--------------------------------------------------------------------------------------------------
MANIFEST FORMAT (CSV) — one row per RNA-seq library (expected: 48 rows; 24 root + 24 shoot):
    library_id,genotype,tissue,treatment,plate_well
    lib01,A68,Root,Flight,a12          # plate_well must match the RSML naming (a12, b2, c11, ...)
    lib02,A68,Root,Flight,b2
    ...
  - genotype in {WT,A68,D130}; tissue in {Root,Shoot}; treatment in {Flight,Ground}
  - plate_well = the lower-case plate+well token used in RSML filenames (see rsml_traits.csv `plant`)
  - If replicates were POOLED across plants, per-plant pairing is impossible — record the pool
    membership instead and stop here (only group-level integration is valid).

USAGE:
  python pair_rnaseq_to_images.py --manifest library_to_plant.csv \
      [--traits ../morphometrics/rsml_traits.csv] [--day 6] [--out results/paired_rnaseq_rsml.csv]

OUTPUT (results/paired_rnaseq_rsml.csv): one row per library, with its plant's RSML traits at --day
joined on {genotype, plate_well, condition}. This is the input for an INDIVIDUAL-level integrated
model (n up to 24 root plants) that breaks the group-level Flight/growth collinearity.
No third-party deps (stdlib).
"""
import argparse, csv, os, sys

REQUIRED_COLS = ["library_id", "genotype", "tissue", "treatment", "plate_well"]
COND = {"Flight": "FL", "Ground": "GC"}


def load_traits(path, day):
    """plant-day RSML traits -> lookup keyed (genotype, plate_well, FL/GC) at the chosen day."""
    lut = {}
    for r in csv.DictReader(open(path, encoding="utf-8")):
        if not r.get("genotype") or str(r.get("day")) != str(day):
            continue
        lut[(r["genotype"], r["plant"], r["condition"])] = r
    return lut


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=False,
                    help="CSV mapping each RNA-seq library to genotype,tissue,treatment,plate_well")
    ap.add_argument("--traits", default=os.path.join(here, "..", "morphometrics", "rsml_traits.csv"))
    ap.add_argument("--day", default="6", help="RSML day to join on (root imaging days 3-6)")
    ap.add_argument("--out", default=os.path.join(here, "results", "paired_rnaseq_rsml.csv"))
    a = ap.parse_args()

    if not a.manifest:
        sys.exit("No --manifest supplied. Pairing needs a library->plant map; see "
                 "PAIRING_WHATS_NEEDED.md for what to ask collaborators for. Nothing written.")
    if not os.path.exists(a.manifest):
        sys.exit(f"Manifest not found: {a.manifest}")

    man = list(csv.DictReader(open(a.manifest, encoding="utf-8")))
    missing_cols = [c for c in REQUIRED_COLS if not man or c not in man[0]]
    if missing_cols:
        sys.exit(f"Manifest missing required column(s): {missing_cols}. Expected: {REQUIRED_COLS}")

    traits = load_traits(a.traits, a.day)
    if not traits:
        sys.exit(f"No RSML traits for day {a.day} in {a.traits}")

    trait_cols = ["n_roots", "n_lateral", "total_length_native", "primary_length_native", "mean_diameter_px"]
    rows, paired, unmatched = [], 0, []
    for m in man:
        if m["tissue"].strip().capitalize() != "Root":
            continue  # RSML is root-only; shoot libraries have no image to pair
        cond = COND.get(m["treatment"].strip().capitalize())
        key = (m["genotype"].strip(), m["plate_well"].strip(), cond)
        t = traits.get(key)
        rec = {"library_id": m["library_id"], "genotype": m["genotype"], "treatment": m["treatment"],
               "plate_well": m["plate_well"], "matched": bool(t)}
        for c in trait_cols:
            rec[c] = t.get(c) if t else ""
        rows.append(rec)
        if t: paired += 1
        else: unmatched.append(f'{m["library_id"]} -> {key}')

    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    with open(a.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["library_id", "genotype", "treatment", "plate_well",
                                           "matched"] + trait_cols)
        w.writeheader(); w.writerows(rows)

    print(f"Root libraries in manifest: {len(rows)} | paired to an imaged plant: {paired} | "
          f"unmatched: {len(rows) - paired}")
    if unmatched:
        print("Unmatched (check plate_well spelling / that the plant was imaged & survived QC):")
        for u in unmatched[:20]:
            print("  ", u)
    print(f"Wrote {a.out}")
    if paired >= 8:
        print("\nNext: fit an individual-level model on this table (n = paired plants) — e.g. within-group "
              "cor(module eigengene, total_length), or mediation Treatment -> expression -> trait — which "
              "separates growth-driving from flight-responding effects.")


if __name__ == "__main__":
    main()
