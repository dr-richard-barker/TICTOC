#!/usr/bin/env python3
"""
Summarise rsml_traits.csv into group means (genotype x condition x day) for Fig 2 / Results R2.
Descriptive only — formal stats (mixed model on the unbalanced design) are the next step (roadmap 4.x).
Usage: python summarise_traits.py [--in rsml_traits.csv] [--out rsml_traits_summary.csv]
"""
import argparse, csv, os, statistics as st
from collections import defaultdict

TRAITS = ['total_length_native', 'primary_length_native', 'n_roots', 'n_lateral', 'mean_diameter_px']


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser()
    ap.add_argument('--in', dest='infile', default=os.path.join(here, 'rsml_traits.csv'))
    ap.add_argument('--out', default=os.path.join(here, 'rsml_traits_summary.csv'))
    a = ap.parse_args()
    rows = list(csv.DictReader(open(a.infile)))
    groups = defaultdict(list)
    for r in rows:
        groups[(r['genotype'], r['condition'], int(r['day']))].append(r)

    def fnum(r, k):
        try: return float(r[k])
        except (ValueError, KeyError): return None

    out = []
    for (g, c, d) in sorted(groups):
        grp = groups[(g, c, d)]
        rec = dict(genotype=g, condition=c, day=d, n_plants=len(grp))
        for t in TRAITS:
            vals = [v for r in grp if (v := fnum(r, t)) is not None]
            rec[f'{t}_mean'] = round(st.mean(vals), 3) if vals else ''
            rec[f'{t}_sd'] = round(st.pstdev(vals), 3) if len(vals) > 1 else ''
        out.append(rec)

    cols = ['genotype', 'condition', 'day', 'n_plants'] + [f'{t}_{s}' for t in TRAITS for s in ('mean', 'sd')]
    with open(a.out, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader(); w.writerows(out)
    print(f"Wrote {len(out)} group rows -> {a.out}")

    # console: day-6 total-length Flight vs Ground per genotype (descriptive headline)
    print("\nDay-6 total_length_native (mean, n) — Flight vs Ground:")
    for g in ['WT', 'A68', 'D130']:
        cell = {c: next((r for r in out if r['genotype'] == g and r['condition'] == c and r['day'] == 6), None)
                for c in ('FL', 'GC')}
        fl, gc = cell['FL'], cell['GC']
        if fl and gc and fl['total_length_native_mean'] != '' and gc['total_length_native_mean'] != '':
            print(f"  {g:5} FL {fl['total_length_native_mean']:>8} (n={fl['n_plants']})   "
                  f"GC {gc['total_length_native_mean']:>8} (n={gc['n_plants']})   "
                  f"FL/GC={fl['total_length_native_mean']/gc['total_length_native_mean']:.2f}")


if __name__ == '__main__':
    main()
