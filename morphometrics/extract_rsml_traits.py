#!/usr/bin/env python3
"""
Extract root-architecture traits from TICTOC SmartRoot RSML tracings.

Reads Data/Final_RSML_format/*.rsml (curated set) and writes one tidy row per
plant-day with per-plant morphometric traits, plus the parsed design factors
(genotype, condition, day) from the filename grammar {Geno}_{plate}{well}_{FL|GC}_{Day}.

Trait sources (see DATA_DICTIONARY.md sec 7):
  - length_native : SmartRoot's own <properties><length> per root (RSML metadata unit; see caveat)
  - px_length     : geometric polyline length in image pixels (independent check)
  - diameter_px   : mean of the per-node <function name='diameter'> samples (pixels)
  - root order    : top-level <root> under <plant> = primary (order 1); nested = lateral (order >=2)

UNIT CAVEAT: RSML <metadata><resolution> (px per unit) and <unit> are recorded per file so lengths
can be calibrated to real units later. SmartRoot's native <length> and the geometric px_length are
reported side by side; if they disagree, calibrate before publishing absolute lengths.

Usage:
  python extract_rsml_traits.py [--in ../Data/Final_RSML_format] [--out rsml_traits.csv]
No third-party deps (stdlib xml + csv).
"""
import argparse, csv, glob, os, re, sys, xml.etree.ElementTree as ET
from math import hypot

SAMPLE_RE = re.compile(r'^([A-Za-z0-9]+)_([a-z]+\d+)_(FL|GC)_(\d+)$')


def local(tag):  # strip XML namespace
    return tag.rsplit('}', 1)[-1]


def polyline_px_length(root_el):
    pts = []
    for geom in root_el:
        if local(geom.tag) != 'geometry':
            continue
        for poly in geom:
            if local(poly.tag) != 'polyline':
                continue
            for p in poly:
                try:
                    pts.append((float(p.get('x')), float(p.get('y'))))
                except (TypeError, ValueError):
                    pass
    return sum(hypot(pts[i+1][0]-pts[i][0], pts[i+1][1]-pts[i][1]) for i in range(len(pts)-1)), len(pts)


def prop_length(root_el):
    for el in root_el:
        if local(el.tag) == 'properties':
            for pr in el:
                if local(pr.tag) == 'length':
                    try:
                        return float(pr.text)
                    except (TypeError, ValueError):
                        return None
    return None


def diameter_mean(root_el):
    vals = []
    for el in root_el:
        if local(el.tag) != 'functions':
            continue
        for fn in el:
            if fn.get('name') == 'diameter':
                for s in fn:
                    try:
                        vals.append(float(s.text))
                    except (TypeError, ValueError):
                        pass
    return sum(vals)/len(vals) if vals else None


def walk_roots(el, order, out):
    """Recurse: collect (order, root_element) for every <root>, nesting = lateral order."""
    for child in el:
        if local(child.tag) == 'root':
            out.append((order, child))
            walk_roots(child, order + 1, out)      # nested roots are laterals
        elif local(child.tag) in ('scene', 'plant'):
            walk_roots(child, order, out)


def parse_file(path):
    try:
        tree = ET.parse(path)
    except ET.ParseError as e:
        sys.stderr.write(f"parse error {path}: {e}\n")
        return None
    root = tree.getroot()
    # metadata
    resolution = unit = None
    for md in root:
        if local(md.tag) == 'metadata':
            for m in md:
                if local(m.tag) == 'resolution':
                    try: resolution = float(m.text)
                    except (TypeError, ValueError): pass
                elif local(m.tag) == 'unit':
                    unit = (m.text or '').strip()
    roots = []
    walk_roots(root, 1, roots)
    if not roots:
        return None
    total_native = sum(v for _, r in roots if (v := prop_length(r)) is not None)
    total_px = 0.0
    diam_vals = []
    primary_native = 0.0
    n_primary = n_lateral = 0
    for order, r in roots:
        pxl, _ = polyline_px_length(r)
        total_px += pxl
        d = diameter_mean(r)
        if d is not None:
            diam_vals.append(d)
        if order == 1:
            n_primary += 1
            pl = prop_length(r)
            if pl is not None:
                primary_native += pl
        else:
            n_lateral += 1
    base = os.path.splitext(os.path.basename(path))[0]
    m = SAMPLE_RE.match(base)
    geno, plant, cond, day = (m.group(1), m.group(2), m.group(3), int(m.group(4))) if m else ('', '', '', '')
    return dict(
        sample=base, genotype=geno, plant=plant, condition=cond, day=day,
        n_roots=len(roots), n_primary=n_primary, n_lateral=n_lateral,
        total_length_native=round(total_native, 4), primary_length_native=round(primary_native, 4),
        total_length_px=round(total_px, 2),
        mean_diameter_px=round(sum(diam_vals)/len(diam_vals), 3) if diam_vals else '',
        resolution_px_per_unit=resolution, unit=unit,
    )


def main():
    ap = argparse.ArgumentParser()
    here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument('--in', dest='indir', default=os.path.join(here, '..', 'Data', 'Final_RSML_format'))
    ap.add_argument('--out', default=os.path.join(here, 'rsml_traits.csv'))
    a = ap.parse_args()
    files = sorted(glob.glob(os.path.join(a.indir, '*.rsml')))
    if not files:
        sys.exit(f"No .rsml files in {a.indir}")
    rows = [r for r in (parse_file(f) for f in files) if r]
    cols = ['sample', 'genotype', 'plant', 'condition', 'day', 'n_roots', 'n_primary', 'n_lateral',
            'total_length_native', 'primary_length_native', 'total_length_px', 'mean_diameter_px',
            'resolution_px_per_unit', 'unit']
    with open(a.out, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    print(f"Wrote {len(rows)} plant-day rows -> {a.out}")


if __name__ == '__main__':
    main()
