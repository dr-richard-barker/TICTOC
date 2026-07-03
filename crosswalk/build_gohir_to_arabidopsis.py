#!/usr/bin/env python3
"""
Build a Gohir (Gossypium hirsutum, Phytozome UTX-TM1 v2.1) -> Arabidopsis thaliana
gene crosswalk from the CottonGen BLASTP best-hit annotation.

WHY THIS EXISTS
---------------
The TICTOC count matrices and DEG tables use Phytozome `Gohir.*` gene IDs (G. hirsutum,
allotetraploid AD1). The repo's older `Ara_vs_Cotton_biomart_export` maps Arabidopsis to
`B456_`/`Gorai` IDs (G. RAIMONDII, the D-genome progenitor) because Ensembl Plants never
carried G. hirsutum -- so it does NOT join to the expression data. This script builds the
correct `Gohir -> AT` bridge needed for GO annotation transfer and PhysioSpace decoding.

SOURCE (authoritative, public, no login):
  CottonGen, G. hirsutum UTX-TM1 v2.1, homology/
  blastp_G.hirsutum_UTX_v2.1_vs_arabidopsis.xlsx.gz
  (BLASTP of every Gohir protein vs TAIR Arabidopsis proteins.)

METHOD
------
1. Download + gunzip the xlsx (4 "Page" worksheets, ~97.6k hit rows) if not present.
2. Collapse transcript-level hits to GENE level (strip `.<isoform>` and the `_UTX-TM1...` suffix).
3. Keep ONE best hit per Gohir gene: min E-value, then max bit-score, then max % identity.
4. Write a tidy TSV + print coverage against the filtered count matrix.

CAVEATS
-------
- A BLASTP best hit is a *proxy* for orthology, not a curated ortholog call. Filter on
  `pid`/`evalue` for stricter sets. Many Gohir genes map to the same AT gene (A/D homoeologs
  + paralog expansion) -- expected for an allotetraploid.
- Join to `Gohir.*` count-matrix IDs is exact; DEG tables use UPPER-CASE `GOHIR.*`, so join
  case-insensitively (an `gohir_gene_upper` column is provided for convenience).

Usage:  python build_gohir_to_arabidopsis.py [--xlsx path/to/file.xlsx] [--counts ../TICTOC_run1_filteredCounts_v3.csv]
Requires: pandas, openpyxl (+ requests OR a pre-downloaded xlsx).
"""
import argparse, gzip, os, re, sys, urllib.request

URL = ("https://www.cottongen.org/cottongen_downloads/Gossypium_hirsutum/"
       "UTX-TM1_v2.1/homology/blastp_G.hirsutum_UTX_v2.1_vs_arabidopsis.xlsx.gz")
UA = "Mozilla/5.0 (compatible; TICTOC-crosswalk-builder)"
HERE = os.path.dirname(os.path.abspath(__file__))

# NB: source cells can carry trailing whitespace -> strip() BEFORE removing the .isoform suffix.
gene_from_query = lambda q: re.sub(r"\.\d+$", "", str(q).strip().split("_")[0])  # Gohir.A03G108500.3_UTX.. -> Gohir.A03G108500
gene_from_match = lambda m: re.sub(r"\.\d+$", "", str(m).strip())                # "AT5G65110.1 " -> AT5G65110


def fetch(xlsx_path):
    gz = xlsx_path + ".gz"
    if os.path.exists(xlsx_path):
        return xlsx_path
    print(f"Downloading {URL}", file=sys.stderr)
    req = urllib.request.Request(URL, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=300) as r, open(gz, "wb") as f:
        f.write(r.read())
    with gzip.open(gz, "rb") as g, open(xlsx_path, "wb") as f:
        f.write(g.read())
    return xlsx_path


def main():
    import pandas as pd
    ap = argparse.ArgumentParser()
    ap.add_argument("--xlsx", default=os.path.join(HERE, "blastp_Gohir_vs_arabidopsis.xlsx"))
    ap.add_argument("--counts", default=os.path.join(HERE, "..", "TICTOC_run1_filteredCounts_v3.csv"))
    ap.add_argument("--out", default=os.path.join(HERE, "gohir_to_arabidopsis.tsv"))
    a = ap.parse_args()

    xlsx = fetch(a.xlsx)
    xl = pd.ExcelFile(xlsx)
    parts = []
    for sh in xl.sheet_names:                    # header note on row 0, blank row 1, labels row 2
        d = xl.parse(sh, header=2)
        parts.append(d[[c for c in d.columns if not str(c).startswith("Unnamed")]])
    df = pd.concat(parts, ignore_index=True).dropna(subset=["Query", "Match"])
    print(f"BLASTP hit rows (all pages): {len(df)}", file=sys.stderr)

    df["gohir_gene"] = df["Query"].map(gene_from_query)
    df["arabidopsis_gene"] = df["Match"].map(gene_from_match)
    for c in ["Exp", "Score", "PID", "Align Length"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    df = df.sort_values(["gohir_gene", "Exp", "Score", "PID"], ascending=[True, True, False, False])
    best = df.groupby("gohir_gene", as_index=False).first()

    out = pd.DataFrame({
        "gohir_gene": best["gohir_gene"],
        "gohir_gene_upper": best["gohir_gene"].str.upper(),
        "arabidopsis_gene": best["arabidopsis_gene"],
        "pid": best["PID"].round(2),
        "evalue": best["Exp"],
        "bit_score": best["Score"],
        "align_length": best["Align Length"],
        "gohir_transcript_hit": best["Query"].str.strip().str.split("_").str[0],
        "arabidopsis_transcript_hit": best["Match"].str.strip(),
        "arabidopsis_description": best["Description"].astype(str).str.strip(),
    }).sort_values("gohir_gene")
    out.to_csv(a.out, sep="\t", index=False)
    print(f"Wrote {len(out)} Gohir->AT best-hit rows to {a.out}", file=sys.stderr)
    print(f"Unique AT genes referenced: {out['arabidopsis_gene'].nunique()}", file=sys.stderr)

    # coverage vs count matrix
    try:
        import csv
        with open(a.counts, encoding="utf-8", errors="replace") as fh:
            r = csv.reader(fh); next(r)
            mat = [row[0] for row in r if row]
        matU = {g.upper() for g in mat}
        cov = len(matU & set(out["gohir_gene_upper"]))
        print(f"Count-matrix coverage: {cov}/{len(mat)} ({100*cov/len(mat):.1f}%)", file=sys.stderr)
    except FileNotFoundError:
        print("(counts file not found; skipped coverage)", file=sys.stderr)


if __name__ == "__main__":
    main()
