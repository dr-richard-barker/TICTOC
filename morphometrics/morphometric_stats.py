#!/usr/bin/env python3
"""
Mixed-model analysis of RSML root-length traits for the unbalanced TICTOC time series.

Model (headline trait = total root length, log-scale for multiplicative growth):
    log(total_length_native) ~ condition * genotype + day + day:condition   [+ random intercept per plant]
  - condition: Flight vs Ground (ref = Ground)
  - genotype : WT / A68 / D130   (ref = WT)
  - day      : 3-6 (continuous)  -> growth rate on log scale
  - random intercept per plant (plant measured repeatedly across days = repeated measures)

Reports the fixed-effect table; the condition:genotype terms test whether the AVP-OX lines' spaceflight
response differs from WT. Writes a coefficients CSV + a text summary.

*** CAVEAT ***: this quantifies the pattern IN THE TRACED DATA. It does NOT resolve the Flight-vs-Ground
image-calibration confound (see README) — a significant `condition` effect could still be an imaging
artifact. Confirm calibration before any biological interpretation.

Usage: python morphometric_stats.py [--in rsml_traits.csv] [--trait total_length_native]
Deps: pandas, numpy, statsmodels.
"""
import argparse, os, sys
import numpy as np, pandas as pd
import statsmodels.formula.api as smf


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser()
    ap.add_argument('--in', dest='infile', default=os.path.join(here, 'rsml_traits.csv'))
    ap.add_argument('--trait', default='total_length_native')
    ap.add_argument('--out-prefix', default=os.path.join(here, 'morphometric_stats'))
    a = ap.parse_args()

    df = pd.read_csv(a.infile)
    df = df[df['genotype'].astype(bool)].copy()
    df['plant_uid'] = df['genotype'] + '_' + df['plant'] + '_' + df['condition']   # unique plant
    df['day'] = df['day'].astype(float)
    df[a.trait] = pd.to_numeric(df[a.trait], errors='coerce')
    df = df.dropna(subset=[a.trait])
    df = df[df[a.trait] > 0]
    df['logtrait'] = np.log(df[a.trait])
    # order factor levels so Ground / WT are the reference
    df['condition'] = pd.Categorical(df['condition'], categories=['GC', 'FL'])
    df['genotype'] = pd.Categorical(df['genotype'], categories=['WT', 'A68', 'D130'])

    n_plants = df['plant_uid'].nunique()
    print(f"{len(df)} obs, {n_plants} plants, trait={a.trait}")
    print(df.groupby(['genotype', 'condition'], observed=True)['plant_uid'].nunique().rename('n_plants'))

    formula = "logtrait ~ C(condition) * C(genotype) + day + day:C(condition)"
    try:
        model = smf.mixedlm(formula, df, groups=df['plant_uid'])
        fit = model.fit(reml=True, method='lbfgs')
    except Exception as e:                                   # fall back to random intercept only / simpler
        sys.stderr.write(f"mixedlm failed ({e}); retrying simpler model\n")
        formula = "logtrait ~ C(condition) * C(genotype) + day"
        fit = smf.mixedlm(formula, df, groups=df['plant_uid']).fit(reml=True, method='lbfgs')

    print("\n=== Mixed model (REML) ===")
    print("formula:", formula)
    print(fit.summary())

    # tidy coefficients
    coef = pd.DataFrame({
        'term': fit.params.index, 'estimate_log': fit.params.values,
        'std_err': fit.bse.reindex(fit.params.index).values,
        'z': fit.tvalues.reindex(fit.params.index).values,
        'p_value': fit.pvalues.reindex(fit.params.index).values,
    })
    coef['fold_change'] = np.where(coef['term'].str.contains('Group Var'), np.nan, np.exp(coef['estimate_log']))
    coef.to_csv(f"{a.out_prefix}_coefficients.csv", index=False)

    with open(f"{a.out_prefix}_summary.txt", 'w') as fh:
        fh.write(f"TICTOC morphometric mixed model\ntrait: {a.trait} (log scale)\nformula: {formula}\n")
        fh.write(f"n_obs={len(df)}, n_plants={n_plants}\n\n")
        fh.write(str(fit.summary()))
        fh.write("\n\nCAVEAT: does not resolve the Flight/Ground image-calibration confound (see README).\n")

    # interpret the key terms
    print("\n=== Key effects (exp(coef) = multiplicative effect on total length) ===")
    for t in coef['term']:
        if 'condition' in t.lower():
            row = coef[coef['term'] == t].iloc[0]
            print(f"  {t:45} x{row['fold_change']:.2f}  p={row['p_value']:.2g}")
    print(f"\nWrote {a.out_prefix}_coefficients.csv and _summary.txt")


if __name__ == '__main__':
    main()
