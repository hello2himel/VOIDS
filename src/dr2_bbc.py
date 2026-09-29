"""VOIDS BBC/selection bound (blinded): is matched Delta driven by Malmquist?
Reads data/dr2_hr.csv (has mB via mu_raw? recompute mB = mu_raw - 0.14*x1 + 3.1*c).
1. mB-vs-z slope (Malmquist): fit mB(z) per env; differential slope x dz-gap bounds bias.
2. z-binned matched Delta (split at median z): consistency check.
3. Magnitude-limited reweight: weight wall to void mB distribution, recompute Delta.
All blinded differences. Pass = |shift| < 0.01 mag vs +0.023 baseline.
Outputs logs/dr2_bbc.log (persisted). No unblinding.
Run: venvs/b3/bin/python src/dr2_bbc.py
"""
import numpy as np
import pandas as pd

q = pd.read_csv("data/dr2_hr.csv")
q["mB"] = q.mu_raw - 0.14 * q.x1 + 3.1 * q["c"]
qv, qw = q[q.env == "void"], q[q.env == "wall"]
print(f"N void={len(qv)} wall={len(qw)}")
# (1) Malmquist proxy: HR-vs-z slope (selection imprints z-trends in residuals;
# mB-vs-z is distance-modulus dominated and unusable).
A = np.column_stack([np.ones(len(qw)), qw.redshift])
bh, *_ = np.linalg.lstsq(A, qw.hr_blind, rcond=None)
print(f"wall HR(z) slope={bh[1]:+.3f} mag/z")
dz = qv.redshift.median() - qw.redshift.median()
print(f"median dz void-wall={dz:+.4f}")
print(f"Malmquist-bound bias ~ slope*dz = {bh[1] * dz:+.4f} mag")
# (2) z-binned matched Delta
med = q.redshift.median()
for tag, s in [("lo-z", q[q.redshift <= med]), ("hi-z", q[q.redshift > med])]:
    sv, sw = s[s.env == "void"], s[s.env == "wall"]
    diffs = []
    for _, r in sv.iterrows():
        if not np.isfinite(r.mass):
            continue
        pool = sw[((sw.mass - r.mass).abs() < 0.25) & ((sw.redshift - r.redshift).abs() < 0.015)]
        if len(pool) == 0:
            continue
        diffs.append(r.hr_blind - pool.hr_blind.mean())
    diffs = np.array(diffs)
    if len(diffs) > 2:
        print(f"{tag}: n={len(diffs)} Delta={diffs.mean():+.4f}+/-{diffs.std(ddof=1)/np.sqrt(len(diffs)):.4f}")
    else:
        print(f"{tag}: too few ({len(diffs)})")
# (3) mB-reweighted wall to void mB distribution
bins = np.quantile(q.mB, [0, 0.25, 0.5, 0.75, 1.0])
qv_h, _ = np.histogram(qv.mB, bins=bins)
qw_h, _ = np.histogram(qw.mB, bins=bins)
w = np.ones(len(qw))
for i in range(4):
    m = (qw.mB >= bins[i]) & (qw.mB <= bins[i + 1])
    if m.sum() and qw_h[i]:
        w[m] = (qv_h[i] / max(qv_h.sum(), 1)) / (qw_h[i] / max(qw_h.sum(), 1))
diffs = []
qw_idx = {idx: k for k, idx in enumerate(qw.index)}
for _, r in qv.iterrows():
    if not np.isfinite(r.mass):
        continue
    pool = qw[((qw.mass - r.mass).abs() < 0.25) & ((qw.redshift - r.redshift).abs() < 0.015)]
    if len(pool) == 0:
        continue
    ww = np.array([w[qw_idx[i]] for i in pool.index])
    diffs.append(r.hr_blind - np.average(pool.hr_blind, weights=ww))
diffs = np.array(diffs)
print(f"mB-reweighted matched Delta={diffs.mean():+.4f} (n={len(diffs)}; baseline +0.0234)")
print("BBC bound done (BLINDED)")
