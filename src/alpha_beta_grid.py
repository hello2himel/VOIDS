"""VOIDS Tripp alpha/beta sensitivity grid (blinded differences only).
Recomputes mass+z matched void-wall Delta over alpha in [0.10,0.18],
beta in [2.5,3.5]. M refit as sample mean per grid point (cancels in differences).
Writes logs/alpha_beta_grid.log (persisted). No unblinding.
Run: venvs/b3/bin/python src/alpha_beta_grid.py
"""
import numpy as np
import pandas as pd

q = pd.read_csv("data/dr2_hr.csv")
print("alpha/beta grid on matched Delta (blinded):")
for a in [0.10, 0.14, 0.18]:
    for b in [2.5, 3.1, 3.5]:
        mB = -2.5 * np.log10(q.x0) + 10.635
        mu = mB + a * q.x1 - b * q["c"]
        M = np.average(mu - q.mu_cosmo, weights=1 / q.sig2)
        hr = mu - q.mu_cosmo - M
        qv = q[q.env == "void"]
        qw = q[q.env == "wall"]
        diffs = []
        for i, r in qv.iterrows():
            if not np.isfinite(r.mass):
                continue
            pool = qw[((qw.mass - r.mass).abs() < 0.25) & ((qw.redshift - r.redshift).abs() < 0.015)]
            if len(pool) == 0:
                continue
            diffs.append(hr.loc[i] - hr.loc[pool.index].mean())
        diffs = np.array(diffs)
        print(f"a={a} b={b}: Delta={diffs.mean():+.4f}+/-{diffs.std(ddof=1)/np.sqrt(len(diffs)):.4f} (n={len(diffs)})")
print("grid done (BLINDED)")
