"""VOIDS DR2 mass-matched void-wall Delta_HR (blinded; blind cancels in differences).
Reads data/dr2_hr.csv. (1) mass+z matched void-wall comparison; (2) joint
HR ~ a + b*(logM-10) + D*I_void fit. All outputs blinded (absolute scale hidden).
Outputs data/dr2_matched.csv + console numbers -> logs/dr2_matched.log (persisted).
Run: venvs/b3/bin/python src/dr2_matched.py
"""
import numpy as np
import pandas as pd

q = pd.read_csv("data/dr2_hr.csv")
qv = q[q.env == "void"].copy()
qw = q[q.env == "wall"].copy()
diffs = []
for _, r in qv.iterrows():
    if not np.isfinite(r.mass):
        continue
    pool = qw[((qw.mass - r.mass).abs() < 0.25) & ((qw.redshift - r.redshift).abs() < 0.015)]
    if len(pool) == 0:
        continue
    diffs.append((r.hr_blind, pool.hr_blind.mean(), len(pool)))
d = pd.DataFrame(diffs, columns=["v", "wmean", "k"])
print("matched voids:", len(d), "/", len(qv))
delta = (d.v - d.wmean).mean()
se = (d.v - d.wmean).std(ddof=1) / np.sqrt(len(d))
print(f"MATCHED Delta_HR BLINDED = {delta:+.4f} se={se:.4f} (mass+z matched, stat only)")
d.to_csv("data/dr2_matched.csv", index=False)

A = np.column_stack([np.ones(len(q)), q.mass.fillna(q.mass.median()) - 10.0,
                     (q.env == "void").astype(float)])
y = q.hr_blind.values
w = 1 / q.sig2.values
import numpy.linalg as la
beta, *_ = la.lstsq(A * np.sqrt(w)[:, None], y * np.sqrt(w), rcond=None)
pred = A @ beta
s2 = ((y - pred) ** 2 * w).sum() / (len(y) - 3)
cov = np.linalg.inv((A * w[:, None]).T @ A) * s2
print(f"JOINT mass-step slope={beta[1]:+.4f}+/-{np.sqrt(cov[1,1]):.4f} "
      f"void-offset={beta[2]:+.4f}+/-{np.sqrt(cov[2,2]):.4f} (blinded)")
print("wrote data/dr2_matched.csv (BLINDED differences only; absolute scale sealed)")
