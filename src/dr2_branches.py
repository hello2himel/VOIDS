"""VOIDS DR2 robustness branches (blinded): mass-step vs linear, z<=0.06 HR slice, zmin grid.
Reads data/dr2_hr.csv. All outputs blinded differences (absolute scale sealed).
Appends human-readable summary to logs/step_zcut.log (persisted).
Run: venvs/b3/bin/python src/dr2_branches.py
"""
import numpy as np
import pandas as pd
import numpy.linalg as la

q = pd.read_csv("data/dr2_hr.csv")
out = []
out.append("=== mass-step branch (step at 10.0) vs linear ===")
for name, sub in [("full", q), ("zmax0.06", q[q.redshift <= 0.06])]:
    y = sub.hr_blind.values
    w = 1 / sub.sig2.values
    A = np.column_stack([np.ones(len(sub)), sub.mass.fillna(sub.mass.median()) - 10.0,
                         (sub.env == "void").astype(float)])
    beta, *_ = la.lstsq(A * np.sqrt(w)[:, None], y * np.sqrt(w), rcond=None)
    pred = A @ beta
    s2 = ((y - pred) ** 2 * w).sum() / (len(y) - 3)
    cov = np.linalg.inv((A * w[:, None]).T @ A) * s2
    B = np.column_stack([np.ones(len(sub)), (sub.mass.fillna(sub.mass.median()) >= 10.0).astype(float),
                         (sub.env == "void").astype(float)])
    g, *_ = la.lstsq(B * np.sqrt(w)[:, None], y * np.sqrt(w), rcond=None)
    pred2 = B @ g
    s2b = ((y - pred2) ** 2 * w).sum() / (len(y) - 3)
    covb = np.linalg.inv((B * w[:, None]).T @ B) * s2b
    nv = int((sub.env == "void").sum())
    out.append(f"{name}: N={len(sub)} Nvoid={nv} linear-void={beta[2]:+.4f}+/-{np.sqrt(cov[2,2]):.4f} "
               f"step-void={g[2]:+.4f}+/-{np.sqrt(covb[2,2]):.4f} step-size={g[1]:+.4f}")
out.append("=== zmin stability (matched) ===")
for zmin in [0.015, 0.02, 0.025, 0.03]:
    s = q[q.redshift >= zmin]
    qv, qw = s[s.env == "void"], s[s.env == "wall"]
    diffs = []
    for _, r in qv.iterrows():
        if not np.isfinite(r.mass):
            continue
        pool = qw[((qw.mass - r.mass).abs() < 0.25) & ((qw.redshift - r.redshift).abs() < 0.015)]
        if len(pool) == 0:
            continue
        diffs.append(r.hr_blind - pool.hr_blind.mean())
    diffs = np.array(diffs)
    if len(diffs) > 2:
        out.append(f"zmin {zmin}: N={len(s)} Nvoid={len(qv)} matched={len(diffs)} "
                   f"Delta={diffs.mean():+.4f}+/-{diffs.std(ddof=1)/np.sqrt(len(diffs)):.4f}")
    else:
        out.append(f"zmin {zmin}: too few ({len(diffs)})")
open("logs/step_zcut.log", "w").write("\n".join(out) + "\n")
print("\n".join(out))
print("wrote logs/step_zcut.log (BLINDED)")
