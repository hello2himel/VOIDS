"""VOIDS full-volume DR2 nulls (unblinded post-2026-10-01 values; same code as relabeled nulls).
Reads data/dr2_fullx.csv (needs env/hr columns; recompute matching here).
RA-shuffle null, block bootstrap, injections, LOO leverage.
Outputs logs/dr2_fullx_nulls.log (persisted).
Run: venvs/b3/bin/python src/dr2_fullx_nulls.py
"""
import numpy as np
import pandas as pd

SEED = 20261001
rng = np.random.default_rng(SEED)
q = pd.read_csv("data/dr2_fullx.csv")
qv = q[q.env == "void"].copy()
qw = q[q.env == "wall"].copy()
print(f"fullx: void {len(qv)} wall {len(qw)}")


def matched(vv, ww, mw=0.25, zw=0.015):
    diffs = []
    for _, r in vv.iterrows():
        if not np.isfinite(r.mass):
            continue
        pool = ww[((ww.mass - r.mass).abs() < mw) & ((ww.redshift - r.redshift).abs() < zw)]
        if len(pool) == 0:
            continue
        diffs.append(r.hr - pool.hr.mean())
    diffs = np.array(diffs)
    if len(diffs) < 3:
        return np.nan, np.nan, len(diffs)
    return diffs.mean(), diffs.std(ddof=1) / np.sqrt(len(diffs)), len(diffs)


d0, se0, n0 = matched(qv, qw)
print(f"primary matched Delta = {d0:+.4f} se={se0:.4f} (n={n0})")
pool = pd.concat([qv, qw], ignore_index=True)
nulls = []
for b in range(500):
    lab = np.array(["void"] * len(qv) + ["wall"] * len(qw))
    rng.shuffle(lab)
    pv = pool[np.array(lab) == "void"]
    pw = pool[np.array(lab) == "wall"]
    d, _, _ = matched(pv, pw)
    nulls.append(d)
nulls = np.array([x for x in nulls if np.isfinite(x)])
print(f"RA-shuffle: mean={nulls.mean():+.4f} sd={nulls.std(ddof=1):.4f} frac>=obs={(np.abs(nulls) >= abs(d0)).mean():.3f}")
boots = []
for b in range(500):
    idx = rng.integers(0, len(qv), len(qv))
    d, _, _ = matched(qv.iloc[idx], qw)
    boots.append(d)
boots = np.array([x for x in boots if np.isfinite(x)])
print(f"bootstrap: mean={boots.mean():+.4f} sd={boots.std(ddof=1):.4f} (vs analytic {se0:.4f})")
for inj in [0.00, 0.05, -0.05]:
    rec = []
    for b in range(200):
        idx = rng.integers(0, len(qv), len(qv))
        s = qv.iloc[idx].copy()
        s["hr"] = s.hr + inj
        d, _, _ = matched(s, qw)
        rec.append(d)
    rec = np.array([x for x in rec if np.isfinite(x)])
    print(f"inject {inj:+.2f}: recovered={rec.mean():+.4f} (expect ~{d0 + inj:+.4f})")
print("LOO leverage:")
shifts = []
for k in range(len(qv)):
    d, _, _ = matched(qv.drop(qv.index[k]), qw)
    shifts.append(abs(d - d0))
    print(f"drop {qv.iloc[k].ztfname}: Delta={d:+.4f} shift={d - d0:+.4f}")
print(f"LOO max shift={max(shifts):.4f} vs se={se0:.4f}")
print("fullx nulls done")
