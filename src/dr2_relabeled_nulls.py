"""VOIDS relabeled nulls + robustness (blinded): DR2-z labels as HR primary.
Reads data/dr2_hr_relabeled.csv (env_drz). (1) RA-shuffle null; (2) block bootstrap;
(3) injections 0.00/0.05; (4) leave-one-out leverage; (5) mass/z window grid.
All blinded differences. No unblinding.
Outputs logs/relabeled_nulls.log (persisted).
Run: venvs/b3/bin/python src/dr2_relabeled_nulls.py
"""
import numpy as np
import pandas as pd

SEED = 20261002
rng = np.random.default_rng(SEED)
q = pd.read_csv("data/dr2_hr_relabeled.csv")
qv = q[q.env_drz == "void"].copy()
qw = q[q.env_drz == "wall"].copy()
print(f"relabeled primary: void {len(qv)} wall {len(qw)}")


def matched(vv, ww, mw=0.25, zw=0.015):
    diffs = []
    for _, r in vv.iterrows():
        if not np.isfinite(r.mass):
            continue
        pool = ww[((ww.mass - r.mass).abs() < mw) & ((ww.redshift - r.redshift).abs() < zw)]
        if len(pool) == 0:
            continue
        diffs.append(r.hr_blind - pool.hr_blind.mean())
    diffs = np.array(diffs)
    if len(diffs) < 3:
        return np.nan, np.nan, len(diffs)
    return diffs.mean(), diffs.std(ddof=1) / np.sqrt(len(diffs)), len(diffs)


d0, se0, n0 = matched(qv, qw)
print(f"primary matched Delta BLINDED = {d0:+.4f} se={se0:.4f} (n={n0})")
# (1) RA-shuffle
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
# (2) block bootstrap
boots = []
for b in range(500):
    idx = rng.integers(0, len(qv), len(qv))
    d, _, _ = matched(qv.iloc[idx], qw)
    boots.append(d)
boots = np.array([x for x in boots if np.isfinite(x)])
print(f"bootstrap: mean={boots.mean():+.4f} sd={boots.std(ddof=1):.4f} (vs analytic {se0:.4f})")
# (3) injections
for inj in [0.00, 0.05, -0.05]:
    rec = []
    for b in range(200):
        idx = rng.integers(0, len(qv), len(qv))
        s = qv.iloc[idx].copy()
        s["hr_blind"] = s.hr_blind + inj
        d, _, _ = matched(s, qw)
        rec.append(d)
    rec = np.array([x for x in rec if np.isfinite(x)])
    print(f"inject {inj:+.2f}: recovered={rec.mean():+.4f} (expect ~{d0 + inj:+.4f})")
# (4) LOO leverage
base, _, _ = matched(qv, qw)
for k in range(len(qv)):
    d, _, _ = matched(qv.drop(qv.index[k]), qw)
    print(f"LOO drop {qv.iloc[k].ztfname}: Delta={d:+.4f} shift={d - base:+.4f}")
# (5) window grid
for mw, zw in [(0.15, 0.01), (0.25, 0.015), (0.35, 0.02)]:
    d, se, n = matched(qv, qw, mw, zw)
    print(f"window m{mw}/z{zw}: Delta={d:+.4f}+/-{se:.4f} (n={n})")
print("relabeled nulls done (BLINDED)")
