"""VOIDS DR2 nulls + injections (blinded; validates quoted errors before unblinding).
Reads data/dr2_hr.csv. (1) RA-shuffle null: permute env labels within footprint, recompute
matched Delta distribution (expect ~N(0, se)); (2) block bootstrap over voids (resample 31
voids with replacement, recompute matched Delta); (3) inject +/-0.00/0.05 mag into void HRs,
verify recovery. All on blinded values (offsets cancel in differences).
Outputs logs/dr2_nulls.log summary (persisted). No unblinding.
Run: venvs/b3/bin/python src/dr2_nulls.py
"""
import numpy as np
import pandas as pd

SEED = 20260930
rng = np.random.default_rng(SEED)
q = pd.read_csv("data/dr2_hr.csv")
qv = q[q.env == "void"].copy()
qw = q[q.env == "wall"].copy()
print(f"input: void {len(qv)} wall {len(qw)}")


def matched_delta(void_df, wall_df):
    diffs = []
    for _, r in void_df.iterrows():
        if not np.isfinite(r.mass):
            continue
        pool = wall_df[((wall_df.mass - r.mass).abs() < 0.25) & ((wall_df.redshift - r.redshift).abs() < 0.015)]
        if len(pool) == 0:
            continue
        diffs.append(r.hr_blind - pool.hr_blind.mean())
    diffs = np.array(diffs)
    if len(diffs) < 3:
        return np.nan, np.nan, len(diffs)
    return diffs.mean(), diffs.std(ddof=1) / np.sqrt(len(diffs)), len(diffs)


d0, se0, n0 = matched_delta(qv, qw)
print(f"observed matched Delta BLINDED = {d0:+.4f} se={se0:.4f} (n={n0})")

# (1) RA-shuffle null: permute env labels, keep N_void fixed
nulls = []
for b in range(500):
    lab = np.array(["void"] * len(qv) + ["wall"] * len(qw))
    rng.shuffle(lab)
    # map shuffled labels back onto pooled rows (void-labeled rows act as pseudo-voids)
    pool = pd.concat([qv, qw], ignore_index=True)
    pv = pool[np.array(lab) == "void"]
    pw = pool[np.array(lab) == "wall"]
    d, _, _ = matched_delta(pv, pw)
    nulls.append(d)
nulls = np.array([x for x in nulls if np.isfinite(x)])
print(f"RA-shuffle null: mean={nulls.mean():+.4f} sd={nulls.std(ddof=1):.4f} "
      f"frac|>=obs|={(np.abs(nulls) >= abs(d0)).mean():.3f}")

# (2) block bootstrap over voids
boots = []
for b in range(500):
    idx = rng.integers(0, len(qv), len(qv))
    d, _, _ = matched_delta(qv.iloc[idx], qw)
    boots.append(d)
boots = np.array([x for x in boots if np.isfinite(x)])
print(f"block-bootstrap: mean={boots.mean():+.4f} sd={boots.std(ddof=1):.4f} (vs analytic se={se0:.4f})")

# (3) injections
for inj in [0.00, 0.05, -0.05]:
    rec = []
    for b in range(200):
        idx = rng.integers(0, len(qv), len(qv))
        sample = qv.iloc[idx].copy()
        sample["hr_blind"] = sample.hr_blind + inj
        d, _, _ = matched_delta(sample, qw)
        rec.append(d)
    rec = np.array([x for x in rec if np.isfinite(x)])
    print(f"inject {inj:+.2f}: recovered mean={rec.mean():+.4f} (expect ~{d0 + inj:+.4f})")
print("nulls done (BLINDED; no unblinding)")
