"""VOIDS DR2-z relabel: void membership with precise DR2 redshifts (not coarse BTS z).
BTS redshifts are rounded/coarse (max |dz|=0.026 vs DR2, ~110 Mpc errors). Recompute
r/Rv for the DR2-HR sample in the recomputed Planck18 frame with DR2 redshifts,
save corrected labels, recompute blinded matched Delta + joint fit.
Reads data/dr2_hr.csv. Writes data/dr2_hr_relabeled.csv + logs/relabel_drz.log.
Blinded differences only. No unblinding.
Run: venvs/b3/bin/python src/dr2_relabel.py
"""
import numpy as np
import pandas as pd
from astropy.cosmology import Planck18 as P18
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy.spatial import cKDTree

H = 0.6766
vc = []
with open("data/vast/VoidFinder-nsa_v1_0_1_Planck2018_comoving_maximal.txt", encoding="utf-8-sig") as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        p = line.split()
        ra, dec, radius = float(p[7]), float(p[8]), float(p[3])
        D = float(p[6]) / H
        lon, lat = np.deg2rad(ra), np.deg2rad(dec)
        vc.append((np.array([D * np.cos(lat) * np.cos(lon), D * np.cos(lat) * np.sin(lon),
                             D * np.sin(lat)]), radius / H))
q = pd.read_csv("data/dr2_hr.csv")
sc = SkyCoord(ra=q.ra.values * u.deg, dec=q.dec.values * u.deg,
              distance=P18.comoving_distance(q.redshift.values))
SX = np.vstack([sc.cartesian.x.value, sc.cartesian.y.value, sc.cartesian.z.value]).T
vx = np.array([v[0] for v in vc])
vr = np.array([v[1] for v in vc])
tree = cKDTree(vx)
mR = vr.max()
rr, hh = [], []
for i in range(len(SX)):
    cand = tree.query_ball_point(SX[i], r=mR * 1.2)
    if not cand:
        rr.append(np.inf)
        hh.append(-1)
        continue
    vals = [(np.linalg.norm(SX[i] - vx[j]) / vr[j], j) for j in cand]
    best = min(vals)
    rr.append(best[0])
    hh.append(best[1])
q["r_Rv_drz"] = rr
q["host_drz"] = hh
q["env_drz"] = ["void" if v < 0.8 else ("shell" if v <= 1.0 else "wall") if np.isfinite(v) else "outside"
                for v in rr]
print(pd.crosstab(q.env, q.env_drz))
print(f"agreement BTS-z vs DR2-z labels: {(q.env == q.env_drz).mean() * 100:.1f}%")
qv = q[q.env_drz == "void"]
print(f"DR2-z void: {len(qv)} (was {(q.env == 'void').sum()} BTS-z)")
q.to_csv("data/dr2_hr_relabeled.csv", index=False)

# Blinded matched Delta on relabeled sample
qw = q[q.env_drz == "wall"]
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
    print(f"RELABELED matched Delta BLINDED = {diffs.mean():+.4f} "
          f"se={diffs.std(ddof=1) / np.sqrt(len(diffs)):.4f} (n={len(diffs)})")
else:
    print(f"too few matched ({len(diffs)}) for Delta")
print("wrote data/dr2_hr_relabeled.csv (BLINDED)")
