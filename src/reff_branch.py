"""VOIDS Reff systematics branch (counts only, blinded-safe): membership with effective
(void-union) radii instead of maximal-sphere cores. Compares Rmax vs Reff branches.
Writes logs/reff_branch.log (persisted, not tmp).
Run: venvs/b3/bin/python src/reff_branch.py
"""
import pickle
import numpy as np
import pandas as pd
from astropy.cosmology import Planck18
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy.spatial import cKDTree

H = 0.6766
with open("data/vast/NSA_main_mask.pickle", "rb") as f:
    mask, _, _ = pickle.load(f)

voids = []
with open("data/vast/VoidFinder-nsa_v1_0_1_Planck2018_comoving_maximal.txt", encoding="utf-8-sig") as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        p = line.split()
        ra, dec, reff = float(p[7]), float(p[8]), float(p[9])
        D = float(p[6]) / H
        lon, lat = np.deg2rad(ra), np.deg2rad(dec)
        X = np.array([D * np.cos(lat) * np.cos(lon), D * np.cos(lat) * np.sin(lon), D * np.sin(lat)])
        voids.append((X, reff / H))
df = pd.read_parquet("data/b3_match_v2.parquet")
sn = df[~df["type"].str.contains("TDE|Ca-rich|^other")].copy()
c = SkyCoord(ra=sn.ra.values * u.deg, dec=sn.dec.values * u.deg,
             distance=Planck18.comoving_distance(sn.z.values))
sx = np.vstack([c.cartesian.x.value, c.cartesian.y.value, c.cartesian.z.value]).T
vx = np.array([v[0] for v in voids])
vr = np.array([v[1] for v in voids])
tree = cKDTree(vx)
mR = vr.max()
rr = np.empty(len(sx))
for i in range(len(sx)):
    cand = tree.query_ball_point(sx[i], r=mR * 1.2)
    rr[i] = np.inf if not cand else min(np.linalg.norm(sx[i] - vx[j]) / vr[j] for j in cand)
lines = []
lines.append(f"Reff branch: void={int((rr < 0.8).sum())} shell={int(((rr >= 0.8) & (rr <= 1.0)).sum())} "
             f"wall={int(((rr > 1.0) & np.isfinite(rr)).sum())} out={int(np.isinf(rr).sum())}")
lines.append("Rmax branch: void=142 shell=187 wall=1996 out=42")
s = sn[sn.z <= 0.06]
ia = s["type"] == "SN Ia"
lines.append(f"z<=0.06 all: N={len(s)} void={int((s.r_Rv < 0.8).sum())} "
             f"wall={int(((s.r_Rv > 1.0) & np.isfinite(s.r_Rv)).sum())}")
lines.append(f"z<=0.06 Ia: N={int(ia.sum())} void={int(((s.r_Rv < 0.8) & ia).sum())} "
             f"wall={int(((s.r_Rv > 1.0) & np.isfinite(s.r_Rv) & ia).sum())}")
open("logs/reff_branch.log", "w").write("\n".join(lines) + "\n")
print("\n".join(lines))
print("wrote logs/reff_branch.log")
