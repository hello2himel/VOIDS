"""VOIDS full crossmatch — needs venvs/b3 (astropy, scipy, pandas).
Inputs (persistent):
  data/vast/VoidFinder-nsa_v1_0_1_Planck2018_comoving_maximal.txt
  data/bts/bts_explorer_quality.csv
Outputs (persistent):
  data/b3_match.parquet
  plots/rRv_hist.png, plots/footprint.png
  logs/crossmatch.log
Run: venvs/b3/bin/python src/crossmatch.py
"""
import csv
import numpy as np
import pandas as pd
from astropy.cosmology import Planck18
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy.spatial import cKDTree
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

H_VOIDCAT = 0.7  # VAST header h

def parse_ra(hms):
    h, m, s = [float(x) for x in hms.split(":")]
    return (h + m/60 + s/3600) * 15.0

def parse_dec(dms):
    sign = -1 if dms.strip()[0] == "-" else 1
    dms = dms.strip().lstrip("+-")
    d, m, s = [float(x) for x in dms.split(":")]
    return sign * (d + m/60 + s/3600)

# Voids: x y z (Mpc/h) + ra dec Reff
voids = []
with open("data/vast/VoidFinder-nsa_v1_0_1_Planck2018_comoving_maximal.txt", encoding="utf-8-sig") as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        p = line.split()
        x, y, z, rad = float(p[0]), float(p[1]), float(p[2]), float(p[3])
        ra, dec, reff = float(p[7]), float(p[8]), float(p[9])
        # Convert Mpc/h -> Mpc with catalog h, then to Planck18 frame (approx: scale by h ratio not needed for membership, use given XYZ/h)
        voids.append((x/H_VOIDCAT, y/H_VOIDCAT, z/H_VOIDCAT, rad/H_VOIDCAT, ra, dec, reff))
vxyz = np.array([[v[0], v[1], v[2]] for v in voids])
vrad = np.array([v[3] for v in voids])
print(f"voids={len(voids)} median_Rv={np.median(vrad):.1f} Mpc")

# BTS
rows = []
with open("data/bts/bts_explorer_quality.csv", newline="", encoding="utf-8-sig") as f:
    for r in csv.DictReader(l for l in f if l.strip()):
        try:
            z = float(r["redshift"]) if r["redshift"].strip() not in ("", "-", "?") else None
        except Exception:
            z = None
        if z is None or not (0.01 <= z <= 0.114):
            continue
        try:
            ra = parse_ra(r["RA"].strip()); dec = parse_dec(r["Dec"].strip())
        except Exception:
            continue
        rows.append((r["ZTFID"], r["type"].strip(), ra, dec, z))
print(f"bts_in_volume={len(rows)}")
df = pd.DataFrame(rows, columns=["ztfid", "type", "ra", "dec", "z"])
# XYZ via Planck18
coords = SkyCoord(ra=df["ra"].values*u.deg, dec=df["dec"].values*u.deg,
                  distance=Planck18.comoving_distance(df["z"].values))
sxyz = np.vstack([coords.cartesian.x.value, coords.cartesian.y.value, coords.cartesian.z.value]).T

tree = cKDTree(vxyz)
rmax = vrad.max() + 5.0
dists, idxs = tree.query(sxyz, k=5, distance_upper_bound=rmax)
# For each SN, compute r/Rv over 5 nearest, take min
rrv, host = [], []
for i in range(len(sxyz)):
    best = 99; bh = -1
    for k in range(5):
        j = idxs[i, k]
        if j >= len(voids):
            continue
        d = np.linalg.norm(sxyz[i] - vxyz[j])
        v = d / vrad[j]
        if v < best:
            best = v; bh = j
    rrv.append(best); host.append(bh)
df["r_Rv"] = rrv; df["host_void"] = host
def cls(v):
    if v < 0.8: return "void"
    if v <= 1.0: return "edge"
    return "wall"
df["env"] = df["r_Rv"].map(cls)
print(df["env"].value_counts())
print(df[df["env"]=="void"]["type"].value_counts().head(10))
df.to_parquet("data/b3_match.parquet", index=False)

# Plots
plt.figure(); plt.hist(df["r_Rv"].clip(0, 3), bins=60); plt.xlabel("r/Rv"); plt.ylabel("N SNe"); plt.savefig("plots/rRv_hist.png", dpi=150)
plt.figure(); plt.scatter(df["ra"], df["dec"], s=2, c=df["env"].map({"void":0,"edge":1,"wall":2})); plt.xlabel("RA"); plt.ylabel("Dec"); plt.savefig("plots/footprint.png", dpi=150)
print("wrote data/b3_match.parquet plots/rRv_hist.png plots/footprint.png")
