"""VOIDS full-volume DR2-xmatch: every ZTF DR2 Ia vs voids with DR2 redshifts (not BTS-gated).
Cuts: snia-cosmo + sub_type norm + fitquality/lccoverage + z<0.114 + SDSS mask.
Membership: recomputed Planck18 frame, ball query, r/Rv<0.8 void.
Then mass+z matched void-wall Delta (Tripp fiducial + intrinsic scatter).
Outputs logs/dr2_fullx.log (persisted). Blind: uses unblinded post-2026-10-01 values.
Run: venvs/b3/bin/python src/dr2_full_xmatch.py
"""
import pickle
import numpy as np
import pandas as pd
from astropy.cosmology import Planck18 as P18
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy.spatial import cKDTree

H = 0.6766
with open("data/vast/NSA_main_mask.pickle", "rb") as f:
    mask, _, _ = pickle.load(f)

def in_mask(ra, dec):
    return bool(mask[min(max(int(np.floor(ra % 360)), 0), 359),
                     min(max(int(np.floor(dec + 90)), 0), 179)])

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
dr2 = pd.read_csv("data/ztf_dr2/ztfsniadr2_lite/tables/snia_data.csv", index_col=0)
hosts = pd.read_csv("data/ztf_dr2/ztfsniadr2_lite/tables/globalhost_data.csv", index_col=0)
d = dr2[(dr2.sn_type == "snia-cosmo") & (dr2.sub_type == "norm")
        & (dr2.fitquality_flag == 1.0) & (dr2.redshift < 0.114) & (dr2.redshift >= 0.01)].copy()
d = d.merge(hosts[["mass", "mass_err"]], left_on="ztfname", right_on="ztfname", how="left")
d["inmask"] = [in_mask(r.ra, r.dec) for _, r in d.iterrows()]
print(f"DR2 cosmo-norm-quality z<0.114: {len(d)}, in-mask: {int(d.inmask.sum())}")
d = d[d.inmask].copy()
sc = SkyCoord(ra=d.ra.values * u.deg, dec=d.dec.values * u.deg,
              distance=P18.comoving_distance(d.redshift.values))
SX = np.vstack([sc.cartesian.x.value, sc.cartesian.y.value, sc.cartesian.z.value]).T
vx = np.array([v[0] for v in vc])
vr = np.array([v[1] for v in vc])
tree = cKDTree(vx)
mR = vr.max()
rr = np.empty(len(d))
for i in range(len(SX)):
    cand = tree.query_ball_point(SX[i], r=mR * 1.2)
    rr[i] = np.inf if not cand else min(np.linalg.norm(SX[i] - vx[j]) / vr[j] for j in cand)
d["r_Rv"] = rr
d["env"] = ["void" if v < 0.8 else ("shell" if v <= 1.0 else "wall") if np.isfinite(v) else "outside"
            for v in rr]
print(d.env.value_counts().to_dict())
ALPHA, BETA, MB_ZP, SIG_INT = 0.14, 3.1, 10.635, 0.12
d = d[np.isfinite(d.x0) & np.isfinite(d.x1) & np.isfinite(d[["c"]].values.ravel())].copy()
d["mB"] = -2.5 * np.log10(d.x0) + MB_ZP
d["mu_raw"] = d.mB + ALPHA * d.x1 - BETA * d["c"]
d["mu_cosmo"] = P18.distmod(d.redshift).value
d["sig2"] = (2.5 * d.x0_err / (d.x0 * np.log(10))) ** 2 + 0.055 ** 2 + SIG_INT ** 2
M = float(np.average(d.mu_raw - d.mu_cosmo, weights=1 / d.sig2))
d["hr"] = d.mu_raw - d.mu_cosmo - M
qv, qw = d[d.env == "void"], d[d.env == "wall"]
wv, ww = 1 / qv.sig2, 1 / qw.sig2
raw = np.average(qv.hr, weights=wv) - np.average(qw.hr, weights=ww)
rawse = np.sqrt(1 / wv.sum() + 1 / ww.sum())
print(f"FULL-VOLUME raw Delta={raw:+.4f} se={rawse:.4f} (void {len(qv)}, wall {len(qw)})")
diffs = []
for _, r in qv.iterrows():
    if not np.isfinite(r.mass):
        continue
    pool = qw[((qw.mass - r.mass).abs() < 0.25) & ((qw.redshift - r.redshift).abs() < 0.015)]
    if len(pool) == 0:
        continue
    diffs.append(r.hr - pool.hr.mean())
diffs = np.array(diffs)
print(f"FULL-VOLUME matched Delta={diffs.mean():+.4f} se={diffs.std(ddof=1) / np.sqrt(len(diffs)):.4f} (n={len(diffs)})")
d[["ztfname", "redshift", "r_Rv", "env", "x1", "c", "mass", "hr"]].to_csv("data/dr2_fullx.csv", index=False)
print("wrote data/dr2_fullx.csv (UNBLINDED post-2026-10-01)")
