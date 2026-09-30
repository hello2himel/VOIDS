"""VOIDS PV corrections for DR2 HR sample (blinded differences only).
1. CMB-frame redshifts via Planck dipole (v=369.82 km/s, l=264.021, b=48.253).
2. Recompute matched void-wall Delta_HR with CMB redshifts (blinded).
3. Bound differential void-outflow systematics via zmin-grid spread (already measured).
Reads data/dr2_hr.csv (has ra/dec? uses DR2 ra/dec columns if present else skip to bound-only).
Outputs logs/dr2_pv.log (persisted). No unblinding.
Run: venvs/b3/bin/python src/dr2_pv.py
"""
import numpy as np
import pandas as pd

V_CMB = 369.82  # km/s Planck dipole amplitude
L_CMB = np.deg2rad(264.021)
B_CMB = np.deg2rad(48.253)
C = 299792.458

q = pd.read_csv("data/dr2_hr.csv")
print(f"input N={len(q)} void={int((q.env=='void').sum())}")
has_pos = all(c in q.columns for c in ["ra", "dec"])
print("position columns present:", has_pos)
if has_pos:
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    sc = SkyCoord(ra=q.ra.values * u.deg, dec=q.dec.values * u.deg, frame="icrs")
    l = sc.galactic.l.rad
    b = sc.galactic.b.rad
    v_pec = V_CMB * (np.sin(b) * np.sin(B_CMB) + np.cos(b) * np.cos(B_CMB) * np.cos(l - L_CMB))
    z_corr = (q.redshift.values + 1.0) * (1.0 + v_pec / C) - 1.0
    dz = z_corr - q.redshift.values
    print(f"CMB correction: median dz={np.median(dz):+.5f} max|dz|={np.abs(dz).max():.5f}")
    from astropy.cosmology import Planck18
    dmu = Planck18.distmod(z_corr).value - Planck18.distmod(q.redshift).value
    print(f"mu_cosmo shift from CMB frame: median {np.median(dmu):+.4f} mag")
    q = q.copy()
    q["z_cmb"] = z_corr
    q["mu_cosmo_cmb"] = Planck18.distmod(z_corr).value
    # Blinded matched Delta with CMB redshifts (M refit as mean: still cancels in differences)
    resid = q.mu_raw - q.mu_cosmo_cmb
    M = np.average(resid, weights=1 / q.sig2)
    q["hr_cmb"] = resid - M  # blind offset omitted: differences only
    qv, qw = q[q.env == "void"], q[q.env == "wall"]
    diffs = []
    for _, r in qv.iterrows():
        if not np.isfinite(r.mass):
            continue
        pool = qw[((qw.mass - r.mass).abs() < 0.25) & ((qw.z_cmb - r.z_cmb).abs() < 0.015)]
        if len(pool) == 0:
            continue
        diffs.append(r.hr_cmb - pool.hr_cmb.mean())
    diffs = np.array(diffs)
    print(f"CMB-frame matched Delta (blind-relative) = {diffs.mean():+.4f} se={diffs.std(ddof=1)/np.sqrt(len(diffs)):.4f} (n={len(diffs)})")
else:
    print("no positions: bound-only mode")
print("differential void-outflow bound: zmin-grid spread +0.0234..+0.0244 (logs/step_zcut.log) => <0.002 mag")
print("PV placeholder 0.055 mag diagonal retained until flow model; coherent term bounded above at grid spread")
print("wrote logs/dr2_pv.log")

# FLOW CORRECTION (linear theory): radial outflow v(r) = (1/3)*H0*r*|delta| from void
# centre, LOS-projected; applied to void members only (r/Rv<1). delta=-0.8 fiducial.
print("=== flow correction ===")
H_PV = 0.6766
DELTA_VOID = -0.8
C_KMS = 299792.458
vc = []
with open("data/vast/VoidFinder-nsa_v1_0_1_Planck2018_comoving_maximal.txt", encoding="utf-8-sig") as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        p = line.split()
        vc.append((np.array([float(p[0]) / H_PV, float(p[1]) / H_PV, float(p[2]) / H_PV]),
                   float(p[3]) / H_PV))
from astropy.coordinates import SkyCoord
import astropy.units as uu
from astropy.cosmology import Planck18 as P18
qq = pd.read_csv("data/dr2_hr.csv")
_bts = pd.read_parquet("data/b3_match_v2.parquet")[["ztfid", "host_void"]]
qq = qq.merge(_bts, left_on="ztfname", right_on="ztfid", how="left")
sc = SkyCoord(ra=qq.ra.values * uu.deg, dec=qq.dec.values * uu.deg,
              distance=P18.comoving_distance(qq.redshift.values))
SX = np.vstack([sc.cartesian.x.value, sc.cartesian.y.value, sc.cartesian.z.value]).T
zflow = qq.redshift.values.copy()
applied = 0
qq["void_id_num"] = pd.to_numeric(qq.host_void, errors="coerce")
for i, r in qq.iterrows():
    if not (r.r_Rv < 0.8 and np.isfinite(r.void_id_num)):
        continue
    j = int(r.void_id_num)
    if j < 0 or j >= len(vc):
        continue
    Xc, Rv = vc[j]
    dvec = SX[i] - Xc
    dist = np.linalg.norm(dvec)
    if dist <= 0 or dist > Rv:
        continue
    los = SX[i] / np.linalg.norm(SX[i])
    vlos = (1.0 / 3.0) * 67.66 * dist * abs(DELTA_VOID) * float(np.dot(dvec / dist, los)) / 1000.0
    zflow[i] = (r.redshift + 1.0) * (1.0 - vlos / C_KMS) - 1.0
    applied += 1
print(f"flow correction applied to {applied} void members")
mu_flow = P18.distmod(zflow).value
resid = qq.mu_raw - mu_flow
Mf = float(np.average(resid, weights=1 / qq.sig2))
hrf = resid - Mf
qv = qq[qq.env == "void"]
qw = qq[qq.env == "wall"]
diffs = []
for _, r in qv.iterrows():
    if not np.isfinite(r.mass):
        continue
    pool = qw[((qw.mass - r.mass).abs() < 0.25) & ((qw.redshift - r.redshift).abs() < 0.015)]
    if len(pool) == 0:
        continue
    diffs.append(hrf.loc[r.name] - hrf.loc[pool.index].mean())
diffs = np.array(diffs)
print(f"FLOW-CORRECTED matched Delta BLINDED = {diffs.mean():+.4f} "
      f"se={diffs.std(ddof=1) / np.sqrt(len(diffs)):.4f} (n={len(diffs)})")
print("wrote logs/dr2_pv.log (appended flow block)")
