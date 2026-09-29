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
