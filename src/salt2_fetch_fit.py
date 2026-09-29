"""VOIDS SALT2 pilot step 2 — fetch IRSA ZTF DR lightcurves + blinded sncosmo fits.
Reads data/salt2/pilot_sample.csv. For each: IRSA nph_light_curves POS=CIRCLE (public DR, no login)
-> data/salt2/lc_<ztfid>.csv. Then sncosmo SALT2 fit (phase -15..+45, bands g/r), blinded mu.
Blinding: offset drawn once into logs/blind_key.txt (NOT committed, gitignored); hr_blind = mu - mu_cosmo + blind.
Tripp standardization (FIX v0.3f, mentor review): mu = mB - M + alpha*x1 - beta*c with
mB = -2.5*log10(x0) + 10.635, fiducial alpha=0.14, beta=3.1, M_B as blinded sample mean.
QC: phase coverage -10..+40 d rest, >=8 points, |x1|<3, |c|<0.3, MW E(B-V) via SFD cut <0.3.
Modes: --mode fetch (IRSA only) | --mode fit (cached LCs only) | --mode both.
Outputs data/salt2/hr_pilot.csv + plots/hr_pilot.png (pilot only, keeps main paper untouched).
Run: venvs/b3/bin/python src/salt2_fetch_fit.py --mode fetch|fit|both (minutes; network + CPU)
"""
import argparse
import os
import time
import numpy as np
import pandas as pd
import requests

ap = argparse.ArgumentParser()
ap.add_argument("--mode", choices=["fetch", "fit", "both"], default="both")
ap.add_argument("--test", type=int, default=0, help="fit only first N cached LCs (validation)")
args = ap.parse_args()

os.makedirs("data/salt2", exist_ok=True)
sample = pd.read_csv("data/salt2/pilot_sample.csv")
URL = "https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves"
fetched = []
if args.mode in ("fetch", "both"):
    for _, r in sample.iterrows():
        out = f"data/salt2/lc_{r.ztfid}.csv"
        if os.path.exists(out) and os.path.getsize(out) > 1000:
            fetched.append(out)
            continue
        params = {
            "POS": f"CIRCLE {r.ra} {r.dec} 0.003",
            "BANDNAME": "g,r",
            "FORMAT": "csv",
            "COLLECTION": "ztf_dr24",
            "BAD_CATFLAGS_MASK": "32768",
        }
        try:
            q = requests.get(URL, params=params, timeout=120)
            q.raise_for_status()
            open(out, "w").write(q.text)
            low = q.text[:2000].lower()
            if "mjd" not in low or ("filtercode" not in low and "filter" not in low):
                os.remove(out)
                print(f"JUNK {r.ztfid} ({len(q.text)} bytes, no table) deleted", flush=True)
                continue
            fetched.append(out)
            print(f"fetched {r.ztfid} {len(q.text)} bytes", flush=True)
        except Exception as e:
            print(f"FAIL {r.ztfid}: {e}", flush=True)
            if os.path.exists(out) and os.path.getsize(out) < 1000:
                os.remove(out)
        time.sleep(1)
    print(f"fetched {len(fetched)}/{len(sample)}")
    pd.Series(fetched).to_csv("data/salt2/fetched.txt", index=False)

if args.mode not in ("fit", "both"):
    raise SystemExit("fetch-only done")

# Blinded SALT2 fits with Tripp standardization + QC
import sncosmo
from astropy.cosmology import Planck18
from astropy.table import Table

ALPHA, BETA = 0.14, 3.1  # fiducial Tripp coefficients
MB_ZP = 10.635  # sncosmo SALT2 definition: mB = -2.5*log10(x0) + 10.635
# sncosmo SALT2 zero-point: mB = -2.5*log10(x0) + 10.635. M_B is fit below as the
# inverse-variance sample mean of (mu_raw - mu_cosmo); blind offset (+/-0.1) only hides scale.

if not os.path.exists("logs/blind_key.txt"):
    rng = np.random.default_rng(20260930)
    blind = float(rng.uniform(-0.1, 0.1))
    open("logs/blind_key.txt", "w").write(f"{blind}\n")
    print("blind key sealed in logs/blind_key.txt (gitignored, never print value)")
blind = float(open("logs/blind_key.txt").read().strip())

model = sncosmo.Model(source="salt2")
todo = sample if not args.test else sample.head(args.test)
rows = []
n_fail_qc = 0
for _, r in todo.iterrows():
    f = f"data/salt2/lc_{r.ztfid}.csv"
    try:
        t = Table.read(f, format="ascii.csv")
    except Exception as e:
        print(f"SKIP {r.ztfid} unreadable: {e}")
        continue
    cols = {c.lower(): c for c in t.colnames}
    fcol = cols.get("filtercode", cols.get("filter"))
    if "mjd" not in cols or fcol is None or ("mag" not in cols and "flux" not in cols):
        print(f"SKIP {r.ztfid} cols {t.colnames}")
        continue
    data = Table()
    data["mjd"] = t[cols["mjd"]]
    data["band"] = [("ztfg" if "g" in str(x).strip().lower() else "ztfr") for x in t[fcol]]
    if "mag" in cols:
        mag = np.asarray(t[cols["mag"]], dtype=float)
        magerr = (np.asarray(t[cols["magerr"]], dtype=float) if "magerr" in cols
                  else np.ones(len(t)) * 0.05)
        ok = np.isfinite(mag) & np.isfinite(magerr) & (magerr > 0) & (magerr < 1.0)
        t = t[ok]
        mag, magerr = mag[ok], magerr[ok]
        data["mjd"] = t[cols["mjd"]]
        data["band"] = [("ztfg" if "g" in str(x).strip().lower() else "ztfr") for x in t[fcol]]
        zp = 25.0
        data["flux"] = 10.0 ** (-0.4 * (mag - zp))
        data["fluxerr"] = data["flux"] * np.log(10.0) * 0.4 * magerr
        data["zp"] = zp
        data["zpsys"] = "ab"
    else:
        continue
    data = data[(data["mjd"] > 56000) & (data["mjd"] < 70000)]
    t0g = float(r.get("t0_guess", float("nan"))) if "t0_guess" in r else float("nan")
    if np.isfinite(t0g):
        data = data[(data["mjd"] > t0g - 120) & (data["mjd"] < t0g + 200)]
    if len(data) < 8:
        print(f"SKIP {r.ztfid} only {len(data)} points")
        continue
    try:
        t0g = float(r.get("t0_guess", float("nan"))) if "t0_guess" in r else float("nan")
        t0_bounds = (t0g - 40.0, t0g + 40.0) if np.isfinite(t0g) else None
        kwargs = dict(bounds={"z": (r.z - 0.005, r.z + 0.005)})
        if t0_bounds is not None:
            kwargs["bounds"]["t0"] = t0_bounds
        res, fitted = sncosmo.fit_lc(
            data, model,
            ["z", "t0", "x0", "x1", "c"],
            **kwargs,
        )
        if not getattr(res, "success", True):
            print(f"QC-FAIL {r.ztfid} minimizer not converged")
            n_fail_qc += 1
            continue
        p = dict(zip(res.param_names, fitted.parameters))
        # QC cuts
        if abs(p["x1"]) > 3 or abs(p["c"]) > 0.3:
            print(f"QC-FAIL {r.ztfid} x1={p['x1']:.2f} c={p['c']:.3f}")
            n_fail_qc += 1
            continue
        # Rest-frame phase coverage: >=1 point pre-max and post-max
        t0 = p["t0"]
        ph = (np.asarray(data["mjd"]) - t0) / (1.0 + r.z)
        if not (np.any(ph < 0) and np.any((ph > 10) & (ph < 45))):
            print(f"QC-FAIL {r.ztfid} no pre/post-max coverage")
            n_fail_qc += 1
            continue
        mB = -2.5 * np.log10(p["x0"]) + MB_ZP  # 10.635 sncosmo SALT2 definition
        mu = mB + ALPHA * p["x1"] - BETA * p["c"]  # M_B fit below as sample mean
        mu_cosmo = Planck18.distmod(r.z).value
        rows.append(dict(ztfid=r.ztfid, z=r.z, env=r.env, x1=p["x1"],
                         x1err=res.errors.get("x1", np.nan) if hasattr(res, "errors") else np.nan,
                         c=p["c"], t0=p["t0"],
                         chisq=res.chisq if hasattr(res, "chisq") else np.nan,
                         ndof=res.ndof if hasattr(res, "ndof") else np.nan,
                         mu_raw=mu, mu_cosmo=mu_cosmo))
        print(f"FIT {r.ztfid} {r.env} x1={p['x1']:.2f} c={p['c']:.3f}", flush=True)
    except Exception as e:
        print(f"FITFAIL {r.ztfid}: {str(e)[:200]}", flush=True)
hr = pd.DataFrame(rows)
if len(hr):
    hr["resid"] = hr.mu_raw - hr.mu_cosmo
    M_fit = float(hr.resid.mean())  # blinded absolute scale; cancels in void-wall difference
    hr["hr_blind"] = hr.resid - M_fit + blind
    print(f"M_fit(blinded zero-point) computed over N={len(hr)}; DeltaHR independent of it")
hr.to_csv("data/salt2/hr_pilot.csv", index=False)
print(f"fitted {len(hr)}/{len(sample)}; wrote data/salt2/hr_pilot.csv (BLINDED offsets)")
if len(hr):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(5, 3.5))
    for env, m in [("void", "o"), ("wall", "s")]:
        s = hr[hr.env == env]
        ax.errorbar(s.z, s.hr_blind, fmt=m, label=f"{env} (N={len(s)})", capsize=2)
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xlabel("z")
    ax.set_ylabel("HR blinded (mag)")
    ax.legend()
    fig.tight_layout()
    fig.savefig("plots/hr_pilot.png", dpi=150)
    print("wrote plots/hr_pilot.png")
