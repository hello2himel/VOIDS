"""VOIDS SALT2 pilot step 2 — fetch IRSA ZTF DR lightcurves + blinded sncosmo fits.
Reads data/salt2/pilot_sample.csv. For each: IRSA nph_light_curves POS=CIRCLE (public DR, no login)
-> data/salt2/lc_<ztfid>.csv. Then sncosmo SALT2 fit (phase -15..+45, bands g/r), blinded mu.
Blinding: offset drawn once into logs/blind_key.txt (NOT committed, gitignored); hr_blind = mu - mu_cosmo + blind.
Outputs data/salt2/hr_pilot.csv + plots/hr_pilot.png (pilot only, keeps main paper untouched).
Run: venvs/b3/bin/python src/salt2_fetch_fit.py (minutes; network + CPU)
"""
import os
import time
import numpy as np
import pandas as pd
import requests

os.makedirs("data/salt2", exist_ok=True)
sample = pd.read_csv("data/salt2/pilot_sample.csv")
URL = "https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves"
fetched = []
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
        fetched.append(out)
        print(f"fetched {r.ztfid} {len(q.text)} bytes", flush=True)
    except Exception as e:
        print(f"FAIL {r.ztfid}: {e}", flush=True)
    time.sleep(1)
print(f"fetched {len(fetched)}/{len(sample)}")
pd.Series(fetched).to_csv("data/salt2/fetched.txt", index=False)

# Blinded SALT2 fits
import sncosmo
from astropy.cosmology import Planck18
from astropy.table import Table

if not os.path.exists("logs/blind_key.txt"):
    rng = np.random.default_rng(20260930)
    blind = float(rng.uniform(-0.1, 0.1))
    open("logs/blind_key.txt", "w").write(f"{blind}\n")
    print("blind key sealed in logs/blind_key.txt (gitignored, never print value)")
blind = float(open("logs/blind_key.txt").read().strip())

model = sncosmo.Model(source="salt2")
rows = []
for _, r in sample.iterrows():
    f = f"data/salt2/lc_{r.ztfid}.csv"
    try:
        t = Table.read(f, format="ascii.csv")
    except Exception as e:
        print(f"SKIP {r.ztfid} unreadable: {e}")
        continue
    cols = {c.lower(): c for c in t.colnames}
    need = ["mjd", "mag", "magerr", "filter"]
    if not all(k in cols for k in ["mjd", "filter"]) or ("mag" not in cols and "flux" not in cols):
        print(f"SKIP {r.ztfid} cols {t.colnames}")
        continue
    data = Table()
    data["mjd"] = t[cols["mjd"]]
    data["band"] = [("ztfg" if str(x).strip().lower().startswith("g") else "ztfr") for x in t[cols["filter"]]]
    if "mag" in cols:
        data["mag"] = t[cols["mag"]]
        data["magerr"] = t[cols.get("magerr", cols["mag"])] if "magerr" in cols else np.ones(len(t)) * 0.05
        data["zp"] = 25.0
        data["zpsys"] = "ab"
    else:
        continue
    data = data[(data["mjd"] > 56000) & (data["mjd"] < 70000)]
    if len(data) < 8:
        print(f"SKIP {r.ztfid} only {len(data)} points")
        continue
    try:
        res, fitted = sncosmo.fit_lc(
            data, model,
            ["z", "t0", "x0", "x1", "c"],
            bounds={"z": (r.z - 0.005, r.z + 0.005)},
        )
        mu = -2.5 * np.log10(fitted.parameters[2]) if fitted.parameters[2] > 0 else np.nan
        mu_cosmo = Planck18.distmod(r.z).value
        rows.append(dict(ztfid=r.ztfid, z=r.z, env=r.env, x1=fitted.parameters[3],
                         c=fitted.parameters[4], mu_raw=mu,
                         hr_blind=(mu - mu_cosmo + blind)))
        print(f"FIT {r.ztfid} {r.env} x1={fitted.parameters[3]:.2f} c={fitted.parameters[4]:.3f}", flush=True)
    except Exception as e:
        print(f"FITFAIL {r.ztfid}: {str(e)[:200]}", flush=True)
hr = pd.DataFrame(rows)
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
