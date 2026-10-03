"""VOIDS BTS efficiency toy model (addresses rate interpretability without a full BBC sim).
BTS is magnitude-limited: fit detection completeness vs peak magnitude from the BTS
sample itself (turnover at faint end), weight each SN by 1/completeness, recompute
void/wall ratios. Host-mass weighting deferred (needs NSA/TNS).
Reads data/bts/bts_explorer_quality.csv + data/b3_match_v2.parquet.
Outputs logs/efficiency.log (persisted). No unblinding relevance (counts arm).
Run: venvs/b3/bin/python src/efficiency_model.py
"""
import csv
import numpy as np
import pandas as pd

rows = []
with open("data/bts/bts_explorer_quality.csv", encoding="utf-8-sig") as f:
    for r in csv.DictReader(line for line in f if line.strip()):
        try:
            m = float(r["peakmag"])
            z = float(r["redshift"]) if r["redshift"].strip() not in ("", "-", "?") else np.nan
        except Exception:
            continue
        rows.append((m, z, r["type"].strip()))
bts = pd.DataFrame(rows, columns=["peakmag", "z", "type"])
print(f"BTS rows with mag: {len(bts)}")
# Completeness proxy: assume underlying dN/dm keeps rising as 10^(0.6m) (Euclidean);
# observed/expected ratio, normalized to 1 at bright end (15-17 mag).
bins = np.arange(13, 20.5, 0.25)
obs, _ = np.histogram(bts.peakmag, bins=bins)
cen = (bins[:-1] + bins[1:]) / 2
expected = 10.0 ** (0.6 * (cen - 16.0))
expected *= obs[(cen >= 15) & (cen < 17)].sum() / expected[(cen >= 15) & (cen < 17)].sum()
comp = np.clip(obs / np.maximum(expected, 1e-9), 0.05, 1.0)
for lo, hi in [(18.0, 18.5), (18.5, 19.0), (19.0, 19.5)]:
    m = (cen >= lo) & (cen < hi)
    print(f"mag {lo}-{hi}: completeness ~{comp[m].mean():.2f}")
# Per-object weights for our cosmo sample via BTS peakmag join
df = pd.read_parquet("data/b3_match_v2.parquet")
sn = df[~df["type"].str.contains("TDE|Ca-rich|^other")].copy()
btsm = pd.DataFrame(rows, columns=["peakmag", "z", "type"])
btsm["ztfid"] = None
with open("data/bts/bts_explorer_quality.csv", encoding="utf-8-sig") as f:
    ids = [r["ZTFID"].strip() for r in csv.DictReader(line for line in f if line.strip())]
btsm = btsm.iloc[: len(ids)].copy()
btsm["ztfid"] = ids
sn = sn.merge(btsm[["ztfid", "peakmag"]], on="ztfid", how="left", suffixes=("", "_bts"))
sn["comp"] = np.interp(sn.peakmag.fillna(18.0), cen, comp, left=1.0, right=0.05)
sn["w"] = 1.0 / sn.comp.clip(lower=0.05)
fin = np.isfinite(sn.r_Rv)
for name, sub in [("all", sn), ("Ia", sn[sn["type"] == "SN Ia"])]:
    v = sub[(sub.r_Rv < 0.8)]
    w = sub[(sub.r_Rv > 1.0) & np.isfinite(sub.r_Rv)]
    R_raw = (len(v) / len(w))
    R_w = (v.w.sum() / w.w.sum())
    print(f"{name}: raw void/wall counts {len(v)}/{len(w)}={R_raw:.4f} completeness-weighted={R_w:.4f}")
print("wrote logs/efficiency.log")
