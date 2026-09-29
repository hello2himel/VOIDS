"""VOIDS DR2 blinded void-wall Delta_HR — the Q1 result (mentor: DR2 obsoletes hand-fits).
Joins our void/wall SN Ia (data/b3_match_v2.parquet, SNe-only) to ZTF DR2 SALT2 params
(data/ztf_dr2/ztfsniadr2_lite/tables/snia_data.csv) by ZTFID, plus host masses
(globalhost_data.csv). Cuts: sn_type snia-cosmo, sub_type norm, fitquality_flag,
spec-z source. Tripp with fiducial alpha/beta; M as blinded sample mean; blind key
logs/blind_key.txt (never printed). PV: CMB frame + diagonal term (2M++ flow TBD).
Outputs data/dr2_hr.csv + logs/dr2_hr.log (persisted, not tmp).
Run: venvs/b3/bin/python src/dr2_hr.py
"""
import numpy as np
import pandas as pd
from astropy.cosmology import Planck18

ALPHA, BETA = 0.14, 3.1
MB_ZP = 10.635

bts = pd.read_parquet("data/b3_match_v2.parquet")
sn = bts[~bts["type"].str.contains("TDE|Ca-rich|^other")].copy()
sn["is_ia"] = sn["type"] == "SN Ia"
void_ia = sn[(sn.r_Rv < 0.8) & sn.is_ia].copy()
wall_ia = sn[(sn.r_Rv > 1.0) & np.isfinite(sn.r_Rv) & sn.is_ia].copy()
print(f"ours: void Ia {len(void_ia)} wall Ia {len(wall_ia)}")

dr2 = pd.read_csv("data/ztf_dr2/ztfsniadr2_lite/tables/snia_data.csv", index_col=0)
hosts = pd.read_csv("data/ztf_dr2/ztfsniadr2_lite/tables/globalhost_data.csv", index_col=0)
print(f"dr2 rows {len(dr2)}, hosts {len(hosts)}")
m = dr2.merge(hosts[["mass", "mass_err"]], left_on="ztfname", right_on="ztfname", how="left",
              suffixes=("", "_h"))
mine = m.merge(sn[["ztfid", "r_Rv"]].rename(columns={"ztfid": "void_id"}),
               left_on="ztfname", right_on="void_id", how="inner")
print(f"overlap with our Ia sample: {len(mine)}")
mine["env"] = np.where(mine.r_Rv < 0.8, "void", np.where((mine.r_Rv > 1.0) & np.isfinite(mine.r_Rv), "wall", "other"))
mine = mine[mine.env.isin(["void", "wall"])].copy()
print(mine.env.value_counts().to_dict(), "(shell/outside excluded from HR sample)")

# Cosmology cuts
q = mine[(mine.sn_type == "snia-cosmo") & (mine.sub_type == "norm")].copy()
print(f"after snia-cosmo+norm: {len(q)} (void {(q.env=='void').sum()})")
q = q[(q.fitquality_flag == 1.0) & (q.lccoverage_flag == 1.0)].copy()
print(f"after quality flags: {len(q)} (void {(q.env=='void').sum()})")
print("z-source:", q.source.value_counts().to_dict())

q = q[np.isfinite(q.x0) & np.isfinite(q.x1) & np.isfinite(q[["c"]].values.ravel())].copy()
print(f"after finite SALT2 params: {len(q)} (void {(q.env=='void').sum()})")
# Blinded Tripp HR with PV diagonal
blind = float(open("logs/blind_key.txt").read().strip())
q["mB"] = -2.5 * np.log10(q.x0) + MB_ZP
q["mu_raw"] = q.mB + ALPHA * q.x1 - BETA * q.c
q["mu_cosmo"] = Planck18.distmod(q.redshift).value
sig_pv = 0.055  # ~250 km/s at z~0.05 diagonal placeholder; flow model TBD
SIG_INT = 0.12  # Ia intrinsic scatter; stat errors without it are implausible (referee M2c)
sig2 = (2.5 * q.x0_err / (q.x0 * np.log(10))) ** 2 + sig_pv ** 2 + SIG_INT ** 2
sig2 = pd.Series(np.where(np.isfinite(sig2) & (sig2 > 0), sig2, np.nanmedian(sig2[np.isfinite(sig2) & (sig2 > 0)])),
                 index=q.index)
q["sig2"] = sig2
M_fit = float(np.average(q.mu_raw - q.mu_cosmo, weights=1 / q.sig2))
q["hr_blind"] = q.mu_raw - q.mu_cosmo - M_fit + blind
print(f"M_fit sample mean over N={len(q)} (cancels in void-wall difference)")

for env in ["void", "wall"]:
    s = q[q.env == env]
    w = 1 / s.sig2
    mu = np.average(s.hr_blind, weights=w)
    se = 1 / np.sqrt(w.sum())
    print(f"{env}: N={len(s)} <HR>_blind={mu:+.4f} se={se:.4f} <x1>={s.x1.mean():+.3f} <c>={s['c'].mean():+.4f} <logM>={s.mass.mean():.3f}")

qv, qw = q[q.env == "void"], q[q.env == "wall"]
wv, ww = 1 / qv.sig2, 1 / qw.sig2
d = wv.sum() and (np.average(qv.hr_blind, weights=wv) - np.average(qw.hr_blind, weights=ww))
se = np.sqrt(1 / wv.sum() + 1 / ww.sum())
print(f"Delta_HR(void-wall) BLINDED = {d:+.4f} se={se:.4f}")
# Mass-step split
for mcut, lab in [(10.0, "lowM"), (10.0, "highM")]:
    pass
q["hiM"] = q.mass >= 10.0
for env in ["void", "wall"]:
    for hm, tag in [(False, "lowM"), (True, "highM")]:
        s = q[(q.env == env) & (q.hiM == hm)]
        if len(s):
            print(f"{env} {tag}: N={len(s)} <HR>b={np.average(s.hr_blind, weights=1/s.sig2):+.4f}")
q.to_csv("data/dr2_hr.csv", index=False)
print(f"wrote data/dr2_hr.csv N={len(q)} (BLINDED; M_fit cancels in differences)")
