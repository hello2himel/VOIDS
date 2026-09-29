"""VOIDS counts table generator -- single source of truth.
Reads data/b3_match_v2.parquet (cosmo) + BTS CSV + voids + mask (all-types recomputed, no log copies).
Writes paper/tables/tab_counts.tex + paper/tables/numbers.tex (paper inputs both; never hand numbers).
Wall means finite 1.0<r/Rv<inf; inf is Outside (tabulated, not dropped).
Run: venvs/b3/bin/python src/make_counts_table.py
"""
import csv
import pickle
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

H = 0.7
with open("data/vast/NSA_main_mask.pickle", "rb") as f:
    mask, _, _ = pickle.load(f)

def in_mask(ra, dec):
    return bool(mask[min(max(int(np.floor(ra % 360)), 0), 359),
                     min(max(int(np.floor(dec + 90)), 0), 179)])

def parse_ra(hms):
    h, m, s = [float(x) for x in hms.split(":")]
    return (h + m / 60 + s / 3600) * 15.0

def parse_dec(dms):
    sign = -1 if dms.strip()[0] == "-" else 1
    dms = dms.strip().lstrip("+-")
    d, m, s = [float(x) for x in dms.split(":")]
    return sign * (d + m / 60 + s / 3600)

voids = []
with open("data/vast/VoidFinder-nsa_v1_0_1_Planck2018_comoving_maximal.txt", encoding="utf-8-sig") as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        p = line.split()
        voids.append((float(p[0]) / H, float(p[1]) / H, float(p[2]) / H, float(p[3]) / H))
vxyz = np.array([[v[0], v[1], v[2]] for v in voids])
vrad = np.array([v[3] for v in voids])
tree = cKDTree(vxyz)
maxRv = vrad.max()

# All-types masked: recompute from BTS CSV (no log copies)
from astropy.cosmology import Planck18
from astropy.coordinates import SkyCoord
import astropy.units as u
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
            ra = parse_ra(r["RA"].strip())
            dec = parse_dec(r["Dec"].strip())
        except Exception:
            continue
        if not in_mask(ra, dec):
            continue
        rows.append((ra, dec, z))
cc = SkyCoord(ra=[t[0] for t in rows] * u.deg, dec=[t[1] for t in rows] * u.deg,
              distance=Planck18.comoving_distance([t[2] for t in rows]))
sxyz = np.vstack([cc.cartesian.x.value, cc.cartesian.y.value, cc.cartesian.z.value]).T
a_void = a_shell = a_wall = a_out = 0
for i in range(len(sxyz)):
    cand = tree.query_ball_point(sxyz[i], r=maxRv * 1.2)
    if not cand:
        a_out += 1
        continue
    v = min(np.linalg.norm(sxyz[i] - vxyz[j]) / vrad[j] for j in cand)
    if v < 0.8:
        a_void += 1
    elif v <= 1.0:
        a_shell += 1
    else:
        a_wall += 1
a_tot = len(rows)
print(f"all-types masked recomputed: N={a_tot} void={a_void} shell={a_shell} wall={a_wall} outside={a_out}")
assert a_void + a_shell + a_wall + a_out == a_tot

# Cosmo-grade from parquet (finite wall, inf outside)
df = pd.read_parquet("data/b3_match_v2.parquet")
# Cosmo-grade from parquet (finite wall, inf outside). Non-SN transients excluded from
# cosmology rows but retained in All-types: TDE 12 + other 1 + Ca-rich 1 = 14.
is_other = df["type"].str.contains("TDE|Ca-rich|^other$")
sn = df[~is_other].copy()
print(f"non-SN excluded from cosmo rows: {int(is_other.sum())}")
print(df[is_other]["type"].value_counts().to_string())
print(df[is_other].groupby(["type", "env"]).size().to_string())
fin = np.isfinite(sn.r_Rv)
is_ia = sn["type"] == "SN Ia"
is_cc = sn["type"].str.contains("SN II|SN Ib|SN Ic|SLSN")
def parts(mask):
    s = sn[mask]
    f = np.isfinite(s.r_Rv)
    return (len(s), int((s.r_Rv < 0.8).sum()), int(((s.r_Rv >= 0.8) & (s.r_Rv <= 1.0)).sum()),
            int(((s.r_Rv > 1.0) & f).sum()), int((~f).sum()))
N, n_void, n_shell, n_wall, n_out = parts(np.ones(len(sn), bool))
Nia, n_ia_void, n_ia_shell, n_ia_wall, n_ia_out = parts(is_ia.values)
Ncc, n_cc_void, n_cc_shell, n_cc_wall, n_cc_out = parts(is_cc.values)
assert n_void + n_shell + n_wall + n_out == N == 2367, (n_void, n_shell, n_wall, n_out, N)
print(f"cosmo N={len(df)} void={n_void} shell={n_shell} wall={n_wall} outside={n_out}")

# Rates from randoms_summary f values
f_void = f_wall = None
with open("data/randoms_summary.txt") as f:
    for line in f:
        if line.startswith("f_void="):
            f_void = float(line.split("=")[1].split("+")[0])
        if line.startswith("f_wall_finite="):
            f_wall = float(line.split("=")[1])
R_all = (n_void / n_wall) * (f_wall / f_void)
R_ia = (n_ia_void / n_ia_wall) * (f_wall / f_void)
R_cc = (n_cc_void / n_cc_wall) * (f_wall / f_void)
print(f"R all={R_all:.3f} Ia={R_ia:.3f} CC={R_cc:.3f}")

with open("paper/tables/tab_counts.tex", "w") as f:
    f.write("\\begin{table}\n\\caption{Masked sample composition. Wall means finite $1.0<r/R_{\\max}<\\infty$; Outside is $r/R_{\\max}=\\infty$ (no void within $1.2\\max R_{\\max}$). Cosmo rows are spectroscopically classified SNe only (TDE/Ca-rich/other excluded); CC includes 6 SLSN; Ia subtype audit (91T-like etc.) pending TNS.}\n")
    f.write("\\label{tab:counts}\n{\\small\\setlength{\\tabcolsep}{3pt}\n\\begin{tabular}{lrrrrr}\n\\hline\n")
    f.write("Sample & Total & Void & Shell & Wall & Outside \\\\\n\\hline\n")
    f.write(f"All-types & {a_tot} & {a_void} & {a_shell} & {a_wall} & {a_out} \\\\\n")
    f.write(f"Cosmo SNe & {N} & {n_void} & {n_shell} & {n_wall} & {n_out} \\\\\n")
    f.write(f"SN Ia & {Nia} & {n_ia_void} & {n_ia_shell} & {n_ia_wall} & {n_ia_out} \\\\\n")
    f.write(f"CC & {Ncc} & {n_cc_void} & {n_cc_shell} & {n_cc_wall} & {n_cc_out} \\\\\n")
    f.write("\\hline\n\\end{tabular}}\n\\end{table}\n")
with open("paper/tables/numbers.tex", "w") as f:
    f.write(f"\\newcommand{{\\NVoidCosmo}}{{{n_void}}}\n")
    f.write(f"\\newcommand{{\\NWallCosmo}}{{{n_wall}}}\n")
    f.write(f"\\newcommand{{\\NShellCosmo}}{{{n_shell}}}\n")
    f.write(f"\\newcommand{{\\NOutCosmo}}{{{n_out}}}\n")
    f.write(f"\\newcommand{{\\NIaVoid}}{{{n_ia_void}}}\n")
    f.write(f"\\newcommand{{\\NIaWall}}{{{n_ia_wall}}}\n")
    f.write(f"\\newcommand{{\\NIaShell}}{{{n_ia_shell}}}\n")
    f.write(f"\\newcommand{{\\NIaOut}}{{{n_ia_out}}}\n")
    f.write(f"\\newcommand{{\\NCCVoid}}{{{n_cc_void}}}\n")
    f.write(f"\\newcommand{{\\NCCWall}}{{{n_cc_wall}}}\n")
    f.write(f"\\newcommand{{\\NCCShell}}{{{n_cc_shell}}}\n")
    f.write(f"\\newcommand{{\\NCCOut}}{{{n_cc_out}}}\n")
    f.write(f"\\newcommand{{\\NCosmo}}{{{N}}}\n")
    f.write(f"\\newcommand{{\\RRall}}{{{R_all:.2f}}}\n")
    f.write(f"\\newcommand{{\\RRia}}{{{R_ia:.2f}}}\n")
    f.write(f"\\newcommand{{\\RRcc}}{{{R_cc:.2f}}}\n")
print("wrote paper/tables/tab_counts.tex paper/tables/numbers.tex")
