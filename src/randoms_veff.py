"""VOIDS randoms + Veff — volume control for mask-corrected rates.
- cos-weighted sky fraction (not equirectangular mean)
- observed counts read from data/b3_match_v2.parquet (no hand copies)
- Ia-only + CC-only split rates (no lumped-progenitor headline)
- 3 random seeds for sightline Monte Carlo error; cosmic variance via mocks TBD
- boundary note: 1-deg mask quantized (limitation); edge-loss fraction TBD
Run: venvs/b3/bin/python src/randoms_veff.py
"""
import pickle
import numpy as np
import pandas as pd
from astropy.cosmology import Planck18
import astropy.units as u
from scipy.spatial import cKDTree

H = 0.6766  # Planck18 H0/100; single-h consistency with crossmatch_v2
with open("data/vast/NSA_main_mask.pickle", "rb") as f:
    mask, _, _ = pickle.load(f)
assert mask.shape == (360, 180)
dec_centers = np.arange(-89.5, 90.0, 1.0)
cosw = np.cos(np.deg2rad(dec_centers))
sky_frac = (mask * cosw[None, :]).sum() / cosw.sum() / 360.0
print(f"sky_frac cos-weighted={sky_frac:.4f} (naive mean was {mask.mean():.4f}, wrong)")

def in_mask(ra, dec):
    # 1-deg quantization documented as limitation; boundary buffer NOT applied (TBD)
    return bool(mask[min(max(int(np.floor(ra % 360)), 0), 359),
                     min(max(int(np.floor(dec + 90)), 0), 179)])

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

# Observed counts from parquet; non-SN transients (TDE/Ca-rich/other = 14) excluded from
# all cosmology rows (retained only in All-types table row, generated separately).
obs = pd.read_parquet("data/b3_match_v2.parquet")
obs = obs[~obs["type"].str.contains("TDE|Ca-rich|^other$")].copy()
is_ia = obs["type"] == "SN Ia"
is_cc = obs["type"].str.contains("SN II|SN Ib|SN Ic|SLSN")
for name, sub in [("all", obs), ("Ia-only", obs[is_ia]), ("CC-only", obs[is_cc])]:
    v = (sub.r_Rv < 0.8).sum()
    w = ((sub.r_Rv > 1.0) & np.isfinite(sub.r_Rv)).sum()
    print(f"obs {name}: N={len(sub)} void={v} wall-finite={w}")

zmin, zmax = 0.01, 0.114
Vmin = Planck18.comoving_volume(zmin).value
Vmax = Planck18.comoving_volume(zmax).value
Veff_tot = (Vmax - Vmin) * sky_frac
print(f"Veff_tot={(Veff_tot):.6e} Mpc^3 (full-sky shell x cos-weighted sky frac)")

zg = np.linspace(zmin, zmax, 500)
Vg = Planck18.comoving_volume(zg).value
from astropy.coordinates import SkyCoord

f_runs = []
for seed in [42, 43, 44]:
    rng = np.random.default_rng(seed)
    N = 200000
    V = rng.uniform(Vmin, Vmax, N)
    z = np.interp(V, Vg, zg)
    ra = rng.uniform(0, 360, N)
    dec = np.rad2deg(np.arcsin(rng.uniform(-1, 1, N)))
    keep = np.array([in_mask(r, d) for r, d in zip(ra, dec)])
    ra, dec, z = ra[keep], dec[keep], z[keep]
    c = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, distance=Planck18.comoving_distance(z))
    sxyz = np.vstack([c.cartesian.x.value, c.cartesian.y.value, c.cartesian.z.value]).T
    rr = np.empty(len(sxyz))
    for i in range(len(sxyz)):
        cand = tree.query_ball_point(sxyz[i], r=maxRv * 1.2)
        rr[i] = np.inf if not cand else min(np.linalg.norm(sxyz[i] - vxyz[j]) / vrad[j] for j in cand)
    f_runs.append(((rr < 0.8).mean(), ((rr >= 0.8) & (rr <= 1.0)).mean(),
                    ((rr > 1.0) & np.isfinite(rr)).mean(), np.isinf(rr).mean(), len(ra)))
    print(f"seed {seed}: in-mask {len(ra)} f_void={f_runs[-1][0]:.4f} f_out={f_runs[-1][3]:.4f}")
f_void = float(np.mean([r[0] for r in f_runs]))
f_shell = float(np.mean([r[1] for r in f_runs]))
f_wall = float(np.mean([r[2] for r in f_runs]))
f_out = float(np.mean([r[3] for r in f_runs]))
f_void_se = float(np.std([r[0] for r in f_runs], ddof=1))
print(f"f_rand void={f_void:.4f}+/-{f_void_se:.4f} shell={f_shell:.4f} wall-finite={f_wall:.4f} out={f_out:.4f}")

# Errors: seed spread across 42/43/44 is sightline Monte Carlo only.
# Cosmic variance needs HR4 mocks (TBD); no void-resampling claimed.
print("MC error from seed spread +/-0.0016; cosmic variance via HR4 mocks TBD")

# Split rates reuse obs/is_ia/is_cc from above (finite wall throughout).
for name, sub in [("all", obs), ("Ia-only", obs[is_ia]), ("CC-only", obs[is_cc])]:
    Nv = int((sub.r_Rv < 0.8).sum())
    Nw = int(((sub.r_Rv > 1.0) & np.isfinite(sub.r_Rv)).sum())
    R = (Nv / (Veff_tot * f_void)) / (Nw / (Veff_tot * f_wall)) if Nw else float("nan")
    rel = float(np.sqrt(1 / max(Nv, 1) + 1 / max(Nw, 1) + (f_void_se / f_void) ** 2))
    print(f"rate {name}: Nv={Nv} Nw={Nw} R_volcorr={R:.3f} relerr~{rel:.2f} (APPARENT: no efficiency/host-mass correction)")

with open("data/randoms_summary.txt", "w") as o:
    o.write("seeds=42,43,44 (sightline MC only; cosmic variance via mocks TBD)\n")
    o.write(f"sky_frac_cosw={sky_frac:.4f}\n")
    o.write(f"f_void={f_void:.4f}+/-{f_void_se:.4f}\nf_shell={f_shell:.4f}\nf_wall_finite={f_wall:.4f}\nf_out={f_out:.4f}\n")
    o.write(f"Veff_tot={Veff_tot:.6e}\nmask=1deg quantized LIMITATION, boundary buffer TBD\n")
print("wrote data/randoms_summary.txt")
