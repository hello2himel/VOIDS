"""VOIDS publication plots — generated from data (no hand numbers in legends).
Reads data/b3_match_v2.parquet (cosmo incl 14 non-SN) + excludes non-SN for SNe panels.
Writes plots/rRv_hist_v2.png + plots/footprint_v2.png (MNRAS single-column 3.3in, serif).
Run: venvs/b3/bin/python src/make_plots.py
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"font.family": "serif", "font.serif": ["Times New Roman", "DejaVu Serif"],
                     "font.size": 9, "axes.labelsize": 9, "xtick.labelsize": 8,
                     "ytick.labelsize": 8, "legend.fontsize": 8, "figure.dpi": 150})
df = pd.read_parquet("data/b3_match_v2.parquet")
sn = df[~df["type"].str.contains("TDE|Ca-rich|^other$")].copy()
print(f"plot sample: SNe N={len(sn)}")

d = sn[np.isfinite(sn.r_Rv) & (sn.r_Rv < 2.6)].copy()
fig, ax = plt.subplots(figsize=(3.3, 2.6))
ax.hist(d["r_Rv"], bins=52, range=(0, 2.6), color="#1f77b4", edgecolor="none")
ax.axvline(0.8, color="red", lw=1.2, label="void $r/R_{max}=0.8$")
ax.axvline(1.0, color="black", lw=1.2, label="wall $r/R_{max}=1.0$")
ax.set_xlabel("$r/R_{max}$")
ax.set_ylabel("N SNe")
ax.set_xlim(0, 2.6)
ax.legend(frameon=True, loc="lower left", fontsize=7)
fig.tight_layout()
fig.savefig("plots/rRv_hist_v2.png", dpi=150)
fig.savefig("paper/figures/rRv_hist_v2.png", dpi=150)
plt.close(fig)

fig, ax = plt.subplots(figsize=(3.3, 2.8))
colors = {"wall": "#2ca02c", "shell": "#ff7f0e", "void": "#d62728"}
import numpy as np
fin = np.isfinite(sn.r_Rv)
for env, sel in [("wall", (sn.r_Rv > 1.0) & fin), ("shell", (sn.r_Rv >= 0.8) & (sn.r_Rv <= 1.0)),
                 ("void", sn.r_Rv < 0.8)]:
    s = sn[sel]
    ax.scatter(s["ra"], s["dec"], s=4, c=colors[env],
               label=f"{env} (N={len(s)})", alpha=0.7, linewidths=0)
ax.set_xlabel("RA (deg)")
ax.set_ylabel("Dec (deg)")
ax.legend(frameon=True, loc="upper right", markerscale=3)
fig.tight_layout()
fig.savefig("plots/footprint_v2.png", dpi=150)
fig.savefig("paper/figures/footprint_v2.png", dpi=150)
plt.close(fig)
print("wrote plots/rRv_hist_v2.png plots/footprint_v2.png (+ copies in paper/figures/)")
