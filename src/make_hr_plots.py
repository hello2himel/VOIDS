"""VOIDS HR diagnostic figures (unblinded values; differences only ever quoted).
Reads data/dr2_hr.csv. Writes plots/hr_vs_rr.png (HR vs r/Rmax by env) and
plots/mass_z.png (mass and redshift split void vs wall) + copies to paper/figures/.
Run: venvs/b3/bin/python src/make_hr_plots.py
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import shutil

plt.rcParams.update({"font.family": "serif", "font.serif": ["Times New Roman", "DejaVu Serif"],
                     "font.size": 9, "axes.labelsize": 9, "xtick.labelsize": 8,
                     "ytick.labelsize": 8, "legend.fontsize": 8, "figure.dpi": 150})
q = pd.read_csv("data/dr2_hr.csv")
qv, qw = q[q.env == "void"], q[q.env == "wall"]
fig, ax = plt.subplots(figsize=(3.3, 2.6))
ax.errorbar(qv.r_Rv, qv.hr_blind, yerr=np.sqrt(qv.sig2), fmt="o", ms=3, label=f"void (N={len(qv)})", alpha=0.7)
ax.errorbar(qw.r_Rv, qw.hr_blind, yerr=np.sqrt(qw.sig2), fmt="s", ms=2, label=f"wall (N={len(qw)})", alpha=0.4)
ax.axhline(0, color="k", lw=0.8)
ax.set_xlabel("$r/R_{\\max}$ (BTS-$z$ labels)")
ax.set_ylabel("HR (mag, unblinded 2026-10-01)")
ax.legend(frameon=True, loc="upper right")
fig.tight_layout()
fig.savefig("plots/hr_vs_rr.png", dpi=150)
shutil.copy("plots/hr_vs_rr.png", "paper/figures/hr_vs_rr.png")
plt.close(fig)
fig, ax = plt.subplots(1, 2, figsize=(6.6, 2.6))
ax[0].hist([qv.redshift, qw.redshift], bins=20, label=[f"void (N={len(qv)})", f"wall (N={len(qw)})"],
           color=["#d62728", "#2ca02c"], alpha=0.6)
ax[0].set_xlabel("$z$ (DR2)")
ax[0].set_ylabel("N")
ax[0].legend(frameon=True)
ax[1].hist([qv.mass.dropna(), qw.mass.dropna()], bins=15, label=["void", "wall"],
           color=["#d62728", "#2ca02c"], alpha=0.6)
ax[1].set_xlabel(r"$\log M_*$")
ax[1].legend(frameon=True)
fig.tight_layout()
fig.savefig("plots/mass_z.png", dpi=150)
shutil.copy("plots/mass_z.png", "paper/figures/mass_z.png")
plt.close(fig)
print("wrote plots/hr_vs_rr.png plots/mass_z.png (+ paper/figures/ copies)")
