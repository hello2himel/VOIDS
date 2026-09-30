"""VOIDS HR diagnostic figures (unblinded values; full-volume DR2-z primary sample).
Reads data/dr2_fullx.csv + DR2 errors for per-SN uncertainties.
Writes plots/hr_vs_rr.png (HR vs r/Rmax by env) and plots/mass_z.png
(mass and redshift split void vs wall) + copies to paper/figures/.
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
q = pd.read_csv("data/dr2_fullx.csv")
dr2 = pd.read_csv("data/ztf_dr2/ztfsniadr2_lite/tables/snia_data.csv", index_col=0)
e = dr2.set_index("ztfname")[["x0_err", "x1_err", "c_err", "cov_x0_x1", "cov_x0_c", "cov_x1_c"]]
q = q.join(e, on="ztfname")
ALPHA, BETA, SIG_INT = 0.14, 3.1, 0.12
sig2 = []
for _, r in q.iterrows():
    row = dr2.loc[dr2.ztfname == r.ztfname].iloc[0]
    x0v = row.x0
    j0v = -2.5 / (x0v * np.log(10))
    v = ((j0v * row.x0_err) ** 2 + (ALPHA * row.x1_err) ** 2 + (BETA * row.c_err) ** 2
         + 2 * ALPHA * j0v * row.cov_x0_x1 - 2 * BETA * j0v * row.cov_x0_c
         - 2 * ALPHA * BETA * row.cov_x1_c + SIG_INT ** 2)
    sig2.append(v if np.isfinite(v) and v > 0 else 0.0225)
q["sig_plot"] = np.sqrt(sig2)
qv, qw = q[q.env == "void"], q[q.env == "wall"]
fig, ax = plt.subplots(figsize=(3.3, 2.6))
ax.errorbar(qv.r_Rv, qv.hr, yerr=qv.sig_plot, fmt="o", ms=3, label=f"void (N={len(qv)})", alpha=0.7)
ax.errorbar(qw.r_Rv, qw.hr, yerr=qw.sig_plot, fmt="s", ms=2, label=f"wall (N={len(qw)})", alpha=0.4)
ax.axhline(0, color="k", lw=0.8)
ax.set_xlabel("$r/R_{\\max}$ (DR2-$z$ labels)")
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
print(f"wrote plots/hr_vs_rr.png plots/mass_z.png (fullx N={len(q)}, void {len(qv)})")
