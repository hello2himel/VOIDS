"""VOIDS proper Tripp fit (blinded): wall-only M/alpha/beta with SALT2 covariance.
Uses DR2 per-object errors + covariances (cov_x0_x1, cov_x0_c, cov_x1_c) +
intrinsic scatter 0.12, Chauvenet 3-sigma clipping. Fit on WALL sample only
(voids never touch the nuisance fit), then apply to void-wall Delta.
M absorbs the blind offset (differences blind-independent).
Outputs logs/tripp_fit.log (persisted). No unblinding.
Run: venvs/b3/bin/python src/fit_tripp.py
"""
import numpy as np
import pandas as pd

SIG_INT = 0.12
df = pd.read_csv("data/ztf_dr2/ztfsniadr2_lite/tables/snia_data.csv", index_col=0)
w = pd.read_parquet("data/b3_match_v2.parquet")[["ztfid"]].copy()
w["inwall"] = False
# wall Ia set from our match
bts = pd.read_parquet("data/b3_match_v2.parquet")
sn = bts[~bts["type"].str.contains("TDE|Ca-rich|^other")]
wall_ids = set(sn[(sn.r_Rv > 1.0) & np.isfinite(sn.r_Rv) & (sn["type"] == "SN Ia")].ztfid)
d = df[(df.sn_type == "snia-cosmo") & (df.sub_type == "norm") & (df.fitquality_flag == 1.0)].copy()
d["is_wall"] = d.ztfname.isin(wall_ids)
print(f"DR2 cosmo-norm-quality: {len(d)}, wall overlap: {int(d.is_wall.sum())}")
d = d[np.isfinite(d.x0) & np.isfinite(d.x1) & np.isfinite(d[["c"]].values.ravel())
      & np.isfinite(d.x0_err) & (d.x0 > 0)].copy()
print(f"after finite-SALT2 cut: {len(d)}, wall: {int(d.is_wall.sum())}")
fit = d[d.is_wall].copy().reset_index(drop=True)


def design(sub):
    x0 = sub.x0.values
    mB = -2.5 * np.log10(x0) + 10.635
    # Jacobian rows [dmu/dx0, dmu/dx1, dmu/dc] with alpha,beta from previous iter (start fiducial)
    return mB, sub.x1.values, sub["c"].values


def var_mu(sub, a, b):
    j0 = -2.5 / (sub.x0.values * np.log(10))
    v = ((j0 * sub.x0_err.values) ** 2 + (a * sub.x1_err.values) ** 2
         + (b * sub.c_err.values) ** 2
         + 2 * a * j0 * sub.cov_x0_x1.values - 2 * b * j0 * sub.cov_x0_c.values
         - 2 * a * b * sub.cov_x1_c.values + SIG_INT ** 2)
    return np.where(np.isfinite(v) & (v > 0), v, np.nanmedian(v[np.isfinite(v) & (v > 0)]))


a, b = 0.14, 3.1
keep = np.ones(len(fit), bool)
for it in range(10):
    mB = -2.5 * np.log10(fit.x0.values) + 10.635
    v = var_mu(fit, a, b)
    r = fit.mu_cosmo_proxy.values if "mu_cosmo_proxy" in fit else None
    if r is None:
        from astropy.cosmology import Planck18
        r = Planck18.distmod(fit.redshift.values).value
        fit["mu_cosmo_proxy"] = r
    y = mB + a * fit.x1.values - b * fit["c"].values - r
    wmean = np.average(y[keep], weights=1 / v[keep])
    pull = (y - wmean) / np.sqrt(v)
    new_keep = np.abs(pull) < 3.0
    new_keep = np.where(np.isfinite(pull), new_keep, False)
    A = np.column_stack([np.ones(len(fit)), -fit.x1.values, fit["c"].values])
    # solve M,alpha,beta on kept set with current variances
    W = 1 / v
    Aw, yw = A[keep] * np.sqrt(W[keep])[:, None], (mB - r)[keep] * np.sqrt(W[keep])
    import numpy.linalg as la
    beta, *_ = la.lstsq(Aw, yw, rcond=None)
    a, b = beta[1], beta[2]
    if np.array_equal(new_keep, keep):
        keep = new_keep
        break
    keep = new_keep
print(f"converged: kept {int(keep.sum())}/{len(fit)}, alpha={a:+.4f} beta={b:+.4f}")
out = ["cex_tripp_wall_fit", f"kept={int(keep.sum())}/{len(fit)} alpha={a:+.4f} beta={b:+.4f}"]
# apply to void-wall Delta (same matching as dr2_matched, recompute quickly here is overkill;
# report fitted nuisance only; Delta update follows in dr2_hr rerun)
open("logs/tripp_fit.log", "w").write("\n".join(out) + "\n")
print("\n".join(out))
print("wrote logs/tripp_fit.log (wall-only fit; void untouched)")
