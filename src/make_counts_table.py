"""VOIDS counts table generator -- single source of truth (fixes hand-copy incoherence).
Reads data/b3_match_v2.parquet, writes tables + macros. Wall means finite 1.0<r/Rv<inf; inf is 'outside' (60).
Run: venvs/b3/bin/python src/make_counts_table.py
"""
import pandas as pd
import numpy as np

df = pd.read_parquet("data/b3_match_v2.parquet")
fin = np.isfinite(df.r_Rv)
# NOTE: b3_match_v2.parquet is masked+cosmo only (2381). All-types masked counts from v2 log:
# MASKED all-types: wall 1950 / shell 271 / void 183 / outside 61 (total 2465). Recompute Ia/CC here:
is_ia = df["type"] == "SN Ia"
is_cc = df["type"].str.contains("SN II|SN Ib|SN Ic|SLSN")
n_void = int((df.r_Rv < 0.8).sum())
n_shell = int(((df.r_Rv >= 0.8) & (df.r_Rv <= 1.0)).sum())
n_wall = int(((df.r_Rv > 1.0) & fin).sum())
n_out = int((~fin).sum())
n_ia_void = int(((df.r_Rv < 0.8) & is_ia).sum())
n_ia_wall = int(((df.r_Rv > 1.0) & fin & is_ia).sum())
n_cc_void = int(((df.r_Rv < 0.8) & is_cc).sum())
n_cc_wall = int(((df.r_Rv > 1.0) & fin & is_cc).sum())
print(f"cosmo N={len(df)} void={n_void} shell={n_shell} wall={n_wall} outside/inf={n_out}")
print(f"Ia void={n_ia_void} wall={n_ia_wall} | CC void={n_cc_void} wall={n_cc_wall}")
assert n_void + n_shell + n_wall + n_out == len(df), "env partition must sum"

with open("paper/tables/tab_counts.tex", "w") as f:
    f.write("\\begin{table}\n\\caption{Masked sample composition. Cosmo-grade excludes 91T/91bg/pec/Iax/CSM/SC.}\n")
    f.write("\\label{tab:counts}\n{\\small\\setlength{\\tabcolsep}{3pt}\n\\begin{tabular}{lrrrr}\n\\hline\n")
    f.write("Sample & Total & Void & Shell & Wall \\\\\n\\hline\n")
    f.write(f"All-types masked & 2465 & 183 & 271 & 1950 \\\\\n")
    f.write(f"Cosmo-grade & {len(df)} & {n_void} & {n_shell} & {n_wall} \\\\\n")
    f.write(f"Normal Ia & {int(is_ia.sum())} & {n_ia_void} & -- & {n_ia_wall} \\\\\n")
    f.write(f"CC & {int(is_cc.sum())} & {n_cc_void} & -- & {n_cc_wall} \\\\\n")
    f.write("\\hline\n\\end{tabular}}\n\\end{table}\n")
with open("paper/tables/numbers.tex", "w") as f:
    f.write(f"\\newcommand{{\\NVoidCosmo}}{{{n_void}}}\n")
    f.write(f"\\newcommand{{\\NWallCosmo}}{{{n_wall}}}\n")
    f.write(f"\\newcommand{{\\NShellCosmo}}{{{n_shell}}}\n")
    f.write(f"\\newcommand{{\\NIaVoid}}{{{n_ia_void}}}\n")
    f.write(f"\\newcommand{{\\NIaWall}}{{{n_ia_wall}}}\n")
    f.write(f"\\newcommand{{\\NCCVoid}}{{{n_cc_void}}}\n")
    f.write(f"\\newcommand{{\\NCCWall}}{{{n_cc_wall}}}\n")
print("wrote paper/tables/tab_counts.tex paper/tables/numbers.tex")
