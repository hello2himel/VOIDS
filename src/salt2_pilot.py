"""VOIDS SALT2 blinded pilot — first Hubble-residual numbers (reviewer top leverage).
Pilot: 20 void normal-Ia (z<=0.06 preferred) + 20 redshift-matched wall Ia.
Steps: IRSA ZTF DR24 forced photometry (public, no login) -> data/salt2/lc_*.csv
       sncosmo SALT2 fit with blinded global offset -> data/salt2/hr_pilot.csv + plots/hr_pilot.png
Blinding: all mu shifted by hidden constant from logs/blind_key.txt (sealed, not in repo text).
Run: venvs/b3/bin/python src/salt2_pilot.py
"""
import os
import numpy as np
import pandas as pd

os.makedirs("data/salt2", exist_ok=True)
df = pd.read_parquet("data/b3_match_v2.parquet")
void_ia = df[(df.r_Rv < 0.8) & (df["type"] == "SN Ia")].sort_values("z")
wall_ia = df[(np.isfinite(df.r_Rv)) & (df.r_Rv > 1.0) & (df["type"] == "SN Ia")].sort_values("z")
pilot_void = void_ia[void_ia.z <= 0.06].head(20)
if len(pilot_void) < 20:
    pilot_void = void_ia.head(20)
# Redshift-matched wall controls: nearest-z wall Ia per void pilot object, WITHOUT replacement
used = set()
picks = []
for _, r in pilot_void.iterrows():
    d = (wall_ia.z - r.z).abs().sort_values()
    for idx in d.index:
        if idx not in used:
            used.add(idx)
            picks.append(idx)
            break
pilot_wall = wall_ia.loc[picks]
pilot = pd.concat([pilot_void.assign(env="void"), pilot_wall.assign(env="wall")])
pilot[["ztfid", "type", "ra", "dec", "z", "env"]].to_csv("data/salt2/pilot_sample.csv", index=False)
print(f"pilot: void {len(pilot_void)} wall {len(pilot_wall)}")
print(pilot[["ztfid", "z", "env"]].to_string())
print("wrote data/salt2/pilot_sample.csv (lightcurve fetch = step 2)")
