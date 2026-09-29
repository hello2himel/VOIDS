"""VOIDS pilot counts — stdlib only, no astropy needed.
Reads VAST maximals (ra,dec,Reff) + BTS CSV, does volume pre-cuts.
Full XYZ KDTree needs venv astropy/scipy (Phase 2).
Output persisted to logs/pilot_prelim.txt
"""
import csv, math

VAST = "data/vast/VoidFinder-nsa_v1_0_1_Planck2018_comoving_maximal.txt"
BTS = "data/bts/bts_explorer_quality.csv"

def parse_ra(hms):
    # "08:13:16.95" -> deg
    h, m, s = [float(x) for x in hms.split(":")]
    return (h + m/60 + s/3600) * 15.0

def parse_dec(dms):
    sign = -1 if dms.strip()[0] == "-" else 1
    dms = dms.strip().lstrip("+-")
    d, m, s = [float(x) for x in dms.split(":")]
    return sign * (d + m/60 + s/3600)

# Voids
voids = []
with open(VAST, encoding="utf-8-sig") as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        p = line.split()
        # x y z radius void edge r ra dec Reff
        ra, dec, reff = float(p[7]), float(p[8]), float(p[9])
        voids.append((ra, dec, reff))
print(f"voids={len(voids)}")

# BTS
sns = []
with open(BTS, newline="", encoding="utf-8-sig") as f:
    for row in csv.DictReader(l for l in f if l.strip()):
        try:
            z = float(row["redshift"]) if row["redshift"].strip() not in ("", "-", "?") else None
        except Exception:
            z = None
        try:
            ra = parse_ra(row["RA"].strip())
            dec = parse_dec(row["Dec"].strip())
        except Exception:
            continue
        sns.append((ra, dec, z, row["type"].strip()))

print(f"bts_total={len(sns)}")
zcut = [s for s in sns if s[2] is not None and 0.01 <= s[2] <= 0.114]
print(f"bts_z0.01-0.114={len(zcut)}")
from collections import Counter
print("types_in_zcut:", Counter(s[3] for s in zcut).most_common(10))
ia = [s for s in zcut if s[3].startswith("SN Ia")]
print(f"SN Ia in volume={len(ia)}")
# Rough sky overlap: SDSS DR7 NGC ~7500 deg2 / 41253 = 18% ; expect ~18% of BTS in footprint before mask
print(f"naive_footprint_estimate_18pct_Ia={int(len(ia)*0.18)}")
with open("logs/pilot_prelim.txt", "w") as out:
    out.write(f"voids={len(voids)}\n")
    out.write(f"bts_total={len(sns)}\n")
    out.write(f"bts_zcut={len(zcut)}\n")
    out.write(f"ia_in_volume={len(ia)}\n")
print("wrote logs/pilot_prelim.txt")
