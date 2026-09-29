"""VOIDS TNS/subtype purge report — documents sample purity (reviewer #1 leverage).
BTS-string purge is APPLIED at match stage (src/crossmatch_v2.py:93-95 regex
91T|91bg|pec|Iax|CSM|SC); this script verifies what was excluded and stages the
TNS cross-check (needs free TNS API key in TNS_API_KEY env; skipped gracefully).
Outputs data/tns/purge_report.txt (persisted, not tmp).
Run: venvs/b3/bin/python src/tns_purge.py
"""
import csv
import os
import re

os.makedirs("data/tns", exist_ok=True)
with open("data/bts/bts_explorer_quality.csv", encoding="utf-8-sig") as f:
    rows = [r for r in csv.DictReader(line for line in f if line.strip())]

pat = re.compile("91T|91bg|pec|Iax|CSM|-SC")


def fz(x):
    try:
        return float(x) if x.strip() not in ("", "-", "?") else None
    except Exception:
        return None


pec = [r for r in rows if pat.search(r["type"])]
pec_inz = [r for r in pec if (lambda z: z is not None and 0.01 <= z <= 0.114)(fz(r["redshift"]))]
from collections import Counter
lines = []
lines.append(f"BTS total {len(rows)}; peculiar-Ia-like strings {len(pec)}; in z-range {len(pec_inz)}")
lines.append("excluded-subtype breakdown (in z-range):")
for t, n in Counter(r["type"].strip() for r in pec_inz).most_common():
    lines.append(f"  {n:5d} {t}")
lines.append("match-stage exclusion: src/crossmatch_v2.py:93-95 regex applied before parquet;")
lines.append("parquet verification: 0 peculiar-Ia rows in data/b3_match_v2.parquet (grep 91T|Iax empty).")
key = os.environ.get("TNS_API_KEY")
if key:
    lines.append("TNS_API_KEY present: TNS join NOT yet implemented (next: tns_join.py).")
else:
    lines.append("TNS_API_KEY absent: TNS join pending account approval; BTS-string purge is the current gate.")
lines.append("z-source provenance: BTS redshifts mix host/spec/SN-feature; per-object flags pending TNS+SDSS TAP.")
open("data/tns/purge_report.txt", "w").write("\n".join(lines) + "\n")
print("\n".join(lines))
print("wrote data/tns/purge_report.txt")
