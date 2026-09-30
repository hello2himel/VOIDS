# VOIDS unblinding protocol (frozen)

## Blinding scheme
All distance moduli carry a single additive blind offset drawn once into
`logs/blind_key.txt` (gitignored, never printed, never committed). Every quoted
void-minus-wall difference is blind-independent; only the absolute scale is sealed.

## Frozen cuts (v0.18+)
- Sample: ZTF DR2 `snia-cosmo` + `sub_type norm` + fitquality/lccoverage flags (N=396, 31 void, finite wall).
- BTS-string purge: 0/396 rows carry peculiar strings (91T/91bg/pec/Iax/CSM/SLSN/TDE);
  log `logs/purge_check.log`.
- Tripp fiducial alpha=0.14, beta=3.1 (grid-stable <0.017); intrinsic scatter 0.12;
  PV diagonal 0.055 + CMB-frame and zmin cross-checks.

## Unblinding thresholds (all must pass)
1. Nulls: RA-shuffle p>0.05, bootstrap==analytic within 10%, injections recover.
2. Closed systematics budget: zero `open` rows in `paper/tables/tab_budget.tex`.
3. Frozen cuts unchanged since blinding; any cut change re-blinds.

## Rotation
At unblinding: generate a fresh key with a new seed, re-run fits, quote unblinded
Delta with the rotation date logged. The old key is deleted, never published.
