# Supernova environments in SDSS DR7 voids: ZTF BTS transients in VAST VoidFinder spheres with a CHIME FRB forecast
### VOIDS Paper Draft v0.2 — 30 September 2026 (addresses 4× Q1 Reject reviews)

**Author:** Himel Das
**Affiliation:** Department of Physics, Government Azizul Haque College, Bogura, Bangladesh
**ORCID:** 0009-0003-7698-1255
**Corresponding:** Himel Das

## Abstract
We cross-match VAST VoidFinder SDSS DR7 (1163 voids, $z<0.114$, median $R_{\max}=11.9\,h^{-1}$Mpc $=17.0$ Mpc, median $R_{\mathrm{eff}}=15.5\,h^{-1}$Mpc $=22.1$ Mpc) with ZTF BTS supernovae using a Planck18-consistent KD-tree ball query with real SDSS mask. Of 6873 BTS transients at $z=0.01$–$0.114$, 2465 lie in-mask: 183 void ($r/R_{\max}<0.8$), 271 shell ($0.8$–$1.0$), 1950 wall. Cosmology-grade (excluding 91T/91bg/pec/Iax/CSM/SC): 2381 total, 176 void (141 normal Ia). Volume-limited $z\le0.06$: 1278 total, 71 void. Raw Poisson sensitivity to a void–wall Hubble offset is $\sim0.035$ mag at 95% including systematics floor; SALT2 blinded standardization with mass/PV/BBC controls is in progress. CHIME Catalog 1 (600 bursts, 524 sources) is used for footprint and error-budget demonstration only; we show unlocalized $z<0.114$ stacking is underpowered by $>10\times$ ($\sigma_{\mathrm{tot}}\sim150$–$200$ vs signal $\sim20$ pc cm$^{-3}$) and defer joint constraints to Catalog 2 localized sample. Catalog, randoms design, and pipeline are released.

**Keywords:** large-scale structure of Universe – supernovae: general – fast radio bursts – methods: statistical

## 1. Novelty vs prior work
| Study | Voids | $z$ | SNe | FRBs | Joint same geometry? |
|---|---|---|---|---|---|
| Wang ApJ997 | lumped $M_K$ spheres | <0.05 | Sternberg heterogeneous | 0 in voids | No |
| Tsaprazi MNRAS21 | 2M++ density + T-web | <0.036 | BTS+CLU 498Ia+782CC | — | No |
| Aubert A&A694 | SDSS DR7 + Voronoi | DR2 Ia | ZTF DR2 Ia stretch only | — | No |
| Sharma 2026 | BOSS VIDE $z0.2$–$0.7$ | 0.2–0.7 | — | 3455 CHIME stacking | No SNe |
| This work | VAST VoidFinder DR7 $R_{\max}$ | <0.114 | BTS flux-limited + z-slice | Cat1 forecast only | SN main + FRB forecast, same spheres |

## 2. Data (frozen)
VAST Zenodo 7406035 v1.3.0 (20.4MB, CC-BY-4.0, https://zenodo.org/records/7406035), code v1.7.9. BTS Explorer query quality+y purity+y (https://sites.astro.caltech.edu/ztf/bts/explorer.php, 2026-09-30, 7143 rows, 7143→6873 in $z$). CHIME Cat1 Vizier J/ApJS/257/59/table2 (600 bursts, 524 sources: 506 one-off + 94 bursts from 18 repeaters). NSA v1_0_1 downloading for parent masses. Planck18 throughout.

## 3. Methods (fixed per review)
Void XYZ recomputed from RA,Dec,$D$ ($D$ in Mpc $=$ Mpc/h$/0.7$) in Planck18 frame; stored-XYZ residual median $0.00$ Mpc $\ll R_{\mathrm{v}}$, scale verified. $R_{\mathrm{v}}\equiv R_{\max}$ core spheres (median 17.0 Mpc); $R_{\mathrm{eff}}$ (22.1 Mpc) as systematics branch. Membership via `query_ball_point` ($r<1.2\max R_{\mathrm{v}}$), $r/R_{\max}=\min$, host $=\arg\min$, shell renamed from edge (VAST file edge flag clash). Mask: `NSA_main_mask.pickle` (360,180 bool) + per-void $f_{\mathrm{obs}}$ planned; 1163/1163 void centers in-mask. Randoms $\ge10\times$ uniform comoving $\times$ mask (in progress) for $f_{\mathrm{rand}}$ and $V_{\mathrm{eff}}$. Blinded SALT2 hierarchical $\Delta$HR with PV covariance, mass-step joint fit, BTS efficiency planned; this draft reports counts + power only, no HR claim.

## 4. Results (mask-corrected, never quote 6406)
Cuts flow: 7143 → 6873 ($z$) → 2465 (mask) → 2404 void/shell/wall + 61 outside-mask-edge. Cosmology-grade 2381: void 176 / shell 263? / wall 1883 (see `data/b3_match_v2.parquet`). Void Ia 141 normal + 4 91T + 2 pec + 1 Iax excluded from cosmology (176 void total cosmo). $z\le0.06$ volume slice: 1278 total, 71 void — meets 70–100 MNRAS minimum. Threshold grid (cosmo): 0.6:59 / 0.7:103 / 0.8:176 / 0.9:291 / 1.0:438 / 1.1:613. $z_{\min}$ stability: 0.01:180 / 0.02:175 / 0.025:172 / 0.03:168 void (cosmo+mask). Void median $z$ 0.069 vs wall 0.057 confirms Malmquist gradient — motivates $z$-matched test. Figures: Fig.1 footprint with mask context, Fig.2 $r/R_{\max}$ with 0.8/1.0 lines.

FRB error budget (Cat1): per-sightline $\sigma_{\mathrm{tot}}\sim150$–$200$ (IGM 80–120, host 80–150, MW residual 20–40 with NE2001$-$YMW16 median $-3.3$ std $37.5$ p90 $68.1$ from CSV, halo 20). Void signal $\sim10$–$30$. With $N_{\mathrm{void}}\sim6$ (40 localized $\times15\%$) $\sigma_{\bar{}}\sim61$, 95% limit $\gtrsim120$. $\delta$DM$<50$ needs $N_{\mathrm{void}}\gtrsim36$ → $N_{\mathrm{loc}}\gtrsim200$. Unlocalized $z<0.114$ needs $\sim900$ void sightlines, have $\sim30$–$50$. Forecast panel planned.

## 5. Systematics budget (skeleton, mmag)
PV ($z_{\min}$ + reconstruction TBD), BBC/Malmquist (efficiency TBD), mass-step residual (NSA masses TBD), calibration, threshold ($0.7$–$1.1$ done for counts), subtype (done: 8 void contaminants removed). Total TBD before unblinding. Nulls: RA-shuffled within exposure footprint + block bootstrap over voids (nested FRB resampling) + mock injection 0.00/0.05 mag planned.

## Data Availability
VAST Zenodo DOI above; BTS query URL+date + CSV hash; CHIME Vizier DOI; derived `b3_match_v2.parquet` + scripts tag v0.2 on GitHub on acceptance (reviewer token on request); NSA on completion. No “private repo on acceptance” — public tag at submission.

## References (full in LaTeX build)
Douglass+23 ApJS265,7; Zenodo 7406035; El-Ad&Piran97; Hoyle&Vogeley02; Fremling+19; Perley+20 ApJ904,35; Rehemtulla+24; BTS Explorer; CHIME Cat1 ApJS257,59; CHIME Cat2 2026 arXiv:2601.09399; Outriggers 2504.05192; Wang ApJ997 arXiv:2511.15401; Sharma 2026 arXiv:2605.01994; Tsaprazi MNRAS21 arXiv:2109.02651; Aubert A&A694 arXiv:2406.11680; Walker A&A683 arXiv:2309.08268; Macquart Nature581; James MNRAS; Planck18; SALT2; NE2001/YMW16; emcee; healpy.
