# Joint Constraints on Supernova Standardization and Fast Radio Burst Dispersion in Cosmic Voids: A ZTF BTS × CHIME × VAST SDSS DR7 Test
### VOIDS Paper Draft v0.1 — 30 September 2026

**Author:** Himel Das
**Affiliation:** Government Azizul Haque College, Bogura, Bangladesh
**ORCID:** 0009-0003-7698-1255

## Abstract
We present a joint test of transient environment in low-redshift cosmic voids using a single geometric void definition. We cross-match VAST VoidFinder SDSS DR7 (1163 voids, $z<0.114$) with ZTF Bright Transient Survey supernovae (7143 objects, 6873 in $z=0.01$–$0.114$) and CHIME FRB Catalog 1 (600 bursts, Catalog 2 exposure in progress). Using comoving KD-tree point-in-sphere matching with $r/R_{\mathrm{v}}<0.8$ void, $0.8$–$1.0$ edge excluded, we find 188 void supernovae (144 Type Ia, 18 Type II), 279 edge, and 6406 wall. This exceeds the minimum for a $0.03$ mag Hubble-residual constraint. Full SALT2 standardization, selection-corrected rates, and FRB dispersion stacking are in progress. Catalog, randoms, and pipeline will be released.

**Keywords:** large-scale structure — supernovae — fast radio bursts — methods: statistical

## 1. Introduction
Voids fill most of the volume. Supernova luminosities trace age/metallicity; FRB dispersion traces ionized baryons. Prior work uses heterogeneous voids or single probes. This work fixes geometry for both probes.

## 2. Data
VAST Zenodo 7406035 v1.3.0 Planck18 maximals/holes + mask. BTS Explorer quality+y purity+y CSV dated 2026-09-30 (7143 rows). CHIME Cat1 Vizier J/ApJS/257/59 (600 rows). NSA v1_0_1 downloading for parent control. Cosmology Planck18.

## 3. Methods
XYZ via Planck18 comoving distance. VAST XYZ converted Mpc/h ($h=0.7$) to Mpc. cKDTree, $r_{\max}=\max R_{\mathrm{v}}+5$ Mpc, $r/R_{\mathrm{v}}=\min$. Volume $z=0.01$–$0.114$, footprint cuts pending full mask. Randoms + NSA parent + $z<0.06$ complete slice + Ia/CC control + CHIME exposure weighting + RA-rotated nulls + void bootstraps planned.

## 4. Preliminary Results
1163 voids, median $R_{\mathrm{v}}=17.0$ Mpc. BTS in volume 6873. Environment: wall 6406, edge 279, void 188. Void composition dominated by Ia (144), consistent with flux-limited selection. Naive SDSS footprint expectation met. Figures: `plots/rRv_hist.png`, `plots/footprint.png`. Table `data/b3_match.parquet` (250K) frozen.

Power: $\sigma_{\mathrm{int}}\sim0.15$ mag, $N_{\mathrm{v,Ia}}=144$, $N_{\mathrm{w}}\sim5000$ gives sensitivity to $\sim0.025$ mag void offset at $2\sigma$. FRB arm requires localized redshifts; Cat1-only gives limit demonstration.

## 5. Systematics and Next Steps
Mask $f_{\mathrm{obs}}>0.8$, BTS efficiency vs $z$, SALT2 refit blinded, hierarchical $\Delta$HR, DM decomposition NE2001/YMW16 both, rate ratios Ia/II/Ibc, threshold grid $0.7$–$1.1$, mock recovery, DESI-Y3 forecast paragraph.

## Data Availability
VAST Zenodo DOI, BTS Explorer query URL+date, CHIME Vizier, code tag v0.1 on private repo (public on acceptance).

## References
Douglass+2023, Perley+2020, Fremling+2019, Rehemtulla+2024, CHIME Cat1, CHIME Cat2 2026, Wang ApJ997, Tsaprazi MNRAS21, Aubert A&A694, Sharma26, Walker A&A683, Macquart20.
