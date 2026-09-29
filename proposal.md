# VOIDS — Void Observations of Intergalactic Dispersion and Supernovae
### Research Proposal

**Author:** Himel Das
**Affiliation:** Government Azizul Haque College, Bogura, Bangladesh
**ORCID:** 0009-0003-7698-1255
**Date:** 30 September 2026
**Contact:** Independent researcher
**Repository:** Private — `hello2himel/VOIDS`, branch `main`

---

## Abstract
I propose a joint test of whether Type Ia supernova distances and fast radio burst (FRB) dispersion measures differ inside versus outside low-redshift cosmic voids. I will use a single geometric void sample throughout — VAST VoidFinder SDSS DR7 (1163 voids, $z < 0.114$) — crossed with ZTF Bright Transient Survey (BTS) supernovae and CHIME FRB Catalog 2. All work uses public archives only and will be executed on a home workstation. The result will be publishable as a detection or as a quantitative limit, accompanied by an open cross-matched catalog and versioned pipeline.

**Keywords:** large-scale structure of Universe — cosmological parameters — supernovae: general — fast radio bursts — methods: statistical

## 1. Background
Cosmic voids fill most of the cosmic volume but host few galaxies. They provide a clean laboratory for progenitor and baryon systematics: supernova rates and standardized brightness trace stellar age and metallicity, while FRB dispersion traces the ionized baryon deficit along void sightlines.

VAST provides the definitive low-redshift geometric catalog for this work: Douglass et al. (2023, ApJS 265, 7) constructed VoidFinder, VIDE, and REVOLVER catalogs in the same SDSS DR7 volume with masks and mocks. ZTF BTS is the largest flux-limited supernova survey (Perley et al. 2020; Fremling et al. 2019; BTSbot, Rehemtulla et al. 2024), now exceeding 11,000 classified transients with a public explorer. CHIME Catalog 2 (2026) provides 4539 bursts from 3641 sources with uniform reprocessing and exposure maps.

Recent work brackets but does not close this test. Wang et al. (ApJ 997, 121) surveyed transients in void galaxies without BTS, CHIME Catalog 2, or VoidFinder. Tsaprazi et al. (MNRAS 2021) crossed BTS with Bayesian density fields at $z < 0.036$ without voids-as-spheres or FRBs. Aubert et al. (A&A 694, A7) studied ZTF DR2 Ia stretch in voids without colour-independent Hubble residuals, core-collapse rates, or FRBs. Sharma et al. (2026) stacked 3455 CHIME sightlines on BOSS VIDE voids at $z = 0.2$–$0.7$ for intergalactic medium physics without supernovae or DR7 geometry. No published paper jointly uses VoidFinder DR7 $+$ BTS $+$ CHIME Catalog 2 in the same low-redshift volumes.

## 2. Research Questions
1. Is the supernova void fraction different from the galaxy void fraction? Does the core-collapse to Type Ia ratio change in voids?
2. Is there a void–wall offset in standardized Type Ia Hubble residuals?
3. Is there a void deficit in FRB extragalactic dispersion as a function of impact parameter?
4. What do the combined limits imply for supernova systematics and the void ionized baryon fraction?

## 3. Data
I will freeze versions at the start and record hashes. Total frozen subset is approximately 2.7 GB. No proprietary data are required.

| Dataset | Version | Access |
|---|---|---|
| VAST SDSS DR7 VoidFinder Planck18 holes, maximals, mask | Zenodo 7406035 v1.3.0; code v1.7.9 | Public, CC-BY-4.0 |
| ZTF BTS Explorer | Quality and purity cuts, mag $\leq 18.5$, query URL and date recorded | Public |
| CHIME Catalog 1 positions and DMs; Catalog 2 exposure | Catalog 1 CSV; Catalog 2 exposure maps | Public; raw intensity/baseband excluded |
| NASA-Sloan Atlas | v1_0_1 | Public |
| Transient Name Server | Names, types, redshifts | Public; free account |

## 4. Methodology
I adopt Planck18 distances and convert void radii to Mpc. Transients and void centers are placed in comoving Cartesian coordinates. I use a KD-tree point-in-sphere match, computing $r/R_{\mathrm{v}}$ for each transient. I define void as $r/R_{\mathrm{v}} < 0.8$, edge $0.8$–$1.0$ excluded from the primary test, and wall as $> 1.0$ for all voids. The primary volume is $z = 0.01$–$0.114$ with SDSS, ZTF, and CHIME footprint cuts and observed volume fraction $> 0.8$ per void.

Selection control is central. I will generate uniform comoving randoms in the SDSS mask to predict the expected void fraction, compare transients to parent NASA-Sloan Atlas galaxies, repeat the supernova test in the highly complete $z < 0.06$ slice, model detection probability versus redshift and extinction, and use Type Ia versus core-collapse as a same-survey control. For FRBs I will weight by published exposure, use declination-matched controls, rotate void positions for null tests, and bootstrap over voids.

For supernovae I take published SALT2 parameters where available, fit the global standardization blind to environment, then test Hubble residuals void versus wall with Welch, Kolmogorov–Smirnov, covariate regression, and a hierarchical Bayesian offset model. With a few hundred void Type Ia events I am sensitive to $\sim 0.02$–$0.03$ mag. For FRBs I decompose dispersion measure into Milky Way, halo, intergalactic, and host terms, stack extragalactic residuals versus transverse impact $b/R_{\mathrm{v}}$ for localized events, and report an upper limit if unlocalized scatter dominates. For rates I compare counts per effective void versus wall volume with exact binomial and Poisson tests, split by Type Ia, Type II, and stripped-envelope.

Required outputs are footprint and $r/R_{\mathrm{v}}$ diagnostic plots, a Hubble-residual void–wall comparison, dispersion versus $b/R_{\mathrm{v}}$ stacking, rate-ratio plots, a systematics grid over thresholds and subsamples, a cuts-flow table, and a statistics summary table. A two-week pilot delivers the cross-match, randoms test, Hubble-residual comparison, and rate ratio. The full analysis adds SALT2 matching, hierarchical limits, exposure-weighted FRB stacking, and the systematics grid.

## 5. Timeline
Weeks 1–2: environment setup, catalog ingest, pilot counts and go/no-go on void numbers. Weeks 3–4: supernova standardization and FRB decomposition. Weeks 5–6: masks, randoms, and exposure controls. Week 7: systematics and mock recovery. Weeks 8–10: draft, reproducibility rerun, Zenodo and GitHub freeze, submission. Minimum for a full article is approximately 70–100 void Type Ia with several hundred field controls, or a joint 80-supernova plus 25-FRB sample. Below that I will target a note or catalog paper.

## 6. Deliverables
1. MNRAS paper (fallbacks A&A or OJAp), submitted with same-day arXiv preprint.
2. Public cross-matched catalog with flags, randoms, and exposure maps, Zenodo DOI, CC-BY-4.0.
3. Versioned pipeline, MIT licensed, with requirements file and rerun notebook.

## References
Douglass et al. 2023, ApJS 265, 7; VAST Zenodo 7406035; El-Ad & Piran 1997; Hoyle & Vogeley 2002; Fremling et al. 2019; Perley et al. 2020; Rehemtulla et al. 2024; BTS Explorer; CHIME Catalog 1, ApJS 257, 59; CHIME Catalog 2, 2026; CHIME Outriggers design; Lanman et al. KKO; Wang et al. ApJ 997, 121; Sharma et al. 2026; Tsaprazi et al. MNRAS 2021; Aubert et al. A&A 694, A7; Walker et al. A&A 683, A71; Macquart et al. 2020 Nature; Rafiei-Ravandi et al. ApJ; Shin et al. ApJ; Pleunis et al. ApJ; Sutter et al. 2012; Chung et al. ApJ; James et al. MNRAS; DESIVAST Rincon et al. ApJ 982, 38.
