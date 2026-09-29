# VOIDS — Void Observations of Intergalactic Dispersion and Supernovae

## Principal Investigator
Independent researcher

## Abstract
I propose a joint test of whether Type Ia supernova distances and fast radio burst dispersion measures differ inside versus outside low-redshift cosmic voids. I will use the same geometric void sample throughout: VAST VoidFinder SDSS DR7 (1163 voids, z < 0.114), crossed with ZTF Bright Transient Survey supernovae and CHIME FRB Catalog 2. All work will use public archives only. The result will be publishable as a detection or as a quantitative limit, with an open catalog and pipeline.

## Background
Cosmic voids fill most of the volume but host few galaxies. They are a clean lab for progenitor and baryon systematics: supernova rates and standardized brightness should track age and metallicity, while FRB dispersion should track the ionized baryon deficit along void sightlines.

VAST provides the definitive low-z geometric catalog I need: Douglass et al. 2023 (ApJS 265, 7) built VoidFinder, VIDE and REVOLVER catalogs in the same SDSS DR7 volume with masks and mocks. ZTF BTS is the largest flux-limited supernova survey (Perley et al. 2020; Fremling et al. 2019; BTSbot, Rehemtulla et al. 2024), now over 11,000 classified transients with a public explorer. CHIME Catalog 2 (2026) gives 4539 bursts from 3641 sources with uniform reprocessing and exposure maps.

Recent work brackets but does not close this test. Wang et al. (ApJ 997, 121) surveyed transients in void galaxies without BTS, CHIME Cat2, or VoidFinder. Tsaprazi et al. (MNRAS 2021) crossed BTS with density fields at z < 0.036 without voids-as-spheres or FRBs. Aubert et al. (A&A 694, A7) studied ZTF DR2 Ia stretch in voids without colors-independent Hubble residuals, core-collapse rates, or FRBs. Sharma et al. (2026) stacked 3455 CHIME sightlines on BOSS VIDE voids at z 0.2–0.7 for IGM physics without supernovae or DR7 geometry. No published paper jointly uses VoidFinder DR7 + BTS + CHIME Cat2 in the same low-z volumes.

## Research Questions
1. Is the supernova void fraction different from the galaxy void fraction? Does the core-collapse to Ia ratio change in voids?
2. Is there a void–wall offset in standardized Ia Hubble residuals?
3. Is there a void deficit in FRB extragalactic dispersion versus impact parameter?
4. What do the combined limits imply for supernova systematics and the void ionized baryon fraction?

## Data
I will freeze versions at the start and record hashes:
- VAST SDSS DR7 VoidFinder Planck18 holes + maximals + mask, Zenodo 7406035 v1.3.0, code v1.7.9.
- ZTF BTS Explorer CSV with quality and purity cuts, mag ≤ 18.5, plus query URL and date.
- CHIME Catalog 1 CSV via the open-data package for positions and DMs; Catalog 2 exposure information for selection control. I will not bulk-download raw intensity or baseband products.
- NASA-Sloan Atlas v1_0_1 for the volume-limited parent galaxy sample.
- TNS for names, types and redshifts where needed.

Total frozen subset is about 2.7 GB. No proprietary data are required.

## Methodology
I adopt Planck18 distances and convert void radii to Mpc. Transients and void centers go to comoving Cartesian coordinates. I use a KD-tree point-in-sphere match, computing r/Rv for each transient. I define void as r/Rv < 0.8, edge 0.8–1.0 excluded from the primary test, wall as > 1.0 for all voids. Primary volume is z 0.01–0.114 with SDSS, ZTF and CHIME footprint cuts and observed volume fraction > 0.8 per void.

Selection control is central. I will generate uniform comoving randoms in the SDSS mask to predict the expected void fraction, compare transients to parent NSA galaxies, repeat the supernova test in the highly complete z < 0.06 slice, model detection probability versus redshift and extinction, and use Ia versus core-collapse as a same-survey control. For FRBs I will weight by published exposure, use declination-matched controls, rotate void positions for null tests, and bootstrap over voids.

For supernovae I take published SALT2 parameters where available, fit the global standardization blind to environment, then test Hubble residuals void versus wall with Welch, KS, covariate regression and a hierarchical Bayesian offset model. With a few hundred void Ia I am sensitive to ~0.02–0.03 mag. For FRBs I decompose DM into Milky Way, halo, IGM and host terms, stack extragalactic DM residuals versus transverse impact b/Rv for localized events, and report an upper limit if unlocalized scatter dominates. For rates I compare counts per effective void versus wall volume with exact binomial and Poisson tests, split by Ia, Type II and stripped-envelope.

Required outputs are footprint and r/Rv diagnostic plots, Hubble-residual void–wall comparison, DM versus b/Rv stacking, rate-ratio plots, a systematics grid over thresholds and subsamples, a cuts-flow table and a statistics summary table.

A two-week pilot will deliver the cross-match, randoms test, Hubble-residual comparison and rate ratio. The full analysis adds SALT2 matching, hierarchical limits, exposure-weighted FRB stacking and the systematics grid.

## Timeline
Weeks 1–2: environment setup, catalog ingest, pilot counts and go/no-go on void numbers. Weeks 3–4: supernova standardization and FRB decomposition. Weeks 5–6: masks, randoms and exposure controls. Week 7: systematics and mock recovery. Weeks 8–10: draft, reproducibility rerun, Zenodo and GitHub freeze, submission.

Minimum for a full article is roughly 70–100 void Ia with several hundred field controls, or a joint 80 supernova plus 25 FRB sample. Below that I will target a note or catalog paper.

## Deliverables
- MNRAS paper (fallbacks A&A or OJAp), submitted with same-day arXiv preprint.
- Public cross-matched catalog with flags, randoms and exposure maps, Zenodo DOI, CC-BY-4.0.
- Versioned pipeline, MIT licensed, with requirements file and rerun notebook.

## References
Douglass et al. 2023, ApJS 265, 7; VAST Zenodo 7406035; El-Ad & Piran 1997; Hoyle & Vogeley 2002; Fremling et al. 2019; Perley et al. 2020; Rehemtulla et al. 2024; BTS Explorer; CHIME Cat1 ApJS 257, 59; CHIME Cat2 2026; CHIME Outriggers design; Lanman et al. KKO; Wang et al. ApJ 997, 121; Sharma et al. 2026; Tsaprazi et al. MNRAS 2021; Aubert et al. A&A 694, A7; Walker et al. A&A 683, A71; Macquart et al. 2020 Nature; Rafiei-Ravandi et al. ApJ; Shin et al. ApJ; Pleunis et al. ApJ; Sutter et al. 2012; Chung et al. ApJ; James et al. MNRAS; DESIVAST Rincon et al. ApJ 982, 38.
