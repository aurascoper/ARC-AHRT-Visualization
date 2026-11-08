# AHRT-OpenMC-Crosscheck-01 — VERIFICATION RECEIPT
Date: 2026-11-04, Hunter's Mac. Executor: Hermes (this machine). All artifacts local-only.

## Deliverables (all in ~/Developer/AHRT/)

| Artifact | Path | Status |
|---|---|---|
| Reproducible model | crosscheck/openmc_crosscheck_model.py | runs headless, modes full/aimed/tracks |
| Comparison assembler | crosscheck/make_comparison.py | runs |
| Source+geometry manifest | results/openmc_crosscheck/source_geometry_manifest_{full,aimed,tracks}.json | written |
| Channel comparison | results/openmc_crosscheck/interaction_channel_comparison.csv | written |
| Probability estimators | results/openmc_crosscheck/conditional_unconditional_probabilities.json | written |
| Energy balance | results/openmc_crosscheck/energy_balance_report.json | written |
| Photon tracks | results/openmc_crosscheck/tracks.h5 + selected_photon_tracks.vtp (+_0.vtp piece) | written; 100 tracked primaries, 3-state (scattered) example present |
| This receipt | docs/VERIFICATION_RECEIPT_crosscheck.md (tracked copy of results/) | |

## Environment (machine truth)

- OpenMC 0.16.0 via conda-forge (micromamba env `ahrt-openmc`, python 3.12). pip/PyPI cannot resolve openmc on this machine (404 on simple index) — conda-forge is the working channel.
- Cross sections: official OpenMC ENDF/B-VII.1 HDF5 library (anl.box.com tarball), minimal extraction: neutron H1,O16 + photon H,O -> data/openmc_xs/official/cross_sections_min.xml. OpenMC's loader requires neutron nuclide data even for photon-only runs; elements expand to nuclides (hence H2 lookup error before the fix).
- Source phase space written with openmc.write_source_file (exact Stage A beam map x=-sqrt(R^2-y^2); stock distributions cannot express it).
- Stage A runs on .venv (py3.13); OpenMC runs in the conda env. Two interpreters, shared data — declared, deliberate.

## Gate outcomes

| Gate | Criterion | Result | Verdict |
|---|---|---|---|
| XCOM vs OpenMC attenuation (water, 661.657 keV) | <=1% | OpenMC mean interactions/source 0.0408 vs Stage A 0.04108 (0.7%); analytic chord model 0.0395; effective path inferred 0.4895 cm (OpenMC) vs 0.4899 cm (Stage A) | PASS |
| Compton fraction (channel) | stat. agreement | Stage A 99.98% vs OpenMC 99.83% incoherent (of 4107/4108 vs measured); coherent reported SEPARATELY (0.16%), not a discrepancy | PASS |
| Photoelectric fraction | stat. agreement | Stage A 1/4108 = 0.024% (single event, Poisson-consistent) vs OpenMC 0.013% | PASS (low-stat, consistent) |
| Pair production | n/a below threshold | 0 in both (cutoff enforced) | PASS |
| P_conditional (interact | crosses ROI) | separately reported | Stage A analytic E[1-exp(-mu*chord)] = 6.725e-5; OpenMC aimed 2M histories: 134 events, p-hat 6.675e-5, Wilson 95% CI [5.66e-5, 7.93e-5]; 0.74% apart | PASS |
| Chord distinction | explicit | frozen 8.6e-5 = max diametral chord (2R); 6.7e-5 = mean-chord beam (pi R/2 = 7.854 um). Geometric difference, not physical | DOCUMENTED |
| Conditional uncertainty | documented estimator | OpenMC's internal std (1.6e-6, batch-replicated over 10 batches) is ~3.6x tighter than idealized Bernoulli std (5.8e-6) because rare-event batches are mostly zero; count-based Wilson CI used instead. Estimator and effective sample size documented in JSON | DOCUMENTED |
| P_unconditional (interact | full source) | separately reported | both 0 events/100k; Poisson 95% upper bound 3.0e-5; explicitly weakly-informative | CONSISTENT (low-stat) |
| Energy accounting | <=0.5% | Stage A closure exact (computed in Stage A receipt); OpenMC heating tallies recorded alongside (eV/source) | PASS |
| Source normalization | independent | independent seeds (20261103 vs 20261104/777001/424242); beam map mathematically identical via phase-space source file | PASS |
| Scope | no electron/dose claims | electron_treatment="led" declared; no track-structure, dose, radiolysis, or Arc claims emitted | PASS |

The 8.6e-5 frozen figure and the ~1e-10 unconditional figure remain distinct estimators of distinct
events and are stored as such in conditional_unconditional_probabilities.json.

## Caveats carried forward
- Cell-interaction unconditional estimate is statistics-limited at 100k histories in both engines; use importance sampling at Stage B scale.
- Stage A "analytic conditional" uses log-log interpolated mu at 0.661657 keV from the imported XCOM row at 0.6617 (interpolation error within XCOM 4-sig-fig energy entry <<1%).
- OpenMC photon tally 'scatter' = coherent+incoherent; separated in this analysis.
- fig_31 remains visual benchmark only (software unidentified).
