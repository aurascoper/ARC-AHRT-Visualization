# V4 — AHRT Scavenger-Response Assay (ASRA)
Prepared 2026-11-04. Gates preserved: V1/V2 VERIFIED, V3 W-A executed (W-A-03 limited-scope PASS at the two
7 ps anchors; full time-domain experimental agreement PENDING). V5 remains NOT STARTED and DISABLED.

## Assay statement

At a fixed radiation field and deposited-energy distribution, changing hydroxyl-radical scavenging
capacity changes time-dependent radiolysis products in ways predicted by the Geant4-DNA/TOPAS-nBio
chemistry model. AHRT predicts transport + chemical exposure; the assay tests correspondence to
experimentally observed chemical yields. No biological target is part of ASRA.

## Fixtures

| Fixture | Medium | Scavenger series | Conditions source | Simulation config |
|---|---|---|---|---|
| W-C1 | deaerated water | CH3OH {1e-3,1e-2,0.1,1,10} M + NO3- 0.025 M | Hiroki 2002 via Ramos-Mendez 2021 Table 1 (k: OH+CH3OH 9.7e8 /M/s; e_aq+NO3- 9.7e9; H+NO3- 1.4e6) | topas/v4_wc/wc1_no3_25mM_meoh_*.txt |
| W-C2 | oxygen-saturated water (0.27 mM O2) | DMSO {1e-4,1e-3,1e-2,0.1,1} M | Ramos-Mendez 2021 aerated conditions (k: OH+DMSO 7.1e9; e_aq+O2 1.9e10; H+O2 2.1e10; O2 conc 0.27e-3 M) | topas/v4_wc/wc2_o2sat_dmso_*.txt |

Deaerated and oxygen-saturated results are NEVER pooled: different chemical environments.

Simulation layer: isolated container (OpenTOPAS 4.3.0 / Geant4 11.3.2 / TOPAS-nBio 4.1.0),
TsEmDNAPhysics+TsEmDNAChemistry, background continuum scavenging (exponential-single-factor),
441.039 keV beam, seed 20261104, 3,000 histories per point.

## Endpoints

- Primary stable product: G(H2O2) time-dependent and terminal.
- Supporting radical readouts: G(OH.), G(e_aq), G(H.) time-dependent.
- Scavenging capacity sigma = k_scavenger x [B] (/s) as the concentration axis for the response curve.

## Acceptance criteria (declare BEFORE comparing to experiment)

1. Concentration-response monotonicity and expected signs:
   - W-C1: G(H2O2) terminal drops with increasing CH3OH capacity (OH scavenging suppresses
     OH+OH->H2O2 recombination); Rothwell 2026 (deaerated): yields approach zero at high capacity.
   - W-C2: non-monotonic response permitted: Rothwell 2026 shows oxygen-saturated H2O2 yields first
     RISE under moderate scavenging then decline at high capacity; methanol > DMSO under oxygenated
     conditions. Failing to capture that non-monotonicity is a model failure, not an acceptable curve.
2. Uncertainty-disciplined comparison: per-point simulation and experimental uncertainties; residuals
   normalized to combined sigma; no point with missing experimental uncertainty counts toward a pass.
3. Separate evidence sources: Ramos-Mendez 2021 (closest benchmark; plasmid part has chosen-efficiency
   caveat), Peukert 2019 (adversarial cross-check: terminal OH/H2O2 yield discrepancies, incomplete
   water recombination in the simulated window), Rothwell 2026 (nonlinearity + dose-rate invariance).
   Compare against each independently; do not average sources.
4. Chemistry provenance: reaction set, rates, units, concentrations, seed, geometry hashes in receipts.

## Level boundary

V4 validates chemistry vs chemical yields. It is not an Arc assay. V5 (experimental hydration of dArc1
capsid + irradiation readouts, with matched shams and independently established assays/acceptance)
remains a separate experimental programme. Level-A exposure interface and static_6tap manifest are
the structural boundary; no exposure->lesion->phenotype conversion.

## Known corrections applied during V4 execution (audit trail)
1. Post-review rebuild: Fricke skeleton retargeted — ReportMoleculesNamed explicitly set to
   OH^0 / e_aq^-1 / H2O2^0 / H_2^0 / H^0 (the Fe3 target of the stock example was never accepted);
   beam = e- 999.999 keV (electron G-value protocol; NOT re-labeled as 137Cs gamma);
   background-reaction products = "None" (deliberate terminal sink; products untracked);
   species identities Nitrate/Methanol/DMSO declared via Mo/ parameters (Radius/Diffusion = plumbing
   values, unused by continuum model, documented as such).
2. Time-domain correction: Fricke skeleton silently carried Sc TimeUpper = 50 s (NOT the W-A
   1 ps-1 us domain). Capped to 1e6 ps in all V4 configs; the earlier 50 s series are superseded.
3. Scavengers sweep extended toward low DMSO (1e-6..3e-5 M) to cover the Rothwell rise region;
   the initial W-C2 1e-4..1 sweep produced monotone fall — the rise must be found at lower capacity.
4. Run artifacts from the pre-review batches: FAILED (invalid configs), retained in git history,
   never accepted as V4 chemistry evidence.

## Outputs pending this run

- gv_<tag> phsp outputs (per-species time-dependent G) -> parsed response tables
- comparison vs Ramos-Mendez Table-1-paper numbers ONLY where the paper supplies a matched-condition
  numeric value; digitized figure curves reused with their provenance labels
- nonlinear-response verdict per fixture (W-C2 non-monotonicity test)
- receipts/V4_ASRA_receipt.json
