# AHRT Stage A — local verification receipt (2026-11-04, Hunter's Mac)

Independent re-verification of the scaffold's claims after install to ~/Developer/AHRT.

## What was re-run locally (not inherited from the archive's report)

| Item | Command | Result |
|---|---|---|
| Full test suite | `.venv/bin/python -m unittest discover -s tests` | 9/9 OK before fix; 11/11 OK after (2 added regressions) |
| Klein-Nishina G1 | `ahrt validate-kn --energy-MeV 0.661657 --samples 150000 --seed 81723` | PASS, bit-identical moments to archive report (deterministic seeded) |
| XCOM G0/G1 | `ahrt validate-xcom --table data/xcom_water.csv --manifest data/xcom_water_manifest.json --energy-MeV 0.661657` | PASS with REAL NIST data |
| Reference transport | `ahrt run --histories 100000 --seed 20261103 ...` | PASS after bug fix |

## Real XCOM acquisition (completed the archive's pending item)

- NIST XCOM compound H2O, standard grid 0.001-1.0 MeV + additional energy 0.6617 MeV (UI hidden-field flow; scripted multipart uploads rejected by NIST CGI, so acquisition done via real browser session).
- 26 rows; raw page text retained at data/original_xcom_output.txt.
- Table point at 0.6617 MeV: mu_total_without_coherent = 8.563E-02 cm2/g; interpolated at 0.661657 -> 0.0856316.
- Manifest sha256: 9c829492a9b1a468887a2d5b7fff7d0cd071d969168f26d02df690c6b6122ac3
- Gate outcome: coefficient closure 0.000% at requested energy (worst raw-row 3.7e-4 relative); Beer-Lambert sampled attenuation 0.14% relative error vs 1% acceptance; mean free path 11.678 cm.

## Bug found and fixed (archive's own tests did not cover it)

- `klein_nishina.rotate_direction` line 119 referenced undefined name `d`; the end-to-end `ahrt run` crashed with NameError on the first Compton scatter. Tests passed because no test exercised rotate_direction through transport.
- Fix: `d = (dx, dy, dz)` before the generator (commit this).
- Regression: tests/test_rotate_direction_regression.py — two cases (axis-aligned, off-axis) asserting unit norm and exact scattering-angle preservation; failed without the fix, pass with it.

## Reference-run tallies (100k histories, Cs-137 661.657 keV, 3mm-radius 6mm-long water cylinder, rho=1.0)

- photon interactions 4108 (4107 Compton, 1 photoelectric); 99,999 escapes — the single PE photon was absorbed in-domain.
- direct_cell_photon_interactions = 0 — quantitatively consistent with the nested-electron-source argument (expected ~1e-10 per photon for cell cross-section x interaction probability).
- Energy accounting closes: 66165.70 incident = 65148.79 escaped + 1016.25 Compton recoil source + 0.66 PE local.
- 4107 weighted Compton recoil-electron phase-space records exported (compton_electron_sources.csv) as the Stage B interface.

## Boundaries unchanged

Stage A remains primary photon transport + Compton source phase space. No electron transport, no dose, no radiolysis, no Arc response, no Rayleigh/pair. The G5 damage gate stays disabled.
