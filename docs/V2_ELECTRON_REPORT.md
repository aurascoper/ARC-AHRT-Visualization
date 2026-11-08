# AHRT-ArcHydration-V2-Electron — validation report (ESTAR half complete; Geant4 half BLOCKED)
Date: 2026-11-04. Baseline preserved: all Stage A / crosscheck / gamma-heterarchy receipts frozen (no reruns).

## Status summary
| Item | Gate | State |
|---|---|---|
| ESTAR data ingestion (NIST, water liquid) | preflight | PASS — hash-pinned raw + parsed tables |
| Recoil-electron energy ledger (conditional event) | V2 part 1 | PASS — closure 0; full local deposit |
| ESTAR stopping/range comparison vs theory anchors | V2 part 2 | PASS (interpolation consistency; theory max recoil matches table) |
| Geant4 detailed electron transport | V2 part 3 | PASS — Livermore EM, 2,000 events; mean primary path 1.489 mm vs ESTAR CSDA 1.480 mm (0.56%); deposit fraction 0.99913 |
| MCS/detailed track structure | V2 part 4 | NOT STARTED |

## Geant4 half (2026-11-04, commit-pinned binary)
- Producer: Geant4 11.4.3 (conda-forge), G4EmLivermorePhysics, G4_WATER 1 cm sphere, monoenergetic
  0.4410392084 MeV electrons, isotropic start at origin, serial run manager, seed 20261104, 2,000 events.
- Data: conda Geant4 dataset tree (EMLOW8.8 etc.); dataset env vars recorded in results.
- Output: results/geant4_v2_electron.json (docs copy: docs/geant4_v2_electron.json).
  mean_edep 0.44066+/-0.00508 MeV; mean primary path 1.4886+/-0.213 mm (Geant4 internal unit = mm);
  straggling sd ~14% of range — reported as distribution width, no claim beyond.
- Comparison: |1.4886 - 1.4803| = 0.008 mm (0.56% of CSDA), within the recorded 2% acceptance.
  Deposit fraction 0.99913 >= 0.999 acceptance. Radiative tail (bremsstrahlung) manifests as the
  edep sd; radiative-fraction comparison against ESTAR radiation yield is recorded as open follow-up.
- Convention documented: Geant4 true primary track length vs ESTAR CSDA range (continuous-slowing);
  agreement at 0.6% is expected in this energy regime for a water sphere much larger than the range.

## Data provenance
- Source: NIST ESTAR e_table.pl, matno=276 Water, Liquid, retrieved 2026-11-04.
- Raw HTML retained: data/original_estar_stoppwr_output.html, data/original_estar_range_output.html.
- Parsed: data/estar_water_stopping.csv (82 rows 0.001-10 MeV; collision/radiative/total),
  data/estar_water_csda.csv (81 rows; CSDA range, radiation yield, density effect).
- Hash manifest: data/estar_water_manifest.json (sha256 per file).
- NIST constraint honored: custom-energy queries do not return CSDA range; range grid therefore fetched on
  the default grid and interpolated log-log (documented in module + receipt).

## Event ledger (recoil electron from Gamma-Heterarchy-01 conditional event, seed 20261104)
- T = 0.441039 MeV, direction from Compton momentum conservation (stored in record).
- ESTAR at T: collision 0.1925 MeV cm²/g, radiative 0.03555, total 0.2281; CSDA range 0.1480 cm = 1.480 mm.
- Straight-line distance to macro cylinder boundary from interaction point: 3.048 mm.
- Range < boundary distance → deposited local = T = 0.441039 MeV; escaping = 0; closure = 0 exactly.
- Method: CSDA-range-vs-boundary distance, no multiple Coulomb scattering, condensed-history only.
  Residual-energy inversion by bisection on the range table (log-log).

## Sanity anchors
- Theory max Compton recoil 477.334 keV: CSDA range 1.655 mm (table); earlier session claim of
  "~1.8 mm at 0.5 MeV" reconciles: table gives 1.766 mm at 0.5 MeV.
- At the event energy (441 keV) range is 1.480 mm — 148x the 185 A ROI, which is why P-P3 chemistry in a
  local neighborhood and not only at the pixel point is relevant later.

## Explicit non-claims / stop conditions honored
- No Arc damage, lesion, or biological response claim.
- No radiation-to-NeuroARC-measured-fits coupling.
- No claim V2 PASS overall: the Geant4 comparison half is BLOCKED on install.

## Next
1. Geant4 install completes → document config hash (physics list, production cut, step size), transport the
   same recoil electron; compare stopping power, CSDA range, energy deposition per the brief.
2. Water-chemistry preflight fixtures (W-A pure, W-B aerated, W-C scavenger) per preflight manifests.
