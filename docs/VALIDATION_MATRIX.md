# AHRT Stage A validation matrix

| Gate | Test | Implementation | Status in this scaffold |
|---|---|---|---|
| G0 | Runtime input provenance | SHA-256 manifest + NIST source and units checks | Implemented; needs real NIST export |
| G0 | Dimensional integrity | CSV column names encode MeV and cm2/g; geometry uses cm | Implemented |
| G1 | XCOM partial-to-total closure | Compare sum of process coefficients to both imported totals | Implemented; needs real NIST export |
| G1 | XCOM attenuation sampling | Compare Monte Carlo free-path interaction probability to Beer-Lambert prediction from imported coefficient | Implemented; needs real NIST export |
| G1 | Klein-Nishina energy-angle | Sampled moments vs numerical quadrature; per-event energy-angle identity | Implemented |
| G2 | Electron stopping/range vs ESTAR | Not yet implemented | Blocked on electron transport stage |
| G3 | Water-radiolysis yields | Not yet implemented | Blocked on Geant4-DNA reference run |
| G4 | Cell ROI absorbed dose | Not yet implemented | Requires transported electrons and energy deposition scoring |
| G5 | Arc structural damage/response | Disabled | Intentionally disabled |

The synthetic CSV under `tests/fixtures/` is only a software test fixture and
must never be used to generate physical results.
