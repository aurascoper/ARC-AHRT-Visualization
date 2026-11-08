# W-A-03 uncertainty report — species/time reconciliation
Date: 2026-11-04. Frozen baselines (Ref A/B, W-A-ARC-441) unchanged; comparison only.

## Experimental anchors (transcribed from article text)
| Species | t (ps) | G (molec/100 eV) | ±1sd | Source | Medium |
|---|---|---|---|---|---|
| OH· | 7 | 4.8 | 0.2 | El Omar 2011 | H2O |
| e⁻_aq | 7 | 4.2 | 0.2 | Wang 2018 | D2O (isotope label preserved per disposition) |

## Anchor comparison (Ref B = 999.999 keV, pure water, frozen run)
| Species | Ref B G @7 ps | sim 1sd | combined σ | residual | verdict |
|---|---|---|---|---|---|
| OH· | 4.8237 | 0.0029 | 0.200022 | +0.119σ | PASS |
| e⁻_aq | 4.1772 | 0.0026 | 0.200017 | −0.114σ | PASS |

Both anchors sit well inside one experimental standard deviation. The paper was tuned to these
anchors (RMS displacement scan driven by Wang 2018), so anchor agreement demonstrates
configuration-equivalent reproduction of the paper's setup — not independent predictive validation.

## Figure-curve consistency
comparison_results.csv holds the full matched grid. Paper curves are digitized figure-data tables
provided with the eScholarship posting. Series identities: panel content implies top-row
pure-water panels contain simulated curves; bottom-row Series 2/3 are UNLABELED and are recorded
with condition UNKNOWN (candidate: scavenger-conditional simulations or experimental series; not
used as experimental evidence).

## H2O2 and H·
Per the paper and the disposition: pure-water curves are NOT expected to reproduce the scavenger
experiments; the relevant comparison is V4 (not started).

## Limitations
- Digitized paper curves carry no per-point uncertainties (curve-to-curve consistency, not
  experimental agreement).
- Paper chemical stage ran to 10 µs; my runs stop at 1 µs; comparisons restricted to the shared
  plotted domain (1 ps – 1 µs); no extrapolation performed.
- Full time-domain agreement vs experiment remains pending until measured tables from the eight
  cited measurement studies are extracted from their original sources.

## Verdict on accepted anchors
W-A-03 (limited scope): PASS at the 7 ps OH· and e⁻_aq anchors. W-A-03 full scope: STILL PENDING
(pointwise experimental tables beyond the two text anchors are not yet available).
