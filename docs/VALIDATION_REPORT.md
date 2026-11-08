# Scaffold validation report

Executed in the build container before packaging.

## Automated tests

`PYTHONPATH=src python -m unittest discover -s tests -v`

Result: **9 tests passed**.

The tests cover cylinder geometry, Klein-Nishina event energy-angle consistency,
Monte Carlo angular moments against numerical quadrature, XCOM CSV parsing and
hash binding, a synthetic-fixture attenuation gate, and a synthetic-fixture
transport energy-balance smoke test.

## Klein-Nishina G1 run

Command:

```bash
PYTHONPATH=src python -m ahrt_stage_a validate-kn \
  --energy-MeV 0.661657 --samples 150000 --seed 81723
```

Result: **PASS**.

- Sample mean cos(theta): 0.3221617139
- Numerical quadrature mean cos(theta): 0.3186931968
- Difference: 2.27 Monte Carlo standard errors
- Sample mean outgoing photon energy: 0.4100467108 MeV
- Numerical quadrature mean: 0.4092017877 MeV
- Difference: 2.20 Monte Carlo standard errors
- Maximum event energy-angle identity error: 0 MeV in sampled run

This validates the sampler against the free-electron Klein-Nishina mathematical
model. It does not validate the free-electron approximation against the
bound-electron incoherent scattering function used in XCOM.

## XCOM G0/G1 status

**Pending real NIST H2O export.** No physical XCOM coefficient table is bundled
or hard-coded. Parser and validation code paths were tested using a deliberately
synthetic CSV fixture marked `NOT PHYSICAL`; that fixture is not suitable for
radiation calculations.

The real-data gate will check manifest provenance/hash, partial-to-total
coefficient closure, and sampled exponential attenuation against Beer-Lambert
using the imported coefficient. This is not an independent re-query of NIST.

## Explicit non-achievements

- No local files under `~/Developer/ArcAdvisor` were changed.
- No Julia implementation was built in this archive; the working reference
  implementation is Python 3.11+ with standard library only.
- No Geant4, Geant4-DNA, OpenMC or electron transport was executed.
- No absorbed dose, radiolysis map, radical concentration, Arc damage or
  biological response was calculated.
