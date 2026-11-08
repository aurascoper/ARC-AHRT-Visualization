# AHRT v0.1 -- FROZEN BASELINE
## Cs-137 / dArc1 structural reference / nested radiation transport

Status: FROZEN 2026-11-04. Source: external review spec, machine-verified before freezing.
Numerical claims below marked [VERIFIED] were computed independently at freeze time.

## 1. Decisions

| Decision                   | Baseline                                                              |
| -------------------------- | --------------------------------------------------------------------- |
| Photon source              | Cs-137 reference line, 661.657 keV                                    |
| Target                     | Arc-containing cellular volume                                        |
| Cell geometry              | Cylindrical, heterogeneous at the biological-response level           |
| Primary photon medium      | Homogeneous water-equivalent medium                                   |
| Molecular structural input | dArc1 6TAP, chain classes kept distinct                               |
| Secondary transport        | Electron transport, then microscopic track structure                  |
| Chemistry                  | Water radiolysis; separately parameterized Arc/RNA reactions          |
| Visualization              | ChimeraX 1.12 (installed)                                             |
| Initial implementation     | Julia/Python prototype; validated transport engine as comparison      |
| Synthetic scenario         | Mammalian Arc packing -> separate `scenario_synthetic_mammalian_arc`; NOT baseline |

## 2. Verified numerics (freeze-time recomputation)

| Claim                         | Spec value      | Recomputed                     | Status     |
| ----------------------------- | --------------- | ------------------------------ | ---------- |
| Compton edge at 661.657 keV   | ~478 keV        | 477.3 keV; backscatter 184.3   | [VERIFIED] |
| Water mu/rho at 662 keV       | ~0.086 cm^2/g   | XCOM canonical ~0.0859-0.0866  | [VERIFIED; import exact table at runtime, do not hard-code] |
| Mean free path (rho=1 g/cm3)  | ~11.6 cm        | 11.6 cm                        | [VERIFIED] |
| P(direct interact, 10 um cell)| "very low"      | 8.6e-5 (0.0086%)               | [VERIFIED] -- secondary-electron argument holds |
| 500 keV e- CSDA range, water  | ~1.8 mm         | consistent (Feather rule ~1.4-1.8 mm); import ESTAR at runtime | [PLAUSIBLE -- verify at G2 gate] |
| e- range vs 5 um cell radius  | >>1             | 360x                           | [VERIFIED] -- nested domains materially necessary, not cosmetic |

## 3. Spatial hierarchy (modeling proposal, not a measured cell)

- Level 1: macro water buffer, 3 mm nominal -- subject to electron-range and leakage convergence tests (G4).
- Level 2: cellular cylinder R=5 um, L=20 um (synthetic geometry).
- Level 3: molecular neighborhood 100 nm ROI around selected Arc assembly.
- Level 4: 6TAP coordinates, chain classes, Zn sites, experimentally grounded contacts.

## 4. Transport architecture (stage-coupled; interfaces fixed)

- Stage A -- primary photon transport: Julia/Python, NIST XCOM partial coefficients + Klein-Nishina Compton sampling; OpenMC as independent photon cross-check. OpenMC does NOT track secondary electrons spatially (local energy deposition / TTB only) -- it cannot produce fig_31-class track geometry. [CONFIRMED against docs.openmc.org description]
- Stage B -- secondary electron transport: Geant4 standard EM (NOT installed on this machine at freeze time).
- Stage C -- microscopic track structure: Geant4-DNA liquid water (NOT installed).
- Stage D -- chemistry + Arc heterarchy response: Julia. 6TAP-derived geometry + typed interaction edges.

Machine reality at freeze: ChimeraX-1.12.app present; Geant4 absent; OpenFOAM case data only. Stage A is locally developable now; B/C enter via container when needed, behind defined phase-space interfaces.

## 5. 6TAP chain-class treatment (binds to chain-class reconciliation record)

Structural input identity: `dArc1_6TAP_full_assembly_chain_class_aware`
- A-D -> 240 CA instances = T=4 capsid shell (T=4 applies ONLY to these).
- E (res 224-250, Arg-rich NC zinc-knuckle fragment) -> 60 instances; separate structural/chemical class; NOT Caspar-Klug lattice subunits.
- F -> Zn ion; zinc coordination represented separately from the polymer graph.
- 300 polymer-chain instances = 240 CA + 60 NC fragment; NOT a 300-subunit T-number.
- 6TAP is insect dArc1; no mammalian Arc capsid claim.

Typed edges: covalent / direct_contact / zinc_coordination / rna_association (independent evidence gate) / literature_association / mechanical_coupling (requires analysis) / hypothetical_packing (synthetic; excluded from baseline).

## 6. Primary/secondary separation (enforced in code architecture)

Primary photon stage: attenuation + Compton ONLY. Never: photon advection, Manning, PB, counterion fields.
Secondary stage: electron track structure, water radiolysis (OH*, e_aq-, H*, H2O2, H2), radical diffusion/reactions, RNA/protein structural response; Manning/PB layers live here and ONLY here, as optional molecular-environment inputs. NaCl and MgCl2 are separate electrolyte species, never one valence dial. A Manning condensed fraction is never a measured radical-scavenging rate.

## 7. Validation gates (release requirements)

- G0 input provenance + dimensional consistency: 100% of active inputs classified.
- G1 water attenuation vs XCOM within 1%; Klein-Nishina energy-angle consistency.
- G2 electron stopping/range vs ESTAR within 5%; weighted global energy balance within 0.5% (CSDA path length != projected range -- tortuous paths).
- G3 water radiolysis yields vs published benchmarks within benchmark uncertainty.
- G4 cell-ROI dose uncertainty <= 5% relative; macro-domain boundary convergence < 2% change.
- G5 Arc damage-response model: DISABLED until independently parameterized.

Smoke test: 1e6 primary histories; then weighted importance splitting around the cell; preserve particle weights across stages. No Gy/s or biological dose until source activity, distance, geometry, duration are real inputs.

## 8. First executable milestone

Cs-137 parallel beam (lateral incidence, unpolarized) -> water-equivalent nested cylinder -> electron source map -> 5 um track region -> 100 nm dArc1 ROI -> energy-deposition field + track history -> ChimeraX export. One synchronized visualization distinguishing simulated trajectories from structural coordinates from conceptual network edges. Every synthetic/assumed/structure-derived quantity recorded in the evidence ledger.

## 9. Honest boundary

Without an experimental Arc-cellular irradiation dataset, AHRT v0.1 establishes radiation transport and calculated local exposure ONLY. It cannot establish that gamma reorganizes the Arc heterarchy, alters RNA packaging, or causes a cellular response. Those are testable hypotheses downstream, not baseline outputs.

## 10. Open items at freeze

- fig_31 software provenance unknown (visual benchmark only until identified).
- Import XCOM water table and ESTAR water data at runtime (data/xcom/, data/estar/).
- cpm-sp-v1 integration: replace prescribed -kM mortality with D(x)=E_dep(x)/m(x) exposure field; first release emits exposure/damage-opportunity maps, NOT cell-death probabilities.
- OpenFOAM and PHREEQC: explicitly out of v0.1.
