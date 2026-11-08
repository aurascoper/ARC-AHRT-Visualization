# W-A: pure-water G-value preflight — AHRT-ArcHydration-V3 / W-A suite
Prepared 2026-11-04. No chemistry results yet. This document freezes the protocol, not any output.

## Environment isolation (VERIFIED constraint)

OpenTOPAS compatibility (verified against opentopas.readthedocs.io Geant4-Versions page 2026-11-04):
- OpenTOPAS 4.3.0: supports Geant4 11.3.2 and 11.4.2 (core, fully regression-tested on both).
- Public TOPAS-nBio chemistry: supported ONLY with Geant4 11.3.2; currently incompatible with 11.4.2.
- Geant4 11.4.3 (our V2 env) is NEWER than the documented compatibility set; do not reuse it for chemistry.

Decision: retain Geant4 11.4.3 (V2 frozen) untouched. Build the chemistry stage separately:
- Preferred route: official container `opentopas/opentopas:v4.3.0-geant4-11.3.2-arm64` (host Docker via OrbStack).
  TOPAS-nBio extension built inside the container via the documented topas-docker extension workflow.
- Fallback: native source build Geant4 11.3.2 + OpenTOPAS 4.3.0 + TOPAS-nBio 4.1.0 (separate prefix).
Both routes produce an isolated runtime; W-A receipts must record which was used, with image tag/digest or build config hash.

## Reference configurations to reproduce (not to silently inherit)

| Ref | Upstream file | Configuration sketch | Role |
|---|---|---|---|
| A | SBSGetGValue/GvalueG4DNADefault.txt | 80 keV electron point source; G4_WATER; g4em-dna + g4em-dna-chemistry; 100 histories; time-dependent G to 1 us | Geant4-DNA default chemistry |
| B | IRTGetGValue/TsIRTGvalue.txt | 999.999 keV electron point source; water-only; TsEmDNAPhysics + TsEmDNAChemistry; ELSEPA elastic model; Meesungnoen solvated-electron thermalization; explicit TOPAS-nBio reactions + diffusion; 3,000 histories | Validated-chemistry reproduction (Ramos-Méndez 2021) |
| C | project-specific W-A-ARC-441 | 441.039 keV recoil electron (frozen conditional event), pure liquid water, validated W-A chemistry config | Transfer test toward Arc target — NOT experimentally validated by W-A passing |

Upstream URLs (github.com/topas-nbio/TOPAS-nBio main): examples/scorers/SBSGetGValue/GvalueG4DNADefault.txt; examples/scorers/IRTGetGValue/TsIRTGvalue.txt; examples/scorers/SBSGetGValue/TOPASDefaultReactions.txt.

## Acceptance matrix

| Test | Acceptance |
|---|---|
| W-A-01 configuration discipline | Physics/chemistry modules named explicitly; no unsupported materials |
| W-A-02 upstream execution | Ref A and Ref B run end-to-end, interpretable time-dependent output |
| W-A-03 species agreement | OH·, e-_aq, H2 within declared uncertainty of cited experiments; H2O2 and H· EXCLUDED from pure-water acceptance (2021 paper: these need scavenger conditions — a miss in W-A is not a chemistry failure) |
| W-A-04 accounting | Units, species identities, reaction stoichiometry correct; deterministic on seed |
| W-A-05 W-A-ARC-441 | Separate documented project-specific result; no validation claim from neighbor energies |
| W-A-06 provenance | Physics/chemistry/seed/geometry/config hashes retained in receipt |

## Declared configuration parameters (values recorded WITH sources, not silently inherited)

- RMS displacement parameter for water molecules after ionization/excitation: the 2021 validation used 1.24 nm (improved agreement vs the older 2 nm default). MUST appear in the Ref B configuration manifest with its source citation; do not accept the code default without noting which value is in effect.
- Time grid, species set, diffusion coefficients, reaction rates: extracted from the upstream .txt files; each numeric entry must carry its source line reference in the manifest.
- Temperature, pH, oxygen condition: declared per fixture (W-A: pure, unaerated).

## Deliverables (pending the isolated runtime)

1. W-A-A + W-A-B execution transcripts and time-dependent G-value outputs.
2. Species-by-species agreement table vs the 2021 experimental comparison set.
3. W-A-ARC-441 run: initial ionization/excitation distribution, initial species, time-dependent G, diffusion/reaction events, statistical uncertainties, deposited-energy accounting.
4. Updated receipts + compatibility-matrix gate status.

## Stop conditions carried over
- No Arc molecular-damage scoring inside W-A (chemistry-first, target later).
- No validation claim for the 441 keV transfer from reference-energy passes alone.
- No claim that a successful TOPAS execution is itself validation.
