# Frozen AHRT v0.1 bindings

This scaffold adopts the user's frozen AHRT v0.1 decisions as implementation
constraints. This file records the supplied freeze state; it does not claim to
have independently inspected the user's local `ArcAdvisor` repository or
claude-mind store.

## Primary photon stage

- Reference source: Cs-137 line, 661.657 keV.
- Primary medium: homogeneous water-equivalent medium.
- Gamma is neutral: no Manning parameter, Poisson-Boltzmann force, electric
  steering, or fluid advection is applied to the primary photon.
- OpenMC is restricted to an independent photon-transport cross-check; no
  electron track-structure claim is made from OpenMC local deposition/TTB.
- XCOM coefficients must be imported as runtime data and provenance-hashed.
- Stage A uses XCOM total attenuation without coherent scattering. Coherent
  scattering is not silently sampled with an approximate angular law.

## Structural reference reserved for downstream stages

- PDB: 6TAP (dArc1).
- Chain classes: A-D are 240 CA shell instances, classified as T=4; E is 60
  NC-fragment instances; F is a separate Zn site.
- The 300 polymer-chain count is 240 CA + 60 NC fragments. It is not a single
  300-subunit Caspar-Klug T-number assignment.
- The structure is not a mammalian Arc capsid structure.
- The actual structural coordinates and assembly transformation file must be
  added with their own accession and checksum before any geometry-derived
  molecular scoring is enabled.

## Image/evidence bindings

- fig_31: electron-track visual benchmark only. Its generating software,
  physics settings, medium and provenance must be verified before treating it
  as a quantitative validation dataset.
- fig_32: tissue/cell context only.
- fig_33: Caspar-Klug didactic geometry context.
- Other images retain the evidence classes in the user's v3 asset ledger;
  this package does not upgrade candidate citations to verified sources.

## Deferred physics and response

- Manning/PB is confined to an optional secondary chemical-environment layer.
- NaCl and MgCl2 remain separate solution conditions in a future chemistry
  stage; they are not represented by one effective valence parameter.
- G5 damage/biological response remains disabled.
- No claim is made that gamma irradiation reorganizes Arc, alters RNA packaging,
  or causes a specified cellular phenotype.
- OpenFOAM and PHREEQC are out of Stage A scope.
