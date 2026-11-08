# AHRT-Fish / AHRT-Sphyrna — workstream specification
Prepared 2026-11-04. Branch of AHRT; V1/V2 baselines and W-A/V3/W-A-03 receipts untouched.
Scope: radiation-transport + radiochemistry model of WHOLE aquatic organisms. No Arc content here —
the Arc track (V5) remains DISABLED. Fish branch shares only the frozen transport stages and the
receipt discipline.

## Work packages and gates

| WP | Deliverable | Status | Gate/acceptance |
|---|---|---|---|
| F0 | Fish anatomy + water-boundary phantom (Medaka reference first) | NEW | Organ segmentation manifest + material table flagged as ICRP-human-borrowed where applicable |
| F1 | External photon/electron transport + organ-dose benchmark | NEW | AHRT Stage A/B transport reused UNCHANGED; organ absorbed fractions vs Medaka voxel-phantom published values |
| F2 | Internal radionuclide source + biokinetic compartments | NEW | Marine-fish compartment model (muscle/bone/organs/gill/digestive uptake) parameterized from literature; separate receipts per nuclide |
| F3 | Fish-tissue chemistry + scavenger validation | NEW | Tissue-specific oxygenation/scavenging; pure-water chemistry is NOT tissue chemistry — new validation suite mandatory |
| F4 | CFD trajectory-to-exposure coupling | NEW | Sphyrnidae cephalofoil/hydro field (2020 family-wide CFD benchmark) -> exposure maps; typed handoff only |
| F5 | Fish biological response validation | DISABLED | Oxidative stress / DNA / gill+liver histology vs measurements; enables ONLY after V4/F3 evidence |
| F6 | SensorimotorMechanics (hydrodynamic-sensory-motor coupling) | SPEC (2026-11-04) | cephalofoil pressure/velocity fields + neuromast positions + body curvature + fin loading + motor feedback synchronized over a swimming cycle; no magnetic/radiation variable silently treated as mechanical stimulus |

## F6 detail (separate channels, never pooled)
- External mechanosensation: lateral line (neuromasts), cutaneous contact receptors
- Internal proprioception: vestibular (acceleration/rotation), muscle/fin proprioceptors, curvature/loading
- Corollary-discharge integration: Skandalis 2021 (lateral-line feedback enabling body-state estimation); Afridi 2025 (EMG-kinematics under laminar vs Karman-vortex flow; koi/carp only)
- Elasmobranch caveat: muscle spindles absent in fishes; corpuscular proprioceptive endings uncommon; two types recognized in elasmobranchs (Bone & Chubb 1976); hammerhead-specific proprioceptive mapping is an EVIDENCE GAP
- PC12 link: Raffa 2018 (magnetic nanoparticles + magnetic-field gradient applied 0.4-0.9 pN stretch -> 0.252 um/h/pN neurite elongation; B-field alone and nanoparticles alone did NOT reproduce; stretch alone insufficient without NGF); magnetic field is the force-delivery means, NOT magnetoreception evidence
- Anchors: Skandalis 2021; Williams 2013 (fin-ray proprioception); Afridi 2025; Webb 2023 (lateral-line review); Raffa 2018; Gardiner 2014 (shark lateral-line behavior); Bone & Chubb 1976

## Hard rules (same logic as the Arc branch)

- External irradiation and internal contamination are SEPARATE exposure problems, never pooled.
- Water is part of the dosimetry geometry (aquarium/container/tank modeled); the >200% dose variation
  reported under differing irradiation conditions (Oxford Academic reference in proposal) justifies
  independent container-geometry dosimetry checks before any organ-dose claim.
- Sphyrna-specific anatomy (cephalofoil, ampullae of Lorenzini, olfactory chambers, cartilaginous
  skeleton) is segmented independently; never obtained by scaling the Medaka phantom.
- CFD results never become radiation-dose or injury claims; chemical fields never become histology claims.
- Material compositions follow the Medaka precedent's provenance flag (borrowed human-ICRP data ->
  explicitly labeled, not validated).

## Visual panels (literature references, provenance-unverified local assets)

A. hydrodynamic field | B. radiation field | C. anatomical segmentation | D. chemical field.
Reference figures downloaded to ~/Downloads (2026-11-04 18:37-38); one already identified as a
cold-plasma tumorsphere/PI + ROS-distance panel (panel D reference). All figure assets are
evidence_class=literature_reference; permitted_claims=none until per-figure provenance audit.
No figure may become a measured fish datum.

## Anchoring literature (titles as given in the proposal)
- Medaka voxel phantom, J Environ Radioact 2023-24 (8-organ absorbed fractions)
- Aquatic-organism gamma dosimetry (water/air, position, rotation; >200% variation) — N. Amer. J. Aquaculture 62(2):95
- Hammerhead cephalofoil family-wide CFD — Sci Rep 2020 (all 8 species + 3 typical sharks)
- Sphyrna tudes nasal-region micro-CT/MRI flow simulation 2013

## Interlocks with existing AHRT state
- Uses the frozen XCOM/ESTAR/KN data + Geant4 11.4.3 V2 binary for transport layers.
- Chemistry layer must reuse the validated TOPAS-nBio path (container 11.3.2), with tissue-scavenger
  configurations as NEW fixtures (do not silently extend W-A).
- AHRT-Fish may not consume ArcAdvisor/NeuroARC artifacts.
