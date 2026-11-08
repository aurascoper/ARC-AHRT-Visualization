# CephaloScope — chromatophore + local sensorimotor MVP (spec)
Prepared 2026-11-04. Independent branches:
- AHRT-Cephalopod: anatomy/seawater/external-vs-internal exposure/dosimetry/radiochemistry (reuses frozen AHRT transport; tissue radiochemistry needs new validation)
- Ceph-Sensorimotor: arm mechanics, sucker sensing, proprioceptive candidates, ganglionic integration, statocysts, central coordination
- CephaloScope (GUI/MVP): skin-patterning observatory; ChromatophoreLab first module

No single undifferentiated model. Every output carries source/species/anatomical location/confidence/validation state.

## Six layers
A. Observation: video frames, animal identity, species, body region, coordinates, illumination, stimulus, timestamps
B. Chromatophore state: segmentation, sacculus expansion, radial-muscle-associated shape, position, local pattern, movement compensation (foundation: CHROMAS pipeline, not recreated from scratch)
C. Neural/receptor hypotheses: glutamate; AMPA/kainate (fast, both regions); NMDA-like (dorsal subpopulation only); calcium fast/slow components; FaRPs slower modulation (Sepia separate species from squid study set)
D. Local feedback/control: stimulus -> neural -> radial muscle -> skin state -> possible sensory feedback
E. Pattern selection/DRT: SVD modes U (spatial), sigma (weights), V (temporal); modes are mathematical decompositions, ONLY neural-circuit claims when independently evidenced; spatial coordinates always preserved per chromatophore (region, not erased)
F. Optional interfaces: Photoreception (octopus arms; sight-independent negative phototaxis 2021 — implement FIRST), CryptochromeRPM (disabled, evidence-gated; bird-centric field, cephalopod compass unvalidated per 2025 review), AHRT (dosimetry handoff only; radiation-response disabled)

## Receptor layer (squid 2003 evidence; spatially typed)
| Component | Evidence | Model |
|---|---|---|
| L-glutamate | evokes contraction dorsal+ventral skin | direct input |
| AMPA/kainate-like | fast response both regions | primary fast excitatory path |
| NMDA-like | dorsal subpopulation only | conditional, location-specific |
| calcium | fast + slow components | separate time constants |
| FaRPs | Sepia (slower modulation) | separate modulatory input, species-labeled |

dc_i/dt = alpha_A g_A(t,i) + alpha_N g_N(t,i) - beta c_i  (proposed abstraction, NOT a validated universal law)
de_i/dt = f(c_i, motor_i, elastic-recoil_i)  (signal / contraction / appearance kept distinct)
Never equate chromatophore expansion with NMDA activation; dorsal vs ventral is mandatory metadata.

## Local-closure acceptance
The boxed criterion: does a MEASURED local perturbation produce a MEASURABLE, causally-dependent local corrective response?
Feedforward-program vs closed-local-loop are competing models evaluated against the same video+stimulus.
Perturbation protocol (proposed, partially pharmacology-grounded in the squid lit): region stimulate; AMPA/kainate block; NMDA block (dorsal population); controlled local deformation; afferent interrupt; afferent restore. Pattern-replica reproduction alone NEVER counts as loop-closure evidence.
Chromatophores are neuromuscular organs; per classic review they apparently lack visual/proprioceptive feedback — treat as the contrasting subsystem to test, not as the assumed closed loop.

## Mandatory identity fields
Organism, anatomical region (dorsal/ventral mantle/head/arm), sample identity, chromatophore identity, normalized+calibrated position, receptor evidence state, measurement time series, provenance (video/pharma/histo/inference), state per AHRT RESULT_STATE taxonomy.

## MVP sequence
MVP-1 CHROMAS ingestion (segmentation, tracking, spatiotemporal measures)
MVP-2 SVD pattern atlas (region/species/stimulus comparisons)
MVP-3 receptor-state overlay (experimental vs simulated display separated)
MVP-4 local-feedback experiment (perturbation protocol above)
MVP-5 photoreception, then cryptochrome (strictly evidence-gated)
MVP-6 AHRT interface (validated exposure histories; response modeling disabled)

## Anchors
CHROMAS (chromatophore tracking); squid glutamate receptor study 2003; octopus-arm neuromuscular control simulation 2211.06767 (simulation != biological closure proof); octopus-arm sight-independent negative phototaxis 2021; Turing-based bimodal population code arXiv 2205.11500; cryptochrome radical-pair review 2508.21350; cephalopod magnetosensation review 2025 (no documented compass; absence of evidence, not proof of absence).

## Naming/hosting
Working names CephaloScope / ChromatophoreLab. No new repo created and no existing repo repurposed until the user identifies the intended host.
