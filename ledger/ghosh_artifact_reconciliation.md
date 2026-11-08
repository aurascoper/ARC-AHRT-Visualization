# Artifact reconciliation: what may go to Dr. Ghosh

Generated 2026-11-03 by ArcAdvisor Release A.1. **Nothing has been sent.**
Grounded in the interaction logs through 2026-11-03 and in defect checks run against
the files themselves today.

## What the logs actually show

Sixteen Ghosh-adjacent turns across 2026-11-01 to 2026-11-03. The decisive one is
`dialectic-turns-2026-11-03.jsonl` line 165, run `324AE82A`, 02:28:04Z, where the
question "should I share the camera racks iso decahedron" was answered by the voice
model with "Yes, send it to Goche". That run contains two turns and ends there.
No email was composed, no draft was written, nothing was sent.

Two naming notes for the ASR glossary: "Santa Neil Goche" and "Centennial Goch" are
both the same collaborator, and "iso decahedron" is not a solid. The measured point
group is icosahedral.

## Defect check run today against the Ghosh-facing files

Each marker below was found still present. These are the defects the 2026-11-03
audit recorded; none has been fixed.

| Marker | Files still containing it | Why it blocks sharing |
|---|---|---|
| `ready_for_ghosh_lab_data` | 3 | A literal `true`. Its test cannot fail, so the readiness it advertises is not measured. |
| `0.942` | 2 | Hard-coded spectral radius. The Julia operator's Perron root is 1.51 to 1.63, so the figure shown disagrees with the code behind it. |
| `9.26` | 3 | Hysteresis loop area as a literal. No loop was integrated. |
| `SEMO-NanoBio-Lab` | 1 | A real laboratory name embedded in a synthetic FCS fixture. Synthetic data carrying a real lab's name is the single worst item here. |
| `0.229232` | 2 | Langevin alignment value. The correct value at the stated eta is 0.229402. |

## Artifact ledger

| Artifact | Size | Status | Decision |
|---|---|---|---|
| `flow-cytometry-model-bridge/workbench/dr_ghosh_discussion_deck_and_synthesis.html` | 163 KB | Carries literals above | **HOLD.** Fix the five markers first. |
| `flow-cytometry-model-bridge/workbench/microfluidic_cytometry_workbench.html` | 135 KB | Carries `0.942` and the synthetic FCS lab name | **HOLD.** |
| `flow-cytometry-model-bridge/docs/track1_pc12_nanocarrier_internalization_brief.md` | 49 KB | Not re-audited today | **HOLD** pending a read. |
| `~/Downloads/PC12 Nanocarrier Working Session.pdf` | 2.5 MB | Honest in text ("Nothing here is lab data") but two layout defects and the stale figures | **HOLD.** |
| `ArcAdvisor/workbench/arc_camera_rack.cxc` | 10 KB | Generated from coordinates; corrections in header; zero em/en dashes | **READY on the science**, still gated by the confidentiality hold. |
| `ArcAdvisor/workbench/arc_advisor.html` | 152 KB | Offline, no network, every number traceable | **READY on the science.** Contains internal critique of your own repositories, so it is a judgement call whether an outside collaborator should see it. |
| `ArcAdvisor/ledger/correction_proposal.md` | 4 KB | Internal correction proposal | **DO NOT SEND.** Internal. |
| `geometry/*.json` | varies | Raw measurements | Optional supporting data if asked. |

## The governing constraint

The confidentiality hold opened by the 2026-11-03 audit is still open. It is not a
science question and nothing in Release A resolves it. Until the supervisor and the
research-security office answer what may be published from personal modelling, the
recommended action is to send prose only, with no attachments, which is how the
accompanying draft is written.

## Recommended sequence

1. Resolve the confidentiality hold.
2. Send the text-only note (the accompanying draft) with the corrected figures.
3. Fix the five markers before any deck or workbench goes out.
4. Rename the synthetic FCS fixture and strip the real laboratory name from it. This
   one is worth doing regardless of whether anything is ever sent.
