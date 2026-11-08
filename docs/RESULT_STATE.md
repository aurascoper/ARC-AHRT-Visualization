# AHRT result-state convention (governing principle, applies to AHRT / AHRT-Fish / Arc branches)
Adopted 2026-11-04 per paste_33. Nulls are typed, not collapsed.

## States

| State | Meaning |
|---|---|
| PENDING | Output exists but run still executing (no result yet) |
| COMPLETED | Completed calculation with a numerical value |
| NOT_REPORTED | Species/quantity not selected by the scorer (no observation) |
| ZERO_ESTIMATE | Completed calculation reports zero (subject to scorer + sampling limits) |
| BELOW_DETECTION | Signal below an experimental detection threshold (interval/bound, not zero) |
| NO_EVENTS | Completed simulation, zero qualifying events (meaningful only after exposure/histories/event-definition checks) |
| UNDEFINED | Quantity mathematically undefined (e.g., zero denominator) — not a physical zero |
| FAILED | Execution terminated abnormally (incl. OOM/SIGKILL) — no valid completed result |
| INVALID | Configuration/model fails acceptance criteria — output cannot support the claim |

Supporting rules:
- `nothing` in anatomical/metadata fields = unknown/unavailable; numerical 0 = represented zero. Never silence one into the other.
- simulation estimates vs experimental measurements are distinguished by `measurement_kind`.
- every result carries: scorer, source, denominator, histories, termination status, evidence/validation reference.

## JSON receipt form (canonical)
{"result_state": "<STATE>", "value": number|null, "uncertainty_1sd": number|null,
 "detection_limit": number|null, "unit": string|null, "measurement_kind": "simulation|experiment",
 "scorer": "...", "source_config_sha256": "...", "histories": int|null, "denominator": "...",
 "termination_status": "...", "evidence_ref": "receipt path or literature ref"}

## Retrofitting rule
Existing receipts keep their historical content. New receipts follow this convention; quantities not
measured carry NOT_REPORTED/PENDING/etc., never null numbers masquerading as zero.

## Current receipts state audit (as of adoption)
- V2 electron ledger: completion states were implicit; effective states = COMPLETED for all tallies.
- W-A outputs: species series COMPLETED; missing species = NOT_REPORTED (none occur, scorer covered 16).
- gamma_heterarchy ordinary run: roi_interactions=0 → classify as NO_EVENTS (histories 200000, event definition verified), not ZERO_ESTIMATE.
- Old V4 batch (killed): all outputs FAILED (invalid configs) — superseded by corrected rerun.
