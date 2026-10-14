# Next hypothesis: frozen ion placement concentrates strain in donor motion

**If the full control comparison also agrees, test whether the unusually large response follows the ion's placement relative to PQQ, rather than adding more shell atoms or suppressing donor motion.** The existing donor-only search fixes both the ion and PQQ. It can relieve an uneven coordination geometry only by moving protein donors. This is a plausible source of structurally sensitive scores, not a demonstrated preparation defect.

This assessment reads saved outputs only. The target common-path differential is DFT +17.135466 versus MACE +16.403481 kcal/mol. Control-La proposal4120 has completed, but the full control/target contrast is not assumed available here. The weaker final convergence evidence for control-La origin and target-Ca proposal remains material; agreement alone does not establish unique electronic states.

## What actually differs

These are previously measured distances in Å, copied from the [saved geometry result](../../../workspaces/pqq_inconclusive_review_20261014/geometry_v1/execution/RESULT.json), SHA256 `b116a1076f14e2d0303b89f31feb42a42311de195e21b985278688e2c30d5ed2`. Both sources use the exact original192-atom state inventoried in [INVENTORY.json](../../../workspaces/exaa_donor_response_20261014/inventory_v1/INVENTORY.json), SHA256 `757a13523549697fd52cf9e94b3ba2ffdd01c463f3969fd8a5d01c1459d28343`.

| Contact | Target two-Ca3 origin → La proposal | Control two-Ca0 origin → La proposal |
|---|---:|---:|
| Metal–Glu213 OE1 |2.099 → 2.493|2.369 → 2.523|
| Metal–Glu213 OE2 |2.975 → 2.693|2.824 → 2.785|
| Metal–Asp350 OD1 |2.774 → 2.701|3.115 → 3.030|
| Metal–Asp350 OD2 |3.566 → 3.229|3.892 → 3.578|
| Metal–Asn300 OD1 |2.371 → unchanged|2.420 → unchanged|
| Metal–PQQ N2 |2.733 → unchanged|2.440 → unchanged|
| Metal–PQQ bound carboxyl O5 |2.845 → unchanged|2.562 → unchanged|
| Metal–PQQ O1 |2.403 → unchanged|2.421 → unchanged|

The target starts with a much more uneven Glu contact and weaker N2/carboxyl PQQ proximity, while its Asp contact is already closer. Its Glu motion redistributes that contact asymmetry substantially. Target Glu χ3 changes +8.709°, versus −14.709° in the control: these are different coupled donor responses, not a uniform radial expansion. Target La reaches the existing0.8Å heavy-displacement bound; control proposals do not. Both retain appreciable **vacuum-MACE** load on frozen metal coordinates, so a nonzero load alone cannot identify the problematic source. No composite-force conclusion follows.

The [1FLG source audit](../../../workspaces/ca_coverage_recovery_20261012/refold/run_v3_h200/prepared/1FLG_AUDIT.json) records chainA Glu179/213 OE1/OE2 at2.839/2.928, Asn266/300 OD1 at2.916 and Asp316/350 OD1 at2.844. Its PQQ donor contacts are2.732,2.792 and2.926 (deposited atom names O5/O7A/N6). Thus the crystal has a more open, balanced donor arrangement, **but long metal–PQQ distances themselves are not a defect**: this admitted chain has only +1.929997kcal/mol operational composite accommodation. Its exact post-proposal donor-distance table has not been computed in this assessment. Redox, dry preparation and construct limitations prevent treating it as a uniquely correct endpoint; chainB remains unavailable.

The older Glu/Trp304/backbone H-bond-switch report concerns **La-conditioned sample0**, not these two-Ca sources. It is not evidence that the same H-bond changes occurred here.

## What completed experiments already narrow

- [Packing expansion](../../exaa_source_packing_20261014/packing/FINAL_REPORT.md) changes target selected differential by only −0.374789kcal/mol and control by −0.246722. Added peptide/aromatic/Leu surroundings penalize both metals similarly. Repeating indiscriminate packing expansion is poorly motivated.
- [Target DFT comparison](../TARGET_RESPONSE_REPORT.md) supports the large same-coordinate vacuum response. A MACE-only exaggeration on this path is not supported. Both methods still share the frozen geometry, dry state and local boundary.
- [The earlier joint-metal pilot](../../adaptive_metal_20260922/FINAL_RESULT.md) relieved metal forces but improved no available class decision on its different four-source panel; one proposal failed at an oversized trial. More movable coordinates are not intrinsically beneficial, and that old experiment does not answer this exact ExaA cause.

## Cheapest discriminating test

Before another optimization, use a **reciprocal metal-position transplant** as an explicitly modeled causal diagnostic:

1. On Slurm, align the exact shared, named24 PQQ heavy atoms of sample0 and3 by one frozen rigid least-squares rule. Transfer only the other source's ion position into each host frame. Keep that host's PQQ, proton inventory, source/cap membership and all protein coordinates unchanged. Do not transplant donor atoms or copy crystal coordinates. This preserves the chemical composition but is a modeled geometry, not another observed structure.
2. Apply each transplanted position to that host's saved origin and `adaptive_La` donor arrangements. Require the existing physical/overlap checks and0.8Å metal-displacement limit before scoring; retain any failed graft explicitly. No nearest-contact target, energy-selected direction or threshold relaxation.
3. Evaluate the eight new vacuum-MACE cells: two hosts × two saved donor arrangements × Ca/La. Reuse all corresponding original cells. This is a fixed-geometry interaction test, with no optimizer or new DFT. Compare individual works and their Ca-minus-La differential under the original versus transplanted position.

Support requires the differential response to change **reciprocally with ion placement**, beyond the existing declared1kcal/mol interpretive resolution, rather than merely producing a desired class sign. If the difference remains attached to the host donor/PQQ arrangement, close the narrow ion-placement explanation. An inadmissible transplant leaves the test unresolved; it does not authorize changing geometry gates.

If supported, evaluate those same admitted cells with the already-qualified frozen-charge solvent components before designing a bounded joint ion/donor proposal. The practical upside would be a small, reusable search that reduces sensitivity to uncertain ion placement without sacrificing the protein's specificity. It earns classifier use only through the existing reference/transfer comparison, not by improving these two consumed cases or by fitting their bands. Keep crystalA as an untouched conditional comparator and all previous failures visible. No new experiment, geometry or submission is executed by this assessment.
