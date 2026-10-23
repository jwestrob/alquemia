# Consumed W8/W9 metal–shell counterpoise diagnostic

## Question

How much does metal–whole-water-shell basis-set superposition change the existing
paired W8→W9 La/Dy reference contrast? The failed model differs by4.338284kcal/mol,
but existing dimer tests do not directly bound this cluster-reference uncertainty.
This is one reference discriminator before further parameter fitting.

## Fixed inputs and cycle

Use the four admitted consumed aquo8_central/aquo9_central La/Dy geometries and
native full-complex outputs linked by the immutable joint-water preparation DATA
record in ROOT_INPUTS.json. Preserve exact atom order, water coordinates, charge
and core/state conventions. Independent saved inspection has verified identical
water coordinates between La/Dy at each size and the same per-metal electronic
recipe across sizes. Preparation must retain these checks and actual parent pins.

Let ΔX = X(Dy,W9)−X(La,W9)−X(Dy,W8)+X(La,W8).
For each of four complexes calculate two electronic single-point fragments:
(1) real ion with every water atom ghosted; (2) real whole water shell with the
metal ghosted. Exactly eight new native evaluations, four full complexes reused.
Define G as the sum of those two fragment electronic energies. The corrected
paired total is C_original,total−ΔG. Own-ion terms cancel across sizes separately
for each metal; own-shell terms cancel across metals at each size. If actual
source checks invalidate either cancellation, stop before molecular calls.

Preserve the original four-complex D4 contribution separately and unchanged.
Fragment basis corrections use electronic energies only; use the established
qualified canonical-water ghost input convention. Pin the final concrete decks,
engine/assets and native settings before launch. Real-ion fragments retain their
original ECP; ghost ions retain orbital/auxiliary basis functions and must not
retain an active ECP. Expected full basis counts391(W8)/434(W9); fragment electron
counts8 for ion,80/90 for neutral shell, singlet throughout under this frozen-core
recipe. Actual native output must confirm ghost roles, charge/electrons, basis,
SCF completion and finite energies; execution success alone is insufficient.

## Resources and limits

Eight independent4MPI-CPU jobs, standard_expand, scheduler-default memory/time,
with memory settings derived from actual grants and the established native runner.
Up to eight jobs in parallel within current user/job limits. Completion watchers
are required. One attempt per declared cell; preserve failures/missing slots.
No new geometries, forces, quantum method sweep, fitting, automatic retries or
production changes. Reuse existing native adapters and qualified ghost conventions;
no new general execution framework. Source/preparation and independent review
precede root binding and submissions. Root owns all submissions.

## Interpretation fixed before execution

Report all eight energies, their exact cancellations, original electronic/D4
contributions, CP shift and unchanged-model discrepancy against this alternative
reference. A shift comparable with4.34kcal/mol would make the reference sensitive
to this basis-sharing correction; a small shift would rule out this particular
explanation at the measured scale. Neither result establishes the exact reference.
This two-fragment diagnostic does not remove intra-shell water–water BSSE, reach
the basis limit, test the frozen-core approximation or qualify forces/free energies.
No automatic replacement of fitting targets, refit, fresh-validation release or
bulk/protein sampling follows. Keep every prior result and failed model intact.
