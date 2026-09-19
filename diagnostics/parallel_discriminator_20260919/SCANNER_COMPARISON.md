# Context preparation by scorer: fixed comparison before new energies

Use the existing masked whole-protein OMOL scorer and existing contextual water-H
proposals for1F6S/6IP9. Preserve whole protein/cofactor/metal/waterO and inventories.
All3GGR and27supportedPQQ entries are dry identity operations, with results reused.
No production/default/reference/band changes or new biological labels.

Input inventory revealed whole-chain MACE's existing water geometry differs by
up to0.0102Å from the normalized DFT water geometry. For the matched preparation
comparison, calculate4normalized bound MACE endpoints in addition to8contextual
bound/detached endpoints. This supplies the originally agreed matched two-by-two
test; do not attribute a water internal-geometry difference to orientation.
Total12new MACE evaluations, zero newDFT. Original archive remains a separate column.

Interaction descriptor is (E_Ca,bound-E_Ca,detached)-(E_La,bound-E_La,detached).
The two-call atom-reference cancellation applies only to identical nonmetal
coordinates; contextual preparations depend on the endpoint so require four calls.
Retain atom-reference-corrected total bound contrast and detached-environment
difference as components, never choose a score definition after inspecting output.

Compare all six alpha/GGR structural pairs (one condition-qualified biological
comparison), retaining all failures and same group interpretation. Positive
alpha-minus-GGR is the supported ordering; do not fit a threshold. DFT and MACE
scales are distinct. PQQ score preservation follows exact unchanged inputs and
the identical frozen scorer, not new quantum computations.

Use existing allocatedGPU runner,1GPU/16CPU; expected minutes from prior receipts.
No overall budget cutoff or new model training. Fast reusable implementation and
actual costs are deliverables; no broad accuracy claim from this consumed pair.
