# Independent mechanics design review

Reviewed 2026-10-19T08:15:06Z by `/root/mex_health_review`.

**The declared fixed-cell intervention is scientifically coherent and appropriately
bounded. Source and actual execution remain unreviewed at this stage.**

Reviewed SPEC.md SHA256:
`07393274ee50b2a4e67fdb60ac999f8a295aacb398e6e6faf466f4171bb52ba9`.
This review reads the saved plan, reports and native adapter source only. It
performs no scientific calculation, submission or source change.

## Why this is the smallest coherent candidate

The complete MM audit14787 found no inventory/copying defect, and the favorable
work spans bonded and nonbonded groups. The tested boundary expansion also leaves
the terminal force discrepancy. These results do not support fitting a terminal
bond, subtracting a selected class, or attaching an empirical force offset.

One whole-system MACE vacuum energy plus the same fixed-charge ALPB is a defined
replacement for the complete partitioned internal-energy model. It includes no
duplicate MM, embedding difference or extra dispersion term. It simultaneously
changes scaffold, boundary, spectator and cross-system interactions, so an
improvement cannot uniquely identify which change mattered.

For LBT this model already exists: qualification12958 passed its frozen component
and total directional tests, and all four broad searches are complete but
nonconverged. Reuse those results. The new material test is the exact Mex mixed-ion
fixture and the separately declared independent whole-peptide reference challenge.

## Subtractive alternatives

A future multilayer alternative can be defined as a whole-system base plus a
matched regional difference, `E_base(full) + E_high(region) - E_base(region)`.
Both regional terms would need the same physical/cap map and a matched physical
state with each native spin convention explicit, complete chain-rule gradients,
a declared environment treatment and fresh
consistency/accuracy checks. It would retain the base treatment of interactions
not replaced by that difference.

Adding the old GFN embedded-minus-vacuum term atop whole MACE has no corresponding
subtraction of the interactions already present in the whole-system model.
Likewise, subtracting selected MM terms does not construct a matched regional
replacement. Neither alternative follows from the inspected decomposition.
The closed Mex external-point-charge DFT derivative failure also prevents using
that route as an already-qualified correction. The declared whole-system model
therefore introduces fewer unsupported choices than a new subtractive scheme.

## Aggregate spin and the eight-cell denominator

The native adapter supplies physical multiplicity through `atoms.info['spin']`
and verifies graph `total_spin`. Its former local La1/Dy6 inputs do not declare
the aggregate state once three Nd spectators enter the model graph.

The frozen La10/Dy15 and La4/Dy9 hypotheses provide a small, matched sensitivity
test. They are aggregate conditioning categories derived from nominal moment
counting. They do not identify localized moments, a particular reversed Nd ion,
actual coupling, spin-orbit physics or a physical ground state. Actual nuclear
electron parity and model categories must pass before each forward; element
coverage alone is insufficient.

One high-spin scout can establish only executability, finite outputs and cost.
The full two-geometry × two-metal × two-hypothesis panel is needed before treating
its response as insensitive to the chosen category. Keep all eight slots and
failures, reuse the admitted scout, and use the same hypothesis at both endpoints.
Do not select a lower-energy hypothesis after inspection. Compare same-model
work and its metal contrast across hypotheses; report force changes as well.
The 1 kcal/mol materiality rule is a reporting convention, not an error bound.
Small sensitivity would leave physical spin and model accuracy unqualified.

## Execution and provenance requirements

The specification retains all 2,049 atoms, waters, target and three Nd spectators,
the exact saved start/terminal arrays, charge +2 and four reused solvent receipts.
The adapter must generalize size/target/state guards while preserving the pinned
checkpoint, float64 forward arithmetic and force convention. Actual full graph
and solvent identities must be checked, not inferred from a task label.

The scout's 600-second model-forward and 90%-of-device peak allocated-memory
limits are explicit continuation gates. Record model loading separately and
define the measured forward interval and peak-memory counter in the adapter.
No scheduler override, smaller system, precision change or energy-only fallback
is implied. Eight gradient calls remain feasibility/sensitivity evidence;
whole-Mex numerical derivative qualification would still precede any later
optimizer. This phase correctly authorizes no optimization.

## Independent evidence and limits

The inspected 7CCO four-cell whole reference is now consumed development and
regression evidence. The older 48-atom references test local chemistry, not the
whole scaffold and termini. Neither establishes prospective scaffold transfer.

Freezing the already-saved 7CCN-La whole-MACE search start/terminal and evaluating
both metals with four new whole-peptide references is a discriminating next test:
it asks whether the candidate's own favorable deformation is supported on the
other source. It needs no new structural search or refitting. This is a new
reference outcome on a consumed source, not a blind protein holdout. Keep the
same-method work definition, reused MACE/ALPB, missing Dy forces, original
native admission, source-state caveats and full failure denominator. Review the
first reference's actual admission/cost before the remaining three.

A contradiction must stop unsupported optimization. Agreement would support
only these conditional fixtures; it would not validate whole Mex, establish
La/Dy affinity, resolve the reference's SCF/force-state limitations or justify
production promotion. No design blocker was found within this declared scope.
