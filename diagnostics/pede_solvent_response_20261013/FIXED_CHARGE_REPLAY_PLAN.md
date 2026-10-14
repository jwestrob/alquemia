# Next diagnostic: replay native ALPB with saved charge distributions

Declared after the eight-cell DFT comparison,14October2026UTC. This is a
development diagnostic, not a fitted classifier or a fresh blind validation.

## Why this test

The two PedE origins differ in their Ca-minus-La solvent correction by
+21.256628kcal/mol with corrected GFN2/ALPB and +10.034216kcal/mol with native
r2SCAN-3c/CPCM. Both show the same direction. The discrepancy does not establish
a unique electronic or solvent cause, and the DFT states have unresolved
stability/reference warnings. Deleting solvent or multiplying it by a fitted
factor is not justified.

We can now separate the native ALPB contribution evaluated on the vacuum charge
distribution from its additional solvent-induced redistribution, using existing
charges and the actual maintained solver. This adds no SCC, DFT, optimization,
new protein structure or changed chemical state.

## Fixed matrix and accounting

Use the exact original sample1/sample3 × Ca/La geometries and their saved
corrected-GFN2 vacuum/water atomic charges, persisted to eight decimal places.
No higher-precision restart is retained. Eight evaluations:
each of the four native ALPB models at its own q_vac and q_water. Atom order,
element/charge/spin, source inventory and parameters must match the archived
endpoint. Do not average charge schemes, import DFT populations or invent shell
charges. Propagate the actual charge-rounding interval rather than inventing
extra precision. Insufficient precision or an unsupported interface is an
explicit blocker.

Let G(q) be native ALPB electrostatic + hydrogen-bond + SASA + reference energy
at that fixed geometry and specified charges. Retain the already computed
intrinsic remainder I = E_water − G(q_water) − E_vac. Then exactly:

`U = G(q_vac) + [G(q_water) − G(q_vac) + I]`.

Report both terms for each endpoint and their source/metal contrasts. The
bracket contains the additional redistribution contribution including its
intrinsic-energy compensation. It is not permissible to drop I, replace U with
G, or declare one isolated component a new solvent correction. G(q_vac) still
contains geometry and vacuum-charge dependence; it does not isolate cavity
geometry by itself.

## Qualification before interpretation

First reproduce all four archived water-charge Gelec/Ghb/Gsasa/Gshift values.
Freeze numerical tolerances in the manifest from persisted charge/output
precision before any replay result is inspected. Retain each residual. A failed
diagonal means the wrapper/model setup is not qualified: diagnose parameters,
units, charge precision or interface before interpreting counterfactuals. Do not
adjust a tolerance to conceal disagreement.

Use pinned xTB6.7.1 native factory and TBorn component routines in an isolated
version-matched build. No new Born/PB engine, manual fitted constants, changed
continuum model or installation over a working environment. Record source,
dependency/compiler, binary and fixture hashes. Native hydrogen bonding is
charge-dependent and must remain separately reported. Preserve the exact
effective water parameters and reference convention from the original outputs.

Source interfaces: [factory](https://raw.githubusercontent.com/grimme-lab/xtb/v6.7.1/src/solv/model.f90)
and [component evaluator](https://raw.githubusercontent.com/grimme-lab/xtb/v6.7.1/src/solv/gbsa.f90).
Installed Python/C APIs do not directly supply internal fixed atomic populations;
an external point-charge field is not a substitute.

## Execution and decision

Root submits a finite isolated build and qualification through standard Slurm,
within current12CPU/job shared policy, no memory/time overrides. All compilation,
tests and component evaluations run on Slurm. Record actual resource receipts
and completion wakes. A build failure may receive a specific technical repair;
do not turn this into an open-ended dependency replacement.

Only after diagonal qualification interpret the eight frozen evaluations.
Dominant additional redistribution directs the next model test toward density
response. Dominant frozen-distribution sensitivity directs attention to cavity,
geometry or the vacuum distribution. Neither result automatically selects a
new scorer. Any changed scorer still needs the designated canonical calibration
and full declared structural test population. Production remains unchanged.

Agent la_parameter_scope owns the wrapper/manifest under fixed_charge_replay/;
root owns execution, interpretation and shared status. The earlier seeded DFT
failure contingency is closed and unrun because the primary converged.
