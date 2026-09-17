# Next candidate: matched metal-coordination contrast

Declared 2026-09-17 before evaluating any separated-metal reference. Covered by
Jacob's renewed blanket authorization toward the active MACE goal. Production
and every completed OMOL/baseline score remain unchanged.

## Hypothesis and energy expression

The GGR readout diagnostic identifies an explicit geometry-independent learned
term `sum(a_element) + N*b(total_charge, spin)`. It contributes -36.19 kcal/mol
to the representation shift. Deleting a convenient term is not sufficient
physical or predictive justification. Instead test a separately defined,
matched *coordination descriptor* that cancels all geometry-independent terms:

```
D_M = E_OMOL(core with metal at source site; Q_M, multiplicity_M)
    - E_OMOL(same atoms with only metal separated; Q_M, multiplicity_M)
R_coord = D_Ca - D_La
```

Convert eV to kcal/mol once. Larger R_coord means greater model stabilization
by La relative to Ca. The geometry-independent element/charge/spin readouts
cancel within each D_M. Other charge-conditioned geometry responses remain.
No chosen atomic subset, fitted coefficient, density substitution or quantum
reference offset is introduced.

This is a finite-range model's frozen coordination response, **not** a physical
ionic dissociation energy or a complete binding free energy. OMOL has a total
charge/spin embedding, not constrained charges on separated fragments; it does
not certify that the distant atom is La(III) or Ca(II). The disconnected inputs
may be outside the training distribution. Failure remains useful evidence.
Do not claim that an exact energetic cancellation validates those states.

## Exact input transformation and numerical qualification

Use the completed benchmark's unchanged 64 bound endpoints: 28 canonical PQQ
pairs and four non-PQQ representations (GGR extended/connected; alpha1F6S/6IP9).
All are consumed cases. Retain every atom, source hydrogen, cap, explicit water,
charge, multiplicity, model, checkpoint and software version.

Construct each separated reference by moving **only atom 0 (the selected metal)**
to `(max(nonmetal_x)+30 A, metal_y, metal_z)`. Its distance from every remaining
atom must exceed the checkpoint's 6 A edge cutoff. This is an explicitly
derived diagnostic geometry, never represented as experimental coordinates.
No ligand or protein coordinates change; no atoms are deleted or introduced.

Stage A: GGR extended/connected separated La/Ca, four primary calls. For the
extended pair also run an exact repeat, a further +10 A metal-x displacement,
and the established rigid rotation about the metal: six more calls, ten total.
Use the existing deterministic rotation, 0.01 kcal/mol energy/contrast tolerance
and 0.001 eV/Angstrom force tolerance. Require zero metal-neighbor edges; a
further displacement must leave energies/forces invariant within these same
tolerances. Native analytic metal force must be within 0.001 eV/Angstrom of zero
when disconnected. Capture native atomic and embedding terms in every call;
their sum must reproduce total energy within 0.01 kcal/mol. Actual input batch
charge/spin and unchanged parameters remain required. If numerical checks fail,
stop this version before expanding it; preserve the failure.

Stage B, conditional on Stage A's numerical checks only: evaluate the remaining
60 separated references, reusing the four exact Stage A primary states. Reuse
all 64 original bound results from jobs 1200797/1200799. Stage B is not gated
on the favorable sign of any GGR result. Total: 70 new native OMOL forwards,
zero DFT, solver, training or optimization calls. No new biological controls
or quantum displacement references are consumed.

## Frozen evaluation and resource scope

Use the existing 25-row canonical calibration rule without alteration:
U=max(Ca R_coord), L=min(La R_coord), require L-U>0.02 kcal/mol and all 25 valid.
Only then release this candidate's own bands and apply them to 1H4I, 4MAE and
1KB0. All three must reach their expected supported regions for the canonical
gate. No inherited OMOL raw-R or baseline bands. Accession overlap and
motif/composition confounding remain explicit; no new independent-validation
or absent-training-overlap claim.

Retain the same non-PQQ direction tests: alpha1F6S minus GGRextended and
alpha6IP9 minus GGRextended must each exceed 0.02 kcal/mol. Repeat both with
GGRconnected; all four must pass for robustness on these consumed comparisons.
Do not choose whichever GGR representation works. Report GGR representation
sensitivity, every missing endpoint, the original OMOL and DFT contrasts, and
all denominators. The two alpha geometries remain one biological group and
the two GGR representations another. These few cases cannot prove broad
affinity discrimination even if they pass.

One A5000, 16 CPUs, 64,474 MiB per job; existing task runner, pinned native
float64 environment and 100M checkpoint. Prior 70-call trial cost 670 allocated
GPU seconds including startup/validation; expect the same order of cost and
measure actual execution. There is no project compute stopping budget. Preserve
scheduler rules, all failures and other jobs. Finite manifests and separate
protocol `mace_omol_matched_coordination_descriptor_v1` prevent substitution for
the baseline or original OMOL score.

After this candidate's fixed evaluation, decide whether it merits additional
independent controls, a physical refinement, or rejection. Do not refit it to
the test results or call the larger active goal complete on this panel alone.
