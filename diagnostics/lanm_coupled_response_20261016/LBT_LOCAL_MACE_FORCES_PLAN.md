# Four fixed matched-core MACE force evaluations

Declared before force evaluation or inspection of unfinished DFT outputs.
Root owns execution; protocol `lbt_local_MACE_forces_v1`.

## Question and fixed denominator

Obtain the missing analytic MACE forces on all four original 48-atom matched
cores: 7CCO/La, 7CCO/Dy, 7CCN/La and 7CCN/Dy. The new observable is the force;
the existing job 12471 energies are preserved and checked during each call.
The previous energy-only implementation explicitly disabled forces, so these
arrays cannot be recovered from its saved output.

Use the exact preparation 12469 XYZ, physical-source and cap-map pins already
used by 12471. Target index is 0, charge −1, physical multiplicity La 1 / Dy 6.
No changed coordinates, protons, caps, electronic method, state, source choice,
water, solvent, field, optimization or DFT-guided direction selection.

## Reused native method and admission

Use the actual job 12500 whole-peptide force executor as the source precedent.
Its warm MACE OMOL float64 model-loading, native batch/state checks and
`compute_force=True` forward call are retained. Change only the declared
representation to the four matched cores, with explicit 48-atom/target-0/
charge-−1 checks and complete task identity against the old energy receipt.
Load one model and evaluate at most four states, in the original task order.
Preserve the exact checkpoint, software versions, physical spin input,
native head, nonperiodic vacuum model, units and unchanged parameters.

The energy replay tolerance is **1e-7 eV absolute**, copied from the accepted
12500 config. Bind that config and result and require the same value before
calls. Each new force-enabled energy must match its own saved 12471 energy
within this existing tolerance. Require finite energy and a finite 48×3 force
array, exact current native input-state metadata, and unchanged model parameters.
Retain each original energy, replay difference and pass/failure. Keep all four
slots; after an execution/admission failure preserve unattempted slots under
the original executor's stop policy. No retry or tolerance adjustment.

## Output and interpretation

Save the full physical-order force array in eV/Å and the target-index-0
Cartesian vector, with the convention force = negative energy derivative.
No direction ranking, fitting, force-error threshold or geometry selection.
Future saved comparison against the same four admitted PBE0 cells is descriptive:
PBE0's saved gradient in kcal/mol/Å converts to force in eV/Å by negating and
dividing by the established 23.06054783061903 kcal/mol per eV. Equivalently,
raw Hartree/bohr gradients use the reference's 627.509474 kcal/mol per Hartree
and 0.529177210903 Å per bohr before that conversion. No force tolerance is
selected. DFT force consistency and physical Dy-state accuracy remain
unqualified; energy replay does not validate either electronic model.

All work runs in one shared GPU allocation, at most 12 CPUs, with threads
bounded by the actual grant and no memory/time overrides. Record model-load
time, four actual call timings, energy replay, GPU/runtime identity and scheduler
cost. No new DFT/GFN calls; the original four-cell DFT matrix is unchanged.
