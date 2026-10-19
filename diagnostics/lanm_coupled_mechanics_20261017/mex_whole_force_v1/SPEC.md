# Whole-Mex initial-geometry directional derivative check

Declared 19 October 2026 UTC, before new displaced evaluations. Root owns
execution under the active mechanics plan, module 1. This is an independent
implementation check while reference14795 runs; it cannot establish accuracy.

## Question and unchanged model

Does the whole2049 MACE vacuum plus fixed-FF-charge ALPB gradient differentiate
its own energy at the saved BEFORE geometry for all four existing state inputs?
The completed fixed panel supplies analytic forces and exact component closure,
but its new whole-system Hamiltonian lacks this directional numerical check.
Original qualification12952 tested a different, partitioned Hamiltonian.

Reuse PANEL.py `dd611ca03982e4f2548a88dabf9477e33aba7592fc2b41d832e60d13edd2add7`
and completed collection14794
`4c0c8b7379a2210f64267f0df94b1df72b54902d7233ec4b9614d4e941349d27`.
Retain all2049 atoms,167waters, EF3 target index2045 and three Nd spectators,
whole charge+2, exact atom order, coordinates, fixed solvent charges and model.
The four states are before_La_aligned (multiplicity10), before_Dy_aligned (15),
before_La_one_Nd_reversed (4), before_Dy_one_Nd_reversed (9).
These are aggregate conditioning hypotheses, not assigned local spin states.

No MM, embedding difference, extra dispersion, cap, SCC or DFT term is added.
Keep native float64 MACE-OMOL weights/kernel/runtime and fixed native ALPB
binary/build/library/parameter lineage unchanged. All numerical stages use Slurm.

## Frozen coordinates and denominator

Use the La vectors and displaced geometries already recorded by12952 in
`workspaces/lanm_coupled_mechanics_20261017/mex_qualification_v1/execution/`.
FD_DIRECTIONS.json SHA
`8972af45b94e6e1eefdc916030b3fa08fc1fc79a1ae581340a95da2f01c9ee47`.
Its La center is the exact current before-La POINT `f6380efc…`.
The old Dy center differs and cannot supply displaced solvent reuse.

Retain all four original sign-independent unit Cartesian directions:

1. `core_backbone`: A/84//C along the original recorded C-to-O vector.
2. `boundary_anchor`: A/83//CA along the recorded perpendicular direction.
   This is a former boundary atom; the whole-system model has no cap boundary.
3. `external_water`: recorded rigid motion of complete water A/301.
4. `spectator`: recorded A/201//ND motion.

Use both original steps h=0.0002 and0.0001 angstrom, and both signs. Bind the
actual archived coordinates, atom IDs, source ordering, direction normalization,
step/sign metadata and origin identity before any new evaluation. The identical
coordinate arrays apply to both target metals and both aggregate inputs.

Denominator: four states × four directions × two steps × two signs =64 displaced
whole-MACE calls. Reuse four central MACE records, two central solvent records
and16 archived displaced La solvent records exactly. Compute only16 missing
Dy solvent records, sharing each across the two spin inputs because fixed-charge
ALPB has no spin dependence. Maximum new calls: **64MACE +16ALPB**. No repeated
center, optimizer, new DFT/GFN/MM, new spin hypothesis or geometry search.

## Admission and derivative gates

Zero-call preflight first re-admits the exact completed centers and all16 La
solvent displacements using original native receipt/coordinate/state/charge
checks. Bind frozen implementation/configuration, manifest and reused inputs.
If exact reuse cannot be admitted, stop before new calls and preserve the reason;
do not silently regenerate centers or relax an equality check.

For each state, direction and step, compare the native central-difference energy
derivative with the saved analytic gradient projection. Subtract native endpoint
energies before conversion. Report MACE, ALPB and their sum separately, with
EV23.06054783061903 and HA627.509474; use the original gradient conversion and
Bohr0.529177210903. A cancellation between components cannot pass a failed term.

Use the unchanged12952 criterion for each term:

`abs(FD - analytic) <= 0.005 + 0.0001 * max(abs(FD), abs(analytic))` kcal/mol/A.

Both steps must pass independently. Preserve all declared unavailable slots if
a native call or gate fails; stop new scheduling at that failure. No adaptive
step refinement, favorable-step selection, precision change, model fallback,
extra center or automatic retry. Record actual attempted calls and partial files.

## Execution and interpretation

Preflight uses1CPU with the pinned scientific Python and no model load. After
actual preflight admission, use one warm native model on1H200/12CPUs and the
existing sequential fixed-charge solvent helper. No scheduler memory/time
override. Reuse the measured panel cost evidence: about1.07s per forward,
6.7s loading and81.7GB peak allocation; this is not a runtime guarantee.
Retain the existing600s per-forward and90% device-memory bounds, default scheduler
limits, actual allocation-derived settings, immutable receipts and completion wake.

A pass establishes only the sampled directional consistency at the BEFORE
geometry for both hypotheses. After geometry, unsampled directions, numerical
accuracy against an independent method, physical electronic states, converged
mechanics and natural-protein transfer remain unqualified. No optimizer, transfer
score, spin selection, failed electronic-base continuation or closed embedded-DFT
retry follows from this check. The full mechanics goal remains active and unmet.
