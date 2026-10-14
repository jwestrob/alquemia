# Reciprocal ion placement: fixed-distribution solvent continuation

Declared after all eight MACE cells completed and before new solvent evaluation.
The frozen reciprocal hypothesis passed in vacuum: control delta_D +3.373337,
target −9.328603 kcal/mol. Jacob's discretionary discriminator-improvement scope
authorizes this contained follow-up. Root owns submissions and integration.

## Question and energy

Does that reciprocal placement effect survive the currently qualified solvent
descriptor? For each metal M and unchanged geometry q:

E(M,q) = E_OMOL,vac(M,q) + G_nativeALPB(q_GFN2,vac(M,q); q).

The second term contains native electrostatic, charge-dependent hydrogen-bond,
SASA and reference components evaluated on the vacuum GFN2 charge distribution.
The solvent reaction field responds; the solute charges remain fixed during this
evaluation. This is the existing perturbative descriptor, not exact self-consistent
solvation or a binding free energy. Charges may differ between metals/geometries;
they are never transplanted from the other source or neutralized by hand.

Use the same definitions of donor work W, D=W_Ca−W_La and delta_D as PLAN.md.
Report MACE and solvent terms separately, then their sum. Keep each host independent;
no absolute comparison between source compositions. No class bands or calibration.

## Fixed denominator and execution

Retain exactly the eight prepared transplanted cells admitted in preparation4129
and evaluated by MACE4130. Reuse all MACE values and the eight exact original
solvent components only after coordinate/state/runtime/parameter compatibility
is demonstrated. Missing or incompatible originals block their comparison and
are not replaced silently. All failed cells remain in the declared denominator.

Reuse the qualified crystal-scoring solvent implementation and native replay:
published-parameter standalone GFN2 with its xTB/tblite vacuum/water component
checks, then native ALPB water-diagonal verification and vacuum-charge replay.
The finite maximum is 32 new SCC calls and 16 native ALPB replay evaluations for
eight new geometry/state cells. Exact compatible caches reduce actual calls.
Do not use ORCA's erroneous default La shell assignments, globally reorder a
parameter file, change weights or remove any established numerical gate. Pin
the actual published/corrected parameter provenance and executable hashes.

Retain identical 192 atoms, source/cap maps, dry oxidized-PQQ3− state, Ca−3/La−2,
singlets, all coordinate files and model pins. No optimization, force evaluation,
new state, folding or DFT. A failed required cell leaves its dependent response
unavailable; no fallback to the vacuum-only result as a successful solvent result.

Use shared standard Slurm within current allocation limits, a finite cell list
and explicit inputs/outputs. Keep independent endpoint work within allocated
threads and memory. Root reviews implementation and input compatibility before
execution, arms the completion wake and records actual calls, wall time, CPU
allocation and observed memory. Production remains unchanged.

## Interpretation fixed before execution

The same descriptive reciprocal rule applies to the composite: target delta_D
below −1 and control above +1 kcal/mol. Report individual works and component
changes even if that rule fails. This scale remains an interpretive resolution,
not a measured error or biological accuracy criterion. A reversal/mixed/null
result limits the mechanism and cannot be hidden by retuning bounds or bands.

If complete and supported, consider a separately declared bounded common search
that allows ion and donor accommodation under the same physical mapping. This
result alone cannot establish such a search helps discrimination. If unsupported,
close this causal branch at its observed limits rather than escalating motions.

Ownership: Heisenberg solvent adapter; Curie original-component/input audit;
Gauss independent review; root plan, execution, reports and next decision.
