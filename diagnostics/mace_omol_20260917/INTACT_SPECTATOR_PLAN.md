# Disconnected charged spectator check

Declared2026-09-17 before preparing or evaluating these states. Jacob's active
MACE goal and blanket research authorization apply. The existing five-case
results are consumed development evidence; no original score, label or criterion
changes. Production and the running canonical experiment remain unchanged.

The native model conditions every atom on the total graph charge. Matched
bound-minus-detached subtraction cancels its additive charge readout, but may
retain charge dependence in the local nonlinear interactions. Ask whether an
essentially noninteracting remote ion changes the coordination descriptor.
This is a representation consistency diagnostic, not another biological test.

## Fixed physical construction

Use exactly all20 primary bound/detached endpoints from the five-case
intact_report_v1. Append one explicitly constructed Na+ spectator to each
endpoint, at x=max(x of all atoms in that endpoint)+10000 Angstrom, with y/z
equal to the original bound metal coordinates. For a paired bound/detached
comparison the spectator position must be identical: use the maximum x over
both states to define it. Preserve every original atom, coordinate, explicit
water, cofactor, protonation, selected metal index and spin multiplicity. Add
one to the exact original total formal charge. Na+ contributes10 electrons;
closed-shell parity and the La/Ca charge difference remain valid. This sodium
coordinate is a declared numerical perturbation, not a deposited ion or a
claim of experimental salt occupancy. Do not change the protein's charge.

Record its zero graph edges and minimum distances. The native6A finite-range
model has no direct interaction with this spectator. In a physical Coulomb
comparison, moving the selected metal by at most300A changes its interaction
with a unit ion at>=10000A by less than0.001kcal/mol. Check the actual geometric
bound rather than assuming this displacement bound. The relevant matched
score cancels the remaining common spectator and protein-only energies.
There is no invented density, energy, label, solvent correction or neutralization.

Na+ and the unchanged protein charge define the intended separated-fragment
state. This checkpoint receives only a total charge and does not predict or
constrain fragment charges. Consequently this test also exposes that unresolved
representation ambiguity; it is not proof of an error in an exact unconstrained
vacuum ground-state energy. It cannot certify the detached Ca/La ionic charges
either. Keep this limitation when interpreting success or failure.

## Method and checks

Same pinned100M OMOL checkpoint, float64, no gradients,1024edge/1024product
batching and existing one-A5000 runner. Keep original full-system Q dependence;
do not remove/rewrite the charge embedding. Total charges after this Na+
addition remain within reported -10..+10 training range for these five cases.
All proteins remain outside reported training sizes.

Twenty new energy-only calls, plus20 explicit earlier baseline endpoint reuses
after actual receipt/geometry/state verification. Compute the same four-endpoint
R_coord per case. Preserve embedding/node component accounting and all raw
energies. Compare each modified R_coord and all three cross-protein contrasts
to the originals. The frozen consistency gate is absolute change<=0.1kcal/mol
for every case and comparison: ten times the existing0.01 numerical tolerance,
well below the smallest9.03kcal/mol development margin. This is a consistency
gate, not a fitted classification criterion. Retain failure without changing
spectator location, charge, model, numerical settings or accepted tolerance.

If this fails, report sensitivity to disconnected charged spectators under
this representation. It does not erase the observed five-case ordering, and
is not evidence of an actual sodium effect on affinity. The charged fragment
allocation problem is distinct from ordinary geometric locality. No threshold
refit or automatic production promotion follows from either outcome.

## Cost and artifacts

Five source proteins, one spectator variant,4calls each=20calls, zero DFT,
training, solver, force or trajectory evaluations. Existing measured endpoint
times suggest several GPU minutes plus preparation/validation overhead. One
A5000/16CPU/64474MiB, expandable-segments allocator, existing scheduler QOS;
no project CPU/time budget. Save immutable input/implementation/task manifests,
execution receipts, reused sources, actual costs and report in a new workspace.
