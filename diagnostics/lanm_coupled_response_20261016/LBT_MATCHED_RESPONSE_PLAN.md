# Matched-solvent LBT response: finite search-objective test

Declared 17 October 2026 UTC before new geometry or molecular evaluation.
Protocol: `lbt3_matched_fixedFF_response_v1`. Root owns submission/integration.

## Question

Does allowing LBT to respond to the energy used for its reported score materially
change the Mex/LBT contrast or its dependence on the two consumed LBT sources?
Mex has been optimized under its fixed-FF solvent model; LBT12720 only rescored
geometries generated with self-consistent GFN2 solvent. This comparison currently
contains a search-objective asymmetry. It is not a new biological validation set.

The second Mex deposit6MI5 cannot supply a matched version of the current8FNS
construct, four-ion occupancy and167-water inventory without modeled additions.
That source-state experiment is deferred; see MEX_SECOND_SOURCE_FEASIBILITY.md.

## Frozen starts and energy

Use exactly the four own-metal ten-coordinate terminal geometries from the
completed12556/12558 search: geometry11/7CCO/La,13/7CCO/Dy,15/7CCN/La,
17/7CCN/Dy in the12720 inventory. Both sources and both metals remain declared.
Whole256-atom LBT, charge-2, target index255, physical multiplicitiesLa1/Dy6,
same proton inventory, no explicit water, same native MACE checkpoint/head/dtype.
The energy is MACE_vacuum + native ALPB(fixed ff19SB charges plus formal+3 target).
No external region means field and MM-ledger terms are exactly absent by partition.
Use the qualified native solvent functional/binary and exact12720 charges.

Reassemble all four starting energies and gradients from saved raw MACE energies/
forces and12720 native ALPB components/gradient.dat. Do not reuse the old POINT
composite gradient, which includes the former solvent response. Verify every
state, physical ID/order, exact coordinates, component and source pin. Starting
gradient = -F_MACE * eV_to_kcal + grad_ALPB * Hartree_to_kcal / Bohr_to_A.
No new molecular call is required to initialize any of the four starts.

## Unchanged finite search and admission

Reuse LBT3_COMMON_BASIS.py's ten physical modes, original-source metric, exact
current Jacobian, covalent/contact checks and original-source coordinate map.
No source rebasing. Source envelope0.40A, target0.20A, total torsions0.20rad,
trial displacement0.05A, unchanged four Armijo multipliers/history/curvature rules.
The existing top-two feasible-coordinate fallback is allowed only after all four
normal trials fail geometry before any energy/cache evaluation. Reset history
because the objective changed. No broader basis, fallback or envelope.

Each start gets at most12 accepted steps and24 new distinct own evaluations.
Before a new molecular call, require a physically admitted predicted own-energy
descent direction under the existing rule. A geometric/no-decrease stop before
any accepted step retains the unchanged admitted start; a later such stop retains
the last accepted point. Both remain in the denominator. A backend failure is
unavailable, preserving valid old evidence separately. Neither triggers a retry
or bound expansion. Actual energy selection uses the same Armijo rule; no score
sign, experimental preference or source gap gates execution. Report all residual
loads, stopping reasons, old/new components, rejected trials and occupied bounds.

Run the four-start zero-call artifact/geometry/force preflight and7CCO/La scout
first. Its execution/physical admission permits the other three unchanged starts.
A geometry-blocked/no-decrease scout is a valid finite result; an interface/backend
failure blocks dependent execution until diagnosed. Keep native MACE warm within
a stage. All numerical work uses Slurm with actual resource grants and completion
wakes; no memory/time overrides, new DFT, GFN SCC or model training.

## Common pool and interpretation

Retain all18 previous geometries/36 matched-model cells from12720. Add each admitted
new terminal, deduplicating exact coordinates with aliases preserved. Each metal
sees every admitted candidate: at most4 new opposite-state evaluations after
the own searches. Missing required cells make that expanded result unavailable;
the earlier pool remains explicitly separate. Rejected trial points do not become
candidates. Maximum100 new MACE and100 native ALPB calls; actual work may be much
less. This is one contained experiment, not a project-wide compute cap.

Report both source-local common pools, global common pool, old18 pool, works and
selected geometries; retain previous origin/prior/control/candidate/union columns.
Use exactly completed Mex common5 (12731) for the balanced exchange; do not repeat
Mex. Positive R=(Dy-La)_Mex-(Dy-La)_LBT favors conditional Mex-La/LBT-Dy.
Keep source hydration, partition, assay, frozen-f/reference and nonstationary
limitations visible. Candidate count is not population; energy is not affinity.

If both source-local and global contrast changes are within the existing1kcal/mol
interpretive resolution and source disagreement remains large, close this finite
search-objective question. A larger change establishes material sampling
sensitivity to the objective, not improved biological accuracy. In either case
do not automatically continue the same optimization or enlarge its bounds.
