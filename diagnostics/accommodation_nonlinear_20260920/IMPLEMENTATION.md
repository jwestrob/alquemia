# Nonlinear pilot implementation/preflight

The frozen PLAN.md is unchanged. Source is the original four-case torsion design,
eight endpoints; no new geometry policy, coefficients or DFT tasks. Ready manifest:
`workspaces/accommodation_nonlinear_20260920/pilot_v2/manifest.json`.
The earlier pilot_v1 preparation is preserved and was never executed.

## Execution design

CPU Python owns scipy L-BFGS-B, existing Kinematics/geometry checks, existing
ORCA manifest runner and the verified native EnGrad parser. A persistent isolated
native float64 OMOL process owns the one allocated GPU; explicit pinned coordinate
requests and native scalar/force receipts cross its pipe. The worker uses one CPU
thread while two native GFN2 media run at eight ranks each through the established
runner. The GPU step uses srun --overlap within this job so retaining the warm
worker does not reserve all32CPU slots or block the job's own MPI tasks. No shared
permissions/runtime changes. Worker initialization failure terminates only its own
child step, retains logs, and writes outer execution failure/accounting.

At every q, the complete Cartesian gradient is projected through J(q), including
link-atom chain rules. The old parser's q0 projection is discarded. MACE force is
negated and converted from eV/A once; parsed GFN2 gradients are already kcal/mol/A,
so ALPB−vacuum is formed without another conversion. The relative objective
subtracts each native/solvent component's q0 value before summing, avoiding the
large native atomic-energy gauge in the stopping test.

Identical requests are cached only within a manifest/endpoint/state/q. Failures
and partial attempts remain visible and are not retried with different chemistry.
An unsupported trial stops that endpoint; independent endpoints can continue.
Final optimizer candidates remain distinct from the stated interior, gradient,
energy and geometry gate, and from the subsequent stable-curvature qualification.
Raw Hessian antisymmetry/refinement are checked before diagnostic symmetrization;
eigenvalues remain unclipped. No entropy or calibrated relaxed score is produced.

Each evaluation stores requested tasks, actual native model-start/complete flags,
ORCA execution receipts, actual native-SCF-start markers and completed analytic
endpoints separately. Nested ORCA timing records are retained for diagnosis; final
allocation cost comes from the outer receipt/scheduler and is not summed with
nested32CPU accounting. Primary-to-Tight q0 changes are retained separately from
accommodation work.

## Actual preflight (zero new molecular calls)

All eight source coordinate/charge pairs and full Ca/La physical maps compare
exactly. Active-axis IDs, moving groups and coordinate units agree. Existing
source/cap Jacobian and bond checks pass at q0 and a fixed joint+0.01radian probe.
All16origin EnGrad inputs pass the existing runtime dry-run. The parent's separate
60corner/axis geometry checks passed, recorded in geometry_preflight_v1.json.

Five focused real-fixture tests pass: paired physical maps/inactive groups, archived
MACE+nativeGFN2 gradient signs/units/current-J projection, actual archived component
works and missing solver states, rejection of nonstationary gradients/unrefined
curvature, and byte-identical supported analytic recipes. These are parser/algebra
and geometry checks, not a fabricated nonlinear scientific integration result.
Parent technical review passed and job1204162 was submitted with pinned pilot_v2.
The warm worker initialized and the first matched analytic GFN2 pair completed;
endpoint optimization is running. Full scientific qualification remains pending.

TESTS_v1.txt retains two test-harness mistakes: demanding bit equality from zero-angle
rotation arithmetic (observed2.22e-16A roundoff), and attempting composite algebra on
an archived failed point. Corrected checks use the existing1e-12A mapping tolerance
and explicitly require the failed composite to remain null. No scientific input,
score, force, selection or acceptance threshold changed. TESTS_v2.txt also retained an overly strict missing-key access in that historical
failed row. TESTS_v3.txt records the corrected five passing tests. Inactive physical coordinates have no finite motion.

After submission, the live collector adds donor-contact distances from each saved
full q through the existing torsion analyzer, and records its own analysis hash.
This is report-only: the pinned running optimizer/energy/gradient implementation
is unchanged. Final manual collection will include these fields; original automatic
outputs remain preserved.
