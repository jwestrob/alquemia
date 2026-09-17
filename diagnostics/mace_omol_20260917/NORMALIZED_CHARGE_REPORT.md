# Normalized QM charges and cap projection pass the potential checks

Eight native CHELPG fits and eight independent quantum-potential evaluations
completed job1200970 without failure. All endpoint and paired potential gates
pass, including the source-graph cap projection. This supports proceeding to a
solvent-energy test; **it is not a solvent-energy or affinity-accuracy result**.

| Case | Ca-minus-La fit relative RMS | Projected relative RMS |
|---|---:|---:|
| Alpha1F6S | 0.4563% | 0.5034% |
| Alpha6IP9 | 0.4081% | 0.4647% |
| GGR connected | 0.3493% | 0.3715% |
| GGR extended | 0.5532% | 0.5852% |

These errors use identical Ca/La exterior probes and the actual saved quantum
potentials. They concern the differential field, not each endpoint's absolute
field. Endpoint fit RMS ranges0.001765–0.002339au; projected endpoint RMS ranges
0.002119–0.002924au. One projected neutral-La endpoint has relative RMS10.0724%,
so it passes the predeclared absolute-error branch rather than the10% branch.
No threshold or charge scheme was changed after examining results.

## Representation and limits

Protocol `normalized_vacuum_CHELPG_source_projection_diagnostic_v1` reuses the
eight normalized vacuum r2SCAN-3c wavefunctions. Native COSMO sampling radii,
0.3A grid,2.8A extent and no dipole constraint are verified from actual utility
output. The charge distribution varies with the metal. Charge-sum error is at
most8e-6e at native printed precision; no renormalization is applied.

Synthetic sigma-cap charge is distributed to its retained and omitted real
heavy anchors using their geometric barycentric weights. This preserves total
charge and dipole while removing synthetic charge sites. All source atoms,
coordinates and covalent connections remain fixed. Maximum conservation errors
are2.44e-15e and4.41e-11eA. Overlapping contributions add through a single mapped
matrix. Higher multipoles are not preserved, and their potential error is
reported. There is no new physical atom or solvent cavity here.

The exterior-potential criterion is the existing RMS<=0.005au OR relative
RMS<=0.10. Meeting it does not establish sub-kcal accuracy after quadratic
reaction-field evaluation. Fixed exterior protein charges and covalent boundary
closure have not yet been assembled. Solvent correction, affinity score and
calibrated class are null. Baseline/default remain unchanged.

## Execution, cost and tests

Job1200970 used **242wall seconds on8CPUs**,1,936allocatedcore-seconds,
897.082reportedCPU-seconds and1,523,816KiB peak batch RSS. No GPU allocation,
DFT recomputation, MACE forward, GB solve, optimization or training.
Summed utility wall time891.622784s: CHELPG856.323842s and vpot35.298943s.
Parallel utilities use one thread each; process startup is included. Source
GBW/densities/index hashes survive unchanged. All copies, logs and attempts
remain under the new manifest. Utility cost is additional to existing quantum
and MACE subtotals, not counted as model calls.

Four distinct real-fixture tests pass. Three preparation/projection tests took
5.094s; actual integration initially skipped, then passed0.027s against the
completed output. Tests cover actual graph overlap, source/cap corruption
rejection, paired probes and receipt-backed charges/potentials. Frozen dry-run
passes; no synthetic successful scientific output.

## Artifacts and next step

[Plan](NORMALIZED_CHARGE_PLAN.md), [commands](NORMALIZED_CHARGE_COMMANDS.md),
[compact result](NORMALIZED_CHARGE_RESULT.json). Full artifacts under
`workspaces/mace_omol_20260917/normalized_charge_v1` and
`normalized_charge_report_v1`; costs in `normalized_charge_cost_v1.json`.
ManifestSHA31caf15249cbfab1ff381f1185c7243e77c69a003c735f52ca80e5ad63bc3fc6.

Next: assemble and test a declared full-protein reaction-field model with these
endpoint charges, identical physical boundaries and a local, explicit covalent
charge ledger. The vacuum hybrid's2/4 ordering failure remains a failure. The
charge pass neither revises that result nor guarantees the next model will win.
