# Full-boundary solvent extension fails despite stable numerics

All48 actual solver calls and46 numerical checks pass. **Partition and all-four
ordering gates fail.** The opt-in OMOL-plus-reaction-field candidate is rejected
for expansion. Baseline/default, prior results and acceptance criteria remain
unchanged. The wider MACE goal remains active.

| Alpha preparation | GGR representation | Alpha minus GGR, kcal-scale |
|---|---|---:|
| 1F6S | Extended | -121.678446972 |
| 1F6S | Connected | -126.161329231 |
| 6IP9 | Extended | -120.104461970 |
| 6IP9 | Connected | -124.587344230 |

All four fail the fixed>0.02 rule. The GGR connected-minus-extended shift is
**+4.482882260**, failing2; the vacuum parent had+0.118496483. The solver term
adds+4.364385776. No modified dielectric, radius, label, threshold or selected
structure is used to rescue the result. These remain two consumed biological
groups, with different assay qualifications; alpha preparations also retain
different frozen water inventories(two/three).

## Components identify an accounting problem worth pursuing

The newly added solvent contrast is-180.639855/-173.931139kcal/mol for the two
alpha preparations, versus-54.631776/-50.267391 for GGRextended/connected.
Much of this disparity comes from the protein–QM reaction cross term. A
separately declared saved-state audit measured its vacuum Coulomb counterpart:

| Preparation | Direct cross R | GB cross R | Their sum |
|---|---:|---:|---:|
| Alpha1F6S | +132.187471 | -135.462497 | -3.275026 |
| Alpha6IP9 | +123.936457 | -128.191930 | -4.255473 |
| GGR extended | -0.233302 | -10.057807 | -10.291109 |
| GGR connected | +9.656199 | -15.356546 | -5.700347 |

These are frozen point-charge cross terms in kcal/mol, not new affinity scores.
The corresponding learned OMOL context contrasts are-58.018728/-59.417523 and
-75.629786/-51.560081. They include inseparable learned effects and cannot be
identified with the direct electrostatic term. **Inference:** adding the full
reaction field to this learned vacuum context lacks a demonstrated compatible
vacuum electrostatic counterpart. This helps explain why charge-fit quality
and stable solver numerics did not produce a useful score. It does not uniquely
assign the failure to a particular learned component or validate the classical
charges as an exact quantum environment.

The direct-coupling audit passes20 accounting/rigid checks, with no new solver
or model call. Radial components show alpha's major direct field originates
within18A. There are no alpha environment atoms beyond36A; GGR's beyond36A term
is about-13.2kcal/mol. Appending only an unquestionably nonlocal tail cannot
repair the approximately120kcal ordering failure. The18A changed-readout bound
and36A shared-readout input bound are distinct; neither proves accurate learned
electrostatics at shorter distances.

## Preparation, checks and cost

Protocol `normalized_QM_projected_ff19SB_full_OBC2_vacuum_hybrid_v1` adds only
full-boundary reaction energy to the matched vacuum DFT+OMOL expression.
No bare Coulomb term or coreCPCM is included. OBC2parameters, physical boundary,
source H/heavy coordinates, donors, protonation, assembly and water inventories
are fixed. Identical Ca/La exterior charges exclude all projected QM atoms.
Explicit residue/formal-charge and bonded-recipient ledgers close; overlapping
source contributions merge. Maximum local charge increment is0.15275e. This is
an approximate, recorded covalent-boundary model, not a unique physical split.

48calls include full/QM/environment components, complete dielectric identity,
rotations/translations, repeats and Reference-versus-CUDA comparisons. Maximum
numerical discrepancy0.000272361kcal, below0.01. All algebra/identity gates pass.
No grid is used; no-grid numerical agreement does not validate OBC against PB
or experiment. Unvalidated common1.8A metal radius and frozen density remain
limitations. Combined gradients, relaxation, entropy and affinity calibration
are unavailable.

Job1200975: **57GPU allocation-seconds,912allocatedcore-seconds,
117.953reportedCPU-seconds**, oneA5000/16CPU/64474MiB. Solver evaluations total
4.379188012seconds; peak batch RSS167,908KiB. GPU memory was not measured, not0.
No newDFT/MACE/chargefit/training. Local preparation/collection/tests and audit
resources are separately retained. Four new real-fixture tests, six existing
native regressions and two legacy solvent tests pass; actual integration
initially skipped, then passed0.004s. The coupling audit's initial preflight
rejected equal source hashes at different live/frozen paths; recovery preserves
all scientific fields and verifies both hashes, with no physics/gate changes.

## Artifacts and next action

[Plan](FULL_BOUNDARY_GB_PLAN.md), [commands](FULL_BOUNDARY_GB_COMMANDS.md),
[compact result](FULL_BOUNDARY_GB_RESULT.json),
[coupling audit](QMFF_COUPLING_AUDIT_RESULT.json).
Full results/receipts: `workspaces/mace_omol_20260917/full_boundary_GB_v1`,
`full_boundary_GB_report_v1/result.json`, `full_boundary_GB_cost_v1.json`,
`QMFF_coupling_audit_v2/result.json`.

**Retain the baseline.** This solver is numerically stable and affordable, but
the assembled candidate is scientifically unsuccessful. The next candidate
must treat direct and reaction electrostatics together and give MACE an
explicitly limited role. The previously tried MACE-POLAR short-range component
is available for that test; its old failures remain recorded, and no new
success is presumed.
