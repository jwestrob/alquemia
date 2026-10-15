# Conditional Asn290 fixed-flip solvent continuation

This implements the continuation already declared in `../ASN_ROTAMER_PLAN.md`.
Root owns execution. The geometry and vacuum stages are immutable inputs; this
adapter changes neither coordinates nor the physical state and performs no
MACE call or optimization.

## Gate and finite scope

Bind only an actual terminal MACE collection with all twenty new cells complete,
all ten sources retained, and the recorded sign-blind mathematical common3→4
paired-change gate of at least 1 kcal/mol passed. Slurm preflight independently
replays that gate from the actual saved MACE receipts before any SCC call. If the
gate closes, leave this continuation unbound and unexecuted.

Use every admitted fixed180° Asn flip, not just the source triggering the gate:
ten conditioning/sample sources × two metals = twenty new204-atom endpoints.
Retain the short physical H contacts represented in the fixed204-atom context;
there is no contact repair, relaxation or further rotamer search. Ca−4 and La−3
remain singlets with the same dry oxidized PQQ3− and neutral Asn state.

Reuse all sixty original common3 cells. Before new calls, revalidate their native
ALPB input identities, component-v2 acceptance, water-diagonal replay, charge
precision, fixed cavity and saved MACE pairing. Reproduce the authoritative
original three-policy pools at the existing1e−8 kcal/mol algebra tolerance.
The recorded states are204 atoms/636 GFN valence electrons for both metal states;
the actual coordinate inventory is rechecked inside Slurm.

## Unchanged energy expression and implementation

Use the existing float64 OMOL vacuum energy plus native ALPB evaluated with the
GFN2 **vacuum** charge distribution. Retain all four native terms: Gelec, Ghb,
Gsasa and Gshift. The qualified standalone xTB6.7.1 parameter file and independent
tblite runtime are pinned through the existing runtime record. Internal native
ALPB uses the already qualified factory/binary/build. The corrected full ORCA
La parameter JSON and its shell-assignment policy are pinned as lineage evidence;
ORCA is not an executed backend here and no uncorrected native ORCA parameter
default is substituted.

The dispatcher is a focused copy of the existing `MATRIX.solvent`: the task
protein metadata becomes `A0ACD6B9F2` instead of its hardcoded ExaA label. Scientific
calls remain `execute_point`, `component` and `run_task` from the pinned existing
implementations. Explicit guards replace equivalent assertions. There is no new
solver, retry policy, tolerance, environment or component calculation.

Each new endpoint uses vacuum/water in both existing SCC backends for component
qualification, followed by native evaluations at saved water/vacuum charges.
Maximum:80 SCC calls and40 native component calls; sixty existing endpoints are
reused. Seven independent one-thread SCC workers run in a real seven-CPU Slurm
allocation. Native replays are serialized. Actual starts, terminal executions,
resource receipts and allocation time are retained separately from these maxima.

## Comparison and failure semantics

The original three geometries and the one flipped origin form a common4 seen by
both metals. Preserve original-first tie order, mathematical minima, the existing
strict−0.1 kcal/mol operational origin-retention rule, and origin as a separate
policy. Report per-metal and Ca−La origin-relative changes, each solvent term,
selected geometries, and full five-source group ranges with matched coverage
explicit when a new cell fails. Reproduce old pools before and after execution.

A failed required new cell makes its entire common4 unavailable. The original3
baseline remains separately named; no successful subset masquerades as common4.
There are no class bands, labels used as admission gates, recalibration, affinity
claims or production changes. A small fixed-flip effect does not exclude all Asn
conformations; a large effect does not validate this crowded, unrelaxed rotamer.

Fresh output directories and workspace locking prevent accidental reruns.
Original absolute-v1 parity failures stay visible beside component-v2 acceptance.
Partial native records and final execution receipts survive collection failures.
All scientific parsing, replay and calculations occur in Slurm; metadata binding
and static source inspection alone are permitted on the login host.
