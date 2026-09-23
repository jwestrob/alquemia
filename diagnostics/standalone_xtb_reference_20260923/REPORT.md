# Standalone static scoring preserves reference fidelity; adaptive scoring fails

**Fixed-origin native OMOL plus standalone GFN2(ALPB−vacuum) retains all25
canonical reference calls and all3 consumed crystal calls.** Its separately
frozen canonical gap is9.759548901model kcal/mol. The same backend combined with
the old adaptive proposals has overlapping canonical classes (gap−1.112975495),
so its decision bands remain null. No threshold adjustment rescues that branch.

All336 logical solvent cells completed:300 fresh standalone calls plus36 exact
pilot reuses. Zero new MACE, DFT, geometry optimization, protonation/water changes
or scientific failures. The old native baseline and all outputs remain intact.

## Separate calibrations, frozen before any structural-fold transfer

Both variants use exactly the same electronic backend/state. Only the offered
nuclear geometries differ: fixedq0 versus the common origin/adaptive_Ca/adaptive_La
pool. Each reference uses only the designated25 canonical members and the same
class-extrema/minimum-gap rule. The three crystals and all noncanonical folds
are excluded from calibration.

| Variant | Ca maximum | La minimum | Gap | Reference status |
|---|---:|---:|---:|---|
|Static origin|−405466.16155897034|−405456.4020100691|+9.75954890123|available|
|Minimal adaptive, operational|−405454.13384922856|−405455.2468247238|−1.11297549523|unsupported separation|
|Minimal adaptive, mathematical|−405454.13384922856|−405455.2468247238|−1.11297549523|unsupported separation|

Values are raw modelkcal/mol, not affinity/free-energy or a universal zero.
A25/25 calibration fit follows the separating extrema rule; it is not25 new
independent tests. The three previously consumed crystal transfers also remain
correct under the static variant's own frozen bands:

| Crystal | Class | Static R | Own-reference decision |
|---|---|---:|---|
|1H4I|Ca|-405501.74774034|Ca-supported|
|4MAE|La|-405440.66787689|La-supported|
|1KB0|Ca|-405476.93515906|Ca-supported|

Adaptive raw scores are still available, but their classifications are
**unavailable because calibration fails**, not because a molecule failed to run.
Old native/reference scores remain alongside every new row in REFERENCE_v1.json.
No old band is passed off as compatible with this backend.

## Why the adaptive branch fails this necessary test

Ca-referenceQ9Z4J7 andLa-referenceC5AXV8 set the overlapping extrema. Both metal
rows select the same adaptive_La candidate in each case; this is not a mismatched
candidate pool or selective metal geometry choice.

| Source | Ca work fromq0 | La work fromq0 | Shift in R |
|---|---:|---:|---:|
|Q9Z4J7|−3.754014|−15.781724|+12.027710|
|C5AXV8|−3.783692|−2.445922|−1.337771|

ForQ9 the nativeMACE work isCa−6.825087/La−15.297848 and standalone solvent work
Ca+3.071073/La−0.483876kcal/mol. Thus the combined representation strongly rewards
La's accommodated pose even in this Ca-class reference. These components explain
the arithmetic failure; they do not identify a unique physical omission or prove
that the proposal itself is a valid protein thermodynamic minimum. We do not
change modes, domain, solvent terms, candidate inventory or labels after seeing it.

The earlier four-source pilot's useful relative ordering did not establish
canonical-wide adaptive fidelity. This full reference test supplies the missing
check. The minimal-adaptive standalone branch stops here; no2480-cell static-plus-
adaptive expansion is launched.

## Backend, execution and tests

Same installedxTB6.7.1,accuracy0.02,pinnedGFN2 parameters,fresh/no-restart,
electronic300K,500SCC maximum,exact charges/zero unpaired electrons. ALPBwater
uses its pinneddefault solvent298.15K/gsolv/P16/GBOBC/230-point surface withH-bond
term and no ionic screening. Both media come from the standalone executable;
all solvent components remain in its printed totals. Actual runs converge in
23–156iterations. The earlier224-cell accuracy/force qualification is supporting
numerical evidence; no extra per-source accuracy sweep or retry was performed.

Job1210476 ran95s on64CPUs/128GiB: **6080 allocated core-seconds, zeroGPU**.
Its concurrent molecular executor took75.1378s (4808.82core-seconds); preflight,
collection and startup contribute the rest. BatchMaxRSS822204KiB is the scheduler
reading. Local setup/analysis and historical reused computations are additional.
These are development throughput measurements, not a cold production speed claim.

Four real-artifact tests pass after completion, zero skips: exact28/336/36/300
scope; designation and source/state reuse; malformed actual-label rejection;
and exact frozen calibration algebra on all real endpoints. The molecular-result
test was explicitly skipped before execution. No successful output was fabricated.

**Next authorized step:** STATIC-only full225 transfer, retaining208physically
prepared q0sources and17preparation exclusions. Actual MACE origins are available
for all416endpoints, including the source whose old nativeGFN failed; that old
solver failure must not prevent this different backend's static test. Use the
new static reference without refitting; keep static/adaptive/native/DFT coverage
and matched comparisons separate. No new MACE or optimization, and no standalone
adaptive decision or default promotion.

Frozen reference: `workspaces/standalone_xtb_reference_20260923/run_v1/REFERENCE_v1.json`.
Actual collection/costs: same directory,`collection_v1.json`,`COSTS.json`.
[Commands](COMMANDS.md) and [plan](PLAN.md). Transfer inventory is read-only;
a separately named static manifest must precede its authorized execution.
