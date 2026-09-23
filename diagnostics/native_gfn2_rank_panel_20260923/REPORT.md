# One-rank native GFN2 qualifies for the fixed union reference pools

All 336 scalar cells and all 28 selected pool contrasts pass the predeclared
numerical gates against their actual eight-rank counterparts. Both metal rows
select the same candidates in every pool. This supports a less expensive scalar
execution profile for subsequent candidate work; it does not qualify native
gradients, change a Hamiltonian, or promote a production default.

## Numerical result

| Check | Actual result | Frozen acceptance |
|---|---:|---:|
| Native vacuum/ALPB cells |336/336 complete; 0 failed|All required|
| Maximum component difference |0.00000780932846 kcal/mol|≤0.1 kcal/mol|
| Maximum fixed-geometry Ca−La difference (84 geometries) |0.00000780931441 kcal/mol|Reported separately|
| Maximum selected-pool Ca−La difference, either policy |0.00000000192085 kcal/mol|≤0.2 kcal/mol|
| Metal-row candidate changes |0/56 for each selection policy|Reported separately|
| SCF cycle-count changes |0/336|Reported separately|

The largest component difference is Q60AR6 origin La/ALPB (100 SCF cycles for
both rank counts). Mathematical minima and the operational 0.1 kcal origin-retention
policy both pass. No favorable component cancellation excuses a failed cell:
all 336 component gates pass independently.

Strict replay of the unchanged reference yields 24 correct / 1 inconclusive among
the 25 calibration sources and 3/3 correct consumed crystals. C5AXV8 is the sole
changed call: its eight-rank score defines the La minimum exactly at
−405456.46881754074; the one-rank score is −405456.4688175414, only
6.40284e−10 kcal/mol lower. Its strict call is therefore inconclusive. This is
retained as an exact-band boundary effect, distinct from numerical qualification.
No threshold tolerance, new calibration or favorable rank selection was added.

## Actual scope and provenance

Exactly 28 frozen union-reference sources, each with origin/adaptive_Ca/adaptive_La,
both metals and both media: 336 fresh single points. All contexts contain 145–202
atoms. Every input/XYZ is copied byte-for-byte from its actual archived eight-rank
execution; all reused MACE coordinates match exactly. The original 284 explicit
MaxIter 500 and 52 default 125 recipes are preserved. Same ORCA 6.1.1, native GFN2,
300 K, NoAutostart, maxcore 2000, exported parameters, charge/spin, geometry and
water inventory. Only runtime rank count changes. No new MACE, DFT, gradients,
optimization or molecular retries.

The eight earlier one-rank pilot cells had no exact geometry/input/state matches;
none was reused. `run_v1` was an unexecuted preparation attempt: the adapter first
expected the newer native worker's `complete` status, then was corrected to
recognize the older successful `computed` status. No energy or state changed.
The executed immutable `run_v2` manifest SHA starts `41a70593`.

## Measured cost and tests

Job 1211024 completed normally on node-224-2t-8gpu-1, 32 one-rank workers / 32 CPUs,
64 GiB requested host RAM and no GPU. Whole allocation: 281 s,
**8,992 allocated core-seconds**. Distinct batch/extern accounting reports
1,732.812 CPU-seconds; parent and child accounting are not double-counted.
Reported batch MaxRSS is 5,833,908 KiB. Scientific executor wall 200.186 s,
median cell 15.887 s, summed cell wall 5,772.495 s; remaining allocation includes
staging/validation/collection. Local preparation/tests/reporting are additional.
This panel is not a matched-hardware speed comparison with historical runs.

Five tests pass, zero skipped in the final run: actual 336 source/input invariants,
actual archived energy/pool replay, deliberately corrupted manifest rejection,
real prelaunch missing-result handling with all 28 denominators, and the actual
completed 336-cell comparison algebra. The scientific-result test was explicitly
skipped before execution. No simulated successful molecular output was used.

## Recommendation and next operation

Use one-rank scalar execution as a separately versioned research profile where
the same native inputs apply. Keep its cache/resource identity explicit and
preserve existing eight-rank histories. This test adds execution evidence, not
independent biological validation or improved affinity accuracy. Production and
the older running experiments were unchanged.

All raw records are under
`workspaces/native_gfn2_rank_panel_20260923/run_v2/`; see `ARTIFACTS.json` for exact
pins and `COMMANDS.md` for the ready report/validation commands. Parent coordinates
the separately declared full 225 precision candidate; this branch launches no
additional chemistry.
