# Two self-continuations settle both tested geometry-sensitive energy jumps

All eight calls succeeded with actual native XTBRESTART confirmed. All four
cells pass the frozen0.1 kcal/mol stage1→2 energy-settling gate; both pairs of
near-identical geometries agree within0.1 kcal/mol after the second pass.
The large cold-start energy jumps therefore have an affordable demonstrated
remedy on these two consumed examples. No classifier score or reference changed.

## Actual result

All values below are kcal/mol; geometry work means new−old at identical metal,
composition, charge and native recipe. Both sources are La/vacuum, charge−1,
multiplicity1, run at one rank.

| Source / geometry pair | Cold geometry work | After pass1 | After pass2 |
|---|---:|---:|---:|
| Q4W6G0 Ca-conditioned sample2 / adaptive_Ca |+10.64877243|−0.00141690|−0.00021469|
| P38539 La-conditioned sample2 / adaptive_La |−11.07967010|−0.00025345|−0.01794691|

| Source / geometry | Pass2−cold | Pass2−pass1 | SCF cycles pass1 / pass2 |
|---|---:|---:|---:|
| Q4 / old |−0.50634900|−0.00151967|48 /21|
| Q4 / new |−11.15533611|−0.00031746|36 /21|
| P385 / old |−11.06221856|+0.01779417|41 /20|
| P385 / new |−0.00049537|+0.00010072|24 /28|

Pass2 is reported uniformly, including the slightly higher P385 energies;
there was no favorable-energy stage selection. The maximum repeated-energy
change is0.01779417 kcal/mol. The archived cold results and original2×2 matrix
remain untouched. No cross-geometry seed, third pass, new MACE, optimization,
DFT, altered threshold or recalculated biological classification was used.

## What this establishes, and what it does not

The preceding2×2 test reproduced each cold result at one and eight ranks.
Together the experiments show that the earlier10–11 kcal/mol changes followed
geometry-sensitive electronic initialization, and that uniform self-continuation
removes those large differences here. This strengthens the practical case for
testing a uniform continued-energy protocol separately from the cold-start
protocol. It does not support replacing only unfavorable cells in an old pool.

**Energy repeatability is the passed criterion.** All eight outputs claim SCF
convergence and satisfy the printed energy-change tolerance, but all eight
exceed their printed MAX-density and RMS-density tolerances. The actual residuals
are retained in the collection. This is not full density convergence, analytic
force qualification, proof of a global electronic minimum, or evidence of
improved biological accuracy. No new full-pool score is available from four
vacuum cells. Prior native derivative failures are not overturned.

## State, execution and cost

Each pass1 reused its own exact cold GBW+xtbw pair, with full XYZ bytes preserved.
Each pass2 used only its confirmed pass1 pair. All receipts confirm one rank;
native mixer, source parameter-export hash, charge, multiplicity and electron
count are unchanged (Q4:452 effective electrons; P385:480). Inputs retain
native GFN2, electronic300K, MaxIter500 and the archived ORCA executable. Only
the verified matching-basename AutoStart activation differs from cold input.
Seed hashes before and after execution are retained independently for each pass.

CPU-only job1211126 completed in **23 seconds on four allocated CPUs =92
allocated core-seconds**,16GiB requested, zeroGPU. Eight actual engine runtimes
span4.593–7.807s and sum46.687s; calls overlap across four workers. The23s job
includes both stages, validation, seed preservation and collection. Local
preparation/test/report work and the reused cold calculations are separate.
There were no molecular failures or retries. Prepared v1 was never submitted;
v2 additionally supports an entirely failed previous stage as four unavailable
rows without silently dropping the denominator.

Five real-artifact tests pass, zero skips in the final run. They check actual
source states/seed pairs, exact coordinates and finite manifests, rejection of
cross-geometry seed metadata, actual two-pass marker/algebra, and missing-cell
accounting using an explicitly corrupted copy of a real collection. Prelaunch
tests explicitly left the unexecuted scientific result untested.

Authoritative result:
`workspaces/precision_geometry_continuation_20260923/run_v2/COMPARISON.json`.
Both stage collections, unrounded energies, outputs, parameters, densities and
execution receipts are pinned there. The bounded diagnostic is complete; no
additional molecular calls are scheduled by this branch.
