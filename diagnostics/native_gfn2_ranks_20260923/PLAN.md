# Matched native GFN2 rank-count experiment — approved

2026-09-23. Root, under Jacob's standing autonomous discriminator-improvement
authorization, explicitly approved: "the24-call matched native-GFN2 rank experiment
exactly as proposed: Q9(145atoms)/Q89(184atoms) actual adaptive_Ca geometries, both
metals/media, fixed1/4/8 ranks ... Gates .1kcal each cell and .2pairedR ... no
per-cell favorable-rank selection or production change."

Question: can the native solvent stage of the stronger adaptive discriminator
use fewer MPI ranks without material numerical changes, and lower its actual
latency/allocated CPU cost? The read-only inventory found512 relevant receipts,
all using eight ranks; no matched1/4/8 experiment. Historical global-protein
memory failures are not a scaling study of these contexts.

Freeze the eight actual precision-pilot adaptive_Ca cells for canonical Q9Z4J7
(145atoms) and Q89GY2 (184atoms): each common geometry evaluated as Ca/La and
vacuum/ALPB. Same source XYZ bytes, charges-3/-2, singlet, ORCA6.1.1 binary,
native GFN2 mixer, NoAutostart, MaxIter500, SmearTemp300, default native parameters,
water ALPB, maxcore2000 and output recipe. Existing native MACE energies are reused
only to report the common composite contrast. No MACE, DFT or geometry call.

Exactly three rank settings1/4/8, eight cells each:24fresh native GFN2 calls.
The runtime renderer changes only PAL rank count. Preserve existing no-binding
MPI launch and one OMP/MKL thread per rank; retain inherited OpenBLAS settings.
Fresh eight-rank cells are the comparator, with old archived energies retained
separately. Do not select a rank's lower/favorable energy as a scientific result.

Three sequential CPU-sharing jobs on node-224-2t-8gpu-1, partitiongpu, no GPU:
eight concurrent workers and8/32/64allocated CPUs,16/64/128GiB host memory. Order
1,4,8 is fixed before output; afterany dependencies preserve all three tests if
one rank setting fails. This avoids self-contention between rank configurations,
but external host load remains a limitation. No custom wall/CPU cutoff. Do not
alter the active225 jobs, their science, thresholds, defaults or shared runner.

Predeclared acceptance: every native vacuum/ALPB energy within0.1kcal/mol of
fresh eight-rank; each source composite Ca-minus-La contrast within0.2kcal/mol.
All cells must converge with matching exported parameters, state and actual
rank/host receipts. Missing/failed cells remain unavailable. Record separate
vacuum, ALPB and transfer contributions, iterations, all raw energies, receipt
wall times, rank-seconds, batch throughput, full allocation cost and memory.
No favorable component cancellation can excuse a failed energy gate. This is
numerical execution qualification, not a new biological validation/calibration.

Report exactly what the24calls show before considering integration. No automatic
retry, extra rank setting, continuation, backend change or production rollout.
