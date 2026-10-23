# PME snapshot exporter: source handoff

The isolated helper and patch are ready for source review. They have **not been compiled or applied** to the frozen periodic prototype. Its existing backend/stability fatal guards remain intact. The [contract](CONTRACT.md) defines the exact initialized-state hook, output, bounds and exclusions.

`PME-NORM-EXPORT 1` requests only the first initialized `udirect1` state per process. The hook is after all actual qfac entries, including the zero mode, have been assigned and before the backward permanent-field FFT. The helper does not touch qgrid or dipoles. Without the keyword it returns before array access/allocation/output. It writes one new file and refuses a collision; later initializations do not reuse a result or receive any stability admission.

The file supplies the actual full qfac, current active atom/alpha mapping and compact spline/igrid/Cartesian-transform metadata. Both native source and the saved replay compute scatter 1-/infinity-/Frobenius norms and the nominal reciprocal norm bound with the mandatory grid-size factor M. Negative or conjugate-unequal multipliers are reported unchanged. Active alpha below the native PCG floor, unsupported stencils/chunks, malformed domains and reporting-capacity violations are rejected. The reporter verifies actual pmetable ownership of every active stencil point.

The independent Python reader enforces record order, atom/spline/grid uniqueness, finite quantities, complete END counters and no trailing data. A complete file is necessary but does not supersede future native exit/error, input/binary/source, initialization-state or numerical qualification joins.

## Artificial checks

All31 checks pass in `controls_v2`: dense versus streamed scatter norms, alpha scaling, a tiny explicit Fourier operator upper bound, even/odd and negative qfac handling, wrapped stencil identity, strict record-failure cases, source hook/patch/absence/collision checks, fixed-form line lengths and preserved periodic guards. Ordinary-double allowance is the prospectively declared1e−12. No factorization backend, physical geometry, native binary or compiler was called.

The original `controls_v1` failed because NumPy2's scalar repr made the artificial text writer emit `np.float64(0.0)`. The parser correctly rejected that token. `TEST_SNAPSHOT_v2.py` changes only that writer to `float(... )` before repr; helper, patch and parser are unchanged. Both source/result pairs are preserved.

Costs, excluding startup/source reading: failed v1 0.01934323 s wall /0.01687695 s CPU; passing v2 0.12694351 s wall /0.12492991 s CPU. Each used one numerical thread with all numerical thread environment variables set to1. No actual native snapshot or bound is available yet.

## Remaining boundary

The norm/error quantities are ordinary-double diagnostics of the intended stored-weight operator. This handoff adds no automatic evenness averaging, interval-proof prerequisite, stability fallback or sampling release. Actual compiler/PME initialization/scatter-gather/FFT parity and error allowances remain to be qualified separately before a static real-plus-reciprocal stability test can consume a snapshot. The exporter is not a periodically refreshed MD certificate.
