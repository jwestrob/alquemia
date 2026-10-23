# Writable sparse-index recovery

The preserved first static attempt failed in the Cython sparse-buffer conversion during symbolic construction: `ValueError: buffer source array is read-only`. It reached zero numerical factorizations and zero solves. This is a representation failure, not a real-split or full-PME instability result.

The sole runner change copies the mmap-backed indices and column pointers into writable C-contiguous arrays before constructing the CSC matrix, just as the data buffer was already copied. Dtypes, values, graph and input payload hashes, nominal shift, factor backend/settings, ten probes, status rules and32 GiB/900 s/10800 CPU-s envelope remain identical. The added int32 copies require approximately65.6 MB; this does not alter the32 GiB cap or promise the full factor fits it. Original source, inputs and failed output remain immutable.

Before any new large attempt, execute one tiny serialized-buffer test in one thread: save a fixed3x3 lower CSC, reload data/indices/indptr as read-only memmaps, apply a fixed artificial diagonal shift0.125, and reproduce the old constructor's read-only error once. Apply only the new writable-index copies, verify all three payloads/dtypes unchanged, then construct one symbolic supernodal factor with the same qualified settings. Require a nonnumeric symbolic state; do not call factorize or solve. Check all file hashes remain unchanged and exact runner/config delta. Preserve any unexpected failure and cost; no retry.

The recovered full-matrix runner remains unreleased. Root may separately version scheduling metadata; this source itself preserves the prior login envelope and cannot run alongside the active13-thread calibration under the current reservation.
