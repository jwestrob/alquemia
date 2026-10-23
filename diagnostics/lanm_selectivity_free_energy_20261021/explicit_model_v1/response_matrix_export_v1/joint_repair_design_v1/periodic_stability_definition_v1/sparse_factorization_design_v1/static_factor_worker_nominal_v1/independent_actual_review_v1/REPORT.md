# Static attempt: technical failure before numerical factorization

The single attempt reached the symbolic constructor, where Cython rejected read-only CSC index buffers. It completed no symbolic object, numerical factorization or solve. No backend factor state or ten-probe result exists; none is inferred. The actual nominal shift and all input/source/release pins match the frozen contract.

The attempt used2.791825 s wall and3.971992 child CPU-s, peak386708 KiB, with12 BLAS/1 OpenMP threads and the declared32 GiB/900 s limits. All outputs are preserved. This failure says nothing about real-split positivity or full PME stability.

A separate source-only recovery may copy the mapped indices and column pointers into writable C-contiguous arrays, preserving all values, hashes, matrix, nominal shift and backend settings. It requires a new release; the active calibration retains its13-thread reservation.
