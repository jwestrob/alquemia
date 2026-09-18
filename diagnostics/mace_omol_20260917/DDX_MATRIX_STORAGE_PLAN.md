# Cache the unchanged ddX operator to test affordable execution

**Status: proposed; not prepared or launched.** The subsequent 2FW0 La coarse
solve failed in the *first* ddPCM dielectric system at 300 iterations. Source
inspection shows `incore` caches the later single-layer system, so it cannot
be assumed to fix that failure. Keep this six-solve proposal unexecuted while
choosing a focused convergence/performance diagnostic.

The first real matrix-free 2FW0 Ca coarse solve takes 214.4703 seconds on
64 CPUs (81 reported iterations). It passes source and energy accounting.
Source inspection of ddX 0.9.0 shows its `incore` option precomputes the sparse
single-layer operator blocks; the alternative rebuilds their action at every
iteration. This is a supported implementation choice for the same equations.

Run one contained performance comparison while the original convergence pilot
continues untouched. Use exactly the same prepared 2FW0 primary physical source
and paired default charges at all three declared basis/grid settings. Switch
only `incore=False` to `incore=True`. Keep PCM, radii, smoothing, direct Coulomb
source evaluation, FMM 12/12, solver tolerance, DIIS and 64-thread execution
unchanged. Three model setups, six forward solves, no added biological cases,
geometry, DFT, MACE, forces or fits. Do not alter the other job or its queue.

Before execution, freeze endpoint and Ca-minus-La numerical equivalence to the
matrix-free outputs at <=1e-6 kcal. Source/energy checks remain unchanged.
Compare every state, including failures. This is an implementation-equivalence
test, not a test of scientific accuracy or a new classification threshold.
The unchanged physical model still needs the original convergence and rigid
checks. Record setup cost as well as solve cost; amortizing setup is legitimate
only for calculations sharing the exact same cavity and numerical settings.

Request 64 CPUs and 128 GiB on standard/memory, with single-thread BLAS and
one model alive at a time. Sparse block storage scales with actual overlapping
sphere pairs; record measured peak memory and CPU model. No resource/time
stopping budget. Preserve the original failed source-evaluation pilot, current
matrix-free pilot and this separate implementation manifest. No default change.
