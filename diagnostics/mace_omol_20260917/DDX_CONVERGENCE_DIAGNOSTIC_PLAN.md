# Log the unchanged failed La dielectric solve

The direct-source 2FW0 coarse Ca solve completes, but La stops in the first
ddPCM dielectric linear system after 300 iterations (232.913 seconds total).
The second single-layer solve is never reached for that state. Matrix storage
accelerates that later system and is not an established fix for this failure.

Run exactly one new forward solve for the same 2FW0 La coarse state. Preserve
all coordinates, source charges, cavity radii, PCM parameters, direct source
potential, FMM orders, thread count, DIIS history and 1e-10 tolerance. Increase
only the numerical iteration ceiling from 300 to 1200 and enable the library's
iteration log. This distinguishes slow convergence from a plateau or unstable
iteration. It does not change the Hamiltonian, loosen a scientific threshold,
select a biological result, or prove a practical production cost.

Keep failed and converged outputs, actual potentials/integrals, solver log,
coefficients when available, timing, CPU and memory. Record that the native
convergence criterion is relative iterate change, not an independently
evaluated equation residual. If it still fails, report failure; no automatic
retry or arbitrary energy. Do not change or stop job1201279.

One model setup and one attempted forward solve, zero DFT/MACE/force calls.
Use the same isolated ddX0.9.0 module and 64-CPU/128-GiB allocation policy.
No project compute/time budget. This is a one-time numerical diagnosis on a
consumed state, not an ordinary proposed scoring operation. Original input
and runtime pins must match before execution. Physical-model changes remain
outside this particular diagnostic.
