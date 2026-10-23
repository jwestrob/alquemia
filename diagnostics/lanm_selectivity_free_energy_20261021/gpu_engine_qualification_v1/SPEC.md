# Tinker9 GPU implementation qualification

Declared 21 October 2026, before GPU execution. Root owns build and submission.

## Question and scope

Does the isolated double-precision CUDA Tinker9 build reproduce its native
periodic/nonperiodic polarization-pair references and the four existing native
CPU La/Dy water/acetate controls, with internally consistent analytic forces?
This is implementation parity. Passing does not establish parameter accuracy,
protein sampling, experimental affinity, or a realistic exchange free energy.

Source is Tinker9 commit `44bcd7c898e68f6827c112df79b8361e443730d2`, with the
recorded Tinker and Catch2 submodule commits. Root owns the isolated build under
`gpu_engine_feasibility_v1/build_v1`; this task changes no build/source files.
The declared build uses CUDA, double precision and compute capability 90. Root
must bind the completed build evidence, CMake cache, compile commands and both
executable hashes in a separate `CONFIG_BOUND.json` before execution. The
configuration template cannot execute GPU work.

The paused draft is preserved under the matching workspace's
`draft_before_resume_v1/`. Build-toolchain recovery may add explicitly reviewed
standard-header includes to Catch2. Root must pin the complete tracked diff for
each repository/submodule as well as its base commit; the adapter rejects any
later difference. This admits documented build compatibility patches without
pretending the source tree is pristine or changing scientific test assertions.

## Frozen sequence and admission gates

1. In an empty working directory, run `all.tests '[polpair]'`, with the native
   Catch2 XML reporter. Admit exactly `PolPair-Ewald`, `PolPair-NonEwald`,
   `PolPair-Ewald-A`, and `PolPair-NonEwald-A`, all successful with nonzero
   assertions and zero failures. Retain the original source assertions and
   fixtures: native energy/gradient tolerances 1e-4 and virial tolerance 1e-3.
2. Only after those tests pass, run four serial commands
   `tinker9 testgrad system.xyz Y Y 0.00001`, in this order:
   `acetate_bi_central_La`, `acetate_bi_central_Dy`, `water_central_La`,
   `water_central_Dy`. Copy the corresponding saved CPU XYZ, key, and parameter
   files byte-for-byte. Retain all La/Dy POLPAIR and VDWPR overrides. No input
   parameter, damping, precision, or solver change is permitted.
3. Require one finite total energy, complete indexed analytic and numerical
   Cartesian gradients for every atom, both total-gradient summary lines, and
   no native fatal/nonconvergence/exception/nonfinite marker. A zero native
   return code alone is insufficient: the source launcher catches exceptions
   and can still return zero.
4. For each control, require absolute GPU-minus-CPU energy difference <=
   **1e-4 kcal/mol**, maximum absolute force-component difference <=
   **1e-3 kcal/mol/Å**, and maximum absolute analytic-minus-numerical gradient
   component difference <= **1e-3 kcal/mol/Å**. Force is minus gradient; every
   Cartesian component of every atom participates, including hydrogens and
   metal. Report component RMS with denominator 3N, alongside the maximum.

Tinker9's testgrad argument 0.00001 Å is the **full finite-difference
separation**: native positions are x ± 0.000005 Å, with denominator 0.00001 Å.
The native printed decimal representation limits this comparison; no precision
is fabricated. The CPU controls are implementation references, not physical
truth. Suite reference tests are NaCl cases; metal-specific coverage comes from
the four CPU comparisons.

## Execution and preservation

One Slurm job, partition `gpu_h200`, QoS `standard_expand`, one H200 GPU, one
task and four allocated CPUs, default memory and walltime (no overrides).
The exact CPU keys retain `openmp-threads 1`; execution environment caps
support-library threads at one. Controls run serially on the allocated GPU.
One suite invocation plus four testgrad invocations includes many native
energy evaluations; it is not reported as five molecular energy calls.

Allow up to 900 s for the native suite and 300 s for each control invocation.
These bounded subprocess limits do not change the scheduler's default walltime.
On timeout, native error, failed gate, changed pin or malformed output, preserve
all output and stop. Remaining rows stay explicitly unrun. No retry, alternate
backend, weakened gate, replacement fixture or production promotion follows.
An explicitly separate root decision is needed for any later recovery.

The adapter records exact commands, input/executable/source hashes, environment,
host, Slurm allocation, GPU inventory, raw outputs, elapsed wall time and child
CPU time. Allocated CPU/GPU seconds are reported separately from measured child
CPU time. Source, inputs and executable pins are rechecked after execution.

## Preparation validation

Saved-only tests exercise the four actual CPU outputs, exact reference parity,
native exception text with otherwise complete output, missing/duplicate rows,
nonfinite values, failed force/energy/FD gates and Catch2 success/failure XML.
They invoke no native binaries or GPU and are not a qualification result.
