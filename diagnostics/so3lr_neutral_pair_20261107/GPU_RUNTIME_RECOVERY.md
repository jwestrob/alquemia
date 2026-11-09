# GPU runtime recovery — 9 November 2026

Job52844 finished two start0 endpoint evaluations and failed on the first
translated frame. The native error states that `jax.pure_callback` could not
find a local CPU device. The launcher had set `JAX_PLATFORMS=cuda`, excluding
the CPU backend required by GLP's neighbor-list update. This is an execution
failure, not an evaluated physical response of Tauri.

Declare one runtime-only recovery before submitting it. Enable both backends
with `JAX_PLATFORMS=cuda,cpu`; assert that the default model backend remains GPU
and that a CPU device is available. Keep model, weights, precision, charge/spin,
head, coordinates, atom order, PME, cutoffs and numerical acceptance criteria
unchanged. Do not alter the potential or neighbor-list algorithm.

Use a fresh `periodic_gpu_v2` directory. Reuse both completed predictions from
v1 after checking raw-array hashes, exact coordinates/elements, energy, dtype
and state records. Complete the remaining28 fixed evaluations and five warm
throughput calls. Preserve v1's complete native output and error. No automatic
retry follows another failure.

Use Slurm standard_expand,12 CPUs/one H200, default scheduler memory/time and
the original900-second internal bound. The first allocation cost62 GPU seconds
and744 allocated CPU seconds, including its failure. Count this recovery and
the necessary host CPUs within the existing one-H200-day envelope.

Apply the same backend availability setting to the still-unsubmitted hydrated
check and MD scout to avoid repeating this known runtime error. This enables
a required host callback; it does not change the primary model or precision.
The previous reference limitations and all scientific stopping criteria remain.
