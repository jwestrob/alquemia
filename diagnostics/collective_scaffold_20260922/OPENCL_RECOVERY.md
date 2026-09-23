# Execution-only recovery before the first molecular evaluation

Job 1210189 finished with all 24 CUDA Context initializations failing with
`CUDA_ERROR_UNSUPPORTED_PTX_VERSION (222)`. There were **zero energy/force
requests and zero evaluations**. Every failure and the original manifest remain
under `workspaces/collective_scaffold_20260922/searches_v1/`.

The existing OpenMM 8.5.1 CUDA plugin links NVRTC 13.2. The current driver cannot
load its generated PTX. OpenMM 8.5.1's
[actual CUDA implementation](https://github.com/openmm/openmm/blob/8.5.1/platforms/cuda/src/CudaContext.cpp)
uses NVRTC; setting a legacy compiler property does not supply a validated remedy.
No driver/library replacement, ABI spoofing or environment installation is made.

This repository already records the identical failure and successful OpenCL/H200
execution in `docs/canonical_amber_1264_relaxation_pilot.md` and
`benchmarks/hans_lanm_amber_1264_v1/run_md_opencl_v2.py`. Reuse that platform route,
with **double precision** retained for this constrained mechanics experiment.
OpenCL has no separate `DeterministicForces` property: do not request a nonexistent
control or claim it was applied. Record all actual properties and require exactly
one GRES-visible NVIDIA H200 from `Platform.getDevices()`, following the existing
device-admission implementation. No automatic fallback exists.

`searches_v2` pins the new execution implementation. It repeats the same 24
declared initializations/searches, with the same parent Systems, inputs, donor
targets, optimizer, numerical guards, objective and admission criteria. No prior
molecular calculation is rerun because v1 reached none. The original CUDA plan
is preserved; this explicit addendum supersedes its platform/toggle request only.
One H200, 32 CPUs and 200000 MiB remain requested. Include the failed allocation
in final total cost; measure the first actual declared force call in v2.
