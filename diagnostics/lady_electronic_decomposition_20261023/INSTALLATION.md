# Isolated CP2K installation — 23 October 2026

**Installed; startup and library loading verified. Scientific compatibility remains untested.** No registration or paid license was needed. The scientific goal remains paused; no molecular calculation followed the build.

## Installed artifacts

All paths below are relative to the Alquemia repository root.

| Artifact | Location / SHA256 |
|---|---|
| Executable | `workspaces/lady_electronic_decomposition_20261023/cp2k_build_v2/build/bin/cp2k.psmp` |
| Executable SHA256 | `19e4c45c2920299863d0d74543b55ded033ab1e145cbba1cd5cfc3d660e23b32` |
| Runtime image | `workspaces/lady_electronic_decomposition_20261023/cp2k_build_v1/cp2k-2026.2-mpich-generic.sif` |
| Image SHA256 | `90c3a7a7368e8025689a3ea13ba773ecf832a94ea5378c984de1b861c9d8deab` |
| Official OCI image reference | `cp2k/cp2k@sha256:28e35f63a5253b95cc0f50e786370669185ad9f58eacacf33170cb1257216f36` |
| Pristine source archive | `workspaces/lady_electronic_decomposition_20261023/cp2k_build_v1/cp2k-2026.2.tar.bz2` |
| Source archive SHA256 | `f9bd86f580f57a53a0768c0045d1417f9f9a1d66d851ed7f662f496200043373` |
| Applied patch | `cp2k_fix/dy55_default_occupation.patch` |
| Patch SHA256 | `7a4fcb3d977f90b108bb1f1f61b45c78e4a2e4236d18b122df9b9f9e96f265e2` |

The official archive's embedded `REVISION` and the binary both report `git:c92cc08`. The previously recorded GitHub release-tag commit `67b5da876dd6a76b8b021d5a04d1c81ba79a4c50` must not be substituted for this banner revision. The publisher-verified archive digest, complete compiler-input comparison and exact patch identify the actual compiled source. See `cp2k_fix/SOURCE_INTEGRITY.json`: all1,585 checked source/CMake files match the pristine archive except the reviewed seven-line reader repair; no files are added or missing.

The repair supplies Dy55's default occupation bookkeeping. It preserves the ECP coefficients, core count, effective charge and existing ALMO equations. The original archive and failed build attempt remain intact.

## Runtime verification

Job21440 successfully ran `--version` and `ldd`; no required library was missing. Root subsequently repeated `--version` using only the runtime image and executable/source bindings, with no build tools mounted. Exact evidence is in `cp2k_build_v2/{BUILD_RESULT.json,version.txt,linked_libraries.txt,RUNTIME_VERIFICATION.json,ROOT_ACK.json}` under the workspace above.

Reported version: CP2K2026.2, GCC15.2.0. Enabled flags: `omp libint fftw3 libxc parallel scalapack mpi_f08 libdftd4 dftd4_v4_2`. The executable was compiled with a Broadwell CPU target; run it on compatible cluster CPUs. The container uses MPICH. Do not combine it with host OpenMPI or treat this as a standalone host binary.

From the repository root, this repeats the metadata-only version check:

```bash
lady_cp2k_root="$PWD/workspaces/lady_electronic_decomposition_20261023"
apptainer exec --cleanenv \
  --bind "$lady_cp2k_root/cp2k_build_v2:/work:ro" \
  --bind "$lady_cp2k_root/cp2k_build_v1/cp2k-2026.2:/source:ro" \
  --pwd /tmp \
  "$lady_cp2k_root/cp2k_build_v1/cp2k-2026.2-mpich-generic.sif" \
  /work/build/bin/cp2k.psmp --version
```

Molecular execution still needs the declared Slurm allocation, explicit rank/thread counts, a writable task directory, exact translated input and the pinned source data path. This command does not launch the pilot or qualify parallel molecular execution.

## Build recovery and costs

| Job | Outcome | Elapsed seconds | Allocated CPUs | Allocated CPU-seconds |
|---|---|---:|---:|---:|
|21396|Image download succeeded; stopped before compilation because the runtime image omitted CMake/pkg-config|1197|24|28728|
|21440|Same image/source reused; build, version and loader checks passed|700|24|16800|
|Total||| |45528|

Both jobs used `standard_expand`. These are allocated resource totals, not measured CPU utilization. The recovery supplied official portable CMake3.28.3, verified against the publisher checksum, and isolated copies of host pkgconf1.8.1 and its library. Quantum dependencies were reused from the immutable image. Exact build scripts, tool hashes and commands are retained; no system installation, production default or unrelated job was changed.

## Remaining scientific gate

The requested CT differential is still unknown. On resumption of the scientific goal, prepare the exact archived two-fragment quartet and run the La-aquo compatibility pilot: charges/electron counts/effective singlets, full-system energy conventions and decomposition closure. Preserve the original PBE0, basis/ECP assets and geometry hashes. Account explicitly for the archived D4 and fragment reference terms. Release the remaining three systems in parallel only if that pilot passes. Installation success is not authorization to expand the experimental scope.
