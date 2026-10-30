# Fast isolated build route for the reviewed Dy55 repair

23 October 2026. Recommendation only: no image layers pulled, source archives downloaded, installation, compilation, submission or molecular execution performed by this reviewer. Public registry manifest/config JSON and official source recipes were read.

## Recommended route

**Use the official CP2K2026.2 generic MPICH container as an immutable compiler/dependency environment; rebuild only the pinned release source plus the reviewed reader patch into a writable project directory.** This avoids compiling Libint, DBCSR and the entire dependency stack. Do not rebuild the Dockerfile: its commands also build dependencies and run broad regressions, outside this bounded task.

Exact official image confirmed through DockerHub and its OCI registry:

- Tag: `cp2k/cp2k:2026.2_mpich_generic_psmp`.
- Immutable reference: `docker://cp2k/cp2k@sha256:28e35f63a5253b95cc0f50e786370669185ad9f58eacacf33170cb1257216f36`.
- Architecture: Linux/amd64; compressed layers 1,432,290,667 bytes.
- Image configuration digest: `sha256:f3d9dc68a19277908fa0656c79c7c3f80179b1e909116ae19266d4cd8d5adb0d`.
- Image history installs GCC/G++/GFortran15 and copies the complete dependency install tree to `/opt/cp2k/spack/spack/opt/spack`, plus the released binary/data under `/opt/cp2k/install`. Its default entrypoint is `/opt/cp2k/install/bin/launch`.

The official [tag listing](https://hub.docker.com/r/cp2k/cp2k/tags) also has an OpenMPI generic image (`sha256:23ad7dc21889b6d0c62579eb4839e43d765158abb4378934cfb37a84d40d57fd`, 1,609,558,705 bytes), but MPICH is smaller and adequate for the intended single-node execution. Do not substitute `latest`, a GPU image or a CPU-specific instruction-set image.

## Why this image is suitable for rebuilding

The released [Spack Dockerfile pattern](https://github.com/cp2k/cp2k/blob/v2026.2/tools/docker/Dockerfile.test_spack_psmp-gcc13) preserves the dependency install tree in the final image. The actual image history confirms this and identifies GCC15 rather than the example Dockerfile's GCC13. The release [dependency specification](https://github.com/cp2k/cp2k/blob/v2026.2/tools/spack/cp2k_deps_p.yaml) explicitly includes CMake, Ninja, GNU make, MPICH, DBCSR, Libint with its Fortran interface, LibXC, FFTW, OpenBLAS and ScaLAPACK. It also includes tblite's D4 dependency. This is stronger evidence than assuming a runtime image contains headers.

The expected merged dependency view is `/opt/cp2k/spack/spack/opt/spack/view`: the release [`make_cp2k.sh`](https://github.com/cp2k/cp2k/blob/v2026.2/make_cp2k.sh) rewrites the YAML's install/view roots under that prefix. Root should verify this view's actual symlinks, compiler versions, CMake/Ninja executables and package-config/module files once after pulling the pinned image. I have not inspected its filesystem, so the recipe below is conditional on that cheap check. Preserve the concrete resolved paths. If a view link is absent, resolve the corresponding existing package prefix inside the same tree; do not silently build or install a replacement dependency.

Use the image's compiler and MPI, not host GCC13 or host OpenMPI4.1.6. This preserves the Fortran `.mod` ABI and MPI library consistency. Run the resulting executable inside the same pinned image.

## Concrete build recipe for root to bind

1. Materialize the pinned image and versioned `cp2k-2026.2.tar.bz2` release asset in the isolated build workspace; record both actual hashes. Preserve a pristine source copy and apply only the reviewed Dy55 reader delta to a separate copy. Confirm release commit `67b5da876dd6a76b8b021d5a04d1c81ba79a4c50` provenance and exact inverse diff.
2. Bind that workspace to `/work` and use **`apptainer exec --cleanenv`**, not the image's default run/launch entrypoint. The latter is for the installed stock executable and must not accidentally select it.
3. With the verified dependency view available, configure out of tree. The following is an unexecuted template; source/build paths and resolved tool paths must be bound by root:

```bash
CP2K_DEPS=/opt/cp2k/spack/spack/opt/spack/view
export PATH="$CP2K_DEPS/bin:/usr/bin:/bin"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export OMP_DYNAMIC=FALSE

cmake -S /work/source-patched -B /work/build -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER=/usr/bin/gcc-15 \
  -DCMAKE_CXX_COMPILER=/usr/bin/g++-15 \
  -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran-15 \
  -DCMAKE_PREFIX_PATH="$CP2K_DEPS" \
  -DMPI_C_COMPILER="$CP2K_DEPS/bin/mpicc" \
  -DMPI_Fortran_COMPILER="$CP2K_DEPS/bin/mpifort" \
  -DBUILD_SHARED_LIBS=OFF \
  -DCP2K_USE_EVERYTHING=OFF \
  -DCP2K_USE_ACCEL=NONE \
  -DCP2K_USE_MPI=ON \
  -DCP2K_USE_MPI_F08=ON \
  -DCP2K_DBCSR_USE_CPU_ONLY=ON \
  -DCP2K_BLAS_VENDOR=OpenBLAS \
  -DCP2K_USE_FFTW3=ON \
  -DCP2K_USE_LIBINT2=ON \
  -DCP2K_USE_LIBXC=ON \
  -DCP2K_USE_DFTD4=ON

cmake --build /work/build --target cp2k-bin --parallel 24
```

`cp2k-bin` is the actual CMake executable target, producing `cp2k.psmp`; `cp2k` is the library target. Build only the requested executable rather than every test/tool. `BUILD_SHARED_LIBS=OFF` makes the CP2K library internal to the executable; external dependency libraries may still be shared. Do not claim a fully static binary. Keep default build-tree runtime paths and execute from the retained build tree inside the container; no installation into `/opt` is needed. Use the source's data directory explicitly for subsequent inputs.

The selected features support CPU PBE0/hybrid exchange and explicit D4 availability without unrelated optional solvers, machine-learning packages or GPU backends. LIBGRPP is bundled and compiled unconditionally in this release. No extra LIBGRPP switch exists. The patch must not change any library, integral or ECP coefficient source. CMake must report all requested dependencies found; missing features are a setup failure, not permission to turn them off or change the Hamiltonian.

## Finite resource and completion plan

One root-submitted build job: `standard` partition, `standard_expand` QoS, `banfield` account, one task and24 CPUs, no GPU, scheduler default memory/time, no requeue. Limit compiler parallelism to24, with library/OpenMP threads1. A prospective two-hour build-process timeout is a conservative finite stop, not an estimated runtime or an instruction to change scheduler defaults. It bounds the attempt at172,800 allocated CPU-seconds; record actual parent accounting and peak RSS. Root may choose a different explicit finite cap before submission.

Save image/source/patch hashes, exact command/environment, complete configure/build logs, `CMakeCache.txt`, dependency paths, binary hash and actual cost. Only version/help/loader inspection follows compilation; no `run_tests`, molecular input, ALMO pilot or broad regression suite is implicit in this build. The archived quartet's scientific verification remains separately controlled.

## If the image lacks a usable development view

Stop this build attempt before compiling replacement dependencies. The documented fallback is the official release's managed toolchain or `make_cp2k.sh` in a new isolated prefix, using a minimal CPU feature set and one bounded package build at a time. That is slower and not the preferred first route. Do not invoke `make_cp2k.sh` inside the production image blindly: it expects a Spack installation/environment beyond the copied install tree and may attempt to recreate the dependency stack. See the official [source-build guide](https://manual.cp2k.org/cp2k-2026_2-branch/getting-started/build-from-source.html) and [container usage](https://github.com/cp2k/cp2k-containers).
