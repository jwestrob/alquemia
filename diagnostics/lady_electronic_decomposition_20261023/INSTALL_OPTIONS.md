# Isolated installation options on biotite

Read-only assessment, 23 October 2026. No software was installed, built, launched, or downloaded as an archive/binary; no environment was changed. Only local tool/version queries and public official documentation/source were read.

## Decision

**Do not install stock CP2K for this quartet now.** Installation is practical, but CP2K 2026.2 rejects the exact archived Dy 55-core ECP. Its ECP reader accepts core counts 0, 2, 10, 18, 28, 36, 46, 54, 60, 68 and 78; 55 reaches `Unknown Core State`. Both inline and file inputs converge on that switch, before the later manual electron-configuration override. La's 46-core case is supported. This is an occupation-initialization/parser restriction, not evidence that the semilocal integral library lacks the required angular channels. See [the exact released source](https://github.com/cp2k/cp2k/blob/67b5da876dd6a76b8b021d5a04d1c81ba79a4c50/src/atom_types.F#L3033) and the separate local `CP2K_COMPATIBILITY.md`/`SOFTWARE_ACCESS.json` assessment.

Do not replace the archived Dy ECP55 with ECP54/60 or GTH, alter electron counts, or claim an ordinary build fixes this. A separately reviewed source extension would be a different software route. This assessment makes no such change. The original quartet and basis/ECP assets remain in `inputs/PROVENANCE.json`, with audit `INPUT_AUDIT.json`.

## Observed local tools

| Component | Observed on current PATH/shell |
|---|---|
| Host | biotite, Linux x86_64 |
| GNU C/C++/Fortran | `/usr/bin/{gcc,g++,gfortran}`, GCC/GFortran 13.3.0 |
| MPI | `/usr/bin/{mpicc,mpifort,mpirun}`, OpenMPI 4.1.6; wrappers select gcc/gfortran |
| Build | `/usr/bin/cmake` 3.28.3; make available; ninja resolves to a pyenv shim, its underlying executable not verified |
| Environment managers | conda 24.7.1 and mamba 1.5.8 under `/shared/software/miniconda3/latest/`; no active conda prefix observed |
| Containers | `/usr/bin/apptainer` 1.5.4; singularity alias also present |
| Confirmed system packages | FFTW3 development 3.3.10, OpenMPI development 4.1.6, LAPACK runtime 3.12.0 |
| Not found on PATH | CP2K executables, micromamba, spack, Intel compiler executables; no module interface in this shell |

The bounded package query did not establish installed DBCSR, Libint2, LibXC, OpenBLAS or ScaLAPACK development stacks. Absence from these queries is not an exhaustive filesystem absence claim.

## Routes if the scientific compatibility blocker is resolved separately

1. **Fastest availability route: official versioned single-node binary.** The [CP2K v2026.2 release](https://github.com/cp2k/cp2k/releases/tag/v2026.2) offers `cp2k-2026.2-Linux-gnu-x86_64.ssmp` (428,713,088 bytes) and `cp2k-2026.2.tar.bz2` (79,898,041 bytes), confirmed through the official GitHub release API. The executable uses the single-node/OpenMP route; it is not an MPI executable. Keep it in an isolated project directory and pin its downloaded hash. Actual loader/CPU compatibility and compiled feature reporting have not been tested. A corresponding data directory or explicitly supplied input assets may also be needed. No binary was requested here.

2. **Controlled source build: versioned 2026.2 release tarball, isolated prefix.** Existing GNU13 and CMake3.28 exceed the release's CMake3.24 requirement. CP2K's [build guide](https://manual.cp2k.org/cp2k-2026_2-branch/getting-started/build-from-source.html) supports a managed toolchain or `make_cp2k.sh`/Spack route. A minimal CPU stack needs DBCSR, BLAS/LAPACK, and suitable hybrid-functional support; MPI adds matching MPI/ScaLAPACK. Use the official release tarball, not an unpinned development checkout. This is feasible but slower than checking the binary because dependencies must be resolved/built.

   For this intended PBE0 calculation, retain Libint2/Hartree–Fock-exchange support and review the XC and D4 choices explicitly. The release exposes `CP2K_USE_LIBINT2`, `CP2K_USE_LIBXC`, and `CP2K_USE_DFTD4`; the Libint2 Fortran interface must exist. **There is no separate LIBGRPP installation/enable switch needed in this release:** [v2026.2 `src/CMakeLists.txt`](https://github.com/cp2k/cp2k/blob/v2026.2/src/CMakeLists.txt#L1461) unconditionally appends the bundled `grpp` C sources and Fortran interface, while [`libgrpp_integrals.F`](https://github.com/cp2k/cp2k/blob/v2026.2/src/libgrpp_integrals.F) implements the semilocal integrals. The [KIND manual](https://manual.cp2k.org/cp2k-2026_1-branch/CP2K_INPUT/FORCE_EVAL/SUBSYS/KIND.html#CP2K_INPUT.FORCE_EVAL.SUBSYS.KIND.ECP_SEMI_LOCAL) distinguishes original semilocal ECPs from the optional nonlocal one-center expansion. Bundled LIBGRPP does not remove the core-count parser restriction above.

3. **Dedicated conda environment or Apptainer image: practical alternatives, not fixes for core55.** CP2K's [distribution guide](https://github.com/cp2k/cp2k/blob/v2026.2/docs/getting-started/distributions.md) documents both conda-forge and official containers. The [official container repository](https://github.com/cp2k/cp2k-containers) documents Apptainer use and notes that 2026.2 container recipes now reside in CP2K itself. Existing local tools can support these routes, but exact package version/build features, image digest, and MPI provenance would need inspection. Do not mix a conda/container MPI library with the system launcher without verifying compatibility. No package solve, image pull, or container run was attempted. Spack is not currently on PATH; introducing it solely for this four-cell task adds setup work.

## Qbics public availability check

The [official download page](https://www.qbics.info/home/qbics-download.html) advertises lean source (~200 MB), source with build environment (~600 MB), and version2.0 executables for Linux, macOS and Windows (release date 1 April 2026). Downloads are presented behind a required registration form collecting email, name, title, institution and research field. **No form was submitted, contact made, archive requested, or download restriction bypassed.**

The public download, installation and feature pages inspected do not state an explicit software license. Therefore source availability is advertised, but license/use/redistribution terms remain unverified; do not call the main Qbics engine open-source based on that alone. The separately advertised Qbics-MolStar toolkit's open-source description is not an engine license.

The [installation manual](https://www.qbics.info/doc/begin/installation.html) lists Linux OpenMP and MPI CPU builds and requires a matching MPI implementation for the latter. The [compilation manual](https://qbics.info/doc/begin/compilation.html) gives compiler >=8 and CMake >=3.16 minima, which the observed local toolchain meets. The [feature page](https://www.qbics.info/home/qbics-features.html) lists hybrid DFT, user-defined pseudopotentials and GKS/TSO EDA. Those declarations are promising availability evidence, **not** proof of the exact requested ALMO constrained state, Dy55 input support, or agreement with the archived ORCA PBE0 energies. Method equivalence and licensing must be resolved before this becomes a usable replacement route.

## Scope and practical conclusion

Local installation infrastructure is adequate. The present CP2K stop is exact-Hamiltonian input compatibility, while Qbics has a required registration step and unresolved public license terms. No installation or scientific calculation is warranted merely to demonstrate that an executable exists. Root retains installation and execution ownership; all archived original inputs and result pointers remain unchanged.
