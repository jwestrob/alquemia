#!/bin/bash
set -euo pipefail
# Executed only inside the pinned official image, under an owned Slurm allocation.
CP2K_DEPS=/opt/cp2k/spack/spack/opt/spack/view
test -d "$CP2K_DEPS"
export PATH="$CP2K_DEPS/bin:/usr/bin:/bin"
export CMAKE_PREFIX_PATH="$CP2K_DEPS"
export PKG_CONFIG_PATH="$CP2K_DEPS/lib/pkgconfig:$CP2K_DEPS/lib64/pkgconfig:$CP2K_DEPS/share/pkgconfig"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_DYNAMIC=FALSE
for program in cmake ninja mpicc mpifort gcc-15 g++-15 gfortran-15; do
    command -v "$program"
    readlink -f "$(command -v "$program")"
done
gcc-15 --version
gfortran-15 --version
mpifort -show
cmake --version
cmake -S /work/cp2k-2026.2 -B /work/build -G Ninja \
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
    -DCP2K_USE_MPI=ON -DCP2K_USE_MPI_F08=ON \
    -DCP2K_DBCSR_USE_CPU_ONLY=ON \
    -DCP2K_BLAS_VENDOR=OpenBLAS \
    -DCP2K_USE_FFTW3=ON -DCP2K_USE_LIBINT2=ON \
    -DCP2K_USE_LIBXC=ON -DCP2K_USE_DFTD4=ON
cmake --build /work/build --target cp2k-bin --parallel 24
test -x /work/build/bin/cp2k.psmp
/work/build/bin/cp2k.psmp --version > /work/version.txt
ldd /work/build/bin/cp2k.psmp > /work/linked_libraries.txt
if grep -q 'not found' /work/linked_libraries.txt; then
    exit 1
fi
sha256sum /work/build/bin/cp2k.psmp > /work/binary.sha256
