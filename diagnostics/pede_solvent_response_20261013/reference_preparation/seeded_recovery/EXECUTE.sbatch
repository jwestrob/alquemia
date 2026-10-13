#!/bin/bash
#SBATCH --job-name=PedE-Ca-seeded-vacuum
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --no-requeue
set -euo pipefail
if [ "$#" -ne 5 ]; then echo 'Usage: EXECUTE.sbatch ABS_MANIFEST SHA ABS_OUTPUT ABS_TEST_RECEIPT TEST_SHA' >&2;exit 2;fi
export PATH=/groups/banfield/users/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg:/usr/bin:/bin:$PATH
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export OMP_THREAD_LIMIT=1 OMP_DYNAMIC=FALSE MKL_DYNAMIC=FALSE METAL_ENV_WORKERS=1
export OMPI_MCA_hwloc_base_binding_policy=none LC_ALL=C LANG=C
exec /usr/bin/python3 "$(dirname "$1")/implementation/SEEDED.py" execute "$1" "$2" "$3" --test-receipt "$4" --test-sha "$5"
