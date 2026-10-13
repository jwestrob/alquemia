#!/bin/bash
#SBATCH --job-name=PedE-r2scan-reference-v2
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --no-requeue
set -euo pipefail
if [ "$#" -lt 8 ]; then
  echo 'Usage: RUN_NATIVE.sbatch ABS_DRIVER ABS_MANIFEST MANIFEST_SHA ABS_REGRESSION_RECEIPT REGRESSION_SHA DRIVER_SHA ABS_OUTPUT TASK_ID [...]' >&2
  exit 2
fi
REFERENCE_DRIVER="$1"
shift
export PATH=/groups/banfield/users/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg:/usr/bin:/bin:$PATH
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export OMP_THREAD_LIMIT=1 OMP_DYNAMIC=FALSE MKL_DYNAMIC=FALSE METAL_ENV_WORKERS=1
export OMPI_MCA_hwloc_base_binding_policy=none
export LC_ALL=C LANG=C
exec /usr/bin/python3 "$REFERENCE_DRIVER" "$@"
