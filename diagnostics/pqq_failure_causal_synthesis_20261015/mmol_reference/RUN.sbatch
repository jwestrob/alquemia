#!/bin/bash
#SBATCH --job-name=MMOL-common3-reference
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --no-requeue
set -euo pipefail
if [ "$#" -lt 3 ]; then
  echo 'Usage: RUN.sbatch ABS_MANIFEST MANIFEST_SHA preflight|execute|collect [TASK_ID ...]' >&2
  exit 2
fi
REFERENCE_MANIFEST="$1"
REFERENCE_SHA="$2"
REFERENCE_ACTION="$3"
shift 3
REFERENCE_WORK=$(dirname "$REFERENCE_MANIFEST")
export PATH=/groups/banfield/users/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg:/usr/bin:/bin:$PATH
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export OMP_THREAD_LIMIT=1 OMP_DYNAMIC=FALSE MKL_DYNAMIC=FALSE METAL_ENV_WORKERS=1
export OMPI_MCA_hwloc_base_binding_policy=none
export LC_ALL=C LANG=C
cd "$REFERENCE_WORK"
exec /usr/bin/python3 "$REFERENCE_WORK/implementation/REFERENCE.py" "$REFERENCE_ACTION" "$REFERENCE_MANIFEST" --manifest-sha "$REFERENCE_SHA" "$@"
