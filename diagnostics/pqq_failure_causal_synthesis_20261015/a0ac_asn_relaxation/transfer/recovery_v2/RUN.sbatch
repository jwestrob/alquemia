#!/bin/bash
#SBATCH --job-name=Asn-transfer-recovery
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --no-requeue
set -euo pipefail
ASN_RECOVERY_MANIFEST="$1"
ASN_RECOVERY_SHA="$2"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1 OPENMM_CPU_THREADS=1 OMP_THREAD_LIMIT=1
export OMP_DYNAMIC=FALSE MKL_DYNAMIC=FALSE
ASN_RECOVERY_WORK=$(dirname "$ASN_RECOVERY_MANIFEST")
ASN_RECOVERY_PYTHON=$(/usr/bin/python3 -c 'import json,sys;print(json.load(open(sys.argv[1]))["python"]["path"])' "$ASN_RECOVERY_MANIFEST")
exec "$ASN_RECOVERY_PYTHON" "$ASN_RECOVERY_WORK/implementation/RECOVER.py" execute --manifest "$ASN_RECOVERY_MANIFEST" --manifest-sha "$ASN_RECOVERY_SHA"
