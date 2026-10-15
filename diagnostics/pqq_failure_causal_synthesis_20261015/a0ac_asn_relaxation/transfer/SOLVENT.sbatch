#!/bin/bash
#SBATCH --job-name=Asn-transfer-solvent
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=7
#SBATCH --no-requeue
set -euo pipefail
ASN_ACTION="$1"
ASN_CONFIG="$2"
ASN_SHA="$3"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1 OPENMM_CPU_THREADS=1 OMP_THREAD_LIMIT=1
export OMP_DYNAMIC=FALSE MKL_DYNAMIC=FALSE
ASN_WORK=$(dirname "$ASN_CONFIG")
ASN_PYTHON=$(/usr/bin/python3 -c 'import json,sys;print(json.load(open(sys.argv[1]))["numerical_python"]["path"])' "$ASN_CONFIG")
exec "$ASN_PYTHON" "$ASN_WORK/SOLVENT.py" "$ASN_ACTION" --config "$ASN_CONFIG" --config-sha "$ASN_SHA"
