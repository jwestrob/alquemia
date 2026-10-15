#!/bin/bash
#SBATCH --job-name=Asn-basin-prep
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --no-requeue
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 OMP_THREAD_LIMIT=1 OPENMM_CPU_THREADS=1
export OMP_DYNAMIC=FALSE MKL_DYNAMIC=FALSE
ASN_CONFIG="$1"
ASN_SHA="$2"
ASN_WORK=$(dirname "$ASN_CONFIG")
ASN_PYTHON=$(/usr/bin/python3 -c 'import json,sys;print(json.load(open(sys.argv[1]))["python"]["path"])' "$ASN_CONFIG")
exec "$ASN_PYTHON" "$ASN_WORK/implementation/PREPARE.py" execute --config "$ASN_CONFIG" --config-sha "$ASN_SHA"
