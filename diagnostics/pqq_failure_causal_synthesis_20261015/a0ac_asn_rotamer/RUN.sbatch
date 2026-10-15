#!/bin/bash
#SBATCH --job-name=Asn-rotamer-prep
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --no-requeue
set -euo pipefail
ASN_MANIFEST="$1"
ASN_SHA="$2"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1 OPENMM_CPU_THREADS=1
export OMP_THREAD_LIMIT=1 OMP_DYNAMIC=FALSE MKL_DYNAMIC=FALSE
ASN_WORK=$(dirname "$ASN_MANIFEST")
ASN_PYTHON=$(/usr/bin/python3 -c 'import json,sys;print(json.load(open(sys.argv[1]))["python"]["path"])' "$ASN_MANIFEST")
exec "$ASN_PYTHON" "$ASN_WORK/implementation/PREPARE.py" execute --manifest "$ASN_MANIFEST" --manifest-sha "$ASN_SHA"
