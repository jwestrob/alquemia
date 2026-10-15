#!/bin/bash
#SBATCH --job-name=Asn-saved-geometry
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --no-requeue
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1
ASN_CONFIG="$1"
ASN_SHA="$2"
ASN_WORK=$(dirname "$ASN_CONFIG")
ASN_PYTHON=$(/usr/bin/python3 -c 'import json,sys;print(json.load(open(sys.argv[1]))["python"]["path"])' "$ASN_CONFIG")
exec "$ASN_PYTHON" "$ASN_WORK/implementation/ANALYZE.py" "$ASN_CONFIG" "$ASN_SHA"
