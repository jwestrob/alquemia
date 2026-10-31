#!/bin/bash
set -euo pipefail
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
workspaces/ln_structure_increment_20261030/chemprop_env/bin/python diagnostics/ln_structure_increment_20261030/enrich_overlap.py
