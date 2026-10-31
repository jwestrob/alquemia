#!/bin/bash
set -euo pipefail
ROOT=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
OUT="${NIKASHA_WORKSPACE:-$ROOT/workspaces/ln_structure_increment_20261030}/historical"
mkdir -p "$OUT"
cp "$ROOT/workspaces/ln_selectivity_calibration_20261027/pairs.tsv" "$OUT/"
cp "$ROOT/workspaces/ln_selectivity_calibration_20261027/smiles.tsv" "$OUT/"
cd "$OUT"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
/groups/banfield/users/jwestrob/envs/lnsel/bin/python "$ROOT/diagnostics/ln_selectivity_calibration_20261027/fit_donors.py"
