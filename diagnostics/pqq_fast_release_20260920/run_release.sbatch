#!/bin/bash
#SBATCH -p gpu
#SBATCH -J PQQ-source-release
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --mem=200000M
#SBATCH --nodelist=node-224-2t-8gpu-1
set -euo pipefail
PQQ_RELEASE_RUN=$1
case "$PQQ_RELEASE_RUN" in /*) ;; *) exit 2;; esac
PQQ_RELEASE_ROOT=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
PQQ_RELEASE_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
PQQ_RELEASE_GPU_PY="$PQQ_RELEASE_ROOT/workspaces/mace_hybrid_20260916/software_v1/venv/bin/python"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
"$PQQ_RELEASE_PY" "$PQQ_RELEASE_RUN/implementation/pqq_fast_prepare.py" prepare \
 --manifest "$PQQ_RELEASE_ROOT/diagnostics/pqq_fast_release_20260920/SOURCES.json" \
 --output "$PQQ_RELEASE_RUN/source_preparation"
"$PQQ_RELEASE_PY" "$PQQ_RELEASE_RUN/implementation/pqq_fast_release.py" prepare-score \
 --preparation "$PQQ_RELEASE_RUN/source_preparation/preparation.json" \
 --comparison "$PQQ_RELEASE_ROOT/workspaces/compact_solvation_20260920/full_v1/comparison_v1.json" \
 --agreement "$PQQ_RELEASE_ROOT/diagnostics/pqq_fast_release_20260920/PLAN.md" \
 --output "$PQQ_RELEASE_RUN/prepared_score" \
 --cpu-python "$PQQ_RELEASE_PY" --gpu-python "$PQQ_RELEASE_GPU_PY"
PQQ_RELEASE_SCORE="$PQQ_RELEASE_RUN/prepared_score/scoring"
"$PQQ_RELEASE_PY" "$PQQ_RELEASE_SCORE/implementation/compact_solvation_scanner.py" execute \
 --manifest "$PQQ_RELEASE_SCORE/manifest.json"
"$PQQ_RELEASE_PY" "$PQQ_RELEASE_RUN/implementation/pqq_fast_release.py" collect \
 --preparation "$PQQ_RELEASE_RUN/source_preparation/preparation.json" \
 --result "$PQQ_RELEASE_SCORE/result_${SLURM_JOB_ID}.json" \
 --output "$PQQ_RELEASE_RUN/RESULT_${SLURM_JOB_ID}.json"
