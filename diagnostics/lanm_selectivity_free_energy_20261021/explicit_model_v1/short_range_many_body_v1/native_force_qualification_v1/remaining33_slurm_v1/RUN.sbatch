#!/bin/bash
#SBATCH --job-name=srmb-force33
#SBATCH --account=banfield
#SBATCH --partition=standard
#SBATCH --qos=standard_expand
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --no-requeue
set -euo pipefail
export OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
unset LD_PRELOAD
# Root submits using this source path; Slurm copies the script, so use pinned argument3.
source_dir="$3"
/home/jwestrob/.pyenv/versions/3.12.0/bin/python "$source_dir/RUN.py" --config "$1" --out "$2"
/home/jwestrob/.pyenv/versions/3.12.0/bin/python "$source_dir/COLLECT.py" --recovery "$2/RESULT.json" --out "$2/collection"
