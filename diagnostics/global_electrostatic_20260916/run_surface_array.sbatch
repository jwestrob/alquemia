#!/bin/bash
#SBATCH -p standard-shared
#SBATCH -J alquemia-tabi-recovery
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --mem=2G
set -euo pipefail
if [ "$#" -ne 1 ]; then
    echo 'Usage: run_surface_array.sbatch ABSOLUTE_EXECUTION_ARRAY_MANIFEST' >&2
    exit 2
fi
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python - "$1" <<'PY'
import os
from pathlib import Path
import subprocess
import sys
sys.path.insert(0, str(Path.cwd() / 'scripts'))
from affordable_common import read_json, verify, write_new
from global_electrostatic import execute_surfaces
manifest = read_json(sys.argv[1])
verify(manifest['parent_campaign'])
index = int(os.environ['SLURM_ARRAY_TASK_ID'])
task = manifest['slices'][index]
assert task['array_index'] == index
path = verify(task['manifest'])
write_new(path.parent / ('allocation_' + os.environ['SLURM_JOB_ID'] + '.json'), {
    'array_job_id': os.environ['SLURM_ARRAY_JOB_ID'], 'array_task_id': index,
    'slurm_job_id': os.environ['SLURM_JOB_ID'],
    'affinity': sorted(os.sched_getaffinity(0)),
    'lscpu': subprocess.check_output(['lscpu'], text=True),
    'cpuinfo_at_start': Path('/proc/cpuinfo').read_text(),
    'scheduler_hint': 'nomultithread; actual topology retained, no frequency or affinity mutation',
})
print(execute_surfaces(path)['status'])
PY
