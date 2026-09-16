#!/bin/bash
#SBATCH -p standard-shared
#SBATCH -J alquemia-tabi-group
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --hint=nomultithread
#SBATCH --mem=24G
set -euo pipefail
if [ "$#" -ne 1 ]; then
    echo 'Usage: run_surface_group.sbatch ABSOLUTE_EXECUTION_SLICE_MANIFEST' >&2
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
from affordable_common import InvalidArtifact, read_json, verify, write_new
from global_electrostatic import execute_surfaces
path = Path(sys.argv[1]).resolve()
manifest = read_json(path)
verify(manifest['parent_campaign'])
affinity = sorted(os.sched_getaffinity(0))
cores = set()
for cpu in affinity:
    root = Path('/sys/devices/system/cpu') / ('cpu' + str(cpu)) / 'topology'
    cores.add((int((root / 'physical_package_id').read_text()), int((root / 'core_id').read_text())))
required = manifest['execution_slice']['required_distinct_physical_cores']
write_new(path.parent / ('allocation_' + os.environ['SLURM_JOB_ID'] + '.json'), {
    'slurm_job_id': os.environ['SLURM_JOB_ID'], 'affinity': affinity,
    'distinct_physical_cores': sorted(cores), 'required_distinct_physical_cores': required,
    'allocated_cpus': int(os.environ['SLURM_CPUS_ON_NODE']),
    'lscpu': subprocess.check_output(['lscpu'], text=True),
    'cpuinfo_at_start': Path('/proc/cpuinfo').read_text(),
})
if len(cores) < required or len(manifest['tasks']) != required or int(os.environ['SLURM_CPUS_ON_NODE']) != required:
    raise InvalidArtifact('allocation does not provide the declared task/core layout')
print(execute_surfaces(path)['status'])
PY
