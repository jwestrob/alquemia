#!/bin/bash
set -euo pipefail
cd "$1"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
ORCADIR=/groups/banfield/users/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg
export LD_LIBRARY_PATH="$ORCADIR:${LD_LIBRARY_PATH:-}"
python - <<'PY'
import os,json,time
from pathlib import Path
assert not Path('job.out').exists() and not Path('job.gbw').exists()
n=int(os.environ['SLURM_NTASKS'])
m=int(os.environ.get('SLURM_MEM_PER_NODE','0')) or int(os.environ['SLURM_MEM_PER_CPU'])*n
maxcore=int(.75*m/n)
Path('job.inp').write_text(Path('prepared.inp').read_text().replace('@NPROCS@',str(n)).replace('@MAXCORE@',str(maxcore)))
Path('RUNTIME.json').write_text(json.dumps(dict(job=os.environ['SLURM_JOB_ID'],nprocs=n,allocated_MB=m,maxcore_MB=maxcore,started=time.time()),indent=2)+'\n')
PY
scontrol show job "$SLURM_JOB_ID" > SCONTROL.txt
set +e
"$ORCADIR/orca" job.inp > job.out 2> job.err < /dev/null
rc=$?
set -e
python - "$rc" <<'PY'
import sys,time,json
from pathlib import Path
Path('EXIT.json').write_text(json.dumps(dict(returncode=int(sys.argv[1]),finished=time.time()),indent=2)+'\n')
PY
exit "$rc"
