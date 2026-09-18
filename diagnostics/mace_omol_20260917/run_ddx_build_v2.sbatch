#!/bin/bash
#SBATCH -p standard,memory
#SBATCH -J alquemia-ddX-build
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=32G
#SBATCH --exclude=node-112-1500g-[1-4],node-224-3t-[1-2]
set -euo pipefail
if [[ $# -ne 1 || $1 != /* ]]; then
    echo 'Usage: run_ddx_build_v2.sbatch ABSOLUTE_BUILD_MANIFEST' >&2
    exit 2
fi
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export CMAKE_BUILD_PARALLEL_LEVEL=64
export FC=/usr/bin/gfortran CXX=/usr/bin/g++
BUILD_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
"$BUILD_PY" - "$1" <<'PY'
import hashlib,json,os,subprocess,sys,time
from pathlib import Path
manifest=Path(sys.argv[1]);m=json.loads(manifest.read_text());root=manifest.parent
for pin in m['pins']:
    if hashlib.sha256(Path(pin['path']).read_bytes()).hexdigest()!=pin['sha256']:raise ValueError('build source pin changed: '+pin['path'])
if int(os.environ['SLURM_CPUS_ON_NODE'])<64:raise ValueError('declared build allocation required')
attempt=root/'attempt_0001';attempt.mkdir(exist_ok=False);start=time.monotonic();commands=[]
def run(args):
    n=len(commands);log=attempt/f'command_{n}.log';resources=attempt/f'command_{n}.resources'
    with log.open('w') as f:r=subprocess.run(['/usr/bin/time','-v','-o',str(resources),*map(str,args)],stdout=f,stderr=subprocess.STDOUT)
    commands.append(dict(command=list(map(str,args)),returncode=r.returncode,log=str(log),resources=str(resources)))
    if r.returncode:raise RuntimeError('build command failed; see '+str(log))
def pin(p):return dict(path=str(p.resolve()),sha256=hashlib.sha256(p.read_bytes()).hexdigest())
status='failed';reason=None
try:
    env=attempt/'venv';run([sys.executable,'-m','venv',env]);py=env/'bin/python'
    run([py,'-m','pip','install','--no-index','--no-deps',*m['dependency_wheels']])
    os.environ['PATH']=str(env/'bin')+os.pathsep+os.environ['PATH']
    run([env/'bin/cmake','--version'])
    wheels=attempt/'wheels';wheels.mkdir()
    run([py,'-m','pip','wheel','--no-index','--no-deps','--no-build-isolation','--wheel-dir',wheels,m['source_archive']])
    built=list(wheels.glob('pyddx*.whl'))
    if len(built)!=1:raise ValueError('expected exactly one ddX wheel')
    run([py,'-m','pip','install','--no-index','--no-deps',built[0]])
    run([py,'-c','import pyddx; print(pyddx.__version__); print(pyddx.banner()); print(pyddx.__file__)'])
    run([py,'-m','pip','freeze','--all'])
    so=list((env/'lib/python3.11/site-packages').glob('pyddx*.so'))
    if len(so)!=1:raise ValueError('expected one native module')
    run(['/usr/bin/ldd',so[0]]);status='built_and_imported'
except Exception as exc:reason=str(exc)
result=dict(status=status,reason=reason,manifest=pin(manifest),slurm_job_id=os.environ['SLURM_JOB_ID'],commands=commands,
            wall_seconds=time.monotonic()-start,scientific_calculations=0,
            output_artifacts=[pin(p) for p in attempt.rglob('*') if p.is_file() and (p.suffix in ('.log','.resources','.whl','.so'))])
(attempt/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ('commands','output_artifacts')}))
if status!='built_and_imported':sys.exit(1)
PY
