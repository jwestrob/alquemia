#!/usr/bin/env python3
"""One immutable 1-MPI x 8-OMP CP2K benchmark attempt; root submits Slurm."""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time
import traceback

HERE = Path(__file__).resolve().parent


def pin(path):
    p = Path(path).resolve()
    h = hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024), b''):
            h.update(b)
    return {'path':str(p), 'sha256':h.hexdigest()}


def save(path, obj):
    tmp = path.with_suffix('.tmp')
    tmp.write_text(json.dumps(obj,indent=2)+'\n')
    tmp.replace(path)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--cutoff',type=int,choices=[400,600],required=True)
    ap.add_argument('--task',type=Path,required=True)
    ap.add_argument('--seconds',type=int,required=True)
    a = ap.parse_args()
    if not os.environ.get('SLURM_JOB_ID') or os.environ.get('SLURM_JOB_QOS') != 'standard_expand':
        raise RuntimeError('Require actual standard_expand Slurm allocation')
    if int(os.environ.get('SLURM_CPUS_PER_TASK','0')) != 8:
        raise RuntimeError('Require exactly 8 allocated CPUs for 1 MPI x 8 OMP')
    if not 1 <= a.seconds <= 3600:
        raise RuntimeError('Root must declare a finite cap of at most 3600 seconds')
    task = a.task.resolve()
    task.mkdir(parents=True, exist_ok=False)
    r = {'status':'starting', 'implementation':pin(__file__), 'cutoff_Ry':a.cutoff,
         'job_id':os.environ['SLURM_JOB_ID'], 'qos':os.environ['SLURM_JOB_QOS'],
         'allocation':{k:v for k,v in os.environ.items() if k in
                       ['SLURM_CPUS_PER_TASK','SLURM_MEM_PER_NODE','SLURM_MEM_PER_CPU','SLURM_JOB_NODELIST']},
         'native_calls':0, 'mpi_ranks':1, 'omp_threads':8, 'error':None,
         'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
    receipt = task/'RUN_RESULT.json'
    start = time.monotonic()
    try:
        manifest = json.loads((HERE/'cp2k/MANIFEST.json').read_text())
        runtime = json.loads((HERE/'cp2k/RUNTIME.json').read_text())
        r['manifest'] = pin(HERE/'cp2k/MANIFEST.json')
        r['runtime'] = pin(HERE/'cp2k/RUNTIME.json')
        for label in ['binary','image']:
            if pin(runtime[label]['path']) != runtime[label]:
                raise RuntimeError(label+' hash differs')
        for name, dest in [(f'water_pbe_{a.cutoff}.inp','input.inp'),
                           ('BASIS_BENCHMARK','BASIS_BENCHMARK'),
                           ('POTENTIAL_BENCHMARK','POTENTIAL_BENCHMARK')]:
            source = HERE/'cp2k'/name
            if pin(source)['sha256'] != manifest['inputs'][name]['sha256']:
                raise RuntimeError('Changed benchmark asset '+name)
            shutil.copyfile(source,task/dest)
        build = Path(runtime['build']).resolve()
        source = Path(runtime['source']).resolve()
        argv = ['/usr/bin/apptainer','exec','--cleanenv',
                '--bind',str(build)+':/work:ro','--bind',str(source)+':/source:ro',
                '--bind',str(task)+':/task','--pwd','/task',
                '--env','OMP_NUM_THREADS=8','--env','OPENBLAS_NUM_THREADS=1',
                '--env','MKL_NUM_THREADS=1','--env','OMP_STACKSIZE=128M',
                '--env','CP2K_DATA_DIR=/source/data',runtime['image']['path'],
                '/opt/cp2k/spack/spack/opt/spack/view/bin/mpiexec',
                '-launcher','fork','-rmk','user','-bind-to','none','-n','1',
                '/work/build/bin/cp2k.psmp','-i','/task/input.inp','-o','/task/endpoint.out']
        r.update(status='running', argv=argv, native_calls=1, cap_seconds=a.seconds,
                 input=pin(task/'input.inp'))
        save(receipt,r)
        t = time.monotonic()
        with (task/'native_launcher.log').open('x') as f:
            proc = subprocess.Popen(argv,stdout=f,stderr=subprocess.STDOUT,start_new_session=True)
            r['native_pid'] = proc.pid
            save(receipt,r)
            try:
                r['returncode'] = proc.wait(timeout=a.seconds)
            except subprocess.TimeoutExpired:
                r['timed_out'] = True
                os.killpg(proc.pid,signal.SIGTERM)
                try:
                    proc.wait(timeout=15)
                except subprocess.TimeoutExpired:
                    os.killpg(proc.pid,signal.SIGKILL)
                    proc.wait()
                r['returncode'] = proc.returncode
                raise RuntimeError('Declared native cap reached; no retry')
            finally:
                r['native_elapsed_seconds'] = time.monotonic()-t
        r['status'] = 'native_complete_needs_comparison' if r['returncode'] == 0 else 'native_failed'
    except Exception:
        r['status'] = 'failed'
        r['error'] = traceback.format_exc()
    finally:
        r['elapsed_seconds'] = time.monotonic()-start
        r['artifacts'] = [pin(p) for p in task.iterdir() if p.is_file() and
                          p.name not in ['RUN_RESULT.json','RUN_RESULT.tmp']]
        save(receipt,r)
    return 0 if r['status'] == 'native_complete_needs_comparison' else 1


if __name__ == '__main__':
    sys.exit(main())
