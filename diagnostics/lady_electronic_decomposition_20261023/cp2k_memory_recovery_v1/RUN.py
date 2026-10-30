#!/usr/bin/env python3
"""One resource-recovery native ALMO call; same inputs, no internal retry."""
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


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def save(path, data):
    tmp = path.with_suffix('.tmp')
    tmp.write_text(json.dumps(data, indent=2)+'\n')
    tmp.replace(path)


def bounded(argv, output, seconds, record):
    record.update(argv=argv, output=str(output), cap_seconds=seconds, started_utc=now())
    started = time.monotonic()
    with output.open('x') as log:
        proc = subprocess.Popen(argv, stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
        record['pid'] = proc.pid
        try:
            record['returncode'] = proc.wait(timeout=seconds)
        except subprocess.TimeoutExpired:
            record['timed_out'] = True
            # Only the child process group created for this declared bounded call.
            os.killpg(proc.pid, signal.SIGTERM)
            try:
                proc.wait(timeout=15)
            except subprocess.TimeoutExpired:
                os.killpg(proc.pid, signal.SIGKILL)
                proc.wait()
            record['returncode'] = proc.returncode
            raise RuntimeError(f'Predeclared process cap reached: {seconds}s')
        finally:
            record['elapsed_seconds'] = time.monotonic()-started
            record['finished_utc'] = now()
    return record['returncode']


def main():
    config_path = Path(sys.argv[1]).resolve()
    task = Path(sys.argv[2]).resolve()
    gid = sys.argv[3]
    config = json.loads(config_path.read_text())
    assert gid in config['systems']
    assert task == Path(config['systems'][gid]['task']).resolve()
    dest = Path(config['diagnostic_results'])/gid
    if dest.exists():
        raise RuntimeError('Existing diagnostic output preserved; refusing overwrite')
    task.mkdir(parents=True, exist_ok=True)
    if (task/'RUN_RESULT.json').exists() or (task/'execution').exists():
        raise RuntimeError('Existing attempt preserved; refusing duplicate or restart')
    execution = task/'execution'
    execution.mkdir()
    result = {'status': 'starting', 'geometry_id': gid, 'started_utc': now(),
              'job_id': os.environ['SLURM_JOB_ID'], 'commands': [],
              'molecular_parent_calls': 0, 'config_sha256': sha(config_path),
              'scope': 'One native molecular guess, block polarization, and full SCF chain.',
              'allocation': {k: v for k,v in os.environ.items() if k in
                    ['SLURM_JOB_ID','SLURM_JOB_NODELIST','SLURM_CPUS_PER_TASK',
                     'SLURM_MEM_PER_CPU','SLURM_MEM_PER_NODE','SLURM_JOB_CPUS_PER_NODE']}}
    receipt = task/'RUN_RESULT.json'
    started = time.monotonic()
    exit_code = 1
    try:
        assert gid in config['systems']
        assert int(os.environ['SLURM_CPUS_PER_TASK']) == config['cpus'] == 24
        assert config['mpi_ranks'] == 12 and config['omp_threads'] == 2
        assert config['native_cap_seconds'] == 7073
        assert int(os.environ.get('SLURM_MEM_PER_NODE', '0')) == 262144
        for pin in config['pins']:
            if sha(pin['path']) != pin['sha256']:
                raise RuntimeError(f"Pin changed: {pin['path']}")
        if gid != config['pilot']:
            admission = Path(config['pilot_admission'])
            a = json.loads(admission.read_text())
            assert a['admitted'] is True and a['config_sha256'] == sha(config_path)
            result['pilot_admission'] = {'path': str(admission), 'sha256': sha(admission)}
        sysrec = config['systems'][gid]
        source_deck = Path(sysrec['deck'])
        assert sha(source_deck) == sysrec['deck_sha256']
        shutil.copyfile(source_deck, execution/'input.inp')
        runtime = config['runtime']
        base = ['/usr/bin/apptainer','exec','--cleanenv',
                '--bind',runtime['build']+':/work:ro',
                '--bind',runtime['source']+':/source:ro',
                '--bind',config['assets']+':/assets:ro',
                '--bind',str(execution)+':/task', '--pwd','/task',
                '--env',f"OMP_NUM_THREADS={config['omp_threads']}",
                '--env','OPENBLAS_NUM_THREADS=1', '--env','MKL_NUM_THREADS=1',
                '--env','OMP_STACKSIZE=128M','--env','CP2K_DATA_DIR=/source/data',
                runtime['image']]
        binary = '/work/build/bin/cp2k.psmp'
        mpi = ['/opt/cp2k/spack/spack/opt/spack/view/bin/mpiexec',
               '-launcher','fork','-rmk','user','-bind-to','none',
               '-n',str(config['mpi_ranks'])]
        # Syntax validation is not a molecular call. One rank avoids needless replicas.
        check = {}
        result['commands'].append(check)
        save(receipt,result)
        rc = bounded(base+[binary,'--check','-i','/task/input.inp'],
                     execution/'syntax.log',180,check)
        save(receipt,result)
        if rc:
            raise RuntimeError(f'Native syntax check failed, exit {rc}; no molecular call')
        run = {}
        result['commands'].append(run)
        result['status'] = 'running_native_almo'
        result['molecular_parent_calls'] = 1
        save(receipt,result)
        rc = bounded(base+mpi+[binary,'-i','/task/input.inp','-o','/task/endpoint.out'],
                     execution/'native_launcher.log',config['native_cap_seconds'],run)
        result['status'] = 'parsing_native_output'
        save(receipt,result)
        parser = [sys.executable,config['parser'],'--output',str(execution/'endpoint.out'),
                  '--geometry-id',gid,'--exit-code',str(rc),'--result',str(task/'RESULT.json')]
        parse = {}
        result['commands'].append(parse)
        prc = bounded(parser,execution/'parse.log',60,parse)
        result['status'] = 'native_finished_review_required' if rc == 0 and prc == 0 else 'native_or_admission_failed'
        result['parser_result'] = str(task/'RESULT.json') if (task/'RESULT.json').exists() else None
        exit_code = 0 if rc == 0 and prc == 0 else 1
    except BaseException as exc:
        result['status'] = 'failed'
        result['error'] = repr(exc)
        result['traceback'] = traceback.format_exc()
        print(result['traceback'], flush=True)
    finally:
        result['elapsed_seconds'] = time.monotonic()-started
        result['allocated_cpu_seconds_worker_elapsed'] = result['elapsed_seconds']*24
        result['cost_caveat'] = 'Allocated capacity, not measured CPU utilization; sacct supplements after termination.'
        result['finished_utc'] = now()
        result['artifacts'] = []
        for path in sorted(execution.glob('*')):
            if path.is_file() and path.suffix in ['.inp','.out','.log']:
                result['artifacts'].append({'path':str(path),'sha256':sha(path),'bytes':path.stat().st_size})
        save(receipt,result)
        dest.mkdir(parents=True,exist_ok=True)
        for name in ['RUN_RESULT.json','RESULT.json']:
            if (task/name).exists():
                shutil.copyfile(task/name,dest/name)
        for artifact in result['artifacts']:
            src = Path(artifact['path'])
            shutil.copyfile(src,dest/src.name)
    return exit_code


if __name__ == '__main__':
    sys.exit(main())
