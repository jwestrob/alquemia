#!/usr/bin/env python3
"""Two explicitly bound operator initializations; no response solve or energy."""
import datetime
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import resource
import shutil
import signal
import socket
import subprocess
import sys
import time


def pin(p):
    p = Path(p).resolve()
    return dict(path=str(p), sha256=hashlib.sha256(p.read_bytes()).hexdigest())


def verify(p):
    if pin(p['path']) != p:
        raise ValueError('Changed pin: ' + p['path'])
    return Path(p['path'])


def read(p):
    return json.loads(verify(p).read_text())


def module(p, name):
    spec = importlib.util.spec_from_file_location(name, verify(p))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def main(config_path, out_path):
    start, cpu = time.monotonic(), time.process_time()
    out = Path(out_path).resolve()
    out.mkdir(parents=True, exist_ok=False)
    result = dict(status='started', error=None, configuration=pin(config_path),
                  implementation=pin(__file__), host=socket.gethostname(), PID=os.getpid(),
                  started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  execution_mode='login', scheduler_allocated_CPU_seconds=0,
                  native_calls=0, link_calls=0, induced_solves=0, energy_calls=0,
                  gradient_calls=0, full_stability_admitted=False, physical_accuracy_qualified=False,
                  cells={m: dict(status='not_attempted') for m in ('La', 'Dy')})

    def save():
        tmp = out/'RESULT.tmp'
        tmp.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
        os.replace(tmp, out/'RESULT.json')

    def command(argv, cwd, name, native=False):
        usage = resource.getrusage(resource.RUSAGE_CHILDREN)
        tick = time.monotonic()
        rec = dict(argv=list(map(str, argv)), cwd=str(cwd), status='starting', wall_limit_seconds=120,
                   address_space_limit_bytes=8*1024**3)
        result[name] = rec
        save()
        def limits():
            resource.setrlimit(resource.RLIMIT_AS, (8*1024**3, 8*1024**3))
        try:
            with (out/(name+'.stdout')).open('w') as so, (out/(name+'.stderr')).open('w') as se:
                proc = subprocess.Popen(rec['argv'], cwd=cwd, stdout=so, stderr=se,
                                        start_new_session=True, preexec_fn=limits)
                rec.update(status='running', PID=proc.pid)
                result['native_calls' if native else 'link_calls'] += 1
                save()
                try:
                    rec['returncode'] = proc.wait(timeout=120)
                except subprocess.TimeoutExpired:
                    os.killpg(proc.pid, signal.SIGKILL)
                    proc.wait()
                    rec['returncode'] = proc.returncode
                    raise RuntimeError(name+' exceeded declared limit')
            if rec['returncode'] != 0:
                raise RuntimeError(name+' returned '+str(rec['returncode']))
            rec['status'] = 'complete'
        except Exception as e:
            rec.update(status='failed', error=repr(e))
            raise
        finally:
            end = resource.getrusage(resource.RUSAGE_CHILDREN)
            rec.update(wall_seconds=time.monotonic()-tick,
                       child_CPU_seconds=end.ru_utime+end.ru_stime-usage.ru_utime-usage.ru_stime,
                       stdout=pin(out/(name+'.stdout')), stderr=pin(out/(name+'.stderr')))
            save()

    try:
        cfg = json.loads(Path(config_path).read_text())
        release = read(cfg['root_release'])
        expected = {k:v for k,v in cfg.items() if k != 'root_release'}
        assert release['configuration_sources'] == expected
        assert release['native_call_cap'] == 2 and release['link_cap'] == 1
        assert cfg['implementation'] == pin(__file__)
        assert not os.environ.get('SLURM_JOB_ID') and not os.environ.get('SLURM_STEP_ID')
        assert socket.gethostname() == cfg['hostname'] and str(out) == cfg['output_directory']
        assert all(os.environ.get(k) == '1' for k in ('OMP_NUM_THREADS', 'OMP_THREAD_LIMIT',
                   'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'BLIS_NUM_THREADS'))
        ready = read(cfg['ready'])
        review = read(cfg['root_source_review'])
        assert review['status'] == 'passed_operator_initialization_source_review'
        assert review['ready'] == cfg['ready'] and review['implementation'] == cfg['implementation']
        for key in ('harness', 'derive', 'inputs', 'specification', 'controls'):
            verify(ready[key])
        inputs = read(ready['inputs'])
        for value in inputs.values():
            if 'path' in value:
                verify(value)
        for value in inputs['native_sources'].values():
            verify(value)
        composition = read(cfg['composition'])
        assert composition['status'] == 'prepared_only_not_executable' and composition['engine_binding'] is None
        assert composition['implementation'] == inputs['composition_source'] and composition['inputs'] == inputs['composition_inputs']
        candidate = composition['candidate']
        fr, ar, selected = read(candidate['result']), read(candidate['actual_review']), read(candidate['selected'])
        assert ar['status'] == 'passed_component_fit_actual_review' and ar['selected_native_admission_passed'] is True
        assert ar['result'] == candidate['result'] and ar['selected'] == candidate['selected'] == fr['selected']
        assert fr['configuration'] == candidate['configuration'] and fr['error'] is None
        read(candidate['configuration'])
        assert selected['status'] == 'complete' and fr['selected_key'] == composition['selected_parameter_key']
        qualification = read(cfg['qualification'])
        assert qualification['status'] == 'initialized_PME_export_qualification_complete' and qualification['error'] is None
        assert qualification['library'] == inputs['read_only_library_authority']
        assert qualification['source_manifest'] == inputs['read_only_source_authority']
        for value in qualification['modules']:
            verify(value)
        derive = module(ready['derive'], 'operator_derive')
        snapshot = module(inputs['snapshot_reader'], 'operator_snapshot')
        exe = out/'init_pme'
        argv = list(qualification['link']['argv'])
        argv[argv.index(qualification['link']['argv'][7])] = ready['harness']['path']
        assert argv[7] == ready['harness']['path']
        argv[-1] = str(exe)
        assert argv[-2] == '-o'
        command(argv, out, 'link')
        result['executable'] = pin(exe)
        result.update(candidate=candidate, composition=cfg['composition'], ready=cfg['ready'],
                      library=inputs['read_only_library_authority'], source_manifest=inputs['read_only_source_authority'])
        for row in composition['files']:
            metal = row['metal']
            assert metal in ('La', 'Dy') and result['cells'][metal]['status'] == 'not_attempted'
            cell = out/metal
            cell.mkdir()
            for key, name in [('xyz', 'system.xyz'), ('parameters', 'amoebabio09.prm')]:
                assert row['prepared'][key]['sha256'] == row['original'][key]['sha256']
                shutil.copyfile(verify(row['prepared'][key]), cell/name)
            key, derivation = derive.derive(verify(row['prepared']['key']).read_text())
            (cell/'system.key').write_text(key)
            result['cells'][metal].update(status='prepared', derivation=derivation,
                inputs={name: pin(cell/name) for name in ('system.xyz', 'system.key', 'amoebabio09.prm')})
            save()
            assert result['native_calls'] < 2
            verify(result['executable'])
            command([exe, 'system.xyz', metal], cell, 'native_'+metal, native=True)
            receipt = (cell/'operator_initialization.dat').read_text()
            assert receipt.splitlines()[-1] == 'END'
            assert 'reciprocal_initializations=1 induced_solves=0 energies=0 gradients=0' in receipt
            with (cell/'pme_norm_export_001.dat').open() as f:
                state = snapshot.read(f)
            assert len(state['sites']) == 11958 and state['shape'] == (64, 64, 64) and state['p'] == 5
            assert state['ewald'] == 0.4 and state['fft'] == 'FFTW'
            replay = snapshot.replay(state)
            assert replay['q_negative_count'] == 0 and replay['q_zero'] == 0
            result['cells'][metal].update(status='complete', operator_norms=replay,
                native_norms=state['native_norms'], export=pin(cell/'pme_norm_export_001.dat'),
                initialization_receipt=pin(cell/'operator_initialization.dat'))
            save()
        assert result['native_calls'] == 2 and all(v['status'] == 'complete' for v in result['cells'].values())
        result['status'] = 'full_box_operator_initialization_complete'
    except Exception as e:
        result.update(status='failed', error=repr(e))
    finally:
        usage = resource.getrusage(resource.RUSAGE_CHILDREN)
        result.update(wall_seconds=time.monotonic()-start, controller_CPU_seconds=time.process_time()-cpu,
                      child_CPU_seconds=usage.ru_utime+usage.ru_stime, peak_child_RSS_KiB=usage.ru_maxrss,
                      threads=1, maximum_native_workers=1)
        save()
    print(json.dumps({k:result[k] for k in ('status', 'error', 'native_calls', 'link_calls', 'wall_seconds')}))
    return 0 if result['status'] == 'full_box_operator_initialization_complete' else 1


if __name__ == '__main__':
    raise SystemExit(main(*sys.argv[1:]))
