#!/usr/bin/env python3
"""Two explicitly bound operator initializations; no response solve or energy."""
import datetime
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import resource
import shutil
import signal
import socket
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
EPS = sys.float_info.epsilon


def require(ok, message):
    if not ok:
        raise ValueError(message)


def norm_comparison(state, replay):
    """Prospective floating consistency check, not an interval error certificate."""
    operations = 3*len(state['sites'])*state['p']**3+128
    n_eps = operations*EPS
    require(n_eps < .01, 'Accumulation domain exceeded')
    # Two implementations accumulate nonnegative terms in different orders.
    # 32 gamma_m is a fixed conservative consistency allowance, not e_app.
    norm_relative = max(1e-12, 32*n_eps/(1-n_eps))
    keys = ('c1','cinf','cF','U_R_nominal','evenness_perturbation_nominal',
            'q_min','q_max','q_abs_max','q_zero','q_conjugate_difference_max')
    raw = list(state['native_norms'])+list(state['q_stats'])
    require(len(raw) == len(keys), 'Complete native summary')
    checks = []
    for i, (key, actual) in enumerate(zip(keys, raw)):
        expected = replay[key]
        require(math.isfinite(actual) and math.isfinite(expected), 'Finite summary '+key)
        relative = norm_relative if i < 5 else 1e-12
        bound = relative*max(1., abs(actual), abs(expected))
        error = abs(actual-expected)
        checks.append(dict(quantity=key,native=actual,replay=expected,error=error,
                           bound=bound,passed=error<=bound))
    return dict(checks=checks,passed=all(x['passed'] for x in checks),
                accumulation_operations=operations,norm_relative_allowance=norm_relative,
                q_stat_relative_allowance=1e-12,
                interpretation='Ordinary-double summary consistency only; not outward rounding or full-box application error.')


def spline5(w):
    """Pinned pmestuf.f bsplgen value/first-derivative recurrence, order five."""
    a = [[0.]*6 for _ in range(6)]
    a[2][1],a[2][2] = 1.-w,w
    a[3][3] = .5*w*a[2][2]
    a[3][2] = .5*((1.+w)*a[2][1]+(2.-w)*a[2][2])
    a[3][1] = .5*(1.-w)*a[2][1]
    for i in (4,5):
        k=i-1;denom=1./k;a[i][i]=denom*w*a[k][k]
        for j in range(1,i-1):
            a[i][i-j]=denom*((w+j)*a[k][i-j-1]+(i-j-w)*a[k][i-j])
        a[i][1]=denom*(1.-w)*a[k][1]
    a[4][5]=a[4][4]
    for i in range(4,1,-1):a[4][i]=a[4][i-1]-a[4][i]
    a[4][1]=-a[4][1]
    return [[a[5][i],a[4][i]] for i in range(1,6)]


def input_sites(xyz, parameters, key):
    rows=[line.split() for line in xyz.splitlines() if line.strip()]
    require(int(rows[0][0]) == 11958 and len(rows) == 11960, 'Exact BIO XYZ denominator')
    cell=list(map(float,rows[1]));require(cell == [49.324]*3+[90.]*3, 'Exact orthogonal cell')
    alpha={}
    for text in (parameters,key):
        for line in text.splitlines():
            fields=line.split()
            if fields and fields[0].lower() == 'polarize':
                kind=int(fields[1]);value=float(fields[2].replace('D','E').replace('d','e'))
                require(kind>0 and math.isfinite(value) and value>=0, 'Only finite positive-type POLARIZE records')
                alpha[kind]=value
    sites=[]
    for i,fields in enumerate(rows[2:],1):
        kind=int(fields[5]);xyz=list(map(float,fields[2:5]))
        require(int(fields[0])==i and kind in (402,403,405,414,900,901) and
                all(math.isfinite(v) for v in xyz) and alpha.get(kind,0)>=1e-8, 'Exact input site/alpha domain')
        sites.append(dict(atom=i,type=kind,xyz=xyz,alpha=alpha[kind]))
    return sites,cell


def check_sites(state, expected, cell):
    require(len(state['sites'])==len(expected)==11958 and state['shape']==(64,64,64) and
            state['p']==5 and state['ewald']==.4 and state['fft']=='FFTW', 'Exact initialized snapshot domain')
    require(state['chunks']==[1,1,1,64,64,64,2,2,4], 'Exact one-thread native chunk convention')
    volume=cell[0]*cell[1]*cell[2]
    require(abs(state['volume']-volume)<=128*EPS*max(1.,abs(volume)), 'Native cell volume')
    for i in range(3):
        for j in range(3):
            target=64./cell[i] if i==j else 0.
            require(abs(float(state['trans'][i,j])-target)<=128*EPS*max(1.,abs(target)), 'Native normalized transform')
    maximum=0.;worst=None
    for site,original in zip(state['sites'],expected):
        require(all(site[k]==original[k] for k in ('atom','type','xyz','alpha')), 'Exact native site identity/type/XYZ/intrinsic alpha')
        for axis in range(3):
            w=original['xyz'][axis]*(1./cell[axis])
            nearest=math.copysign(math.floor(abs(w)+.5),w)
            fr=64.*(w-nearest+.5);ifr=int(fr-1e-8)
            require(site['igrid'][axis]==ifr-5, 'Exact native spline grid index')
            coefficients=spline5(fr-ifr)
            for slot in range(5):
                for derivative in range(2):
                    actual=float(site['theta'][axis,slot,derivative]);target=coefficients[slot][derivative]
                    error=abs(actual-target);bound=256*EPS*max(1.,abs(actual),abs(target))
                    require(error<=bound, 'Native B-spline value/first derivative')
                    if error>maximum:
                        maximum=error;worst=dict(atom=site['atom'],axis=axis+1,slot=slot+1,derivative=derivative,native=actual,replay=target,error=error,bound=bound)
    return dict(passed=True,sites=len(expected),spline_scalars=30*len(expected),
                maximum_spline_difference=maximum,worst_spline_difference=worst,
                coordinate_alpha_comparison='exact parsed doubles',spline_relative_allowance=256*EPS)


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
                  resource_limits=dict(child_address_space_bytes=8*1024**3,child_wall_seconds=120,
                      controller_replay_wall_cap=None,controller_replay_domain='two fixed11958-site64^3 snapshots, one thread; no native application'),
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
        proc = None
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
            if proc is not None and proc.poll() is None:
                os.killpg(proc.pid, signal.SIGKILL)
                proc.wait()
                rec['returncode'] = proc.returncode
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
        assert release['status'] == 'released_two_full_box_operator_initializations'
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
        assert fr['status'] == 'bounded_component_fit_complete' and fr['candidate_frozen'] is True
        assert ar['configuration'] == fr['configuration'] and ar['implementation'] == fr['implementation']
        read(candidate['configuration'])
        assert selected['status'] == 'complete' and fr['selected_key'] == composition['selected_parameter_key']
        assert selected['parameters'] == composition['selected_parameters']
        assert hashlib.sha256(json.dumps(selected['parameters'],sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest() == fr['selected_key']
        assert [row['metal'] for row in composition['files']] == ['La','Dy']
        composition_review = read(composition['source_review'])
        assert composition_review['status'] == 'passed_periodic_model_composition_source_review'
        assert composition_review['implementation'] == composition['implementation']
        assert read(composition['configuration'])['candidate'] == candidate
        qualification = read(cfg['qualification'])
        assert qualification['status'] == 'initialized_PME_export_qualification_complete' and qualification['error'] is None
        assert qualification['library'] == inputs['read_only_library_authority']
        assert qualification['source_manifest'] == inputs['read_only_source_authority']
        actual = read(inputs['finite_export_root_actual_review'])
        assert actual['status'] == 'admitted_finite_initialized_PME_export_operator_qualification'
        assert actual['execution'] == cfg['qualification'] and actual['source_manifest'] == qualification['source_manifest']
        assert actual['build_result'] == qualification['build_result']
        build = read(actual['build_result'])
        assert build['status'] == 'complete' and build['manifest'] == qualification['source_manifest']
        assert build['library'] == qualification['library']
        verify(actual['saved_replay'])
        assert qualification['link']['status'] == 'complete' and qualification['link']['returncode'] == 0
        assert cfg['link_dependencies'] == pin(HERE/'LINK_DEPENDENCIES.json')
        dependencies = read(cfg['link_dependencies'])
        assert dependencies['qualified_execution'] == cfg['qualification']
        assert dependencies['qualified_argv'] == qualification['link']['argv']
        for value in dependencies['files']:
            verify(value)
        for value in qualification['modules']:
            verify(value)
        derive = module(ready['derive'], 'operator_derive')
        snapshot = module(inputs['snapshot_reader'], 'operator_snapshot')
        exe = out/'init_pme'
        argv = list(qualification['link']['argv'])
        assert len(argv) == 15 and argv[8] == qualification['library']['path'] and argv[6].startswith('-I')
        assert {str(Path(value['path']).parent) for value in qualification['modules']} == {argv[6][2:]}
        assert [value['path'] for value in dependencies['files']] == [str(Path(argv[i]).resolve()) for i in (0,9,10,11,12)]
        argv[7] = ready['harness']['path']
        argv[-1] = str(exe)
        assert argv[-2] == '-o'
        command(argv, out, 'link')
        for value in dependencies['files']:
            verify(value)
        verify(inputs['read_only_library_authority'])
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
                verify(row['original'][key])
                shutil.copyfile(verify(row['prepared'][key]), cell/name)
            key, derivation = derive.derive(verify(row['prepared']['key']).read_text())
            (cell/'system.key').write_text(key)
            result['cells'][metal].update(status='prepared', derivation=derivation,
                inputs={name: pin(cell/name) for name in ('system.xyz', 'system.key', 'amoebabio09.prm')})
            save()
            expected_sites, expected_cell = input_sites((cell/'system.xyz').read_text(),
                (cell/'amoebabio09.prm').read_text(), key)
            assert result['native_calls'] < 2
            verify(result['executable'])
            command([exe, 'system.xyz', metal], cell, 'native_'+metal, native=True)
            receipt = (cell/'operator_initialization.dat').read_text()
            records = [line.split() for line in receipt.splitlines()]
            assert records == [['PME_OPERATOR_ONLY_INITIALIZATION_V1'],['METAL',metal],
                ['SITES','11958','11958'],['COUNTS','reciprocal_initializations=1','induced_solves=0','energies=0','gradients=0'],
                ['ENERGY_MODEL_ADMISSION','FALSE'],['FULL_STABILITY_ADMISSION','FALSE'],['END']]
            assert list(cell.glob('pme_norm_export_*.dat')) == [cell/'pme_norm_export_001.dat']
            result['cells'][metal].update(export=pin(cell/'pme_norm_export_001.dat'),
                initialization_receipt=pin(cell/'operator_initialization.dat'))
            save()
            with (cell/'pme_norm_export_001.dat').open() as f:
                state = snapshot.read(f)
            assert len(state['sites']) == 11958 and state['shape'] == (64, 64, 64) and state['p'] == 5
            assert state['ewald'] == 0.4 and state['fft'] == 'FFTW'
            site_check = check_sites(state, expected_sites, expected_cell)
            replay = snapshot.replay(state)
            comparison = norm_comparison(state, replay)
            result['cells'][metal].update(site_checks=site_check,summary_comparison=comparison,
                operator_norms=replay,native_norms=state['native_norms'],native_q_stats=state['q_stats'])
            save()
            assert comparison['passed'], 'Native/replay summary mismatch'
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
