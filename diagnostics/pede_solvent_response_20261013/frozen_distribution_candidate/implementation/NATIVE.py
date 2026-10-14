#!/usr/bin/env python3
"""Finite archived-charge native replay; no SCC, geometry search, or calibration."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import lru_cache
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import signal
import socket
import subprocess
import threading
import time

PROTOCOL = 'native_frozen_distribution_components_v1'
TERMS = ('Gelec', 'Ghb', 'Gsasa', 'Gshift')
VALENCE = {'H': 1, 'C': 4, 'N': 5, 'O': 6, 'S': 6, 'Ca': 2, 'La': 3}
STOP = threading.Event()


def read(p):
    return json.loads(Path(p).read_text())


def pin(p):
    p = Path(p).resolve()
    h = hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda: f.read(1024*1024), b''):
            h.update(b)
    return {'path': str(p), 'sha256': h.hexdigest()}


@lru_cache(None)
def verified_immutable(path, sha):
    # These inputs are immutable; avoid rehashing large shared manifests per cell.
    if pin(path) != {'path':path,'sha256':sha}:
        raise ValueError('changed pin: '+path)
    return Path(path)


def verify(p):
    return verified_immutable(p['path'],p['sha256'])


@lru_cache(None)
def checked_immutable(path, sha):
    return read(verified_immutable(path,sha))


def checked(p):
    return checked_immutable(p['path'],p['sha256'])


def dump(p, obj):
    with Path(p).open('x') as f:
        json.dump(obj, f, indent=2, allow_nan=False)
        f.write('\n')


def require(ok, text):
    if not ok:
        raise ValueError(text)


def scalar(text, label):
    a = re.findall(r'^\s*::\s*(?:->\s*)?'+re.escape(label)+r'\s+([-+0-9.Ee]+)\s+Eh\s*::', text, re.M)
    require(len(a) == 1, 'expected one archived summary component: '+label)
    v = float(a[0])
    require(math.isfinite(v), 'nonfinite archived component')
    return v


def command(args, cwd, stem, env, timer):
    """GNU time supplies per-child accounting even with parallel endpoints."""
    if STOP.is_set():
        raise RuntimeError('interrupted before dispatch')
    start = time.monotonic()
    measured = [timer, '-f', '%e %U %S %M', '-o', str(cwd/(stem+'.time')), '--', *args]
    with (cwd/(stem+'.out')).open('x') as out, (cwd/(stem+'.err')).open('x') as err:
        proc = subprocess.Popen(measured, cwd=cwd, env=env, stdout=out, stderr=err, start_new_session=True)
        stopped = False
        while proc.poll() is None:
            if STOP.is_set() and not stopped:
                try:
                    os.killpg(proc.pid, signal.SIGTERM)
                except ProcessLookupError:
                    pass
                stopped = True
            time.sleep(.02)
    row = {'command': args, 'measured_command': measured, 'cwd': str(cwd),
           'returncode': proc.returncode, 'wall_seconds': time.monotonic()-start,
           'stdout': pin(cwd/(stem+'.out')), 'stderr': pin(cwd/(stem+'.err')),
           'resource_receipt': pin(cwd/(stem+'.time'))}
    lines = (cwd/(stem+'.time')).read_text().splitlines()
    if lines:
        parts = lines[-1].split()
        if len(parts) == 4:
            try:
                elapsed, user, system, rss = map(float, parts)
                row.update(measured_wall_seconds=elapsed, child_CPU_seconds=user+system, peak_RSS_KiB=int(rss))
            except ValueError:
                pass
    dump(cwd/(stem+'_EXECUTION.json'), row)
    require(proc.returncode == 0, 'command failed; see '+str(cwd/(stem+'_EXECUTION.json')))
    return row


def common_preflight(m, mp):
    require(m['protocol'] == PROTOCOL, 'wrong protocol')
    require(pin(__file__) == m['implementation']['NATIVE.py'], 'runner differs from frozen snapshot')
    for p in m['implementation'].values():
        verify(p)
    for p in m['provenance'].values():
        verify(p)
    for p in m['build']['pins']:
        verify(p)
    old = read(verify(m['provenance']['qualified_manifest']))
    require(m['tolerances'] == old['tolerances'] and m['native_model'] == old['native_model'], 'native physics or precision changed')
    require(read(verify(m['provenance']['qualified_collection']))['status'] == 'complete', 'prior replay unavailable')
    inventory = read(verify(m['inventory']))
    require(inventory['summary']['unique_native_inputs_by_stage'] == m['task_counts'], 'finite scope differs')
    require(len({t['task_id'] for t in inventory['native_inputs']}) == len(inventory['native_inputs']), 'duplicate task ID')
    reuse = read(verify(m['reuse']))
    require(reuse['inventory'] == m['inventory'] and reuse['native_model'] == m['native_model'], 'reuse identity differs')
    return inventory, {r['native_task_id']: r for r in reuse['rows']}


def task_preflight(m, t):
    xyz = verify(t['xyz'])
    lines = xyz.read_text().splitlines()
    n = t['atom_count']
    require(int(lines[0]) == n and len(lines[2:]) == n, 'exact atom count differs')
    atoms = [line.split() for line in lines[2:]]
    require(all(len(a) == 4 and all(math.isfinite(float(v)) for v in a[1:]) for a in atoms), 'invalid coordinates')
    require(atoms[0][0] == t['metal'] and all(a[0] in ('H','C','N','O','S',t['metal']) for a in atoms), 'atom identity/order differs')
    require(sum(VALENCE[a[0]] for a in atoms)-t['charge'] == t['valence_electrons'], 'electron state differs')
    require(t['multiplicity'] == 1 and t['valence_electrons'] % 2 == 0, 'unsupported electronic state')
    state = read(verify(t['state_signature']))
    require(state['endpoint_charges'][t['metal']] == t['charge'], 'source state differs')
    verify(t['source_preparation'])
    for p in t['maps'].values():
        verify(p)
    runtime = read(verify(t['scientific_input_key']['runtime']))
    verify(runtime['xtb'])
    verify(runtime['xtb_parameter_file'])
    for use in t['uses']:
        r = checked(use['MACE'])
        require(r['status'] in ('complete','computed') and r['energy_eV'] == use['MACE_eV'], 'MACE receipt differs')
        s = r.get('state', r.get('input_state_check'))
        require(s['atoms'] == n and s['charge'] == t['charge'] and s['spin_multiplicity'] == 1 and s['head'] == 'omol', 'MACE physical state differs')
        if r.get('request'):
            request = checked(r['request'])
            model = r['model']
        else:
            manifest = checked(r['manifest'])
            request = next(q for q in manifest['tasks'] if q['task_id'] == r['task_id'])
            model = manifest['model']
        require(model == m['MACE_model'], 'MACE model differs')
        require(request['xyz']['sha256'] == t['xyz']['sha256'] and request['charge'] == t['charge'] and request.get('multiplicity',request.get('spin_multiplicity')) == 1, 'MACE exact geometry/state differs')
    archived = None
    for medium in ('water', 'vacuum'):
        e = t['endpoints'][medium]
        point = read(verify(e['point']))
        receipt = read(verify(e['execution']))
        text = verify(e['stdout']).read_text()
        verify(e['point_manifest'])
        require(point['manifest'] == e['point_manifest'] and point['task'] == e['task'], 'point provenance differs')
        require(point['backends']['xtb']['execution'] == e['execution'] and receipt['stdout'] == e['stdout'], 'native execution evidence differs')
        require(verify(e['executed_xyz']).read_bytes() == xyz.read_bytes() and receipt['returncode'] == 0, 'executed geometry/completion differs')
        qs = verify(e['charges']).read_text().split()
        require(len(qs) == n and all(re.fullmatch(r'[-+]?\d+\.\d{8}', q) for q in qs), 'saved charge precision/count differs')
        q = [float(v) for v in qs]
        require(q == point['backends']['xtb']['charges_e'] and point['backends']['xtb']['status'] == 'available', 'charge/result mismatch')
        require(abs(math.fsum(q)-t['charge']) <= n*m['tolerances']['charge_rounding_e']+1e-12, 'rounded charge closure failed')
        require(all(point['task'][k] == t[k] for k in ('metal','charge','multiplicity','valence_electrons','atom_count','xyz')), 'endpoint state differs')
        require(point['task']['medium'] == medium, 'endpoint medium differs')
        expected = [runtime['xtb']['path'], 'core.xyz', '--gfn', '2', '--chrg', str(t['charge']), '--uhf', '0', '--etemp', '300', '--acc', '0.0001', '--iterations', '500', '--norestart', '--json', '--pop']
        if medium == 'water':
            expected += ['--alpb', 'water', 'gsolv']
            for line in m['required_setup_lines']:
                require(line in text, 'native solvent setup differs: '+line)
            archived = {k: scalar(text, k) for k in TERMS}
            require(abs(math.fsum(archived.values())-scalar(text,'Gsolv')) <= m['tolerances']['native_component_Eh'], 'archived solvent closure differs')
        require(receipt['command'] == expected, 'native endpoint command differs')
    return archived


def compile_wrapper(m, mp, out, env):
    b = m['build']
    args = [b['compiler']['path'], '-O3', '-DNDEBUG', '-ffree-line-length-none']
    args += ['-I'+p for p in b['include_directories']]
    args += [m['implementation']['REPLAY.f90']['path'], '-o', str(out/'replay'), '-Wl,-rpath,'+b['library_directory'],
             b['xtb_static']['path'], b['blas']['path'], b['mctc_static']['path']]
    receipt = command(args,out,'compile_wrapper',env,b['timer']['path'])
    ld = command(['/usr/bin/ldd',str(out/'replay')],out,'libraries',env,b['timer']['path'])
    linked = []
    for line in verify(ld['stdout']).read_text().splitlines():
        found = re.search(r'(?:=>\s+)?(/\S+)\s+\(',line)
        if found:
            linked.append(pin(found.group(1)))
    require({p['path']: p['sha256'] for p in linked} == {p['path']: p['sha256'] for p in b['linked_libraries']}, 'linked runtime differs from qualified v3')
    row = {'status':'built','manifest':pin(mp),'binary':pin(out/'replay'), 'linked_libraries':linked,
           'compile':receipt,'engine_rebuilt':False,'qualified_engine_build':m['provenance']['qualified_build']}
    dump(out/'BUILD.json',row)
    return pin(out/'BUILD.json'), row


def evaluate(m, t, medium, binary, directory, env):
    directory.mkdir()
    for name, p in (('core.xyz',t['xyz']),('charges',t['endpoints'][medium]['charges'])):
        shutil.copyfile(verify(p),directory/name)
    local_env = dict(env, XTBPATH=str(directory), XTBHOME=str(directory))
    receipt = command([binary['path'],str(directory/'core.xyz'),str(directory/'charges'),str(t['atom_count'])], directory,'native',local_env,m['build']['timer']['path'])
    text = verify(receipt['stdout']).read_text()
    values = re.findall(r'^COMPONENTS\s+(.+)$',text,re.M)
    bounds = re.findall(r'^ROUNDING_BOUNDS\s+(.+)$',text,re.M)
    require(len(values) == len(bounds) == 1, 'missing unique native result')
    energy = dict(zip(TERMS,map(float,values[0].split())))
    rounding = dict(zip(('Gelec','Ghb'),map(float,bounds[0].split())))
    require(len(energy) == 4 and len(rounding) == 2 and all(math.isfinite(v) for v in [*energy.values(),*rounding.values()]), 'invalid native result')
    require('internal GFN2-xTB/ALPB' in text, 'internal native factory unproven')
    row = {'charge_source':medium,'components_Eh':energy,'rounding_bounds_Eh':rounding,
           'G_Eh':math.fsum(energy.values()),'execution':pin(directory/'native_EXECUTION.json'),
           'xyz':pin(directory/'core.xyz'),'charges':pin(directory/'charges'),
           'cavity':pin(directory/'cavity.tsv'),'born_matrix_with_HB':pin(directory/'born_matrix.dat')}
    dump(directory/'RESULT.json',row)
    return row | {'result':pin(directory/'RESULT.json')}


def checks(m, water, archived):
    rows = []
    for term in TERMS:
        bound = water['rounding_bounds_Eh'].get(term,0.0)
        tolerance = bound+m['tolerances']['native_component_Eh']
        residual = water['components_Eh'][term]-archived[term]
        informative = tolerance*m['HA_TO_KCAL'] <= m['tolerances']['maximum_informative_bound_kcal_mol']
        rows.append({'component':term,'residual_Eh':residual,'rounding_bound_Eh':bound,
                     'tolerance_Eh':tolerance,'informative':informative,'passed':informative and abs(residual)<=tolerance})
    return rows


def reused(m, t, reuse):
    oldm = read(verify(reuse['prior_manifest']))
    old = read(verify(reuse['prior_collection']))
    require(old['status'] == 'complete' and old['manifest'] == reuse['prior_manifest'], 'reuse completion differs')
    require(oldm['native_model'] == m['native_model'] and oldm['tolerances'] == m['tolerances'], 'reuse physics differs')
    oldcase = next(c for c in oldm['cases'] if c['case_id'] == reuse['prior_case_id'])
    require(all(oldcase[k] == t[k] for k in ('metal','charge','multiplicity','atom_count','valence_electrons','xyz','state_signature')), 'reuse source state differs')
    rows = []
    for medium in ('water','vacuum'):
        p = reuse[medium+'_replay']
        r = read(verify(p))
        require(r['charge_source'] == medium, 'reuse charge medium differs')
        for field in ('xyz','charges','execution','cavity','born_matrix_with_HB'):
            verify(r[field])
        require(r['xyz']['sha256'] == t['xyz']['sha256'] and r['charges']['sha256'] == t['endpoints'][medium]['charges']['sha256'], 'reuse exact input differs')
        rows.append(r | {'result':p})
    return rows


def run_task(m, mp, t, binary, out, env, reuse):
    directory = out/'tasks'/t['task_id']
    directory.mkdir(parents=True)
    row = {'task_id':t['task_id'],'status':'unavailable','G_qvac_hartree':None,'reason':None,
           'manifest':pin(mp),'input':t,'native_component_dispatch_attempts':0,'reuse':reuse}
    start = time.monotonic()
    try:
        if STOP.is_set():
            raise RuntimeError('interrupted before task')
        archived = task_preflight(m,t)
        row['archived_components_Eh'] = archived
        if reuse:
            water, vacuum = reused(m,t,reuse)
        else:
            row['native_component_dispatch_attempts'] += 1
            water = evaluate(m,t,'water',binary,directory/'water',env)
        row['water'] = water
        row['checks'] = checks(m,water,archived)
        require(all(c['passed'] for c in row['checks']), 'water component gate failed or uninformative; vacuum replay blocked')
        if not reuse:
            row['native_component_dispatch_attempts'] += 1
            vacuum = evaluate(m,t,'vacuum',binary,directory/'vacuum',env)
        row['vacuum'] = vacuum
        require(all(water[k]['sha256'] == vacuum[k]['sha256'] for k in ('cavity','born_matrix_with_HB')), 'cavity changed between distributions')
        row['vacuum_precision'] = [{'component':k,
            'rounding_bound_Eh':vacuum['rounding_bounds_Eh'].get(k,0.0),
            'tolerance_Eh':vacuum['rounding_bounds_Eh'].get(k,0.0)+m['tolerances']['native_component_Eh'],
            'informative':(vacuum['rounding_bounds_Eh'].get(k,0.0)+m['tolerances']['native_component_Eh'])*m['HA_TO_KCAL'] <= m['tolerances']['maximum_informative_bound_kcal_mol']} for k in TERMS]
        require(all(c['informative'] for c in row['vacuum_precision']), 'vacuum stored-charge uncertainty uninformative')
        row.update(status='available',G_qvac_hartree=vacuum['G_Eh'])
    except Exception as exc:
        row['reason'] = type(exc).__name__+': '+str(exc)
    row['wall_seconds'] = time.monotonic()-start
    row['native_component_calls'] = sum((directory/medium/'native_EXECUTION.json').exists() for medium in ('water','vacuum'))
    row['successful_component_evaluations'] = sum((directory/medium/'RESULT.json').exists() for medium in ('water','vacuum'))
    dump(directory/'RESULT.json',row)
    return {k:row[k] for k in ('status','G_qvac_hartree','reason','native_component_calls')} | {'result':pin(directory/'RESULT.json')}


def transfer_gate(m, mp, inventory, cp):
    cal = read(verify(cp))
    require(cal['protocol'] == 'frozen_distribution_candidate_v1' and cal['inventory'] == m['inventory'] and cal['native_manifest'] == pin(mp), 'calibration model/input identity differs')
    require(cal['calibrations']['operational']['status'] == 'available', 'canonical operational calibration failed; transfer blocked')
    prior = read(verify(cal['canonical_native_collection']))
    require(prior['status'] == 'complete' and prior['stage'] == 'canonical' and prior['manifest'] == pin(mp) and prior['inventory'] == m['inventory'], 'canonical receipt differs')
    expected = {t['task_id'] for t in inventory['native_inputs'] if t['stage'] == 'canonical'}
    require(set(prior['tasks']) == expected, 'canonical task completeness differs')
    for task_id, r in prior['tasks'].items():
        rr = read(verify(r['result']))
        require(rr['manifest'] == pin(mp) and rr['task_id'] == task_id, 'canonical task identity differs')
        require(r['status'] == rr['status'] == 'available' and r['G_qvac_hartree'] == rr['G_qvac_hartree'], 'unavailable canonical task')
    bp = prior['build']
    build = read(verify(bp))
    require(build['manifest'] == pin(mp), 'canonical build identity differs')
    verify(build['binary'])
    for p in build['linked_libraries']:
        verify(p)
    return bp, build


def main():
    p = argparse.ArgumentParser()
    p.add_argument('stage',choices=('canonical','transfer'))
    p.add_argument('manifest')
    p.add_argument('--sha',required=True)
    p.add_argument('--output',required=True)
    p.add_argument('--calibration')
    p.add_argument('--calibration-sha')
    a = p.parse_args()
    require(bool(os.environ.get('SLURM_JOB_ID')), 'actual Slurm allocation required; no files created')
    workers = int(os.environ.get('SLURM_NTASKS','0'))
    require(1 <= workers <= 12 and len(os.sched_getaffinity(0)) >= workers, 'requires one to twelve allocated CPUs')
    require(int(os.environ.get('SLURM_JOB_NUM_NODES',os.environ.get('SLURM_NNODES','0'))) == 1, 'all replay workers require one allocated node')
    require(int(os.environ.get('SLURM_CPUS_ON_NODE','0').split('(')[0]) >= workers, 'insufficient CPUs allocated on batch host')
    mp = Path(a.manifest).resolve()
    require(pin(mp)['sha256'] == a.sha, 'manifest pin differs')
    m = read(mp)
    out = Path(a.output).resolve()
    out.mkdir(parents=True,exist_ok=False)
    for sig in (signal.SIGTERM,signal.SIGINT,signal.SIGHUP):
        signal.signal(sig,lambda *_:STOP.set())
    result = {'protocol':PROTOCOL,'status':'failed','terminal':True,'stage':a.stage,'manifest':pin(mp),
              'inventory':m['inventory'],'build':None,'tasks':{},'SCC_calls':0,'DFT_calls':0,
              'slurm_job_id':os.environ['SLURM_JOB_ID'],'hostname':socket.gethostname(),
              'affinity':sorted(os.sched_getaffinity(0)),'workers':workers,'numerical_threads_per_worker':1}
    env = {k:v for k,v in os.environ.items() if k not in ('LD_PRELOAD','PYTHONPATH','CPATH','LIBRARY_PATH','CMAKE_PREFIX_PATH')}
    env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',LD_LIBRARY_PATH=m['build']['library_directory'])
    started = time.monotonic()
    try:
        inventory, reuse = common_preflight(m,mp)
        if a.stage == 'transfer':
            require(a.calibration and a.calibration_sha, 'transfer requires frozen canonical calibration')
            cp = {'path':str(Path(a.calibration).resolve()),'sha256':a.calibration_sha}
            bp, build = transfer_gate(m,mp,inventory,cp)
            result['calibration'] = cp
        else:
            require(not a.calibration and not a.calibration_sha, 'canonical stage cannot import calibration')
            bp, build = compile_wrapper(m,mp,out,env)
        result['build'] = bp
        tasks = [t for t in inventory['native_inputs'] if t['stage'] == a.stage]
        dump(out/'PREFLIGHT.json',{'status':'passed','manifest':pin(mp),'inventory':m['inventory'],
                                 'stage':a.stage,'task_count':len(tasks),'build':bp})
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = {pool.submit(run_task,m,mp,t,build['binary'],out,env,reuse.get(t['task_id'])):t['task_id'] for t in tasks}
            for future in as_completed(futures):
                result['tasks'][futures[future]] = future.result()
        result['status'] = 'interrupted' if STOP.is_set() else 'complete'
    except Exception as exc:
        result['reason'] = type(exc).__name__+': '+str(exc)
    finally:
        result['wall_seconds'] = time.monotonic()-started
        result['available'] = sum(r['status']=='available' for r in result['tasks'].values())
        result['unavailable'] = sum(r['status']=='unavailable' for r in result['tasks'].values())
        result['native_component_calls'] = sum(r['native_component_calls'] for r in result['tasks'].values())
        dump(out/'FINAL_COLLECTION.json',result)
    raise SystemExit(0 if result['status']=='complete' else 1)


if __name__ == '__main__':
    main()
