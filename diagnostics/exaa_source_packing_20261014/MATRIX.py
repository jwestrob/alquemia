#!/usr/bin/env python3
"""Finite fixed-geometry bridge to the existing warm OMOL and native ALPB kernels."""
import argparse
import copy
from concurrent.futures import ThreadPoolExecutor
import fcntl
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import signal
import socket
import subprocess
import sys
import time
import traceback

ROOT = next(p for p in Path(__file__).resolve().parents if (p/'scripts/context_factorial_gfn2.py').is_file())
HERE = ROOT/'diagnostics/exaa_source_packing_20261014'
ORDER = ('origin', 'adaptive_Ca', 'adaptive_La')
METALS = ('Ca', 'La')
THREADS = ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS')
TEMPLATE = ROOT/'workspaces/ca_coverage_recovery_20261012/two_ca_scoring_q9_v1/lane0/worker_manifest.json'
NATIVE = ROOT/'workspaces/pede_solvent_response_20261013/frozen_distribution_candidate_v1/manifest.json'
BUILD = NATIVE.parent/'canonical_v1/BUILD.json'
RUNTIME = ROOT/'workspaces/orca_parameter_audit_20261012/context_gfn2_v1/RUNTIME.json'
COMPONENT = ROOT/'diagnostics/orca_parameter_audit_20261012/expanded_transfer_inventory/REEVALUATE.py'
FACTORIAL = ROOT/'scripts/context_factorial_gfn2.py'
EV = 23.06054783061903
HA = 627.509474


def pin(path):
    p = Path(path).resolve(); h = hashlib.sha256()
    with p.open('rb') as f:
        for chunk in iter(lambda: f.read(1048576), b''): h.update(chunk)
    return {'path': str(p), 'sha256': h.hexdigest()}


def verify(p):
    assert pin(p['path']) == p, 'changed artifact: '+p['path']
    return Path(p['path'])


def read(p):
    return json.loads((verify(p) if isinstance(p, dict) else Path(p)).read_text())


def write(p, data):
    p = Path(p); p.parent.mkdir(parents=True, exist_ok=True)
    with p.open('x') as f: json.dump(data, f, indent=2, allow_nan=False); f.write('\n')


def load(name, p):
    spec = importlib.util.spec_from_file_location(name, verify(p)); m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m); return m


def prepare(input_path, output):
    """Metadata-only; molecular admission is performed inside Slurm."""
    prepared = read(input_path)
    assert prepared['status'] == 'prepared' and len(prepared['rows']) == 2
    old = read(TEMPLATE); native = read(NATIVE); build = read(BUILD)
    base = Path(output).resolve(); base.mkdir(parents=True, exist_ok=False)
    rows = prepared['rows']; tasks = []
    for row in rows:
        assert row['status'] == 'prepared'
        for metal in METALS:
            assert set(row['matrix'][metal]) == set(ORDER)
            for candidate in ORDER:
                cell = row['matrix'][metal][candidate]
                tasks.append({'task_id': row['case_id']+'__'+metal+'__'+candidate,
                    'case_id': row['case_id'], 'metal': metal, 'candidate': candidate,
                    **{k: cell[k] for k in ('xyz','charge','multiplicity','atom_count','valence_electrons')},
                    **{k: row[k] for k in ('maps','state_signature','source_preparation')}})
    assert len(tasks) == 12 and len({t['task_id'] for t in tasks}) == 12
    fixtures = [next(t for t in old['tasks'] if t['task_id'] == f['task_id']) for f in old['fixtures']]
    worker = copy.deepcopy(old)
    worker.update(protocol_id='Nikasha_ExaA_fixed_pool_packing_v1', tasks=tasks+fixtures,
        cases=[], origin_task_count=len(tasks), lanes=1, preparation=pin(input_path),
        authorization='Jacob approved ExaA fixed-pool packing test; Slurm only',
        finite_scope={'fixed_geometry_cells':12, 'qualification_cells':2, 'optimizer_starts':0})
    write(base/'worker_manifest.json', worker)
    runtime_source = read(ROOT/'workspaces/orca_parameter_audit_20261012/context_gfn2_v1/manifest.json')
    config = {'protocol_id':'Nikasha_ExaA_fixed_pool_packing_v1', 'implementation':pin(__file__),
        'plan':pin(HERE/'PLAN.md'), 'input':pin(input_path), 'template':pin(TEMPLATE),
        'worker_manifest':pin(base/'worker_manifest.json'), 'component':pin(COMPONENT),
        'factorial':pin(FACTORIAL), 'runtime':pin(RUNTIME), 'settings':runtime_source['settings'],
        'native_manifest':pin(NATIVE), 'native_build':pin(BUILD), 'native_binary':build['binary'],
        'native_adapter':native['implementation']['NATIVE.py'], 'MACE_model':old['model'],
        'cells':tasks, 'rows':rows, 'CPUs':7, 'new_optimization':False, 'new_calibration':False,
        'production_changed':False}
    write(base/'CONFIG.json', config)
    print(json.dumps({'config':pin(base/'CONFIG.json'),'cells':12,'maximum_SCC_calls':48,'MACE_calls':14}))


def placement():
    job = os.environ.get('SLURM_JOB_ID'); assert job and job.isdigit(), 'Actual Slurm required'
    raw = subprocess.run(['scontrol','show','job','-o',job],capture_output=True,text=True,check=True).stdout
    fields = dict(v.split('=',1) for v in raw.split() if '=' in v)
    assert fields['JobState'] == 'RUNNING' and int(fields['NumCPUs']) == 7 and int(fields['NumNodes']) == 1
    assert int(os.environ['SLURM_CPUS_ON_NODE']) == 7 and len(os.sched_getaffinity(0)) >= 7
    assert all(os.environ.get(k) == '1' for k in THREADS)
    return {'job_id':job,'scontrol':raw,'hostname':socket.gethostname(),'allocated_CPUs':7,
            'affinity':sorted(os.sched_getaffinity(0))}


def preflight(c, f, b):
    m = read(c['worker_manifest']); old = read(c['template']); native = read(c['native_manifest'])
    for k in ('model','settings','pool_settings','optimizer_software','source_modules','calculator_source',
              'model_python','optimizer_python','model_executable','optimizer_executable','fixtures','gates'):
        assert m[k] == old[k], 'original kernel changed: '+k
    for p in m['source_modules'].values(): verify(p)
    for k in ('model_executable','optimizer_executable','calculator_source','implementation'): verify(m[k])
    verify(m['model']['checkpoint']); assert m['model'] == native['MACE_model'] == c['MACE_model']
    runtime = read(c['runtime'])
    for name in ('xtb','xtb_parameter_file','tblite_python'): verify(runtime[name])
    for p in runtime['xtb_dynamic']['libraries']+runtime['tblite_distribution']['distribution_files']: verify(p)
    for d in runtime['tblite_dynamic'].values():
        for p in d['libraries']: verify(p)
    for p in native['build']['pins']: verify(p)
    build = read(c['native_build']); verify(build['binary'])
    assert build['manifest'] == c['native_manifest'] and build['binary'] == c['native_binary']
    for p in build['linked_libraries']: verify(p)
    b.modules(m)
    rows = []
    for row in c['rows']:
        state = read(row['state_signature']); verify(row['source_preparation'])
        for p in row['maps'].values(): verify(p)
        reference = {}
        for candidate in ORDER:
            coords = {}
            for z in METALS:
                t = row['matrix'][z][candidate]; symbols, positions = f.coordinates(t['xyz'])
                assert symbols[0] == z and all(s in f.VALENCE for s in symbols)
                assert len(symbols) == t['atom_count'] and t['charge'] == state['endpoint_charges'][z]
                assert sum(f.VALENCE[s] for s in symbols)-t['charge'] == t['valence_electrons']
                assert t['multiplicity'] == 1 and t['valence_electrons'] % 2 == 0
                assert all(math.isfinite(v) for xyz in positions for v in xyz)
                reference.setdefault(z, symbols); assert reference[z] == symbols
                coords[z] = (symbols, positions, t['charge'])
            ca, la = coords['Ca'], coords['La']
            assert ca[0][1:] == la[0][1:] and ca[1] == la[1] and la[2]-ca[2] == 1
        rows.append({'case_id':row['case_id'],'paired_coordinates_exact':True,'composition_fixed':True})
    return {'status':'passed','rows':rows,'model_and_native_runtime_pins':'verified','new_model_calls':0}


def mace(c, base, b):
    m = read(c['worker_manifest']); mp = verify(c['worker_manifest'])
    attempt = base/'mace'; (attempt/'lane0').mkdir(parents=True)
    transport = None; old_handlers = b.handlers(); results = {}; error = None
    try:
        import numpy as np
        transport = b.WarmCPU(mp, m, attempt, 0)
        checks = []
        for fixture in m['fixtures']:
            t = next(t for t in m['tasks'] if t['task_id'] == fixture['task_id'])
            request = {'manifest':c['worker_manifest'],'task_id':t['task_id'],
                       **{k:t[k] for k in ('xyz','charge','multiplicity')}}
            r, rp = b.cell(attempt/'qualification'/fixture['id'],request,transport,m['model'])
            check = {'fixture':fixture['id'],'receipt':rp,'pass':False}
            if r['status'] == 'complete':
                de = (r['energy_eV']-fixture['energy_eV'])*EV
                df = float(np.max(np.abs(np.load(verify(r['forces']),allow_pickle=False)-np.load(verify(fixture['forces']),allow_pickle=False))))
                check.update(energy_residual_kcal_mol=de, force_residual_eV_A=df,
                             **{'pass':abs(de)<=m['gates']['energy_kcal_mol'] and df<=m['gates']['force_max_eV_A']})
            checks.append(check)
        write(attempt/'QUALIFICATION.json',{'checks':checks,'pass':all(r['pass'] for r in checks)})
        assert all(r['pass'] for r in checks), 'MACE runtime qualification failed'
        for t in c['cells']:
            request = {'manifest':c['worker_manifest'],'task_id':t['task_id'],
                       **{k:t[k] for k in ('xyz','charge','multiplicity')}}
            r, rp = b.cell(attempt/'cells'/t['task_id'],request,transport,m['model'])
            results[t['task_id']] = {'receipt':rp,'status':r['status'],'energy_eV':r['energy_eV']}
    except BaseException as exc:
        error = str(exc); raise
    finally:
        if transport:
            if error: transport.abort()
            else: transport.close()
        for sig, handler in old_handlers.items(): signal.signal(sig,handler)
        write(attempt/'COLLECTION.json',{'rows':results,'reason':error,'expected_cells':12})
    return results


def solvent(c, base, f, r, n, mr):
    directory = base/'solvent'; directory.mkdir()
    tasks = []
    for cell in c['cells']:
        for medium in ('vacuum','water'):
            tasks.append({'id':cell['task_id']+'__'+medium,'protein':'Q9Z4J7','context':cell['case_id'],
                'geometry':cell['candidate'],'medium':medium,'uhf':0,'reference':None,'scout':False,'reuse':None,
                **{k:cell[k] for k in ('xyz','charge','multiplicity','valence_electrons','atom_count','metal')}})
    m = {'protocol_id':c['protocol_id']+'_GFN2', 'implementation':c['factorial'],
         'runtime':c['runtime'], 'settings':c['settings'], 'tasks':tasks, 'config':pin(base.parent/'CONFIG.json')}
    write(directory/'manifest.json',m); runtime = read(c['runtime'])
    with ThreadPoolExecutor(max_workers=7) as pool:
        points = list(pool.map(lambda t:f.execute_point(t,m,directory,runtime),tasks))
    native = read(c['native_manifest']); out = directory/'replay'; out.mkdir(); rows = []
    env = f.clean_env(runtime['xtb_parameter_file']['path'])
    for cell in c['cells']:
        tid = cell['task_id']; record = {'task_id':tid,'status':'unavailable','reason':None,'G_qvac_hartree':None}
        try:
            v,w = [next(t for t in tasks if t['id'] == tid+'__'+medium) for medium in ('vacuum','water')]
            component = r.component(directory,v,w); record['component_parity'] = component
            assert component['accepted'], component.get('reason','component qualification failed')
            assert mr[tid]['status'] == 'complete', 'MACE unavailable'
            t = {k:cell[k] for k in ('xyz','atom_count','metal','charge','multiplicity','valence_electrons','state_signature','source_preparation','maps')}
            t.update(task_id=tid,scientific_input_key={'runtime':c['runtime']},
                uses=[{'MACE':mr[tid]['receipt'],'MACE_eV':mr[tid]['energy_eV']}],endpoints={})
            for medium in ('vacuum','water'):
                p = directory/'points'/(tid+'__'+medium); point = read(p/'RESULT.json')
                ep = p/'xtb'; execution = read(ep/'EXECUTION.json')
                t['endpoints'][medium] = {'point':pin(p/'RESULT.json'),'point_manifest':point['manifest'],
                    'task':point['task'],'execution':pin(ep/'EXECUTION.json'),'stdout':execution['stdout'],
                    'charges':pin(ep/'charges'),'executed_xyz':pin(ep/'core.xyz')}
            replay = n.run_task(native,verify(c['native_manifest']),t,c['native_binary'],out,env,None)
            record.update(replay)
        except Exception as exc: record['reason'] = type(exc).__name__+': '+str(exc)
        rows.append(record)
    write(directory/'COLLECTION.json',{'rows':rows,'raw_inputs':len(points),'maximum_backend_calls':2*len(tasks)})
    return rows


def summarize(c, mr, sr):
    solvent_rows = {r['task_id']:r for r in sr}; rows = []
    for source in c['rows']:
        matrix = {z:{} for z in METALS}; failures = []
        for z in METALS:
            for q in ORDER:
                tid = source['case_id']+'__'+z+'__'+q; mm = mr.get(tid); ss = solvent_rows.get(tid)
                if not mm or mm['status'] != 'complete' or not ss or ss['status'] != 'available':
                    matrix[z][q] = None; failures.append({'metal':z,'candidate':q,'MACE':mm,'solvent':ss}); continue
                matrix[z][q] = {'MACE_eV':mm['energy_eV'],'G_qvac_hartree':ss['G_qvac_hartree'],
                    'MACE':mm['receipt'],'solvent':ss['result']}
        row = {'case_id':source['case_id'],'status':'unavailable' if failures else 'available',
               'matrix':matrix,'failures':failures,'works':None,'selected':None,'classification':None,
               'classification_reason':'changed_context_has_no_calibration'}
        if not failures:
            works = {z:{} for z in METALS}
            for z in METALS:
                zero = matrix[z]['origin']
                for q in ORDER:
                    v = matrix[z][q]; dw = (v['MACE_eV']-zero['MACE_eV'])*EV
                    ds = (v['G_qvac_hartree']-zero['G_qvac_hartree'])*HA
                    works[z][q] = {'MACE_kcal_mol':dw,'solvent_kcal_mol':ds,'total_kcal_mol':dw+ds}
            selected = {}
            for policy in ('mathematical','operational'):
                selection = {}
                for z in METALS:
                    best = min(ORDER,key=lambda q:works[z][q]['total_kcal_mol'])
                    q = best if policy == 'mathematical' or works[z][best]['total_kcal_mol'] < -0.1 else 'origin'
                    selection[z] = {'candidate':q,'work':works[z][q]}
                selection['differential_work_kcal_mol'] = selection['Ca']['work']['total_kcal_mol']-selection['La']['work']['total_kcal_mol']
                selected[policy] = selection
            row.update(works=works,selected=selected,
                same_geometry_differential_works_kcal_mol={q:works['Ca'][q]['total_kcal_mol']-works['La'][q]['total_kcal_mol'] for q in ORDER})
        rows.append(row)
    return {'protocol_id':c['protocol_id'],'denominator':len(rows),'rows':rows,'production_changed':False,
            'new_calibration':False,'new_optimization':False,'new_DFT':False}


def run(config):
    place = placement(); c = read(config); assert c['implementation'] == pin(__file__)
    for key in ('plan','input','component','factorial','native_adapter'): verify(c[key])
    base = Path(config).parent/'execution'
    with (Path(config).parent/'EXECUTION.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB); base.mkdir(exist_ok=False)
        start = time.monotonic(); receipt = {'config':pin(config),'placement':place,'status':'started'}
        write(base/'START.json',receipt)
        try:
            sys.path.insert(0,str(ROOT/'scripts'))
            f = load('packing_factorial',c['factorial']); r = load('packing_component',c['component'])
            n = load('packing_native',c['native_adapter']); b = load('packing_cpu',read(c['worker_manifest'])['implementation'])
            write(base/'PREFLIGHT.json',preflight(c,f,b))
            mr = mace(c,base,b); sr = solvent(c,base,f,r,n,mr)
            result = summarize(c,mr,sr); result['config'] = pin(config)
            write(base/'FINAL_COLLECTION.json',result)
            receipt.update(status='complete',collection=pin(base/'FINAL_COLLECTION.json'))
        except BaseException as exc:
            receipt.update(status='failed',reason=str(exc),traceback=traceback.format_exc()); raise
        finally:
            receipt.update(elapsed_seconds=time.monotonic()-start,allocated_CPU_seconds=7*(time.monotonic()-start))
            write(base/'EXECUTION.json',receipt)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__); sub = p.add_subparsers(dest='action',required=True)
    a = sub.add_parser('prepare'); a.add_argument('--input',required=True); a.add_argument('--output',required=True)
    a = sub.add_parser('run'); a.add_argument('--config',required=True)
    args = p.parse_args()
    if args.action == 'prepare': prepare(args.input,args.output)
    else: run(args.config)
