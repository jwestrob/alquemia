#!/usr/bin/env python3
"""Two real crystal sources through the existing four-mode CPU proposal kernel."""
import argparse
import copy
import fcntl
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import signal
import sys
import time
import traceback

ROOT = next(p for p in Path(__file__).resolve().parents if (p/'scripts/mace_omol.py').is_file())
HERE = Path(__file__).resolve().parent
TEMPLATE = ROOT/'workspaces/ca_coverage_recovery_20261012/two_ca_scoring_q9_v1/lane0/worker_manifest.json'
HELPER = ROOT/'diagnostics/ca_coverage_recovery_20261012/disulfide/NATIVE.py'
PROTOCOL = 'Nikasha_ExaA_1FLG_two_chain_original_four_mode_CPU_v1'
ORDER = ('origin', 'adaptive_Ca', 'adaptive_La')
METALS = ('Ca', 'La')
VALENCE = {'H':1, 'C':4, 'N':5, 'O':6, 'S':6, 'Ca':2, 'La':3}


def pin(path):
    p = Path(path).resolve(); h = hashlib.sha256()
    with p.open('rb') as f:
        for chunk in iter(lambda:f.read(1048576), b''): h.update(chunk)
    return {'path':str(p), 'sha256':h.hexdigest()}


def verify(p):
    assert pin(p['path']) == p, ('changed pin', p)
    return Path(p['path'])


def read(p):
    return json.loads((verify(p) if isinstance(p,dict) else Path(p)).read_text())


def write(p, value):
    p = Path(p); p.parent.mkdir(parents=True, exist_ok=True)
    with p.open('x') as f:
        json.dump(value, f, indent=2, allow_nan=False); f.write('\n')


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def admission(row, path):
    """Only recorded, independently prepared chemistry can enter execution."""
    assert read(path) == row
    if row['status'] != 'prepared': return False
    assert row['mapping_eligible'] and row['paired_maps_equal']
    assert row['local_state_identical_to_Q9'] and row['construct_compatibility']['mature_subset_exact']
    spectator = read(row['spectator_admission'])
    assert spectator['passed'] and spectator['status'] == 'passed'
    if row.get('terminal_completion'):
        terminal = read(row['terminal_completion'])
        assert terminal['status'] == 'passed' and terminal['terminal_absent_from_context_and_caps']
    assert read(row['state_signature']) == read(row['reference_state_signature'])
    for k in ('source_preparation', 'protein_identity', 'raw_atom_map'): verify(row[k])
    for p in row['maps'].values(): verify(p)
    verify(row['representations']['context']['preparation'])
    for z in METALS: verify(row['representations']['context']['endpoints'][z]['xyz'])
    return True


def prepare(preparation, output):
    """Metadata only. Coordinate/kinematic and runtime checks run on Slurm."""
    pp = Path(preparation).resolve(); prepared = read(pp); old = read(TEMPLATE)
    assert prepared['terminal'] and prepared['denominator'] == 2 and len(prepared['rows']) == 2
    assert [r['chain'] for r in prepared['rows']] == ['A','B']
    assert len({r['case_id'] for r in prepared['rows']}) == 2
    base = Path(output).resolve(); base.mkdir(parents=True, exist_ok=False)
    declared = []; tasks = []; cases = []
    for row in prepared['rows']:
        ap = pin(pp.parent/('chain_'+row['chain'])/'ADMISSION.json')
        eligible = admission(row, ap)
        declared.append({'chain':row['chain'], 'source_case_id':row['case_id'], 'admission':ap,
                         'eligible':eligible, 'preparation_status':row['status'], 'reason':row.get('reason')})
        if not eligible: continue
        cid = row['case_id']+'__'+row['selection_id'].rsplit('__',1)[1]
        for z in METALS:
            ep = row['representations']['context']['endpoints'][z]
            tasks.append({'task_id':cid+'__context__'+z, 'case_id':cid, 'source_case_id':row['case_id'],
                'selection_id':row['selection_id'], 'metal':z,
                **{k:ep[k] for k in ('xyz','charge','multiplicity')}, 'mapping':row['maps'][z],
                'state_signature':row['state_signature'], 'context_preparation':row['representations']['context']['preparation'],
                'prepared_admission':ap})
        cases.append({'case_id':cid, 'source_case_id':row['case_id'], 'chain':row['chain'], 'lane':0,
            'selection_id':row['selection_id'], 'prepared_admission':ap, 'state_signature':row['state_signature'],
            'source_preparation':row['representations']['context']['preparation']})
    fixtures = [next(t for t in old['tasks'] if t['task_id'] == f['task_id']) for f in old['fixtures']]
    assert {t['case_id'] for t in fixtures}.isdisjoint(c['case_id'] for c in cases)
    worker = copy.deepcopy(old)
    for k in ('scientific_scope','origin_inputs'): worker.pop(k, None)
    worker.update(protocol_id=PROTOCOL, parent_manifest=pin(TEMPLATE), preparation=pin(pp),
        plan=pin(HERE/'PLAN.md'), tasks=tasks+fixtures, origin_task_count=len(tasks), cases=cases,
        lanes=1, maximum_active_numerical_threads=7,
        finite_scope={'declared_sources':2, 'admitted_sources':len(cases), 'origin_calls':len(tasks),
                      'optimizer_starts':len(tasks), 'maximum_cross_calls':len(tasks), 'qualification_calls':2,
                      'new_GFN2_calls':0, 'new_DFT_calls':0},
        authorization='Jacob authorized finite 1FLG A/B current-method source test; actual Slurm only')
    write(base/'worker_manifest.json', worker)
    config = {'protocol':PROTOCOL, 'implementation':pin(__file__), 'plan':pin(HERE/'PLAN.md'),
        'wrapper':pin(HERE/'RUN.sbatch'), 'helper':pin(HELPER), 'template':pin(TEMPLATE),
        'preparation':pin(pp), 'worker_manifest':pin(base/'worker_manifest.json'), 'declared':declared,
        'scope':worker['finite_scope'], 'state_interpretation':prepared.get('conditional_state'),
        'new_model':False, 'new_calibration':False, 'production_changed':False}
    write(base/'CONFIG.json', config)
    write(base/'READY.json', {'status':'metadata_prepared_scientific_preflight_UNRUN' if cases else 'blocked_no_admitted_source',
                            'config':pin(base/'CONFIG.json'), 'new_molecular_calls':0})
    print(json.dumps({'config':pin(base/'CONFIG.json'), 'admitted_sources':len(cases), 'declared_sources':2}))


def runtime_preflight(c, m, b):
    old = read(c['template']); verify(c['preparation']); verify(c['plan'])
    for k in ('implementation','source_modules','source_snapshot','model','settings','pool_settings',
              'optimizer_software','optimizer_executable','model_executable','optimizer_python',
              'model_python','calculator_source','fixtures','gates','threads_per_model'):
        assert m[k] == old[k], 'original scientific/runtime field differs: '+k
    for k in ('implementation','calculator_source','optimizer_executable','model_executable'): verify(m[k])
    for p in m['source_modules'].values(): verify(p)
    verify(m['model']['checkpoint']); assert pin(sys.executable) == m['optimizer_executable']
    assert m['lanes'] == 1 and m['maximum_active_numerical_threads'] == 7
    assert 1 <= len(m['cases']) <= 2 and m['origin_task_count'] == 2*len(m['cases'])
    assert len(m['tasks']) == m['origin_task_count']+2
    b.modules(m)
    import scipy
    import scipy.optimize._slsqplib as kernel
    import scipy.optimize._slsqp_py as wrapper
    import slsqp_precision
    assert scipy.__version__ == m['optimizer_software']['version']
    assert pin(kernel.__file__) == m['optimizer_software']['kernel'] and pin(wrapper.__file__) == m['optimizer_software']['wrapper']
    assert slsqp_precision.engine()[0].completion.SETTINGS == m['settings']
    for fixture in m['fixtures']:
        for k in ('xyz','forces','native_receipt'): verify(fixture[k])
    return {'status':'passed', 'method_changes':[], 'new_model_calls':0}


def source_preflight(m, case, b):
    import numpy as np
    from mace_site_kinematics import Kinematics
    from accommodation_torsion_profiles import geometry_check
    from affordable_common import paired
    prepared = read(case['prepared_admission']); assert admission(prepared, case['prepared_admission'])
    ts = [t for t in m['tasks'][:m['origin_task_count']] if t['case_id'] == case['case_id']]
    assert [t['metal'] for t in ts] == list(METALS)
    state = read(prepared['state_signature']); checks = []
    for t in ts:
        ep = prepared['representations']['context']['endpoints'][t['metal']]
        assert all(t[k] == ep[k] for k in ('xyz','charge','multiplicity'))
        assert t['mapping'] == prepared['maps'][t['metal']] and t['prepared_admission'] == case['prepared_admission']
        assert t['state_signature'] == prepared['state_signature'] and t['charge'] == state['endpoint_charges'][t['metal']]
        for k in ('mapping','state_signature','context_preparation'): verify(t[k])
        a = b.xyz(verify(t['xyz'])); kin = Kinematics(read(t['mapping'])['context']); q = np.zeros(len(kin.modes))
        residual = float(np.max(np.abs(kin.evaluate(q)[1]-np.asarray([v[1:] for v in a]))))
        assert residual <= 1e-12 and geometry_check(kin,q,[v[0] for v in a])['pass']
        assert t['multiplicity'] == 1 and (sum(VALENCE[v[0]] for v in a)-t['charge']) % 2 == 0
        checks.append({'metal':t['metal'], 'q0_max_coordinate_error_A':residual, 'physical_map_pass':True})
    paired(verify(ts[1]['xyz']),verify(ts[0]['xyz']),ts[1]['charge'],ts[0]['charge'])
    assert read(ts[0]['mapping']) == read(ts[1]['mapping'])
    return {'case_id':case['case_id'], 'passed':True, 'checks':checks}


def export(c, base, output):
    """Read actual native pools; never evaluate or manufacture an absent endpoint."""
    m = read(c['worker_manifest']); rows = []
    for declared in c['declared']:
        prep = read(declared['admission']); case = next((x for x in m['cases'] if x['source_case_id'] == declared['source_case_id']), None)
        result_path = base/'sources'/case['case_id']/'SOURCE_RESULT.json' if case else None
        result = read(result_path) if result_path and result_path.exists() else None
        if result: assert result['batch_manifest'] == c['worker_manifest'] and result['source_case_id'] == declared['source_case_id']
        row = {'case_id':declared['source_case_id'], 'chain':declared['chain'], 'status':'unavailable',
               'admission':declared['admission'], 'source':prep.get('source'), 'reason':declared.get('reason'),
               'native_source_result':pin(result_path) if result else None,
               'matrix':{z:{} for z in METALS}, 'failures':[], 'conditional_source_state':True}
        if case:
            row.update(selection_id=case['selection_id'], maps=prep['maps'], state_signature=prep['state_signature'],
                       source_preparation=prep['representations']['context']['preparation'],
                       reference_state_signature=prep['reference_state_signature'],
                       terminal_completion=prep.get('terminal_completion'))
        for z in METALS:
            for q in ORDER:
                cell = result['matrix'][z].get(q) if result else None
                if not cell or cell['status'] != 'complete':
                    row['matrix'][z][q] = None
                    row['failures'].append({'metal':z,'candidate':q,'reason':(cell or {}).get('reason',declared.get('reason') or 'missing_terminal_native_cell')})
                    continue
                mm = read(cell['MACE']); request = read(mm['request'])
                assert mm['status'] == 'complete' and mm['model'] == m['model'] and mm['energy_eV'] == cell['energy_eV']
                assert request['xyz'] == cell['xyz'] and request['charge'] == cell['charge'] and request['multiplicity'] == cell['multiplicity']
                lines = verify(cell['xyz']).read_text().splitlines(); symbols = [line.split()[0] for line in lines[2:] if line.strip()]
                assert int(lines[0]) == len(symbols) and symbols[0] == z
                electrons = sum(VALENCE[s] for s in symbols)-cell['charge']
                assert cell['multiplicity'] == 1 and electrons % 2 == 0
                verify(cell['forces'])
                row['matrix'][z][q] = {k:cell[k] for k in ('xyz','charge','multiplicity','MACE','forces')}
                row['matrix'][z][q].update(status='MACE_complete_solvent_uncomputed', MACE_eV=cell['energy_eV'],
                    atom_count=len(symbols), valence_electrons=electrons, native_alias_of=cell.get('alias_of'),
                    original_MACE_request=mm['request'])
        if not row['failures']: row.update(status='prepared', reason=None)
        rows.append(row)
    value = {'protocol':PROTOCOL+'_existing_MACE_pool_export', 'status':'prepared' if all(r['status']=='prepared' for r in rows) else 'partial_or_unavailable',
        'rows':rows, 'denominator':2, 'config':pin(base/'CONFIG.json'), 'worker_manifest':c['worker_manifest'],
        'new_MACE_calls_in_export':0, 'new_GFN2_calls':0, 'new_calibration':False, 'production_changed':False,
        'intended_use':'Use these actual MACE receipts in frozen-charge solvent scoring; do not rerun MACE.'}
    write(output, value); return value


def execute(path):
    assert os.environ.get('SLURM_JOB_ID'), 'Actual Slurm required before scientific imports'
    c = read(path); assert c['implementation'] == pin(__file__)
    verify(c['wrapper']); helper = load('crystal_allocation_helper',verify(c['helper'])); placement = helper.allocation()
    mp = verify(c['worker_manifest']); m = read(mp); base = Path(path).resolve().parent
    with (base/'RUN.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX|fcntl.LOCK_NB)
        attempt = base/'attempts'/'001'; attempt.mkdir(parents=True,exist_ok=False); (attempt/'lane0').mkdir()
        write(attempt/'PLACEMENT.json',placement)
        b = load('crystal_existing_cpu',verify(m['implementation']))
        handlers = b.handlers(); transport = None; error = None; started = time.monotonic(); preflights = []
        try:
            write(attempt/'RUNTIME_PREFLIGHT.json',runtime_preflight(c,m,b))
            admissible = []
            for case in m['cases']:
                try: preflights.append(source_preflight(m,case,b)); admissible.append(case)
                except Exception as exc: preflights.append({'case_id':case['case_id'],'passed':False,'reason':type(exc).__name__+': '+str(exc)})
            write(attempt/'SOURCE_PREFLIGHT.json',{'rows':preflights,'new_model_calls':0})
            assert admissible, 'no source passed the unchanged physical preflight'
            transport = b.WarmCPU(mp,m,attempt,0)
            import numpy as np
            checks = []
            for fixture in m['fixtures']:
                task = next(t for t in m['tasks'] if t['task_id'] == fixture['task_id'])
                request = {'manifest':pin(mp),'task_id':task['task_id'],'xyz':fixture['xyz'],
                           'charge':task['charge'],'multiplicity':task['multiplicity']}
                r,rp = b.cell(attempt/'lane0'/'qualification'/fixture['id'],request,transport,m['model'])
                check = {'fixture':fixture['id'],'receipt':rp,'pass':False}
                if r['status'] == 'complete':
                    de = (r['energy_eV']-fixture['energy_eV'])*b.EV_TO_KCAL
                    df = float(np.max(np.abs(np.load(verify(r['forces']),allow_pickle=False)-np.load(verify(fixture['forces']),allow_pickle=False))))
                    check.update(energy_error_kcal_mol=de, force_max_error_eV_A=df,
                                 **{'pass':abs(de)<=m['gates']['energy_kcal_mol'] and df<=m['gates']['force_max_eV_A']})
                checks.append(check)
            write(attempt/'QUALIFICATION.json',{'checks':checks,'pass':all(x['pass'] for x in checks)})
            assert all(x['pass'] for x in checks), 'unchanged CPU force qualification failed'
            for case in admissible: b.run_source(mp,m,case,transport)
        except BaseException as exc:
            error = type(exc).__name__+': '+str(exc); traceback.print_exc()
        finally:
            if transport:
                try:
                    if error: transport.abort()
                    else: transport.close()
                except BaseException as exc:
                    error = (error+'; ' if error else '')+str(exc); transport.abort()
            for sig,handler in handlers.items(): signal.signal(sig,handler)
            exported = None
            try: export(c,base,base/'MATRIX_INPUT.json'); exported = pin(base/'MATRIX_INPUT.json')
            except BaseException as exc: error = (error+'; ' if error else '')+'export: '+type(exc).__name__+': '+str(exc)
            calls = [read(p) for p in sorted(base.glob('**/mace_result.json'))]
            elapsed = time.monotonic()-started
            write(base/'FINAL_COLLECTION.json', {'protocol':PROTOCOL,'config':pin(path),'worker_manifest':c['worker_manifest'],
                'terminal':True,'denominator':2,'error':error,'source_preflights':preflights,'matrix_input':exported,
                'costs':{'elapsed_seconds_through_export':elapsed,'allocated_CPU_seconds_through_export':7*elapsed,
                    'new_recorded_model_calls':len(calls),'actual_forwards':sum(r['actual_model_forwards'] for r in calls),
                    'recorded_process_CPU_seconds':sum(r['process_cpu_seconds'] for r in calls),'placement':placement,
                    'started_without_terminal_receipt':[str(p) for p in base.glob('**/CPU_CALL_STARTED.json') if not (p.parent/'mace_result.json').exists()]},
                'new_GFN2_calls':0,'new_DFT_calls':0,'new_calibration':False,'production_changed':False})
        if error: raise RuntimeError(error)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__); s = p.add_subparsers(dest='command',required=True)
    a = s.add_parser('prepare'); a.add_argument('--preparation',required=True); a.add_argument('--output',required=True)
    a = s.add_parser('execute'); a.add_argument('--config',required=True)
    a = s.add_parser('export'); a.add_argument('--config',required=True); a.add_argument('--output',required=True)
    a = p.parse_args()
    if a.command == 'prepare': prepare(a.preparation,a.output)
    elif a.command == 'execute': execute(a.config)
    else:
        assert os.environ.get('SLURM_JOB_ID'), 'actual Slurm required for saved-source export'
        c = read(a.config); assert c['implementation'] == pin(__file__)
        base = Path(a.config).resolve().parent
        with (base/'RUN.lock').open('a') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
            export(c,base,a.output)
