#!/usr/bin/env python3
"""Slurm-only negative controls for the reviewed response adapter; no engines."""
import argparse
import contextlib
import copy
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path
import sys
import time
import traceback
import types

ADAPTER_SHA = '3acf08878f5d9b11275145ce65355e34ad01d792118c1635dda437db09ee1324'
TEMPLATE_SHA = 'a35bde861f7a22cb575e76c1670ec4fa9cf8c6030cbafa1c9bf79fbd132d9c5a'
PROTOCOL = 'mex_response_failure_gate_regression_v2'
COUNTERS = ('new_DFT_calls_requested', 'new_DFT_calls_attempted', 'new_MM_calls',
            'new_ALPB_calls', 'new_MACE_calls', 'new_GFN_calls', 'optimizer_steps')


def require(ok, message):
    if not ok:
        raise ValueError(message)


def pin(path):
    p = Path(path).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def read(record):
    require(pin(record['path']) == record, 'Changed artifact: '+str(record))
    return json.loads(Path(record['path']).read_text())


def write(path, value):
    with Path(path).open('x') as f:
        json.dump(value, f, indent=2, allow_nan=False)
        f.write('\n')
    return pin(path)


def load(name, record):
    require(pin(record['path']) == record, 'Changed source')
    spec = importlib.util.spec_from_file_location(name, record['path'])
    obj = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


def private_functions(mod):
    """Exact code objects in a test-local namespace; original globals are untouched."""
    ns = dict(vars(mod))
    for name, obj in vars(mod).items():
        if isinstance(obj, types.FunctionType) and obj.__module__ == mod.__name__:
            cloned = types.FunctionType(obj.__code__, ns, name, obj.__defaults__, obj.__closure__)
            cloned.__kwdefaults__ = obj.__kwdefaults__
            ns[name] = cloned
    return ns


def run(args, out, result):
    template_pin = pin(args.template)
    require(template_pin['sha256'] == TEMPLATE_SHA, 'Reviewed template differs')
    template = read(template_pin)
    prep_pin, flight_pin = pin(args.preparation_result), pin(args.preflight_result)
    prep, flight = read(prep_pin), read(flight_pin)
    cfg = read(flight['config'])
    require(cfg['implementation']['sha256'] == template['implementation']['sha256'] == ADAPTER_SHA,
            'Reviewed adapter differs; a future repair needs a separately reviewed regression revision')
    mod = load('mex_response_under_test', cfg['implementation'])
    require(prep['status'] == 'prepared' and prep['stage'] == 'prepare'
            and flight['status'] == 'complete' and flight['stage'] == 'preflight'
            and prep['input_admission_passed'] and flight['input_admission_passed']
            and prep['protocol'] == flight['protocol'] == mod.PROTOCOL
            and flight['preparation_result'] == cfg['preparation_result'] == prep_pin
            and all(prep[k] == flight[k] == 0 for k in COUNTERS),
            'Successful real zero-call preparation and preflight are mandatory')
    require(mod.binding(read(prep['config'])) == mod.binding(cfg), 'Actual preparation binding differs')
    for key in mod.binding(template):
        if key in ('implementation', 'specification'):
            require(pin(cfg[key]['path']) == cfg[key]
                    and cfg[key]['sha256'] == template[key]['sha256'], 'Source/specification differs')
        else:
            require(cfg[key] == template[key], 'Frozen scientific binding differs: '+key)
    qc = read(cfg['recipe_config'])
    require(pin(sys.executable)['sha256'] == qc['python']['sha256'], 'Use the pinned recipe Python')
    inv = read(cfg['inventory']); full = read(prep['full_manifest'])
    tasks = {read(p)['task_id']: read(p) for p in prep['tasks']}
    require(tuple(tasks) == mod.IDS and list(tasks.values()) == full['tasks'], 'Real eight-task fixture differs')
    shared = {tid: read(prep['shared_components'][tid]) for tid in mod.IDS}
    replay = read(flight['source_replay'])
    qualifier = load('mex_response_qualifier_under_test', cfg['qualifier'])
    trips = []

    def forbidden(*unused, **unused_kw):
        trips.append('forbidden engine or full qualification replay reached')
        raise AssertionError(trips[-1])

    # Only unavailable/incompatible evidence is supplied. Real admission must
    # reject it before collector/runner entry; these tripwires never calculate.
    qstub = types.SimpleNamespace(PROTOCOL=qualifier.PROTOCOL, IDS=qualifier.IDS,
                                  SETTINGS=qualifier.SETTINGS, collect=forbidden)
    context = (inv, qc, qstub, None, None,
               {'run_orca_task_manifest': types.SimpleNamespace(run_manifest=forbidden)})
    cfg.update(preflight_result=flight_pin, qualification_result=None)
    for key in ('first_scout_decision', 'response_scout_result', 'remaining_decision'):
        cfg[key] = None
    cfg['endpoint_results'] = {tid: None for tid in mod.IDS}
    mod.preflight_gate(cfg)  # The real, unchanged admission of the supplied fixture.
    result.update(template=template_pin, preparation_result=prep_pin, preflight_result=flight_pin,
                  adapter=cfg['implementation'], fixture_config=flight['config'], cases=[])

    def case(name, stage, qualification, expected_status, message=None, before=None, fault=None):
        d = out/name; d.mkdir()
        c = copy.deepcopy(cfg); c['qualification_result'] = qualification
        cp = write(d/'CONFIG.json', c)
        ns = private_functions(mod)
        # These four hooks replay already admitted preparation data. They do not
        # replace qualification_gate, preflight_gate, endpoint, main or force_report.
        ns['inputs'] = lambda unused: context
        ns['prepared'] = lambda unused, unused_context: (prep, full, tasks)
        ns['source_replay'] = lambda unused, unused_out, unused_diagnosis: (None, copy.deepcopy(replay))
        ns['shared_cell'] = lambda cell, *unused: copy.deepcopy(shared[cell['task_id']])
        if fault is not None:
            ns['shared_cell'] = fault
        argv = sys.argv
        try:
            sys.argv = [cfg['implementation']['path'], '--stage', stage, '--config', cp['path'],
                        '--output', str(d/'execution')]
            if stage == 'endpoint':
                sys.argv += ['--task', mod.SCOUT]
            with contextlib.redirect_stdout(io.StringIO()):
                code = ns['main']()
        finally:
            sys.argv = argv
        rp = pin(d/'execution'/'RESULT.json'); r = read(rp)
        require(r['status'] == expected_status and code == (0 if expected_status == 'complete' else 1),
                'Unexpected status/exit code: '+name)
        require(all(r[k] == 0 for k in COUNTERS) and not trips, 'A forbidden call was reached: '+name)
        if message is not None:
            require(message in r['failure'], 'Wrong failure was exercised: '+name)
        if before is not None:
            require(r['status_before_failure'] == before, 'Exception status did not replace provisional success')
        result['cases'].append({'name': name, 'passed': True, 'result': rp,
            'expected_status': expected_status, 'observed_status': r['status'],
            'status_before_failure': r.get('status_before_failure'), 'failure': r.get('failure'),
            'calls': {k:r[k] for k in COUNTERS}})

    case('saved_preflight_control', 'preflight', None, 'complete')
    case('missing_qualification', 'endpoint', None, 'execution_failure',
         'Missing/changed pinned artifact: None')
    case('real_preparation_is_not_qualification', 'endpoint', prep_pin, 'execution_failure',
         'Actual full derivative qualification required')
    case('preflight_late_gate_failure', 'preflight', prep_pin, 'execution_failure',
         'Actual full derivative qualification required', before='complete')

    # The finite-array regression uses the real, admitted origin assembly and
    # its same-geometry approximate array, not fabricated molecular observations.
    # A single in-memory NaN is explicitly a test fault, never a saved observation.
    import numpy as np
    origin_result = read(qc['origin_result'])
    origin = next(r for r in origin_result['rows'] if r['metal'] == 'La')
    require(origin_result['status'] == 'complete' and origin['status'] == 'complete', 'Real origin fixture unavailable')
    approx_pin = origin['reference_evidence']['assembled']; approx = read(approx_pin)
    require(approx['source_ids'] == origin['source_ids']
            and np.asarray(approx['gradient_kcal_mol_A']).shape == (2049, 3)
            and np.isfinite(approx['gradient_kcal_mol_A']).all(), 'Original approximate array is not finite')
    corrupted = copy.deepcopy(approx); corrupted['gradient_kcal_mol_A'][0][0] = float('nan')
    force_ns = private_functions(mod)
    force_ns['read'] = lambda p: copy.deepcopy(corrupted) if p == approx_pin else mod.read(p)
    dft = {'mapped_physical_gradient': {'gradient_kcal_mol_A': origin['component_gradients_kcal_mol_A']['PBE0_embedded']}}
    same_shared = {'source_ids': origin['source_ids'], 'approximate_assembled_gradient': approx_pin,
                  'component_gradients_kcal_mol_A': {k:origin['component_gradients_kcal_mol_A'][k]
                                                   for k in ('MM_ledger', 'fixed_FF_charge_ALPB')}}

    def nonfinite_fault(*unused):
        force_ns['force_report']({}, dft, same_shared, {'directions': {}}, np)
        raise AssertionError('Nonfinite approximate gradient was accepted')

    case('nonfinite_approximate_gradient', 'preflight', None, 'execution_failure',
         'Full physical force unavailable', fault=nonfinite_fault)
    result['nonfinite_fixture'] = {'origin_result': qc['origin_result'], 'approximate_assembled': approx_pin,
        'test_fault_only': 'gradient_kcal_mol_A[0][0] replaced in memory by NaN',
        'original_artifacts_modified': False, 'new_molecular_observation': False}
    require(pin(cfg['implementation']['path']) == cfg['implementation'] and not trips, 'Source changed or tripwire reached')
    result.update(status='passed', forbidden_entries=trips, molecular_calls=0,
        scope='Test-local saved-data hooks; real gate/main/force-report code; no preparation or engine rerun')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ('template', 'preparation-result', 'preflight-result', 'output'):
        p.add_argument('--'+key, type=Path, required=True)
    args = p.parse_args()
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'This regression runner requires Slurm')
    out = args.output.resolve(); out.mkdir(parents=True, exist_ok=False); start = time.monotonic()
    result = {'protocol': PROTOCOL, 'status': 'execution_failure', 'implementation': pin(__file__),
              'slurm_job_id': os.environ['SLURM_JOB_ID'], 'molecular_calls': 0}
    try:
        run(args, out, result)
    except Exception as exc:
        result.update(status='execution_failure', failure=repr(exc), traceback=traceback.format_exc())
    result['elapsed_seconds'] = time.monotonic()-start
    write(out/'RESULT.json', result)
    print(json.dumps({'status': result['status'], 'result': pin(out/'RESULT.json')}))
    return 0 if result['status'] == 'passed' else 1


if __name__ == '__main__':
    raise SystemExit(main())
