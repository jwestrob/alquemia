#!/usr/bin/env python3
"""Saved-only diagnosis of the frozen common7 raw POINT replay failure."""
import argparse
import copy
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys
import time
import traceback

PROTOCOL = 'mex_response_common7_point_pin_diagnosis_v1'
ADAPTER_SHA = 'bd29908d44cd975eb37b220ec583706c32346fad62272240d0b6ee28795bb608'


def require(ok, message):
    if not ok:
        raise ValueError(message)


def pin(path):
    p = Path(path).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def read(record):
    require(pin(record['path']) == record, 'Changed pinned file: '+str(record))
    return json.loads(Path(record['path']).read_text())


def differences(left, right, path='$'):
    """Exact typed JSON-leaf comparison; no numerical tolerance or normalization."""
    if type(left) is not type(right):
        return [{'path': path, 'archived': left, 'replayed': right, 'kind': 'type'}]
    if isinstance(left, dict):
        rows = []
        for key in sorted(set(left) | set(right)):
            if key not in left or key not in right:
                rows.append({'path': path+'.'+key, 'archived_present': key in left,
                             'replayed_present': key in right, 'kind': 'key'})
            else:
                rows.extend(differences(left[key], right[key], path+'.'+key))
        return rows
    if isinstance(left, list):
        if len(left) != len(right):
            return [{'path': path, 'archived_length': len(left), 'replayed_length': len(right), 'kind': 'length'}]
        return [row for i, (a, b) in enumerate(zip(left, right))
                for row in differences(a, b, path+'['+str(i)+']')]
    return [] if left == right else [{'path': path, 'archived': left, 'replayed': right, 'kind': 'value'}]


def forbidden(*args, **kwargs):
    raise RuntimeError('Diagnostic forbids evaluator calls, including cache requests')


def diagnose(cfg, out, result):
    require(cfg['protocol'] == PROTOCOL and pin(__file__)['sha256'] == cfg['implementation']['sha256'],
            'Diagnostic source/protocol differs')
    for name in ('implementation', 'specification'):
        require(pin(cfg[name]['path']) == cfg[name], 'Changed diagnostic artifact: '+name)
    source_cfg = read(cfg['source_config']); failure = read(cfg['failed_result'])
    require(failure['status'] == 'execution_failure' and failure['slurm_job_id'] == '13063'
            and failure['failure'] == "ValueError('Common7 raw cell replay differs')"
            and failure['config'] == cfg['source_config']
            and source_cfg['implementation'] == failure['implementation']
            and source_cfg['implementation']['sha256'] == ADAPTER_SHA,
            'Exact preserved preparation failure required')
    require(pin(source_cfg['implementation']['path']) == source_cfg['implementation'], 'Changed frozen adapter')
    spec = importlib.util.spec_from_file_location('frozen_mex_response_reference', source_cfg['implementation']['path'])
    ref = importlib.util.module_from_spec(spec); spec.loader.exec_module(ref)
    inv, _, _, _, _, _ = ref.inputs(source_cfg)
    old = ref.read(inv['pool_result']); pcfg = ref.read(old['config'])
    require(old['status'] == 'complete' and old['slurm_job_id'] == '12981'
            and old['declared_cells'] == 14 and old['declared_geometries'] == 7, 'Actual common7 required')
    pool = ref.module('mex_response_saved_pool', old['implementation'])
    require(old['protocol_id'] == pool.PROTOCOL and pcfg['implementation']['sha256'] == old['implementation']['sha256'],
            'Pool implementation/protocol differs')
    ecfg = ref.read(pcfg['evaluator_config'])
    em = ref.module('mex_response_saved_evaluator', ecfg['implementation'])
    sm = ref.module('mex_response_saved_search', pcfg['search_implementation'])
    opt = ref.module('mex_response_saved_optimizer', pcfg['optimizer'])
    ev = em.MexEvaluator(ecfg, out/'saved_evaluator')
    result['_ev'] = ev
    ev.evaluate = forbidden; ev.old.evaluate = forbidden
    guard, _ = sm.build_guard(ev, opt, out)
    oldpool = ref.module('mex_response_old_common5', ev.cfg['pool_implementation'])
    matrix, by_coordinate, _, baseline_pin = pool.baseline(ev, oldpool, guard, out)
    slots, failures = pool.admit_searches(pcfg, ev, sm, opt, guard, matrix, by_coordinate)
    require(all(s['status'] == 'admitted' for s in slots), 'Both actual searches must be admitted')
    for actual in old['matrix']:
        for metal, cell in actual['cells'].items():
            if cell.get('cartesian_request'):
                pool.import_request(ev, metal, cell['cartesian_request'])
    result.update(baseline_replay=baseline_pin, terminal_slots=slots, failed_exact_requests=failures)
    for row in matrix:
        actual = next(a for a in old['matrix'] if a['coordinate_sha256'] == row['coordinate_sha256'])
        require(row['aliases'] == actual['aliases'], 'Common7 candidate identity differs')
        for metal in ref.METALS:
            cached = pool.cached_cell(ev, metal, row['_x']); archived = actual['cells'][metal]
            a = ref.read(archived['point']); b = ref.read(cached['point'])
            physical_a = ref.read(a['physical']); physical_b = ref.read(b['physical'])
            gradient_a = ref.read(a['assembled_gradient']); gradient_b = ref.read(b['assembled_gradient'])
            candidate = copy.deepcopy(b); candidate['physical'] = a['physical']
            exact = cached['point'] == archived['point'] and cached['native_energies'] == archived['native_energies']
            result['rows'].append({'aliases': row['aliases'], 'metal': metal,
                'coordinate_sha256': row['coordinate_sha256'], 'cache_key': ev.key(metal, row['_x']),
                'cache_record': ev.old.cache[ev.key(metal, row['_x'])],
                'archived_cell': archived, 'replayed_cell': cached, 'frozen_raw_cell_gate_passes': exact,
                'POINT_differences': differences(a, b),
                'native_energies_exact': cached['native_energies'] == archived['native_energies'],
                'physical_byte_hash_exact': a['physical']['sha256'] == b['physical']['sha256'],
                'physical_JSON_exact': physical_a == physical_b,
                'physical_ids_exact': [x['id'] for x in physical_a] == [x['id'] for x in physical_b] == ev.ids,
                'assembled_gradient_pin_exact': a['assembled_gradient'] == b['assembled_gradient'],
                'assembled_gradient_JSON_exact': gradient_a == gradient_b,
                'POINT_exact_if_only_physical_pin_substituted': candidate == a,
                'substitution_is_diagnostic_only': True})
    require(len(result['rows']) == 14 and len(matrix) == 7, 'Full saved 7/14 denominator differs')
    require(ev.new_evaluations == 0 and all(v == 0 for v in ev.counters.values()) and ev.old.native is None,
            'Saved diagnosis must not initialize a model or evaluate components')
    import numpy, gemmi
    result.update(imported_environment={'python': sys.executable, 'numpy_version': numpy.__version__,
        'numpy_file': numpy.__file__, 'gemmi_version': gemmi.__version__, 'gemmi_file': gemmi.__file__},
        mismatch_slots=[{'aliases': r['aliases'], 'metal': r['metal']} for r in result['rows']
                        if not r['frozen_raw_cell_gate_passes']],
        diagnostic_complete=True, replay_admitted=False)
    result['status'] = 'diagnostic_mismatch' if result['mismatch_slots'] else 'no_mismatch_observed'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True); parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(); require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Slurm required')
    out = args.output.resolve(); out.mkdir(parents=True, exist_ok=False); started = time.monotonic()
    cfg = json.loads(args.config.read_text())
    result = {'protocol': PROTOCOL, 'status': 'execution_failure', 'implementation': pin(__file__),
        'config': pin(args.config), 'slurm_job_id': os.environ['SLURM_JOB_ID'], 'declared_cells': 14,
        'rows': [], 'diagnostic_complete': False, 'replay_admitted': False, 'new_molecular_calls': 0}
    try:
        diagnose(cfg, out, result)
        require(pin(args.config) == result['config'], 'Diagnostic config changed')
    except Exception as exc:
        result.update(status='execution_failure', failure=repr(exc), traceback=traceback.format_exc())
    ev = result.pop('_ev', None)
    if ev is not None:
        result.update(counters=ev.counters, model_initialized=ev.old.native is not None,
                      new_composite_evaluations=ev.new_evaluations, requests=ev.requests)
    result['elapsed_seconds'] = time.monotonic()-started
    ref = out/'RESULT.json'
    with ref.open('x') as handle:
        json.dump(result, handle, indent=2, allow_nan=False); handle.write('\n')
    print(json.dumps({'status': result['status'], 'result': pin(ref)}))
    return 0 if result['status'] == 'no_mismatch_observed' else 1


if __name__ == '__main__':
    raise SystemExit(main())
