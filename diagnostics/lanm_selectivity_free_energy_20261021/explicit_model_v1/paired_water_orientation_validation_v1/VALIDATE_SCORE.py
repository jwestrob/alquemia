#!/usr/bin/env python3
"""New subset/missing-data controls, using only metadata and synthetic outcomes."""
import copy
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]

def pin(p):
    p = Path(p).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}

def module(p, name):
    s = importlib.util.spec_from_file_location(name, p)
    m = importlib.util.module_from_spec(s)
    s.loader.exec_module(m)
    return m

def metadata():
    stage_path = HERE.parent / 'paired_donor_validation_stageA_v1/SCORE.py'
    body_path = HERE.parent / 'matched_dry_v1/OBJECTIVE.py'
    work = ROOT / 'workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1'
    geometry_path = work / 'matched_dry_v1/reserved_validation_v1/RESERVED_PANEL.json'
    training_path = work / 'matched_dry_v1/training_v1/RESULT.json'
    assert pin(stage_path)['sha256'] == '1c4bb0fa4b899991d3ca064480a10f199f2e8668469146389a1be403c257f104'
    assert pin(body_path)['sha256'] == 'bf65a20a379811f8d35e2cb5ae78365a5d34022691b332bae5df223196078819'
    assert pin(geometry_path)['sha256'] == 'eed99b0fb6dbe7b49857d68853fe809481f6b657c928a7e75197d77c61f3432b'
    assert pin(training_path)['sha256'] == '99ea27f86871b505dc622c1391a6321381f7430d3ad4a86a72e34ab9a629f6bb'
    stage = module(stage_path, 'fixed_stage')
    body = module(body_path, 'fixed_body')
    panel, old = (json.loads(p.read_text()) for p in (geometry_path, training_path))
    rows = {}
    for g in panel['geometries']:
        gid = g['geometry_id']
        if gid not in stage.CELLS:
            continue
        family = next(f for f in stage.FAMILIES if gid.startswith('validation_' + f + '_orientation_'))
        parent = family + '_central_' + g['metal']
        row = {'geometry': copy.deepcopy(g), 'input': {'bonds_zero_based': copy.deepcopy(old['rows'][parent]['input']['bonds_zero_based'])}, 'reference': None}
        row['bodies'] = body.bodies(row)
        rows[gid] = row
    return {'rows': rows}, stage, body

def main(out):
    start, cpu = time.monotonic(), time.process_time()
    assert all(os.environ.get(k) == '1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'))
    data, stage, body = metadata()
    score = module(HERE / 'SCORE_WATER.py', 'water_score')
    model = {gid: None for gid in stage.CELLS}
    water = [gid for gid in stage.CELLS if gid.startswith('validation_aquo8_')]
    for gid in water:
        value = {'energy': 0., 'forces': [[0., 0., 0.] for _ in range(25)]}
        model[gid] = value
        data['rows'][gid]['reference'] = copy.deepcopy(value)
    anchors = {k: {m: {'energy': 0.} for m in stage.METALS} for k in ('model', 'reference')}
    checks = {}
    def check(k, v):
        checks[k] = bool(v)
        assert v, k
    def run(d=data, m=model, a=anchors):
        return score.evaluate(d, m, body, a, stage)
    z = run()
    check('water_zero_passes', z['water_subset_passed'])
    check('full44_retained_partial', z['full_stage_A']['status'] == 'partial' and not z['full_stage_A']['validation_passed'] and sum(v['declared'] for v in z['full_stage_A']['counts'].values()) == 44)
    check('eleven_available', sum(v['available'] for v in z['counts'].values()) == 11)
    check('no_sampling_release', not z['protein_sampling_released'])
    for side in ('model', 'reference'):
        for gid in water:
            d, m = copy.deepcopy(data), copy.deepcopy(model)
            if side == 'model': m[gid] = None
            else: d['rows'][gid]['reference'] = None
            r = run(d=d, m=m)
            check(side + gid + '_missing_no_pass', not r['water_subset_passed'] and r['status'] == 'water_subset_partial' and not r['water_failure_observed'])
    m = copy.deepcopy(model)
    m[water[0]]['energy'] = 2.
    check('single_energy_failure_retained', run(m=m)['water_failure_observed'] and not run(m=m)['water_subset_passed'])
    m = copy.deepcopy(model)
    m[water[0]]['forces'][0][0] = 10.
    check('single_force_failure_retained', run(m=m)['water_failure_observed'] and not run(m=m)['water_subset_passed'])
    a = copy.deepcopy(anchors); a['reference']['Dy'] = None
    r = run(a=a)
    check('missing_anchor_cannot_pass', not r['water_subset_passed'] and r['counts']['balanced_exchange']['available'] == 0)
    donor = next(gid for gid in stage.CELLS if gid not in water)
    for side in ('model', 'reference'):
        d, m = copy.deepcopy(data), copy.deepcopy(model)
        if side == 'model': m[donor] = {'status': 'unavailable'}
        else: d['rows'][donor]['reference'] = {'status': 'unavailable'}
        try: run(d=d, m=m)
        except ValueError: check(side + '_donor_record_rejected', True)
        else: raise AssertionError('Donor outcome accepted')
    result = {'status': 'passed', 'source': pin(HERE / 'SCORE_WATER.py'), 'validation_source': pin(__file__), 'checks': checks,
              'check_count': len(checks), 'reserved_outcomes_opened': False, 'all_numerical_outcomes_synthetic': True,
              'native_calls': 0, 'local_threads': 1, 'wall_seconds': time.monotonic()-start, 'CPU_seconds': time.process_time()-cpu}
    out = Path(out); out.parent.mkdir(parents=True, exist_ok=True)
    with out.open('x') as f: json.dump(result, f, indent=2); f.write('\n')
    print(json.dumps({'result': pin(out), 'checks': len(checks)}))

if __name__ == '__main__': main(sys.argv[1])
