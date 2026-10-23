#!/usr/bin/env python3
"""Synthetic scoring controls; reserved geometry metadata only, no native outputs."""
import argparse
import copy
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
WORK = ROOT / 'workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1'
THREADS = ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')


def pin(path):
    p = Path(path).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def main(output):
    start = time.monotonic()
    cpu = time.process_time()
    assert all(os.environ.get(k) == '1' for k in THREADS)
    score = module(HERE / 'SCORE.py', 'stage_a_synthetic_score')
    body_path = HERE.parent / 'matched_dry_v1/OBJECTIVE.py'
    assert pin(body_path)['sha256'] == 'bf65a20a379811f8d35e2cb5ae78365a5d34022691b332bae5df223196078819'
    body = module(body_path, 'stage_a_fixed_body')
    geometry_path = WORK / 'matched_dry_v1/reserved_validation_v1/RESERVED_PANEL.json'
    training_path = WORK / 'matched_dry_v1/training_v1/RESULT.json'
    assert pin(geometry_path)['sha256'] == 'eed99b0fb6dbe7b49857d68853fe809481f6b657c928a7e75197d77c61f3432b'
    assert pin(training_path)['sha256'] == '99ea27f86871b505dc622c1391a6321381f7430d3ad4a86a72e34ab9a629f6bb'
    panel = json.loads(geometry_path.read_text())
    old = json.loads(training_path.read_text())  # Consumed graph metadata only below.
    rows, model = {}, {}
    for g in panel['geometries']:
        gid = g['geometry_id']
        if gid not in score.CELLS:
            continue
        family = next(f for f in score.FAMILIES if gid.startswith('validation_' + f + '_orientation_'))
        parent = family + '_central_' + g['metal']
        row = {'geometry': copy.deepcopy(g), 'input': copy.deepcopy(old['rows'][parent]['input'])}
        row['bodies'] = body.bodies(row)
        row['reference'] = {'energy': 0.0, 'forces': [[0.0]*3 for _ in range(g['atoms'])]}
        rows[gid] = row
        model[gid] = copy.deepcopy(row['reference'])
    data = {'rows': rows}
    anchors = {k: {m: {'energy': 0.0} for m in score.METALS} for k in ('model', 'reference')}
    checks = {}
    def check(name, value):
        checks[name] = bool(value)
        assert value, name
    def run(d=data, m=model, a=anchors):
        return score.evaluate(d, m, body, a)
    def rejected(name, d=data, m=model, a=anchors):
        try:
            run(d, m, a)
        except (ValueError, KeyError, AssertionError):
            check(name, True)
        else:
            raise AssertionError(name)

    zero = run()
    check('synthetic_zero_pass', zero['validation_passed'] and zero['status'] == 'complete')
    check('fixed_denominators', zero['denominators'] == {'cells':16, 'configurations':8, 'balanced_exchange':8, 'orientation_work':12, 'rigid_response':24, 'total_comparisons':44})
    check('zero_all44', sum(v['passed'] for v in zero['counts'].values()) == 44)
    check('no_fit_objective', not zero['fit_objective_defined'] and 'objective' not in zero)
    check('no_physical_or_B_release', not zero['physical_accuracy_qualified'] and not zero['stage_B_released'])
    shifted = copy.deepcopy(model)
    offset_anchors = copy.deepcopy(anchors)
    for gid, v in shifted.items():
        v['energy'] += 10 if gid.endswith('_La') else 20
    offset_anchors['model']['La']['energy'] = 10
    offset_anchors['model']['Dy']['energy'] = 20
    check('metal_offset_cancellation', run(m=shifted, a=offset_anchors)['validation_passed'])
    for family in score.FAMILIES:
        q = 'validation_' + family + '_orientation_plus'
        m = copy.deepcopy(model)
        m[q + '_La']['energy'] = 2.5
        s = run(m=m)
        check(family + '_contrast_sign', s['balanced_exchange'][q]['error_kcal_mol'] == -2.5)
        check(family + '_La_work_sign', s['orientation_work'][family + '_La']['error_kcal_mol'] == 2.5)
        check(family + '_paired_work_sign', s['orientation_work'][family + '_paired']['error_kcal_mol'] == -2.5)
        check(family + '_unaffected_Dy_work', s['orientation_work'][family + '_Dy']['error_kcal_mol'] == 0)
        check(family + '_failure_not_averaged', not s['validation_passed'] and s['calculated_failure_observed'])
    q = score.CONFIGURATIONS[0]
    m = copy.deepcopy(model)
    m[q + '_Dy']['energy'] = 1.0
    check('inclusive_energy_limit', run(m=m)['validation_passed'])
    a = copy.deepcopy(anchors)
    a['model']['La']['energy'] = 2.0
    s = run(a=a)
    check('anchor_sign_all8', all(v['error_kcal_mol'] == 2.0 for v in s['balanced_exchange'].values()))
    check('anchor_does_not_change_work_or_forces', s['counts']['orientation_work']['passed'] == 12 and s['counts']['rigid_response']['passed'] == 24)
    m = copy.deepcopy(model)
    for metal in score.METALS:
        m[q + '_' + metal]['forces'][0][0] = 12
    s = run(m=m)
    check('individual_shared_force_fails', not s['per_fixture']['La'][q]['passed'] and not s['per_fixture']['Dy'][q]['passed'])
    check('paired_shared_force_cancels', s['per_fixture']['paired'][q]['passed'] and s['per_fixture']['paired'][q]['rigid_RMS_kcal_mol_A'] == 0)
    check('individual_failures_required', not s['validation_passed'] and s['counts']['rigid_response']['failed'] == 2)
    check('atomic_diagnostic_only', all(not v['used_as_pass_criterion'] for v in s['atomic_diagnostics'].values()))
    for side in ('model', 'reference'):
        d, m = copy.deepcopy(data), copy.deepcopy(model)
        gid = q + '_La'
        if side == 'model':
            m[gid] = None
        else:
            d['rows'][gid]['reference'] = {'status':'failed', 'reason':'synthetic unavailable'}
        s = run(d=d, m=m)
        check(side + '_missing_not_pass_or_zero', not s['validation_passed'] and not s['calculated_failure_observed'] and s['status'] == 'partial')
        check(side + '_all44_slots_retained', sum(v['declared'] for v in s['counts'].values()) == 44)
        check(side + '_exact_dependency_counts', [s['counts'][k]['available'] for k in ('balanced_exchange','orientation_work','rigid_response')] == [7,10,22])
        check(side + '_missing_pass_unknown', s['balanced_exchange'][q]['passed'] is None and s['per_fixture']['La'][q]['passed'] is None)
    a = copy.deepcopy(anchors)
    a['reference']['Dy'] = None
    s = run(a=a)
    check('missing_anchor_balanced_unavailable', s['counts']['balanced_exchange']['available'] == 0)
    check('missing_anchor_work_force_retained', s['counts']['orientation_work']['passed'] == 12 and s['counts']['rigid_response']['passed'] == 24)
    check('missing_anchor_overall_incomplete', not s['validation_passed'] and s['status'] == 'partial')
    m = copy.deepcopy(model);m.pop(score.CELLS[0]);rejected('omitted_model_slot_rejected', m=m)
    d = copy.deepcopy(data);d['rows'].pop(score.CELLS[0]);rejected('omitted_geometry_slot_rejected', d=d)
    m = copy.deepcopy(model);m[score.CELLS[0]]['energy'] = float('nan');rejected('nonfinite_energy_rejected', m=m)
    m = copy.deepcopy(model);m[score.CELLS[0]]['forces'][0][0] = float('inf');rejected('nonfinite_force_rejected', m=m)
    m = copy.deepcopy(model);m[score.CELLS[0]]['forces'].pop();rejected('force_shape_rejected', m=m)
    m = copy.deepcopy(model);m[score.CELLS[0]]['status'] = 'running';rejected('ambiguous_status_rejected', m=m)
    d = copy.deepcopy(data);d['rows'][q + '_Dy']['geometry']['positions_A'][0][0] += .01;rejected('paired_geometry_rejected', d=d)
    d = copy.deepcopy(data);d['rows'][q + '_Dy']['geometry']['physical_multiplicity'] = 1;rejected('physical_Dy_label_rejected', d=d)
    d = copy.deepcopy(data);d['rows'][q + '_La']['bodies'].pop();rejected('body_denominator_rejected', d=d)
    a = copy.deepcopy(anchors);a['model'].pop('La');rejected('omitted_anchor_slot_rejected', a=a)
    a = copy.deepcopy(anchors);a['model']['La']['energy'] = True;rejected('boolean_anchor_rejected', a=a)
    check('source_metadata_unchanged', pin(geometry_path)['sha256'] == 'eed99b0fb6dbe7b49857d68853fe809481f6b657c928a7e75197d77c61f3432b'
          and pin(training_path)['sha256'] == '99ea27f86871b505dc622c1391a6321381f7430d3ad4a86a72e34ab9a629f6bb')
    result = {'status':'passed', 'implementation':pin(HERE / 'SCORE.py'), 'specification':pin(HERE / 'SCORE_SPEC.md'),
              'validation_source':pin(__file__), 'geometry_metadata':pin(geometry_path), 'consumed_graph_metadata':pin(training_path),
              'body_helper':pin(body_path), 'checks':checks, 'check_count':len(checks), 'all_energy_force_values_synthetic':True,
              'reserved_energy_force_outputs_opened':False, 'candidate_bound':False, 'native_calls':0,
              'local_numerical_threads':1, 'elapsed_seconds':time.monotonic()-start, 'CPU_seconds':time.process_time()-cpu}
    out = Path(output);out.mkdir(parents=True, exist_ok=False)
    path = out / 'RESULT.json';path.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'result':pin(path), 'checks':len(checks)}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser();parser.add_argument('--output',required=True)
    main(parser.parse_args().output)
