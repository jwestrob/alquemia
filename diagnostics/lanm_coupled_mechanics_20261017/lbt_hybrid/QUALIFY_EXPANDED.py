#!/usr/bin/env python3
"""Finite saved-center LBT hybrid qualification; no optimizer or new engine."""
import argparse
import copy
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import resource
import time
import traceback

PROTOCOL = 'lbt_hybrid_remaining_7CCO_derivatives_v1'
EVALUATOR_SHA = 'a482559ca4febd240a40ec7d6042f4fa24ec4358e4f5ba05770bcfaacea4b784'
CASES = [('7CCO', 'baseline', 'Dy'), ('7CCO', 'expanded', 'La'), ('7CCO', 'expanded', 'Dy')]
FD = {'steps_A': [0.0002, 0.0001], 'absolute_kcal_mol_A': 0.005, 'relative': 0.0001,
      'directions': ['boundary_anchor', 'external_charge'], 'scout_composites': 9}
ANCHORS = {'baseline': ('A/2//C', 'A/2//CA'), 'expanded': ('A/1//C', 'A/1//CA')}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def pin(path):
    p = Path(path).resolve()
    with p.open('rb') as f:
        h = hashlib.file_digest(f, 'sha256').hexdigest()
    return {'path': str(p), 'sha256': h}


def verify(p):
    require(pin(p['path'])['sha256'] == p['sha256'], 'Pin changed: ' + p['path'])
    return Path(p['path'])


def read(p):
    return json.loads(verify(p).read_text())


def write(path, value):
    with Path(path).open('x') as f:
        json.dump(value, f, indent=2, allow_nan=False)
        f.write('\n')
    return pin(path)


def same_bytes(a, b):
    verify(a)
    verify(b)
    return a['sha256'] == b['sha256']


def restore(ev, record):
    """Existing full-point replay plus its already admitted global ALPB cache."""
    p = ev.restore(record)
    sid, region, metal = p['source_id'], p['region'], p['metal']
    raw = read(p['provenance']['component_receipts']['ALPB'])
    for value in raw['artifacts'].values():
        verify(value)
    start = read(raw['artifacts']['START.json'])
    require(start['command'] == [str(verify(ev.cfg['ALPB_binary'])), str(verify(start['xyz'])), str(verify(start['charges']))],
            'Saved ALPB command/runtime differs')
    sy, x = ev.admit.xyz(start['xyz'])
    q = read(ev.tasks[(sid, region, metal)]['global_solvent_charges'])['charges_e']
    require(sy == [a['element'] for a in ev.partitions[(sid, region)][1]['physical_atoms'][:-1]] + [metal]
            and ev.np.array_equal(x, p['positions_A'])
            and [float(v) for v in verify(start['charges']).read_text().split()] == q
            and raw['rounding_bounds_hartree'] == [0., 0.] and abs(raw['charge_e'] + 2) <= 1e-5,
            'Saved ALPB exact coordinates/charges differ')
    grad = ev.np.loadtxt(verify(raw['artifacts']['gradient.dat']))
    value = {'energy_hartree': raw['total_hartree'], 'gradient_hartree_per_bohr': grad.tolist(),
             'components_hartree': raw['components_hartree'], 'record': p['provenance']['component_receipts']['ALPB'],
             'source': 'exact_admitted_parent_cell'}
    key = ev.solvent_key(sid, metal, x, q)
    if key in ev.solvent_cache:
        require(all(ev.solvent_cache[key][k] == value[k] for k in
                    ('energy_hartree', 'gradient_hartree_per_bohr', 'components_hartree')),
                'Conflicting exact saved ALPB cache entry')
    else:
        ev.solvent_cache[key] = value
    return p


def checks(base, ev, center, vector, pairs):
    np = ev.np
    steps = []
    for h, plus, minus in pairs:
        delta = base.difference(plus, minus)
        rows = []
        for k in (*base.COMPONENTS, 'composite_kcal_mol'):
            g = center['gradient_kcal_mol_A'] if k == 'composite_kcal_mol' else center['component_gradients_kcal_mol_A'][k]
            analytic = float(np.sum(np.asarray(g) * vector))
            finite = delta[k] / (2 * h)
            residual = finite - analytic
            tol = FD['absolute_kcal_mol_A'] + FD['relative'] * max(abs(finite), abs(analytic))
            require(all(math.isfinite(v) for v in (analytic, finite, residual, tol)), 'Nonfinite derivative')
            rows.append({'component': k, 'analytic_kcal_mol_A': analytic, 'FD_kcal_mol_A': finite,
                         'residual_kcal_mol_A': residual, 'tolerance_kcal_mol_A': tol, 'pass': abs(residual) <= tol})
        steps.append({'h_A': h, 'plus': plus['provenance']['point'], 'minus': minus['provenance']['point'],
                      'components': rows, 'pass': all(r['pass'] for r in rows)})
    return {'steps': steps, 'pass': all(r['pass'] for r in steps),
            'step_halving_change_kcal_mol_A': {a['component']: b['FD_kcal_mol_A'] - a['FD_kcal_mol_A']
                for a, b in zip(steps[0]['components'], steps[1]['components'])}}


def directions(ev):
    np = ev.np
    old = ev.directions
    require(old['source'] == '7CCO' and old['region'] == 'baseline' and old['metal'] == 'La'
            and old['settings'] == FD and old['frozen_before_calls'], 'Frozen scout directions differ')
    x = ev.origins[('7CCO', 'baseline')]
    result = {}
    for region in ('baseline', 'expanded'):
        part = ev.partitions[('7CCO', region)][1]
        ids = [a['id'] for a in part['physical_atoms']]
        require(ids == old['physical_ids'] and len(set(ids)) == 256
                and np.array_equal(x, ev.origins[('7CCO', region)]), 'Physical source differs between partitions')
        ix = {s: i for i, s in enumerate(ids)}
        cap = part['caps'][0]
        retained, omitted = ANCHORS[region]
        require((cap['retained_source_id'], cap['omitted_source_id']) == (retained, omitted)
                and ix[retained] in part['qm_source_indices'] and ix[omitted] not in part['qm_source_indices'],
                'Declared first cap boundary differs')
        v = x[ix[omitted]] - x[ix[retained]]
        v = v / np.linalg.norm(v)
        axis = np.eye(3)[int(np.argmin(np.abs(v)))]
        v = np.cross(v, axis)
        v = v / np.linalg.norm(v)
        column = np.zeros((256, 3))
        column[ix[omitted]] = v
        anchor = {'physical_id': omitted, 'index': ix[omitted], 'vector': column.tolist()}
        if region == 'baseline':
            require(anchor == old['directions']['boundary_anchor'], 'Baseline boundary direction replay differs')
        else:
            require(ix[ANCHORS['baseline'][1]] in part['qm_source_indices'], 'Expected expanded boundary change absent')
        external = copy.deepcopy(old['directions']['external_charge'])
        require(external['physical_id'] == 'A/1//N' and external['index'] == ix['A/1//N'], 'External direction ID differs')
        field = [a for a in part['field'] if a['id'] == external['physical_id']]
        require(len(field) == 1 and field[0]['source_index'] == external['index']
                and abs(field[0]['charge_e']) > 1e-12 and external['index'] not in part['qm_source_indices'],
                'Same external atom is not valid in this partition')
        items = {'boundary_anchor': anchor, 'external_charge': external}
        for data in items.values():
            vec = np.asarray(data['vector'])
            require(vec.shape == (256, 3) and np.isfinite(vec).all()
                    and np.array_equal(np.flatnonzero(np.any(vec != 0, axis=1)), [data['index']])
                    and abs(float(np.linalg.norm(vec)) - 1) <= 8 * np.finfo(float).eps, 'Direction is not one physical unit displacement')
        result[region] = {'source_id': '7CCO', 'region': region, 'physical_ids': ids,
                          'partition': ev.partitions[('7CCO', region)][0]['partition'],
                          'boundary': cap, 'external_shifted_charge_e': field[0]['charge_e'], 'directions': items}
    return {'protocol': PROTOCOL, 'frozen_before_new_calls': True, 'FD': FD, 'regions': result}


def load_inputs(cfg, out, result):
    require(cfg['protocol'] == PROTOCOL and same_bytes(cfg['implementation'], pin(__file__)), 'Qualification source differs')
    verify(cfg['specification'])
    verify(cfg['plan'])
    require(cfg['evaluator_implementation']['sha256'] == EVALUATOR_SHA, 'Frozen evaluator changed')
    ecfg = read(cfg['evaluator_config'])
    require(same_bytes(ecfg['implementation'], cfg['evaluator_implementation'])
            and same_bytes(ecfg['plan'], cfg['plan']) and ecfg['FD'] == FD, 'Evaluator configuration differs')
    spec = importlib.util.spec_from_file_location('lbt_hybrid_qualification_parent', verify(ecfg['implementation']))
    base = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(base)
    require(base.FD == FD, 'Frozen derivative settings differ')
    scout, origins = read(cfg['scout_result']), read(cfg['origins_result'])
    for parent in (scout, origins):
        require(parent['status'] == 'complete' and parent['protocol'] == base.PROTOCOL
                and same_bytes(parent['implementation'], ecfg['implementation']), 'Parent execution unavailable or changed')
        pcfg = read(parent['config'])
        require({k: v for k, v in pcfg.items() if k != 'scout_result'} ==
                {k: v for k, v in ecfg.items() if k != 'scout_result'}, 'Parent molecular configuration differs')
        require(parent['new_DFT_calls'] == 0 and parent['optimizer_steps'] == 0, 'Unexpected parent molecular scope')
    ocfg = read(origins['config'])
    require(same_bytes(ocfg['scout_result'], cfg['scout_result'])
            and same_bytes(origins['config'], cfg['evaluator_config']), 'Origins are not the configured admitted lineage')
    require(scout['stage'] == 'scout' and scout['qualification_passed']
            and scout['actual_cost']['composites_attempted'] == 9 and len(scout['declared_cells']) == 9,
            'Nine-cell scout gate unavailable')
    require(origins['stage'] == 'origins' and len(origins['rows']) == 8
            and origins['actual_cost']['composites_attempted'] == 7, 'Complete eight-origin stage required')
    ev = base.HybridEvaluator(ecfg, out / 'evaluator')
    # Publish the evaluator before later admission checks, retaining zero-call costs on failure.
    result['_evaluator'] = ev
    require(scout['identity_sha256'] == origins['identity_sha256'] == ev.identity_hash, 'Parent model identity differs')
    saved_directions = read(read(scout['preflight'])['FD_directions'])
    require(saved_directions == ev.directions, 'Actual scout direction receipt differs')
    center = restore(ev, scout['center'])
    require((center['source_id'], center['region'], center['metal']) == ('7CCO', 'baseline', 'La')
            and ev.np.array_equal(center['positions_A'], ev.origins[('7CCO', 'baseline')]), 'Scout center source differs')
    scout_slots = {r['cell']: r for r in scout['declared_cells']}
    require(len(scout_slots) == 9 and all(r['status'] == 'complete' for r in scout_slots.values())
            and same_bytes(scout_slots['center']['point'], scout['center']), 'Scout slots incomplete')
    replay = []
    for name in FD['directions']:
        v = ev.np.asarray(saved_directions['directions'][name]['vector'])
        pairs = []
        for h in FD['steps_A']:
            values = []
            for sign, factor in (('plus', 1), ('minus', -1)):
                slot = scout_slots[name + ':' + format(h, '.17g') + ':' + sign]
                p = restore(ev, slot['point'])
                require((p['source_id'], p['region'], p['metal']) == ('7CCO', 'baseline', 'La')
                        and ev.np.array_equal(p['positions_A'], ev.origins[('7CCO', 'baseline')] + factor * h * v),
                        'Scout displaced source differs')
                values.append(p)
            pairs.append((h, *values))
        check = checks(base, ev, center, v, pairs)
        prior = next(r for r in scout['checks'] if r['direction'] == name)
        require(prior['pass'] and check['pass'] and check['steps'] == prior['steps'], 'Actual scout derivatives do not replay')
        replay.append(dict(check, direction=name))
    centers = {}
    for row in origins['rows']:
        state = (row['source_id'], row['region'], row['metal'])
        require(state in base.STATES and state not in centers and row['status'] == 'complete', 'Origin state missing/duplicated')
        p = restore(ev, row['point'])
        require(state == (p['source_id'], p['region'], p['metal'])
                and ev.np.array_equal(p['positions_A'], ev.origins[state[:2]])
                and all(row[k] == p[k] for k in ('native_energies', 'component_energies_kcal_mol', 'total_energy_kcal_mol')),
                'Origin row/native source differs')
        centers[state] = p
    require(set(centers) == set(base.STATES)
            and same_bytes(centers[('7CCO', 'baseline', 'La')]['provenance']['point'], scout['center']), 'Origin denominator or reused scout differs')
    frozen = directions(ev)
    dp = write(out / 'FD_DIRECTIONS.json', frozen)
    require(ev.cost['composites_attempted'] == ev.cost['MACE_calls_attempted'] == ev.cost['GFN_calls_attempted']
            == ev.cost['MM_configurations_attempted'] == ev.cost['ALPB_calls_attempted'] == 0 and ev.native is None,
            'Preflight made a new molecular call')
    preflight = {'status': 'passed', 'scout_result': cfg['scout_result'], 'origins_result': cfg['origins_result'],
                 'identity_sha256': ev.identity_hash, 'replayed_scout_cells': 9, 'replayed_origin_rows': 8,
                 'scout_checks': replay, 'directions': dp, 'new_molecular_calls': 0,
                 'restored_composite_cache_entries': len(ev.cache),
                 'restored_global_ALPB_cache_entries': len(ev.solvent_cache),
                 'centers': [{'source_id': s, 'region': r, 'metal': m, 'point': centers[(s, r, m)]['provenance']['point']}
                             for s, r, m in CASES]}
    result.update(preflight=write(out / 'PREFLIGHT.json', preflight), directions=dp, identity_sha256=ev.identity_hash,
                  parent_costs={'scout': scout['actual_cost'], 'origins': origins['actual_cost']})
    return base, ev, centers, frozen


def run(base, ev, centers, frozen, result):
    result['checks'] = []
    slots = result['declared_displaced_cells']
    for sid, region, metal in CASES:
        case = {'source_id': sid, 'region': region, 'metal': metal, 'status': 'unavailable',
                'center': centers[(sid, region, metal)]['provenance']['point'], 'directions': []}
        result['checks'].append(case)
        x = ev.origins[(sid, region)]
        for name in FD['directions']:
            data = frozen['regions'][region]['directions'][name]
            v = ev.np.asarray(data['vector'])
            pairs = []
            for h in FD['steps_A']:
                values = []
                for sign, factor in (('plus', 1), ('minus', -1)):
                    slot = next(a for a in slots if (a['source_id'], a['region'], a['metal'], a['direction'], a['h_A'], a['sign'])
                                == (sid, region, metal, name, h, sign))
                    slot['status'] = 'unavailable_attempted'
                    before = copy.deepcopy(ev.cost)
                    request_count = len(ev.requests)
                    try:
                        p = ev.evaluate(sid, region, metal, x + factor * h * v)
                        slot.update(status='complete', point=p['provenance']['point'], cache_key=p['provenance']['cache_key'])
                        values.append(p)
                    except Exception as exc:
                        slot['failure'] = repr(exc)
                        raise
                    finally:
                        slot['cost_delta'] = {k: ev.cost[k] - before.get(k, 0) for k in ev.cost}
                        slot['attempted_point'] = ev.requests[-1] if len(ev.requests) > request_count else None
                pairs.append((h, *values))
            case['directions'].append(dict(checks(base, ev, centers[(sid, region, metal)], v, pairs),
                                           direction=name, physical_id=data['physical_id']))
        case['pass'] = all(r['pass'] for r in case['directions'])
        case['status'] = 'complete' if case['pass'] else 'failed_derivative_gate'
    require(ev.cost['composites_attempted'] == ev.cost['MACE_calls_attempted'] == ev.cost['MM_configurations_attempted'] == 24
            and ev.cost['GFN_calls_attempted'] == 48 and ev.cost['MM_group_queries_attempted'] == 96
            and ev.cost['ALPB_calls_attempted'] <= 24, 'Finite call denominator differs')
    result['qualification_passed'] = all(r['pass'] for r in result['checks'])
    result['status'] = 'complete' if result['qualification_passed'] else 'failed_derivative_gate'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage', choices=('preflight', 'run'), required=True)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    a = parser.parse_args()
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Slurm required')
    out = a.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    tick = time.monotonic()
    usage0 = resource.getrusage(resource.RUSAGE_SELF)
    result = {'status': 'execution_failure', 'stage': a.stage, 'protocol': PROTOCOL,
              'implementation': pin(__file__), 'config': pin(a.config), 'slurm_job_id': os.environ['SLURM_JOB_ID'],
              'qualification_passed': False, 'declared_reused_centers': 3, 'declared_new_composites': 24,
              'declared_displaced_cells': [dict(source_id=s, region=r, metal=m, direction=d, h_A=h, sign=z,
                  status='unavailable_not_attempted') for s, r, m in CASES for d in FD['directions']
                  for h in FD['steps_A'] for z in ('plus', 'minus')],
              'new_DFT_calls': 0, 'optimizer_steps': 0, 'optimization_qualified': False,
              'interpretation': 'Two-direction model derivative diagnostic at three 7CCO states; not electronic-state or global force validation.'}
    try:
        cfg = json.loads(a.config.read_text())
        base, ev, centers, frozen = load_inputs(cfg, out, result)
        if a.stage == 'preflight':
            result['status'] = 'complete'
        else:
            run(base, ev, centers, frozen, result)
    except Exception as exc:
        result.update(failure=repr(exc), traceback=traceback.format_exc())
    ev = result.pop('_evaluator', None)
    if ev is not None:
        result.update(actual_cost=ev.cost, requests=ev.requests, points=ev.points)
        if ev.native is not None:
            result.update(native_MACE_calls=ev.native.calls, native_MACE_device=ev.native.device,
                          native_MACE_versions=ev.native.versions,
                          peak_CUDA_allocated_bytes=ev.native.torch.cuda.max_memory_allocated())
    result['elapsed_seconds'] = time.monotonic() - tick
    usage = resource.getrusage(resource.RUSAGE_SELF)
    result.update(peak_host_RSS_KiB=usage.ru_maxrss,
                  observed_parent_CPU_seconds=usage.ru_utime + usage.ru_stime - usage0.ru_utime - usage0.ru_stime,
                  application_allocated_CPU_seconds=result['elapsed_seconds'] * int(os.environ.get('SLURM_CPUS_PER_TASK', '1')))
    write(out / 'RESULT.json', result)
    print(json.dumps({'status': result['status'], 'result': pin(out / 'RESULT.json')}))
    return 0 if result['status'] == 'complete' else 1


if __name__ == '__main__':
    raise SystemExit(main())
