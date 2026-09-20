"""Thin native EnGrad validation of the frozen nonlinear donor candidates.

No optimization, rescore, numerical gradient or automatic cluster submission.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import re
import shutil
import numpy as np
from affordable_common import HA_TO_KCAL, InvalidArtifact, read_json, record, verify, write_new, xyz
from affordable_response import extract
from hydration_square import endpoint
from mace_site_kinematics import Kinematics

PROTOCOL = 'native_r2scan3c_nonlinear_donor_validation_v1'
BASE_HEADER = '! r2SCAN-3c NoAutostart CPCM(Water) DefGrid3'
HEADER = BASE_HEADER + ' EnGrad'
CASES = ('1H4I', '4MAE', 'PQQSEQ_83440678cbbd658047c9', 'PQQSEQ_07ab500e3df76b30d71c')
PRECISION_A = 1e-12  # Existing physical-map replay policy; never rewrite coordinates.


def recipe(charge, multiplicity=1):
    if multiplicity != 1:
        raise InvalidArtifact('unsupported native spin')
    return HEADER + f'\n* xyzfile {charge} {multiplicity} core.xyz\n'


def check_recipe(path, charge, multiplicity, gradient):
    lines = [s.strip() for s in Path(path).read_text().splitlines()
             if s.strip() and not s.lstrip().startswith('#')]
    wanted = HEADER if gradient else BASE_HEADER
    if ([s for s in lines if s.startswith('!')] != [wanted]
            or [s for s in lines if s.startswith('*')] != [f'* xyzfile {charge} {multiplicity} core.xyz']
            or any(not s.startswith(('!', '*', '%maxcore ')) for s in lines)):
        raise InvalidArtifact('native recipe differs; only EnGrad may be added')


def scf_tolerances(text):
    result = {}
    for name in ('TolE', 'TolG', 'TolErr', 'TolRMSP', 'TolMaxP'):
        values = re.findall(r'\b' + name + r'\s+\.{4}\s+([-+0-9.eE]+)', text)
        if len(values) > 1:
            raise InvalidArtifact('ambiguous SCF tolerance ' + name)
        result[name] = float(values[0]) if values else None
    if any(result[k] is None for k in ('TolE', 'TolG')):
        raise InvalidArtifact('missing actual native SCF tolerances')
    return result


def analytic_drivers(text, metal):
    patterns = {
        'SCF': r'^\s*ORCA SCF GRADIENT CALCULATION\s*$',
        'XC': r'^XC gradient\s+\.{3}\s+done',
        'CPCM': r'^CPCM gradient\s+\.{3}\s+done',
        'dispersion': r'^DISPERSION GRADIENT\s*$',
        'gCP': r'^gCP correction\s+\.{3}\s+done',
    }
    if metal == 'La':
        patterns['ECP'] = r'^ECP gradient\s+\(SHARK\)\s+\.{3}\s+done'
    if re.search(r'numerical gradient|numerical differentiation', text, re.I):
        raise InvalidArtifact('numerical gradient fallback is unsupported')
    found = {k: bool(re.search(v, text, re.M)) for k, v in patterns.items()}
    if not all(found.values()):
        raise InvalidArtifact('native analytic drivers missing: ' + ','.join(k for k, v in found.items() if not v))
    return found


def coordinate_difference(actual, expected):
    if [a[0] for a in actual] != [a[0] for a in expected]:
        raise InvalidArtifact('coordinate atom identity/order differs')
    error = float(np.max(np.abs(np.asarray([a[1:] for a in actual]) - np.asarray([a[1:] for a in expected]))))
    if error > PRECISION_A:
        raise InvalidArtifact('physical coordinates differ beyond existing mapping precision')
    return error


def origins(source, controls_manifest, torsion_manifest):
    """Read actual expanded q0 receipts; absent/failed origins are not zero."""
    design = read_json(source)
    if tuple(c['case_id'] for c in design['cases']) != CASES:
        raise InvalidArtifact('original four-case scope differs')
    control = read_json(controls_manifest); torsion = read_json(torsion_manifest)
    result = []
    for case in design['cases']:
        for metal in ('Ca', 'La'):
            cid = case['case_id']; original = case['origins'][metal]
            if cid in CASES[:2]:
                manifest = controls_manifest
                t = next(t for t in control['tasks'] if t['task_id'] == cid + '__expanded__' + metal)
            else:
                manifest = torsion_manifest
                t = next(t for t in torsion['all_tasks'] if t['point_id'] == cid + '__origin__' + metal)
            if (t['charge'], t['multiplicity']) != (original['charge'], original['multiplicity']):
                raise InvalidArtifact('original native charge/spin differs')
            check_recipe(verify(t['input']), t['charge'], t['multiplicity'], False)
            delta = coordinate_difference(xyz(verify(t['xyz'])), xyz(verify(original['xyz'])))
            row = {'task_id': cid + '__' + metal, 'case_id': cid, 'metal': metal,
                   'source_manifest': record(manifest), 'source_task': t, 'mapping': case['maps'][metal],
                   'source_origin': original, 'executed_vs_source_max_abs_A': delta,
                   'status': 'unavailable', 'result': None, 'SCF_tolerances': None}
            out = Path(t['output_path']); receipt = Path(str(out) + '.execution.json')
            try:
                if not out.exists() or not receipt.exists():
                    raise InvalidArtifact('compatible native origin not complete')
                row.update(result=endpoint(record(out), record(receipt), t['xyz'], t['input']),
                           SCF_tolerances=scf_tolerances(out.read_text()), status='complete')
            except (OSError, ValueError, KeyError) as exc:
                row['reason'] = str(exc)
            result.append(row)
    return result


def pilot_source(collection, execution):
    c = read_json(collection); x = read_json(execution); m = read_json(verify(c['manifest']))
    if c['protocol_id'] != 'bounded_composite_donor_accommodation_v1' or x['manifest'] != c['manifest']:
        raise InvalidArtifact('pilot execution/collection protocol differs')
    wanted = [(cid, z) for cid in CASES for z in ('Ca', 'La')]
    if ([(t['case_id'], t['metal']) for t in m['tasks']] != wanted
            or [(r['case_id'], r['metal']) for r in c['endpoints']] != wanted):
        raise InvalidArtifact('all eight pilot statuses required')
    if 'wall_seconds' not in x or 'error' not in x or 'endpoints' not in x:
        raise InvalidArtifact('terminal pilot execution receipt required')
    recorded = {p['path']: p for p in x['endpoints']}
    for row in c['endpoints']:
        if 'receipt' in row:
            p = row['receipt']
            # Final reporting adds these two descriptive fields; the original
            # endpoint/candidate record still has to match its actual receipt.
            report_only = {'receipt', 'reported_donor_contacts_A', 'reported_curvature_components'}
            if recorded.get(p['path']) != p or {k: v for k, v in row.items() if k not in report_only} != read_json(verify(p)):
                raise InvalidArtifact('pilot endpoint differs from terminal execution')
        elif row.get('candidate') is not None or x['error'] is None:
            raise InvalidArtifact('missing terminal pilot endpoint receipt')
    qualified = [cid for cid in CASES if all(next(r for r in c['endpoints'] if
                 (r['case_id'], r['metal']) == (cid, z)).get('stationary_minimum_qualified', False)
                 for z in ('Ca', 'La'))]
    return c, m, qualified


def candidate_coordinates(task, candidate):
    kin = Kinematics(read_json(verify(task['mapping']))['context'])
    q = np.asarray(candidate['full_q'], dtype=float)
    active = task['active_indices']; inactive = [i for i in range(len(kin.modes)) if i not in active]
    if (q.shape != (len(kin.modes),) or not np.isfinite(q).all()
            or np.any(q[inactive] != 0) or np.max(np.abs(q[active])) > .8 + 1e-12
            or not np.array_equal(q[active], candidate['active_q_radian'])):
        raise InvalidArtifact('candidate has changed physical degrees of freedom')
    if candidate['status'] != 'complete' or not candidate['geometry_checks']['pass']:
        raise InvalidArtifact('candidate scientific/geometry status is unsupported')
    rows = xyz(verify(candidate['coordinate'])); original = xyz(verify(task['xyz']))
    expected = [(a[0], *p) for a, p in zip(original, kin.evaluate(q)[1])]
    delta = coordinate_difference(rows, expected)
    if not kin.check(q)['pass']:
        raise InvalidArtifact('candidate physical mapping is invalid')
    return delta


def stage_candidate(source_task, state, output, ordinal):
    """Keep every runner-owned file inside its own shard, including outputs."""
    candidate = state['candidate']; delta = candidate_coordinates(source_task, candidate)
    td = Path(output).resolve()/f'shard_{ordinal//4}'/'tasks'/source_task['task_id']
    td.mkdir(parents=True, exist_ok=False)
    xp = td/'core.xyz'; shutil.copyfile(verify(candidate['coordinate']), xp)
    ip = td/'endpoint.inp'; ip.write_text(recipe(source_task['charge'], source_task['multiplicity']))
    task = {k: source_task[k] for k in ('task_id', 'case_id', 'metal', 'charge', 'multiplicity', 'mapping', 'active_indices', 'active_mode_ids')}
    task.update(case=task['case_id'], xyz=record(xp), input=record(ip), output_path=str(td/'endpoint.out'),
                engrad_path=str(td/'endpoint.engrad'), candidate=candidate,
                candidate_mapping_roundoff_A=delta, origin=state['origin'], cheap_work=state['cheap_work'])
    return task


def prepare(collection, execution, agreement, controls_manifest, torsion_manifest, output):
    c, pilot, qualified = pilot_source(collection, execution)
    originals = origins(verify(pilot['source']), controls_manifest, torsion_manifest)
    by_id = {r['task_id']: r for r in originals}
    out = Path(output).resolve(); out.mkdir(parents=True, exist_ok=False)
    impl = out/'implementation'; impl.mkdir(); pins = {}
    for p in Path(__file__).parent.glob('*.py'):
        dest = impl/p.name; shutil.copyfile(p, dest); pins[p.name] = record(dest)
    tasks = []; states = []
    for source_task, row in zip(pilot['tasks'], c['endpoints']):
        state = {'task_id': source_task['task_id'], 'case_id': source_task['case_id'], 'metal': source_task['metal'],
                 'cheap_minimum_qualified': row.get('stationary_minimum_qualified', False),
                 'cheap_stationarity': row.get('stationarity'), 'origin': by_id[source_task['task_id']],
                 'status': 'unrun_missing_candidate', 'candidate': row.get('candidate'),
                 'cheap_work': row.get('accommodation_work'), 'source_task': source_task}
        if row.get('candidate') is None:
            state['reason'] = row.get('reason', 'optimizer final candidate unavailable')
        else:
            tasks.append(stage_candidate(source_task, state, out, len(tasks)))
            state['status'] = 'prepared_candidate'
        states.append(state)
    selection = {'protocol_id': PROTOCOL, 'collection': record(collection), 'execution': record(execution),
                 'agreement': record(agreement), 'source': pilot['source'], 'states': states,
                 'qualified_pairs': qualified, 'execution_permitted_by_physical_gate': bool(qualified),
                 'endpoint_denominator': 8, 'prepared_endpoints': len(tasks), 'baseline_changed': False}
    write_new(out/'selection.json', selection)
    shared = {'protocol_id': PROTOCOL, 'selection': record(out/'selection.json'), 'agreement': record(agreement),
              'orca': pilot['orca'], 'implementation': pins,
              'execution_policy': {'task_runner': pins['run_orca_task_manifest.py'], 'runtime_renderer': pins['render_orca_runtime_input.py']},
              'execution_resources': {'mpi_ranks': 16, 'concurrent_tasks': 4}, 'new_DFT_endpoints': len(tasks),
              'allocation_request': {'cpus': 64, 'memory_MiB': 614400},
              'baseline_changed': False, 'execution_permitted_by_physical_gate': bool(qualified)}
    shards = []
    for i in range(0, len(tasks), 4):
        p = out/f'shard_{i//4}'/'manifest.json'
        write_new(p, {**shared, 'tasks': tasks[i:i+4]}); shards.append(record(p))
    write_new(out/'manifest.json', {**shared, 'tasks': tasks, 'shards': shards})
    return dry_run(out/'manifest.json')


def dry_run(manifest):
    from affordable_workflow import dry_run as existing
    m = read_json(manifest); s = read_json(verify(m['selection']))
    if m['protocol_id'] != PROTOCOL or m['execution_permitted_by_physical_gate'] != bool(s['qualified_pairs']):
        raise InvalidArtifact('validation protocol or physical gate differs')
    for p in m['implementation'].values(): verify(p)
    for key in ('collection', 'execution', 'agreement', 'source'): verify(s[key])
    _, pilot, qualified = pilot_source(verify(s['collection']), verify(s['execution']))
    if s['qualified_pairs'] != qualified:
        raise InvalidArtifact('physical execution gate changed')
    eligible = {r['task_id'] for r in s['states'] if r['candidate'] is not None}
    task_ids = [t['task_id'] for t in m['tasks']]
    if len(set(task_ids)) != len(task_ids) or not set(task_ids) <= eligible:
        raise InvalidArtifact('validation task selection changed')
    if 'shards' in m and set(task_ids) != eligible:
        raise InvalidArtifact('combined inventory must retain all available candidates')
    for t in m['tasks']:
        state = next(r for r in s['states'] if r['task_id'] == t['task_id'])
        source_task = next(r for r in pilot['tasks'] if r['task_id'] == t['task_id'])
        if any(t[k] != source_task[k] for k in ('case_id', 'metal', 'charge', 'multiplicity', 'mapping', 'active_indices', 'active_mode_ids')):
            raise InvalidArtifact('candidate chemical state or physical mapping changed')
        if t['candidate'] != state['candidate'] or t['xyz']['sha256'] != state['candidate']['coordinate']['sha256']:
            raise InvalidArtifact('candidate selection or coordinates changed')
        candidate_coordinates(state['source_task'], t['candidate'])
        check_recipe(verify(t['input']), t['charge'], t['multiplicity'], True)
        origin = t['origin']
        if origin != state['origin']:
            raise InvalidArtifact('native origin selection changed')
        if origin['result'] is not None:
            old = origin['result']; st = origin['source_task']
            if endpoint(old['output'], old['receipt'], st['xyz'], st['input']) != old:
                raise InvalidArtifact('archived native origin no longer matches')
    checks = [existing(verify(pin)) for pin in m.get('shards', [])] if 'shards' in m else [existing(manifest)]
    return {'status': 'dry_run_pass', 'manifest': record(manifest), 'tasks': len(m['tasks']),
            'execution_permitted_by_physical_gate': bool(s['qualified_pairs']), 'runner_checks': checks}


def execute(manifest):
    from affordable_workflow import execute as existing
    dry_run(manifest); m = read_json(manifest)
    if 'shards' in m:
        raise InvalidArtifact('execute one of the fixed shard manifests, not the combined inventory')
    if not m['execution_permitted_by_physical_gate']:
        raise InvalidArtifact('no case has two qualified cheap minima; execution disabled')
    return existing(manifest)


def parse_native(task, result):
    text = verify(result['output']).read_text()
    check_recipe(verify(task['input']), task['charge'], task['multiplicity'], True)
    drivers = analytic_drivers(text, task['metal'])
    receipt = read_json(verify(result['receipt'])); gp = verify(receipt['artifacts']['engrad'])
    parsed = extract(gp, verify(result['output']), verify(task['input']), verify(task['xyz']))
    gradient = np.asarray(parsed['gradient_kcal_mol_per_A'])
    return {'engrad': record(gp), 'analytic_drivers': drivers, 'SCF_tolerances': scf_tolerances(text),
            'gradient_kcal_mol_A': gradient.tolist(), 'quantity': 'gradient_not_force',
            **project_gradient(task['mapping'], task['candidate']['full_q'], gradient, task['active_indices']),
            'native_stationarity_qualified': None, 'numerical_gradient_used': False}


def project_gradient(mapping, q, gradient, active):
    kin = Kinematics(read_json(verify(mapping))['context']); gradient = np.asarray(gradient)
    if gradient.shape != kin.core.shape or not np.isfinite(gradient).all():
        raise InvalidArtifact('Cartesian gradient layout differs from candidate mapping')
    if any(kin.modes[i]['unit'] != 'radian' for i in active):
        raise InvalidArtifact('active donor mode is not angular')
    all_modes = np.einsum('mij,ij->m', kin.evaluate(q)[3], gradient)
    return {'active_gradient_kcal_mol_rad': all_modes[active].tolist(),
            'all_mode_derivatives': [{'id': mode['id'], 'coordinate_unit': mode['unit'],
                                     'derivative_kcal_mol_per_unit': float(all_modes[i]), 'active': i in active}
                                    for i, mode in enumerate(kin.modes)]}


def work_comparison(candidate_energy, original_energy, cheap):
    observed = (candidate_energy-original_energy)*HA_TO_KCAL if original_energy is not None else None
    return {'native_work_kcal_mol': observed, 'cheap_work': cheap,
            'model_minus_native_work_kcal_mol': {k: v-observed for k, v in cheap.items() if k in
                ('native_MACE_kcal_mol', 'composite_kcal_mol', 'MACE', 'composite')} if observed is not None and cheap else None}


def collect(manifest, output):
    from compact_solvation import completed
    dry_run(manifest); m = read_json(manifest); s = read_json(verify(m['selection']))
    if 'shards' not in m:
        raise InvalidArtifact('collect the combined inventory to retain all eight states')
    lookup = {t['task_id']: (pin, t) for pin in m['shards'] for t in read_json(verify(pin))['tasks']}
    rows = []
    for state in s['states']:
        row = {'task_id': state['task_id'], 'case_id': state['case_id'], 'metal': state['metal'],
               'status': 'unavailable', 'native': None, 'gradient': None, 'native_work_kcal_mol': None,
               'cheap_minimum_qualified': state['cheap_minimum_qualified'], 'cheap_stationarity': state['cheap_stationarity'],
               'cheap_work': state['cheap_work'], 'origin': state['origin'],
               'cheap_residual_torque_kcal_mol_rad': state['candidate']['gradient_kcal_mol_rad'] if state['candidate'] else None}
        if state['task_id'] not in lookup:
            row.update(status='unrun_missing_candidate', reason=state['reason']); rows.append(row); continue
        pin, t = lookup[state['task_id']]
        try:
            result = completed(verify(pin), t['task_id'])
            if result is None: raise InvalidArtifact('native endpoint not complete')
            native = endpoint(result['output'], result['receipt'], t['xyz'], t['input'])
            gradient = parse_native(t, native)
            old = state['origin']['result']; row.update(status='complete', native=native, gradient=gradient)
            row.update(work_comparison(native['energy_hartree'], old['energy_hartree'] if old else None, state['cheap_work']))
            row['SCF_tolerance_change_from_origin'] = ({k: {'origin': state['origin']['SCF_tolerances'][k], 'candidate': v}
                for k, v in gradient['SCF_tolerances'].items() if state['origin']['SCF_tolerances'][k] != v}
                if state['origin']['SCF_tolerances'] else None)
        except (OSError, ValueError, KeyError) as exc:
            row['reason'] = str(exc)
        rows.append(row)
    cases = []
    for cid in CASES:
        pair = {r['metal']: r for r in rows if r['case_id'] == cid}
        native = (pair['Ca']['native_work_kcal_mol'] - pair['La']['native_work_kcal_mol']
                  if all(r['native_work_kcal_mol'] is not None for r in pair.values()) else None)
        cheap = ({k: pair['Ca']['cheap_work'][k]-pair['La']['cheap_work'][k] for k in pair['Ca']['cheap_work']}
                 if all(r['cheap_work'] is not None for r in pair.values()) else None)
        cases.append({'case_id': cid, 'native_delta_R_kcal_mol': native, 'cheap_delta_R_kcal_mol': cheap,
                      'model_minus_native_delta_R_kcal_mol': {k: v-native for k, v in cheap.items() if k in
                          ('native_MACE_kcal_mol', 'composite_kcal_mol')} if cheap is not None and native is not None else None})
    result = {'protocol_id': PROTOCOL, 'manifest': record(manifest), 'rows': rows, 'cases': cases,
              'completed': sum(r['status'] == 'complete' for r in rows), 'endpoint_denominator': 8,
              'numerical_accuracy_threshold': None, 'biological_accuracy': None, 'entropy': None, 'baseline_changed': False}
    write_new(output, result)
    return {'completed': result['completed'], 'endpoint_denominator': 8, 'output': record(output)}


def main():
    p = argparse.ArgumentParser(description=__doc__); sub = p.add_subparsers(dest='operation', required=True)
    q = sub.add_parser('origins')
    for name in ('source', 'controls-manifest', 'torsion-manifest'): q.add_argument('--'+name, required=True)
    q = sub.add_parser('prepare')
    for name in ('collection', 'execution', 'agreement', 'controls-manifest', 'torsion-manifest', 'output'): q.add_argument('--'+name, required=True)
    for name in ('dry-run', 'execute', 'collect'):
        q = sub.add_parser(name); q.add_argument('--manifest', required=True)
        if name == 'collect': q.add_argument('--output', required=True)
    args = vars(p.parse_args()); op = args.pop('operation').replace('-', '_')
    print(json.dumps(globals()[op](**args), indent=2))


if __name__ == '__main__': main()
