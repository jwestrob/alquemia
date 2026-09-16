"""Collect and assess the approved global electrostatic physical pilot.

Read-only with respect to scientific inputs and executions. Missing, failed,
and superseded attempts remain visible; no quantum or solver fallback exists.
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np

from affordable_common import InvalidArtifact, HA_TO_KCAL, read_json, record, verify, write_new
from affordable_environment import physical_boundary_key
import affordable_tabi as tabi
from global_electrostatic import (PROTOCOL, CASES, TABI_COULOMB_KCAL_A,
                                 collect_endpoints, endpoint_components, paired_components)

STATE_KEYS = tuple(f'{case}_{metal}' for case in CASES for metal in ('La', 'Ca'))
THRESHOLDS = {'numerical_kcal_mol': .5, 'partition_kcal_mol': 2., 'reduction_kcal_mol': .01}
ERRORS = (ValueError, OSError, KeyError, TypeError)


def expected_tasks():
    tasks = {}
    for state in STATE_KEYS:
        for level in ('primary', 'refined', 'tree'):
            tasks[state + '/' + level] = (state, level, 'identity', 'total')
        tasks[state + '/core_only'] = (state, 'primary', 'identity', 'core')
        if state.endswith('_La'):
            tasks[state + '/environment_only'] = (state, 'primary', 'identity', 'environment')
    for metal in ('La', 'Ca'):
        state = '1h4i_qm33_' + metal
        for transform in ('translated', 'rotated'):
            tasks[state + '/' + transform] = (state, 'primary', transform, 'total')
        tasks[state + '/isolated'] = (state, 'primary', 'identity', 'total')
    tasks['1h4i_qm33_La/repeat'] = ('1h4i_qm33_La', 'primary', 'identity', 'total')
    return tasks


def direct_coulomb(state, transform):
    """All recorded core/environment pairs, including caps; no MM exclusions."""
    core = tabi.transformed(state['core_atoms'], transform)
    environment = tabi.transformed(state['environment_atoms'], transform)
    if not environment:
        return 0.
    coords = np.array([a['xyz_A'] for a in environment])
    charge = np.array([a['charge_e'] for a in environment])
    values = []
    for atom in core:
        distances = np.linalg.norm(coords - atom['xyz_A'], axis=1)
        if np.any(distances < 1.):
            raise InvalidArtifact('core/environment charge overlap in executed transform')
        values.append(float(np.sum(TABI_COULOMB_KCAL_A * atom['charge_e'] * charge / distances)))
    result = math.fsum(values)
    if not math.isfinite(result):
        raise InvalidArtifact('nonfinite direct Coulomb interaction')
    return result


def _check(name, status, **details):
    return {'name': name, 'status': status, **details}


def _limit(name, value, tolerance, **details):
    if value is None:
        return _check(name, 'unavailable', tolerance_kcal_mol=tolerance, **details)
    return _check(name, 'passed' if abs(value) <= tolerance else 'failed',
                  signed_difference_kcal_mol=value, absolute_difference_kcal_mol=abs(value),
                  tolerance_kcal_mol=tolerance, **details)


def _status(checks):
    if any(c['status'] == 'failed' for c in checks):
        return 'failed'
    if not checks or any(c['status'] != 'passed' for c in checks):
        return 'unavailable'
    return 'passed'


def _state_pair_checks(states):
    checks = []
    try:
        for case in CASES:
            la, ca = (states[case + '_' + m] for m in ('La', 'Ca'))
            for key in ('source', 'assembly', 'microstate', 'explicit_waters',
                        'environment_atoms', 'expected_environment_charge_e'):
                if la[key] != ca[key]:
                    raise InvalidArtifact(f'{case}: paired {key} differs')
            if la['core_total_charge_e'] - ca['core_total_charge_e'] != 1:
                raise InvalidArtifact('paired charge difference is not one')
            if len(la['core_atoms']) != len(ca['core_atoms']):
                raise InvalidArtifact('paired core atom count differs')
            for a, b in zip(la['core_atoms'], ca['core_atoms']):
                if any(a[k] != b[k] for k in ('id', 'xyz_A', 'radius_A', 'kind')):
                    raise InvalidArtifact('paired core geometry/mapping differs')
                if a['id'] != 'metal' and a['element'] != b['element']:
                    raise InvalidArtifact('paired nonmetal element differs')
        boundaries = {physical_boundary_key(s['physical_atoms']) for s in states.values()}
        sources = {(s['source']['path'], s['source']['sha256']) for s in states.values()}
        if len(boundaries) != 1 or len(sources) != 1:
            raise InvalidArtifact('partitions do not share physical source/cavity')
        for metal in ('La', 'Ca'):
            totals = [states[c + '_' + metal]['core_total_charge_e'] +
                      states[c + '_' + metal]['expected_environment_charge_e'] for c in CASES]
            if abs(totals[0] - totals[1]) > 1e-6:
                raise InvalidArtifact('partition changes whole-system charge')
        checks.append(_check('state_and_partition_invariants', 'passed', physical_boundary_hash=next(iter(boundaries))))
    except ERRORS as exc:
        checks.append(_check('state_and_partition_invariants', 'failed', reason=str(exc)))
    return checks


def _serialized_coordinates(actual, expected, coordinate_columns):
    """Allow one last-place serialization rounding unit, never atom/charge drift."""
    if actual == expected:
        return 0.
    arows, erows = actual.splitlines(), expected.splitlines()
    if len(arows) != len(erows):
        raise InvalidArtifact('serialized atom count differs')
    maximum = 0.
    for ar, er in zip(arows, erows):
        aa, ee = ar.split(), er.split()
        if len(aa) != len(ee):
            raise InvalidArtifact('serialized atom record differs')
        for i, (a, e) in enumerate(zip(aa, ee)):
            if i not in coordinate_columns:
                if a != e:
                    raise InvalidArtifact('serialized atom/order/charge/radius differs')
            else:
                delta = abs(float(a) - float(e))
                if not math.isfinite(delta) or delta > 1.01e-10:
                    raise InvalidArtifact('serialized coordinates differ beyond one tenth-decimal rounding unit')
                maximum = max(maximum, delta)
    return maximum


def _validate_task_state(task, manifest, base, base_record):
    actual = read_json(verify(manifest['state']))
    isolated = task['task_id'].endswith('/isolated')
    if isolated:
        if actual.get('cavity_role') != 'isolated_reduction' or actual['environment_atoms'] or actual['expected_environment_charge_e'] != 0:
            raise InvalidArtifact('isolated state does not remove environment charges')
        for key in ('source', 'assembly', 'microstate', 'explicit_waters', 'boundary_mapping',
                    'core_atoms', 'core_total_charge_e', 'charge_quality', 'settings', 'endpoint_output'):
            if actual[key] != base[key]:
                raise InvalidArtifact(f'isolated state changed endpoint chemistry: {key}')
        if physical_boundary_key(actual['physical_atoms']) != physical_boundary_key(base['core_atoms']):
            raise InvalidArtifact('isolated cavity is not the unchanged core-plus-cap geometry')
    elif manifest['state'] != base_record:
        raise InvalidArtifact('task state differs from campaign base state')
    if manifest['level'] != task['level'] or manifest['component'] != task['component']:
        raise InvalidArtifact('task numerical level/component differs from campaign')
    if manifest['transform'] != tabi.transform_definition(task['transform']):
        raise InvalidArtifact('task transform differs from campaign')
    expected_settings = {**tabi.PHYSICS, **tabi.NUMERICAL_CONSTANTS, **tabi.LEVELS[task['level']]}
    if manifest['settings'] != expected_settings:
        raise InvalidArtifact('task solver settings differ from frozen numerical level')
    # Every base state was fully validated once by collect(), and its exact
    # immutable record was just matched above. Only isolated derivatives need
    # another complete ownership walk here.
    if isolated:
        tabi.validate_state(actual)
    definition = manifest['transform']
    physical = tabi.transformed(actual['physical_atoms'], definition)
    physical_roundoff = _serialized_coordinates(verify(manifest['physical_xyzr']).read_text(),
                                                tabi.physical_xyzr(physical), {0, 1, 2})
    core = tabi.transformed(actual['core_atoms'], definition)
    environment = tabi.transformed(actual['environment_atoms'], definition)
    charges = core + environment
    if task['component'] != 'total':
        active = {a['id'] for a in (core if task['component'] == 'core' else environment)}
        charges = [dict(a, charge_e=a['charge_e'] if a['id'] in active else 0.) for a in charges]
    charge_roundoff = _serialized_coordinates(verify(manifest['charges']).read_text(),
                                              tabi.pqr_lines(charges), {5, 6, 7})
    return actual, {'physical_xyzr_max_reconstruction_roundoff_A': physical_roundoff,
                    'charges_pqr_max_reconstruction_roundoff_A': charge_roundoff}


def _attempts(manifest_path):
    attempts = []
    chosen = None
    for directory in sorted(manifest_path.parent.glob('attempt_*')):
        receipt = directory / 'execution.json'
        row = {'attempt': str(directory), 'status': 'partial', 'execution_receipt': None}
        if receipt.exists():
            row['execution_receipt'] = record(receipt)
            try:
                value = tabi.collect(manifest_path, receipt)
                row['status'] = 'valid'
                chosen = value
            except ERRORS as exc:
                row.update(status='failed', reason=str(exc))
        else:
            row['reason'] = 'attempt has no completed execution receipt'
        attempts.append(row)
    if chosen is not None:
        for row in attempts:
            row['selected'] = row['execution_receipt'] == chosen['execution_receipt']
    return attempts, chosen


def collect(campaign_path):
    campaign = read_json(campaign_path)
    if campaign['protocol_id'] != PROTOCOL or campaign['thresholds'] != THRESHOLDS:
        raise InvalidArtifact('campaign model or frozen tolerances differ')
    for key in ('states_manifest', 'endpoint_manifest', 'numerics'):
        verify(campaign[key])
    expected = expected_tasks()
    tasks = campaign['tasks']
    if len(tasks) != len(expected) or {t['task_id'] for t in tasks} != set(expected):
        raise InvalidArtifact('campaign does not contain the frozen 25-task schedule')
    for t in tasks:
        if tuple(t[k] for k in ('state_key', 'level', 'transform', 'component')) != expected[t['task_id']]:
            raise InvalidArtifact('campaign task metadata differs from frozen schedule')
    sm = read_json(verify(campaign['states_manifest']))
    if sm['protocol_id'] != PROTOCOL or sm['endpoint_manifest'] != campaign['endpoint_manifest'] or set(sm['states']) != set(STATE_KEYS):
        raise InvalidArtifact('campaign states and quantum manifest disagree')
    quantum = collect_endpoints(verify(campaign['endpoint_manifest']))
    qrows = {row['task_id']: row for row in quantum['rows']}
    states, state_checks, summaries = {}, [], {}
    for name, rec in sm['states'].items():
        try:
            state = read_json(verify(rec)); tabi.validate_state(state)
            q = qrows[name]
            if q['status'] != 'vacuum_endpoint_and_mbis_available':
                raise InvalidArtifact('matching completed gas-phase endpoint unavailable')
            if state['endpoint_output'] != q['output'] or [a['charge_e'] for a in state['core_atoms']] != q['charges']['charge_e']:
                raise InvalidArtifact('state charge distribution differs from gas-phase endpoint')
            quality = read_json(verify(state['charge_quality']['receipt']))
            # Python-version summation can differ in the last binary bit even
            # when all printed atomic charges and source bytes are identical.
            charge_keys = ('charge_e', 'model', 'ecp_convention', 'source_output', 'source_xyz')
            if (any(quality['charges'][k] != q['charges'][k] for k in charge_keys)
                    or abs(quality['charges']['sum_e'] - q['charges']['sum_e']) > 1e-12
                    or quality['status'] != 'passed'):
                raise InvalidArtifact('ESP receipt is not for this gas-phase charge distribution')
            states[name] = state
            summaries[name] = {'state': rec, 'physical_boundary_hash': physical_boundary_key(state['physical_atoms']),
                               'core_total_charge_e': state['core_total_charge_e'],
                               'environment_total_charge_e': state['expected_environment_charge_e'],
                               'core_atom_count': len(state['core_atoms']), 'environment_atom_count': len(state['environment_atoms'])}
            state_checks.append(_check('state/' + name, 'passed'))
        except ERRORS as exc:
            state_checks.append(_check('state/' + name, 'failed', reason=str(exc)))
    if len(states) == len(STATE_KEYS):
        state_checks.extend(_state_pair_checks(states))
    else:
        state_checks.append(_check('state_and_partition_invariants', 'unavailable'))
    rows, direct_cache = [], {}
    for task in tasks:
        row = {**{k: task[k] for k in ('task_id', 'state_key', 'level', 'transform', 'component', 'purpose')},
               'manifest': task['manifest'], 'status': 'unavailable', 'attempts': [], 'solver': None, 'endpoint': None}
        try:
            mp = verify(task['manifest']); tm, _ = tabi._verified_manifest(mp)
            if tm['numerics'] != campaign['numerics']:
                raise InvalidArtifact('task numerical record differs from campaign')
            actual, coordinate_audit = _validate_task_state(task, tm, states[task['state_key']], sm['states'][task['state_key']])
            row['coordinate_serialization_audit'] = coordinate_audit
            row['physical_boundary_hash'] = tm['physical_boundary_hash']
            direct_key = (tm['state']['sha256'], task['transform'])
            if direct_key not in direct_cache:
                direct_cache[direct_key] = direct_coulomb(actual, tm['transform'])
            row['direct_coulomb_kcal_mol'] = direct_cache[direct_key]
            row['attempts'], solver = _attempts(mp)
            if solver is None:
                row['status'] = 'failed' if any(a['status'] == 'failed' for a in row['attempts']) else 'unavailable'
                row['reason'] = 'no successfully collected solver attempt'
            else:
                row.update(status='computed', solver=solver)
                if task['component'] == 'total':
                    row['endpoint'] = endpoint_components(qrows[task['state_key']]['energy_hartree'],
                                                         row['direct_coulomb_kcal_mol'], solver['reaction_field_kJ_mol'])
        except ERRORS as exc:
            row.update(status='failed', reason=str(exc))
        rows.append(row)
    return {'schema_version': 'alquemia.global_electrostatic_collection.v1', 'protocol_id': PROTOCOL,
            'campaign_manifest': record(campaign_path), 'thresholds': campaign['thresholds'],
            'states': summaries, 'state_checks': state_checks, 'quantum': quantum, 'rows': rows,
            'direct_coulomb_constant_kcal_A_mol_e2': TABI_COULOMB_KCAL_A,
            'counts': {s: sum(r['status'] == s for r in rows) for s in ('computed', 'failed', 'unavailable')},
            'collector': record(__file__), 'adapter_collector': record(tabi.__file__),
            'scientific_claim': 'physical feasibility only; no calibrated reference or accuracy trial'}


def assess(collection):
    if collection['protocol_id'] != PROTOCOL or collection['thresholds'] != THRESHOLDS:
        raise InvalidArtifact('collection protocol/tolerances differ')
    if len(collection['rows']) != len(expected_tasks()):
        raise InvalidArtifact('collection does not contain exactly the frozen 25 task rows')
    rows = {r['task_id']: r for r in collection['rows']}
    if len(rows) != 25 or set(rows) != set(expected_tasks()):
        raise InvalidArtifact('collection does not cover frozen task schedule')
    checks = list(collection['state_checks'])
    task_checks = [_check('task/' + name, 'passed' if r['status'] == 'computed' else r['status'])
                   for name, r in rows.items()]
    checks.extend(task_checks)
    pairs, endpoint_changes = {}, {}
    for case in CASES:
        variants = ['primary', 'refined', 'tree'] + (['translated', 'rotated', 'isolated'] if case == CASES[0] else [])
        for variant in variants:
            selected = [rows[case + '_' + m + '/' + variant] for m in ('Ca', 'La')]
            pair = paired_components(*(r['endpoint'] or {'status': 'unavailable'} for r in selected))
            pair['task_ids'] = [r['task_id'] for r in selected]
            pairs[case + '/' + variant] = pair
        for variant in variants:
            if variant in ('primary', 'isolated'):
                continue
            primary, changed = pairs[case + '/primary'], pairs[case + '/' + variant]
            delta = None
            if all(p['status'] == 'computed_uncalibrated_descriptor' for p in (primary, changed)):
                delta = changed['R_global_kcal_mol'] - primary['R_global_kcal_mol']
            checks.append(_limit('paired_numerical/' + case + '/' + variant, delta, THRESHOLDS['numerical_kcal_mol']))
            for metal in ('La', 'Ca'):
                a, b = (rows[case + '_' + metal + '/' + v] for v in ('primary', variant))
                if a['endpoint'] is not None and b['endpoint'] is not None:
                    endpoint_changes[case + '_' + metal + '/' + variant] = {
                        k: b['endpoint'][k] - a['endpoint'][k]
                        for k in ('vacuum_energy_hartree', 'direct_coulomb_kcal_mol', 'reaction_field_kcal_mol', 'total_kcal_mol')}

    # Same physical geometry and meshing resolution must yield identical actual meshes.
    mesh_groups = {}
    for name, row in rows.items():
        if name.endswith('/isolated'):
            mesh_groups.setdefault('isolated_qm33/identity', []).append(row)
            continue
        level = 'primary' if row['level'] == 'tree' else row['level']
        mesh_groups.setdefault(level + '/' + row['transform'], []).append(row)
    for name, group in mesh_groups.items():
        if not all(r['status'] == 'computed' for r in group):
            checks.append(_check('common_mesh/' + name, 'unavailable'))
            continue
        hashes = {r['solver']['mesh_geometry_sha256'] for r in group}
        sources = {r['solver']['physical_xyzr_sha256'] for r in group}
        checks.append(_check('common_mesh/' + name, 'passed' if len(hashes) == len(sources) == 1 else 'failed',
                             mesh_geometry_hashes=sorted(hashes), physical_source_hashes=sorted(sources), task_count=len(group)))

    anchor = rows['1h4i_qm33_La/primary']; repeated = rows['1h4i_qm33_La/repeat']
    repeat_delta = None
    if anchor['solver'] is not None and repeated['solver'] is not None:
        repeat_delta = (repeated['solver']['reaction_field_kJ_mol'] - anchor['solver']['reaction_field_kJ_mol']) / 4.184
    checks.append(_limit('repeat_reaction_field', repeat_delta, THRESHOLDS['reduction_kcal_mol']))

    components = {}
    for state in STATE_KEYS:
        case = state.rsplit('_', 1)[0]
        total, core, env = (rows[x] for x in (state + '/primary', state + '/core_only', case + '_La/environment_only'))
        if any(r['solver'] is None for r in (total, core, env)):
            components[state] = {'status': 'unavailable'}
            checks.append(_check('component_accounting/' + state, 'unavailable'))
            continue
        tr, cr, er = (r['solver']['reaction_field_kJ_mol'] / 4.184 for r in (total, core, env))
        native_cross = (total['solver']['coulomb_kJ_mol'] - core['solver']['coulomb_kJ_mol'] - env['solver']['coulomb_kJ_mol']) / 4.184
        error = native_cross - total['direct_coulomb_kcal_mol']
        components[state] = {'status': 'computed', 'reaction_field_total_kcal_mol': tr,
                             'reaction_field_core_kcal_mol': cr, 'reaction_field_environment_kcal_mol': er,
                             'reaction_field_cross_kcal_mol': tr - cr - er,
                             'native_coulomb_cross_kcal_mol': native_cross,
                             'analytic_direct_coulomb_kcal_mol': total['direct_coulomb_kcal_mol'],
                             'coulomb_cross_discrepancy_kcal_mol': error,
                             'environment_component_source': env['task_id'],
                             'production_expression': 'gas_QM + analytic_direct + RF_total; no RF subtraction'}
        checks.append(_limit('component_accounting/' + state, error, THRESHOLDS['reduction_kcal_mol']))

    reduction = {}
    for metal in ('La', 'Ca'):
        key = '1h4i_qm33_' + metal + '/isolated'; row = rows[key]
        error = None
        if row['endpoint'] is not None:
            e = row['endpoint']
            reduced = e['vacuum_energy_hartree'] * HA_TO_KCAL + e['reaction_field_kJ_mol'] / 4.184
            error = e['total_kcal_mol'] - reduced
            reduction[key] = {'status': 'computed', 'without_environment_kcal_mol': reduced,
                              'full_formula_kcal_mol': e['total_kcal_mol'], 'direct_coulomb_kcal_mol': row['direct_coulomb_kcal_mol'],
                              'reaction_field_kcal_mol': e['reaction_field_kcal_mol'],
                              'interpretation': 'algebraic reduction to gas_QM+isolated_RF; neither term is required to vanish'}
        else:
            reduction[key] = {'status': 'unavailable'}
        checks.append(_limit('isolated_reduction/' + metal, error, THRESHOLDS['reduction_kcal_mol'], check_type='algebra'))

    partitions = {}
    for level in ('primary', 'refined', 'tree'):
        small, large = (pairs[case + '/' + level] for case in CASES)
        difference = None
        if all(p['status'] == 'computed_uncalibrated_descriptor' for p in (small, large)):
            difference = large['R_global_kcal_mol'] - small['R_global_kcal_mol']
        partitions[level] = _limit('partition/' + level, difference, THRESHOLDS['partition_kcal_mol'])
    # Frozen before execution: each resolution must pass, not a favorable one.
    checks.extend(partitions.values())
    gate = _status(checks)
    return {'schema_version': 'alquemia.global_electrostatic_assessment.v1', 'protocol_id': PROTOCOL,
            'campaign_manifest': collection['campaign_manifest'], 'thresholds': THRESHOLDS,
            'gate': gate, 'stage2_eligible': gate == 'passed', 'physical_checks': checks,
            'paired_contrasts': pairs, 'partition_differences': partitions,
            'endpoint_numerical_changes': endpoint_changes, 'reaction_field_components': components,
            'isolated_reduction': reduction, 'task_counts': collection['counts'],
            'task_rows': collection['rows'], 'quantum': collection['quantum'],
            'reference_status': 'unavailable', 'calibrated_decision': None,
            'baseline_default_changed': False, 'accuracy_trial_executed_by_assessor': False,
            'claim': 'Physical/numerical feasibility only; a passed gate does not demonstrate predictive usefulness.',
            'assessor': record(__file__)}


def report(assessment):
    def number(value):
        return 'unavailable' if value is None else f'{value:.6f}'
    lines = ['# Global electrostatic physical pilot', '',
             f"**Physical gate: {assessment['gate']}.**", '',
             f"Solver tasks: {assessment['task_counts']['computed']} computed, "
             f"{assessment['task_counts']['failed']} failed, {assessment['task_counts']['unavailable']} unavailable.", '',
             '| Case / numerical condition | R global (kcal/mol) |', '|---|---:|']
    for key, pair in assessment['paired_contrasts'].items():
        lines.append(f"| {key} | {number(pair['R_global_kcal_mol'])} |")
    lines.extend(['', '| Partition condition | qm36 minus qm33 (kcal/mol) | Status |', '|---|---:|---|'])
    for level, check in assessment['partition_differences'].items():
        lines.append(f"| {level} | {number(check.get('signed_difference_kcal_mol'))} | {check['status']} |")
    lines.extend(['', '## Checks requiring attention', ''])
    attention = [c for c in assessment['physical_checks'] if c['status'] != 'passed']
    if not attention:
        lines.append('All required physical checks passed.')
    unavailable = [c for c in attention if c['status'] == 'unavailable']
    failed_tasks = [c for c in attention if c['status'] == 'failed' and c['name'].startswith('task/')]
    if unavailable:
        lines.append(f'- {len(unavailable)} required checks await valid results; their exact statuses are in the JSON.')
    if failed_tasks:
        lines.append(f'- {len(failed_tasks)} task collections failed; per-attempt reasons are retained in the JSON.')
    for check in (c for c in attention if c['status'] == 'failed' and not c['name'].startswith('task/')):
        detail = check.get('reason', '')
        if 'signed_difference_kcal_mol' in check:
            detail = f"{check['signed_difference_kcal_mol']:.6f} kcal/mol; tolerance {check['tolerance_kcal_mol']:g}"
        lines.append(f"- {check['name']}: {check['status']}{'; ' + detail if detail else ''}.")
    lines.extend(['', 'Full reaction-field, direct Coulomb, endpoint, mesh, and attempt records are retained in the JSON.',
                  'Isolated reduction checks gas QM + isolated reaction energy; it is not a zero-transfer test.',
                  'No old aquo reference, zero threshold, or calibrated classification is applied. Baseline remains default.',
                  'This report establishes no accuracy result and launches no conditional Stage 2 jobs.', ''])
    return '\n'.join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='operation', required=True)
    cmd = sub.add_parser('collect'); cmd.add_argument('--campaign', type=Path, required=True); cmd.add_argument('--output', type=Path, required=True)
    cmd = sub.add_parser('assess'); cmd.add_argument('--collection', type=Path, required=True); cmd.add_argument('--output', type=Path, required=True)
    cmd = sub.add_parser('report'); cmd.add_argument('--assessment', type=Path, required=True); cmd.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.operation == 'collect':
        value = collect(args.campaign); write_new(args.output, value); print(value['counts'])
    elif args.operation == 'assess':
        value = assess(read_json(args.collection)); value['collection'] = record(args.collection)
        write_new(args.output, value); print(value['gate'])
    else:
        text = report(read_json(args.assessment))
        with args.output.open('x') as stream:
            stream.write(text)
        print(args.output)


if __name__ == '__main__':
    main()
