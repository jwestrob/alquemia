"""Admit actual collective protein proposals to the unchanged shared scorer."""
from __future__ import annotations
import argparse
import json
import numpy as np
from affordable_common import InvalidArtifact, read_json, record, verify, write_new, xyz
from mace_hybrid import check_atoms
import nikasha_pool as pool
import scaffold_accommodation as mechanics

TARGET_IDS = {'origin': 'scaffold_origin', 'Ca_adaptive': 'scaffold_Ca', 'La_adaptive': 'scaffold_La'}


def context_coordinates(case, positions, metal):
    """Independent source/cap replay, retaining the archived cap offset."""
    mapping = pool.pinned(case['context_parent_mapping'])
    original = np.array([a[1:] for a in xyz(verify(case['origins'][metal]['xyz']))])
    for item in mapping['source_context_atoms']:
        original[item['context_index']] = positions[item['parent_index']]
    for item in mapping['caps']:
        a = positions[item['retained_parent_index']]
        b = positions[item['omitted_parent_index']]
        distance = np.linalg.norm(b-a)
        if distance < 1e-6: raise InvalidArtifact('collapsed actual cap boundary')
        original[item['context_index']] += a + item['length_A']*(b-a)/distance - np.array(item['q0_xyz_A'])
    return original


def physical_candidate(candidate, tasks, settings):
    if settings != {'maximum_heavy_displacement_A': .8}:
        raise InvalidArtifact('changed collective proposal domain')
    receipt = pool.pinned(candidate['mechanics_receipt'])
    manifest = pool.pinned(receipt['manifest'])
    if receipt['protocol_id'] != mechanics.PROTOCOL or manifest['protocol_id'] != mechanics.PROTOCOL:
        raise InvalidArtifact('unsupported collective parent protocol')
    for pin in manifest['implementation'].values(): verify(pin)
    if manifest['implementation']['scaffold_accommodation.py']['sha256'] != record(mechanics.__file__)['sha256']:
        raise InvalidArtifact('use the recorded mechanics implementation for geometry checks')
    if manifest['settings'] != mechanics.SETTINGS:
        raise InvalidArtifact('collective physical settings differ')
    for key in ('agreement', 'target_mapping_agreement', 'execution_plan', 'technical_addendum', 'inputs'): verify(manifest[key])
    if receipt['status'] != 'complete' or receipt['candidate_admitted'] is not True:
        raise InvalidArtifact('unavailable or rejected actual mechanics candidate')
    if receipt['forcefield_energy_added_to_score'] is not False:
        raise InvalidArtifact('protein force field is a proposal only')
    matches = [t for t in manifest['tasks'] if t['task_id'] == receipt['task_id']]
    if len(matches) != 1: raise InvalidArtifact('missing unique mechanics task')
    task = matches[0]
    cid = tasks['Ca']['case_id']; target = candidate['target']
    if (receipt['case_id'], task['case_id'], receipt['target'], task['target']) != (cid, cid, target, target):
        raise InvalidArtifact('source or target mismatch')
    if target not in TARGET_IDS or candidate['id'] != TARGET_IDS[target]:
        raise InvalidArtifact('undeclared collective target identity')
    if receipt['parent_case'] != task['parent_case'] or receipt['parent_manifest'] != manifest['parent_manifest']:
        raise InvalidArtifact('changed parent provenance')
    case = pool.pinned(task['parent_case']); parents = pool.pinned(manifest['parent_manifest'])
    if [c for c in parents['cases'] if c['case_id'] == cid] != [case]:
        raise InvalidArtifact('case differs from prepared parent manifest')
    for key, parent_key in (('source_atom_ids', 'source_atom_ids'), ('parent_system', 'system')):
        if receipt[key] != case['parent'][parent_key]: raise InvalidArtifact('changed parent atom order or system')
        verify(receipt[key])
    initial_pin = case['targets'][target]['parent_positions_A']
    if receipt['initial_parent_positions_A'] != initial_pin or task['initial_parent_positions_A'] != initial_pin:
        raise InvalidArtifact('changed starting target')
    _, initial, _, geometry = mechanics.target_problem(case, target)
    final = np.array(pool.pinned(receipt['final_parent_positions_A']), dtype=float)
    if final.shape != initial.shape or not np.isfinite(final).all():
        raise InvalidArtifact('invalid final parent coordinates')
    checks = geometry.check(final)
    if not checks['pass'] or checks != receipt['checks']:
        raise InvalidArtifact('collective geometry no longer passes its physical guards')
    energies = [receipt[k] for k in ('initial_parent_energy_kcal_mol', 'final_parent_energy_kcal_mol', 'parent_work_kcal_mol')]
    if not np.isfinite(energies).all() or abs(energies[1]-energies[0]-energies[2]) > 1e-8 or energies[2] > 1e-6:
        raise InvalidArtifact('parent work is nonfinite, inconsistent or uphill')
    for label, positions, energy in (('initial', initial, energies[0]), ('final', final, energies[1])):
        evaluation = receipt[label+'_evaluation']
        if evaluation['status'] != 'complete' or evaluation['energy_kcal_mol'] != energy:
            raise InvalidArtifact('missing actual parent evaluation')
        with np.load(verify(evaluation['coordinates_gradient']), allow_pickle=False) as saved:
            if not np.array_equal(saved['positions_A'], positions): raise InvalidArtifact('parent energy geometry mismatch')
            if not np.isfinite(saved['gradient_kcal_mol_A']).all(): raise InvalidArtifact('nonfinite actual gradient')
    rows = {}
    for metal, endpoint in tasks.items():
        origin = case['origins'][metal]
        if (origin['charge'], origin['multiplicity']) != (endpoint['charge'], endpoint['multiplicity']):
            raise InvalidArtifact('changed endpoint electronic state')
        if case['source_preparation'] != endpoint['source_preparation'] or case['context_maps'][metal] != endpoint['mapping']:
            raise InvalidArtifact('changed endpoint preparation or map')
        if xyz(verify(origin['xyz'])) != xyz(verify(endpoint['xyz'])):
            raise InvalidArtifact('changed source origin')
        actual = xyz(verify(receipt['context_coordinates'][metal]))
        symbols = [a[0] for a in xyz(verify(origin['xyz']))]
        if [a[0] for a in actual] != symbols or not np.allclose([a[1:] for a in actual], context_coordinates(case, final, metal), atol=1e-12, rtol=0):
            raise InvalidArtifact('context source, fixed atoms or cap map changed')
        check_atoms(actual, endpoint['charge']); rows[metal] = actual
    if not pool.same_geometry(rows['Ca'], rows['La']): raise InvalidArtifact('paired collective coordinates differ')
    return rows, {'parent_geometry': checks, 'stationary': receipt['stationary'],
                  'optimizer_status': receipt['optimizer_status'], 'parent_work_kcal_mol': energies[2],
                  'forcefield_energy_added_to_score': False}


def specification(collection, inputs, agreement, output):
    collected = read_json(collection); manifest = pool.pinned(collected['manifest'])
    registry = read_json(inputs)
    if manifest['inputs'] != record(inputs) or manifest['agreement'] != record(agreement):
        raise InvalidArtifact('mechanics collection has different sources/agreement')
    expected = [(c['case_id'], t) for c in registry['cases'] for t in TARGET_IDS]
    if [(r['case_id'], r['target']) for r in collected['rows']] != expected:
        raise InvalidArtifact('incomplete declared target denominator')
    if collected['terminal_searches'] != len(expected): raise InvalidArtifact('mechanics work is not terminal')
    cases = []
    for source in registry['cases']:
        rows = [r for r in collected['rows'] if r['case_id'] == source['case_id']]
        complete = all(r['candidate_admitted'] and r['status'] == 'complete' for r in rows)
        case = {'case_id': source['case_id'], 'status': 'prepared' if complete else 'unavailable',
                'reason': None if complete else 'one_or_more_required_matched_targets_unavailable',
                'mechanics_targets': rows, 'candidates': []}
        if complete:
            case['candidates'] = [{'id': TARGET_IDS[r['target']], 'target': r['target'], 'mechanics_receipt': r['receipt']} for r in rows]
        cases.append(case)
    result = {'branch': 'collective_scaffold', 'inputs': record(inputs), 'agreement': record(agreement),
              'mechanics_collection': record(collection), 'coordinate_limits': {'maximum_heavy_displacement_A': .8},
              'maximum_candidates_per_case': 3, 'cases': cases}
    write_new(output, result)
    return {'specification': record(output), 'prepared': sum(c['status']=='prepared' for c in cases), 'denominator': len(cases)}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('collection', 'inputs', 'agreement', 'output'): parser.add_argument('--'+name, required=True)
    print(json.dumps(specification(**vars(parser.parse_args()))))
