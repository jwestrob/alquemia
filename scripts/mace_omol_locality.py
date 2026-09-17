"""Explain saved intact OMOL scores by source residue and fixed radial shells."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import resource
import time

import numpy as np

from affordable_common import InvalidArtifact, read_json, record, verify, write_new, xyz
from mace_file_checks import cached_file_checks
from mace_global_prepare import CASES
from mace_hybrid import accepted_attempt, EV_TO_KCAL
from mace_omol import TOL
from mace_omol_intact import PROTOCOL
from mace_omol_readout import checked_arrays

SHELLS = ('metal', '(0,6]A', '(6,12]A', '(12,18]A', '>18A')


def endpoint(result):
    """Check the selected actual receipt without replaying unrelated campaigns."""
    mp = verify(result['manifest'])
    manifest = read_json(mp)
    for pin in manifest['implementation'].values():
        verify(pin)
    verify(manifest['software'])
    verify(manifest['model']['checkpoint'])
    selected = [t for t in manifest['tasks'] if t['task_id'] == result['task_id']]
    if len(selected) != 1:
        raise InvalidArtifact('missing or duplicate endpoint task')
    task = selected[0]
    if not any(accepted_attempt(a, task, mp) == result for a in
               sorted((mp.parent/'execution'/task['task_id']).glob('attempt_*'))):
        raise InvalidArtifact('saved endpoint differs from actual execution')
    if not result['energy_only'] or not task.get('capture_native_readout'):
        raise InvalidArtifact('energy-only native readout required')
    layers = result['execution_adapter']['layers']
    if len(layers) != 3 or manifest['model']['r_max_A'] != 6.0:
        raise InvalidArtifact('declared three-block/6A locality differs')
    return task, checked_arrays(result, task)


def group_rows(atoms, values, distances, metal_index):
    grouped = {'shell': {}, 'kind': {}, 'residue': {}}
    detail = []
    for i, (atom, value, distance) in enumerate(zip(atoms, values, distances)):
        shell = SHELLS[0] if i == metal_index else SHELLS[
            1 if distance <= 6 else 2 if distance <= 12 else 3 if distance <= 18 else 4]
        residue = atom['id'].rsplit('/', 1)[0] if '/' in atom['id'] else atom['id']
        row = {'atom_index': i, 'id': atom['id'], 'element': atom['element'],
               'kind': atom['kind'], 'residue': residue, 'shell': shell,
               'distance_A': float(distance), 'R_node_kcal_mol': float(value)}
        detail.append(row)
        for key in grouped:
            name = row[key]
            dest = grouped[key].setdefault(name, {'atoms': 0, 'sum_kcal_mol': 0.,
                                                  'sum_abs_kcal_mol': 0., 'max_abs_kcal_mol': 0.})
            dest['atoms'] += 1
            dest['sum_kcal_mol'] += float(value)
            dest['sum_abs_kcal_mol'] += abs(float(value))
            dest['max_abs_kcal_mol'] = max(dest['max_abs_kcal_mol'], abs(float(value)))
    return grouped, detail


@cached_file_checks
def report(source_report, agreement, output):
    start = time.monotonic()
    original = read_json(source_report)
    if (original['protocol_id'] != PROTOCOL or original['status'] != 'complete'
            or not original['numerical_gate_pass'] or set(original['scores']) != set(CASES)):
        raise InvalidArtifact('completed five-case intact report required')
    for field in ('collection', 'physical_preparation', 'agreement', 'report_implementation'):
        verify(original[field])
    scores = {}; checks = []
    for name in CASES:
        source = original['scores'][name]
        preparation = read_json(verify(source['preparation']))
        atoms = preparation['physical_atoms']; index = source['metal_index']
        if index != len(atoms)-1 or atoms[index]['element'] != 'M':
            raise InvalidArtifact('selected metal/source map differs')
        positions = np.array([a['xyz_A'] for a in atoms])
        distances = np.linalg.norm(positions-positions[index], axis=1)
        changes = {}; saved_endpoints = {}; terms = {}
        for metal in ('La', 'Ca'):
            pair = source['endpoints'][metal]; tasks = {}; arrays = {}
            for position in ('bound', 'detached'):
                result = pair[position]
                tasks[position], arrays[position] = endpoint(result)
                task = tasks[position]
                if task['metal'] != metal or task['position'] != position or task['metal_index'] != index:
                    raise InvalidArtifact('endpoint identity/order differs')
                coords = xyz(verify(task['xyz']))
                expected = [(metal if a['element']=='M' else a['element'], *a['xyz_A']) for a in atoms]
                if len(coords) != len(expected) or coords[:index] != expected[:index]:
                    raise InvalidArtifact('nonmetal source coordinates changed')
                if position == 'bound' and coords != expected:
                    raise InvalidArtifact('bound source coordinates changed')
                if position == 'detached':
                    if result['input_state_check']['metal_neighbor_edge_count'] != 0:
                        raise InvalidArtifact('detached metal is not isolated')
                    if min(np.linalg.norm(np.array(coords[index][1:])-positions[:index], axis=1)) <= 18.:
                        raise InvalidArtifact('detached metal is inside declared locality radius')
                saved_endpoints[metal+'_'+position] = {
                    'manifest': result['manifest'], 'task_id': result['task_id'],
                    'readout': result['native_readout'], 'charge': task['charge']}
            if tasks['bound']['charge'] != tasks['detached']['charge']:
                raise InvalidArtifact('charge changed upon detachment')
            expected_charge = preparation['endpoints'][metal]['charge']
            if tasks['bound']['charge'] != expected_charge:
                raise InvalidArtifact('physical source charge changed')
            changes[metal] = {k: arrays['bound'][k]-arrays['detached'][k] for k in arrays['bound']}
            delta = changes[metal]
            terms[metal] = {k: float(v.sum()) for k, v in delta.items()}
            error = (terms[metal]['node_energy_eV']+terms[metal]['embedding_energy_eV']
                     - pair['bound_minus_detached_eV'])*EV_TO_KCAL
            checks.append({'name':name+'_'+metal+'_endpoint_sum', 'error_kcal_mol':float(error),
                           'pass':abs(error)<=TOL['energy_kcal_mol']})
        delta = {k:(changes['Ca'][k]-changes['La'][k])*EV_TO_KCAL for k in changes['Ca']}
        groups, detail = group_rows(atoms, delta['node_energy_eV'], distances, index)
        embedding = float(delta['embedding_energy_eV'].sum())
        total = float(delta['node_energy_eV'].sum())+embedding
        error = total-source['R_coord_kcal_mol']
        checks.append({'name':name+'_score_sum','error_kcal_mol':error,
                       'pass':abs(error)<=TOL['energy_kcal_mol']})
        beyond = distances > 18.
        remote = {metal:{'atoms':int(beyond.sum()),
                         'sum_abs_kcal_mol':float(np.abs(changes[metal]['node_energy_eV'][beyond]).sum()*EV_TO_KCAL),
                         'max_abs_kcal_mol':float(np.max(np.abs(changes[metal]['node_energy_eV'][beyond]),initial=0)*EV_TO_KCAL)}
                  for metal in ('La','Ca')}
        scores[name] = {'source_score_kcal_mol':source['R_coord_kcal_mol'],
                        'summed_score_kcal_mol':total, 'groups':groups, 'per_atom':detail,
                        'endpoint_components_eV':terms, 'endpoint_readouts':saved_endpoints,
                        'preparation':source['preparation'], 'embedding_score_kcal_mol':embedding,
                        'atomic_reference_score_kcal_mol':float(delta['atomic_reference_eV'].sum()),
                        'remote_endpoint_changes':remote}
    contrasts = []
    for c in original['contrasts']:
        positive, negative = (scores[c[k]] for k in ('positive_case','negative_case'))
        value = positive['summed_score_kcal_mol']-negative['summed_score_kcal_mol']
        error = value-c['delta_R_kcal_mol']
        checks.append({'name':c['positive_case']+'_minus_'+c['negative_case']+'_sum',
                       'error_kcal_mol':error, 'pass':abs(error)<=TOL['energy_kcal_mol']})
        contrasts.append({**c,'decomposed_delta_R_kcal_mol':value})
    result = {'status':'complete','source_report':record(source_report),'agreement':record(agreement),
              'implementation':{p.name:record(p) for p in Path(__file__).parent.glob('*.py')},
              'scores':scores,'contrasts':contrasts,'checks':checks,
              'all_checks_pass':all(c['pass'] for c in checks), 'shells':SHELLS,
              'new_energy_evaluations':0, 'score_changed':False, 'unique_physical_cause_established':False,
              'interpretation':'native learned atomic bookkeeping; not observable atomic energies',
              'wall_seconds':time.monotonic()-start,
              'process_CPU_seconds':resource.getrusage(resource.RUSAGE_SELF).ru_utime+resource.getrusage(resource.RUSAGE_SELF).ru_stime}
    out = Path(output).resolve(); out.mkdir(parents=True, exist_ok=False)
    write_new(out/'result.json', result)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('source-report','agreement','output'):
        parser.add_argument('--'+key, required=True)
    result = report(**vars(parser.parse_args()))
    print(json.dumps({'all_checks_pass':result['all_checks_pass'],
                      'shells':{k:v['groups']['shell'] for k,v in result['scores'].items()}}, indent=2))
