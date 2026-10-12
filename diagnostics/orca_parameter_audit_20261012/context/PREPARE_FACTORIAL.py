"""Replay archived physical donor motions in two exact finite contexts; no energies.

Run under the currently authorized login thread allowance or a Slurm allocation.
Inputs and outputs are explicit. Historical artifacts are read only.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys

for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[name] = '1'
import numpy as np


def pin(path):
    path = Path(path).resolve()
    return {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}


def read(p):
    if isinstance(p, dict):
        assert pin(p['path']) == p, ('pin mismatch', p)
        p = p['path']
    return json.loads(Path(p).read_text())


def xyz(p):
    if isinstance(p, dict):
        assert pin(p['path']) == p, ('coordinate pin mismatch', p)
        p = p['path']
    lines = Path(p).read_text().splitlines()
    rows = [line.split() for line in lines[2:] if line.strip()]
    assert int(lines[0]) == len(rows)
    return [r[0] for r in rows], np.array([[float(x) for x in r[1:4]] for r in rows])


def atom_id(atom):
    return tuple(atom.get(k, '') for k in ('chain', 'resnum', 'insertion_code', 'atom', 'altloc', 'element'))


def mapping_id(row):
    if row['kind'] == 'source':
        return ('source',) + atom_id(row['source'])
    if row['kind'] == 'opaque_cofactor':
        return tuple(row.get(k, '') for k in ('kind', 'fragment', 'name', 'element'))
    assert row['kind'] == 'sigma_link_H', row['kind']
    return ('cap', atom_id(row['retained']), atom_id(row['omitted']))


def source_ids(prep):
    rows = prep['mapping']['source_to_qm']
    result = {mapping_id(a): a['qm_index'] for a in rows}
    assert len(result) == len(rows), 'nonunique atom identity'
    result[('metal',)] = 0
    return result


def context(root, manifest, case_id):
    manifest = root / manifest
    parent = read(manifest)
    case = next(c for c in parent['cases'] if c['case_id'] == case_id)
    union = case['source']['union']
    pp = union['representations']['context']['preparation']
    prep = read(pp)
    candidates = {c['id']: c['xyz'] for c in case['candidates']}
    symbols, origin = xyz(candidates['origin'])
    assert len(origin) == union['atom_count']
    proposals = {}
    for candidate in ('adaptive_Ca', 'adaptive_La'):
        receipt = case['aliases'][candidate]['proposal_receipt']
        result = read(receipt)
        assert result['status'] == 'proposal_available'
        assert result['final_geometry']['physical_feasible']
        assert result['final_geometry']['mapping_checks']['pass']
        search = read(result['manifest'])
        task = next(t for t in search['tasks'] if t['task_id'] == result['task_id'])
        mapping = read(task['mapping'])
        assert task['charge'] == (-3 if task['metal'] == 'Ca' else -2)
        assert task['multiplicity'] == 1
        proposals[candidate] = {'receipt': receipt, 'result': result,
                                'mapping_pin': task['mapping'], 'mapping': mapping,
                                'task': task, 'search_manifest': result['manifest']}
    assert proposals['adaptive_Ca']['mapping']['context'] == proposals['adaptive_La']['mapping']['context']
    return {'case_id': case_id, 'manifest': pin(manifest), 'preparation': pp,
            'prep': prep, 'union': union, 'symbols': symbols, 'origin': origin,
            'ids': source_ids(prep), 'candidates': candidates, 'proposals': proposals,
            'selection': case['source']['selection']}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    root, out = Path(args.root).resolve(), Path(args.output).resolve()
    out.mkdir(parents=True, exist_ok=False)
    sys.path.insert(0, str(root / 'scripts'))
    from mace_site_kinematics import Kinematics
    specs = [('q9z4j7-pqq-la_model', 0, '5eba0c9b39ead1ee', 'slsqp_precision_20260923'),
             ('c5axv8-pqq-la_model', 1, 'aa36ec84c4cbdc62', 'slsqp_precision_expansion_20260923')]
    result = {'protocol_id': 'two_boundary_context_motion_factorial_v1',
              'implementation': pin(__file__), 'kinematics': pin(root / 'scripts/mace_site_kinematics.py'),
              'new_energy_calls': 0, 'new_optimizations': 0,
              'coordinate_copy_tolerance_A': 1e-7, 'physical_bond_tolerance_A': 1e-9,
              'overlap_policy': 'existing geometry_check thresholds: new H-containing <0.55 A, other <1.0 A',
              'cases': [], 'rows': []}
    for protein, sample, envelope, original_run in specs:
        source_case = f'{protein}__conditioned_La__seed-1_sample-{sample}'
        contexts = {
            'original': context(root, f'workspaces/{original_run}/pool_v1/manifest.json', protein),
            'expanded': context(root, 'workspaces/motion_envelope_20260923/pilot34_pool_v1/manifest.json',
                                f'{source_case}__envelope_{envelope}')}
        o, e = contexts['original'], contexts['expanded']
        assert o['prep']['source'] == e['prep']['source'], 'different protonated parent'
        shared = o['ids'].keys() & e['ids'].keys()
        for k in shared:
            assert np.array_equal(o['origin'][o['ids'][k]], e['origin'][e['ids'][k]]), ('origin mismatch', k)
        only_o = o['ids'].keys() - e['ids'].keys()
        assert all(k[0] == 'cap' for k in only_o), 'expanded context lost physical atoms'
        inventory = {'protein': protein, 'source_case': source_case, 'source': o['prep']['source'],
                     'shared_coordinate_identity': True, 'shared_atom_count_including_metal': len(shared),
                     'original_only_atom_ids': sorted(only_o, key=str),
                     'expanded_only_atom_ids': sorted(e['ids'].keys() - o['ids'].keys(), key=str),
                     'contexts': {}}
        for name, data in contexts.items():
            inventory['contexts'][name] = {
                'manifest': data['manifest'], 'preparation': data['preparation'],
                'atom_count': len(data['origin']), 'cap_count': len(data['prep']['mapping']['cut_bonds_and_caps']),
                'formal_charges': {'Ca': -3, 'La': -2}, 'multiplicity': 1,
                'origin': data['candidates']['origin'], 'candidate_pins': data['candidates'],
                'source_id_to_qm': [{'id': k, 'index': v} for k, v in sorted(data['ids'].items(), key=lambda x: str(x[0]))],
                'proposal_receipts': {k: v['receipt'] for k, v in data['proposals'].items()},
                'mapping_pins': {k: v['mapping_pin'] for k, v in data['proposals'].items()},
                'selection': data['selection'], 'added_fragments': data['prep'].get('added_fragments'),
                'charge_ledger': data['prep'].get('source_formal_charge_ledger')}
        for target_name, target in contexts.items():
            kin = Kinematics(target['proposals']['adaptive_Ca']['mapping']['context'])
            target_modes = {m['id']: i for i, m in enumerate(kin.modes)}
            assert len(target_modes) == len(kin.modes)
            assert np.max(np.abs(kin.core - target['origin'])) <= 1e-7
            for generator_name, generator in contexts.items():
                for candidate in ('adaptive_Ca', 'adaptive_La'):
                    entry = generator['proposals'][candidate]
                    source_kin = Kinematics(entry['mapping']['context'])
                    sq = np.array(entry['result']['proposal']['full_q'])
                    ss, sx = xyz(generator['candidates'][candidate])
                    assert np.max(np.abs(source_kin.evaluate(sq)[1] - sx)) <= 1e-7
                    q = np.zeros(len(kin.modes))
                    for mode, angle in zip(source_kin.modes, sq):
                        if angle:
                            assert mode['id'] in target_modes, ('missing mode', mode['id'])
                            q[target_modes[mode['id']]] = angle
                    pp, xx, _, _ = kin.evaluate(q)
                    checks = kin.check(q)
                    assert checks['pass'], checks
                    assert kin.displacement(q) <= .8 + 1e-7
                    assert np.max(np.abs(q)) <= .8 + 1e-12
                    moved = {k for k, i in generator['ids'].items() if np.max(np.abs(sx[i] - generator['origin'][i])) > 1e-8}
                    assert moved <= shared, 'context-exclusive source atom moved'
                    assert not any(k[0] in ('cap', 'opaque_cofactor', 'metal') for k in moved), 'unexpected non-donor motion'
                    shared_error = max(float(np.max(np.abs(xx[target['ids'][k]] - sx[generator['ids'][k]]))) for k in shared)
                    assert shared_error <= 1e-7
                    exclusive = [i for k, i in target['ids'].items() if k not in shared]
                    exclusive_error = float(np.max(np.abs(xx[exclusive] - target['origin'][exclusive]))) if exclusive else 0.
                    assert exclusive_error <= 1e-12, 'exclusive atoms/caps moved'
                    overlaps = []
                    for i in range(len(xx)):
                        for j in range(i):
                            limit = .55 if 'H' in (target['symbols'][i], target['symbols'][j]) else 1.
                            old, new = np.linalg.norm(target['origin'][i]-target['origin'][j]), np.linalg.norm(xx[i]-xx[j])
                            if new < limit and old >= limit:
                                overlaps.append({'i': i, 'j': j, 'old_A': float(old), 'new_A': float(new)})
                    assert not overlaps, overlaps
                    case_dir = out / protein / target_name / f'{generator_name}_{candidate}'
                    case_dir.mkdir(parents=True)
                    endpoints = {}
                    for metal, charge in [('Ca', -3), ('La', -2)]:
                        symbols = target['symbols'].copy(); symbols[0] = metal
                        path = case_dir / f'{metal}.xyz'
                        path.write_text(str(len(xx))+'\nSource-mapped donor motion; no energetic evaluation\n'+''.join(
                            f'{s} {v[0]:.16f} {v[1]:.16f} {v[2]:.16f}\n' for s, v in zip(symbols, xx)))
                        endpoints[metal] = {'xyz': pin(path), 'charge': charge, 'multiplicity': 1}
                    result['rows'].append({'protein': protein, 'target_context': target_name,
                        'geometry_id': f'{generator_name}_{candidate}', 'candidate_source': generator['candidates'][candidate],
                        'source_receipt': entry['receipt'], 'target_mapping': target['proposals']['adaptive_Ca']['mapping_pin'],
                        'physical_q': q.tolist(), 'mode_ids': [m['id'] for m in kin.modes],
                        'geometry_checks': checks, 'shared_max_error_A': shared_error,
                        'exclusive_max_error_A': exclusive_error, 'new_severe_overlaps': overlaps,
                        'moved_shared_ids': sorted(moved, key=str), 'maximum_heavy_displacement_A': kin.displacement(q),
                        'previously_scored_same_context': target_name == generator_name,
                        'endpoints': endpoints})
            result['rows'].append({'protein': protein, 'target_context': target_name, 'geometry_id': 'origin',
                'previously_scored_same_context': True, 'endpoints': target['union']['representations']['context']['endpoints']})
        result['cases'].append(inventory)
    result['status'] = 'prepared_no_energy_calls'
    result['logical_geometry_rows'] = len(result['rows'])
    result['cross_context_geometry_rows'] = sum(not x['previously_scored_same_context'] for x in result['rows'])
    (out / 'PREPARATION.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'status': result['status'], 'rows': len(result['rows']),
                      'cross_context_geometry_rows': result['cross_context_geometry_rows'], 'output': str(out)}))


if __name__ == '__main__':
    main()
