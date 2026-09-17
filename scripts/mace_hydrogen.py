"""Versioned protein-H bond preparation and full-system MACE comparisons."""
import argparse
import copy
import csv
import json
from pathlib import Path
import shutil
import numpy as np
from affordable_common import InvalidArtifact, cache_key, read_json, record, verify, write_new, xyz
from mace_response_trace import TRACE_ID

SCHEMA = 'alquemia.mace_hydrogen.v1'
RECIPE = 'ff19sb_protein_H_radial_bond_projection_v1'


def repair(physical, bonds):
    new = copy.deepcopy(physical)
    byid = {a['id']: i for i, a in enumerate(new)}
    if len(byid) != len(new):
        raise InvalidArtifact('duplicate physical identity')
    expected = {a['id'] for a in new if a['element'] == 'H' and a.get('source_index') is not None}
    moves = []
    seen = set()
    for b in bonds:
        left, right = byid[b['atom_a_id']], byid[b['atom_b_id']]
        pair = [left, right]
        hydrogens = [i for i in pair if new[i]['element'] == 'H']
        if not hydrogens:
            continue
        if len(hydrogens) != 1:
            raise InvalidArtifact('unsupported H-H bond')
        hi = hydrogens[0]
        heavy = right if hi == left else left
        h, x = new[hi], new[heavy]
        if h['id'] in seen or h['id'] not in expected:
            raise InvalidArtifact('multiple/unexpected hydrogen attachment')
        seen.add(h['id'])
        origin = np.array(x['xyz_A'])
        vector = np.array(h['xyz_A']) - origin
        length = float(np.linalg.norm(vector))
        target = float(b['forcefield_equilibrium_A'])
        if not .5 < length < 2. or not .5 < target < 2.:
            raise InvalidArtifact('unsupported hydrogen bond length')
        before = h['xyz_A'][:]
        h['xyz_A'] = (origin + vector * (target / length)).tolist()
        moves.append({'physical_index': hi, 'source_id': h['id'], 'heavy_index': heavy,
                      'heavy_source_id': x['id'], 'before_A': before, 'after_A': h['xyz_A'],
                      'old_bond_length_A': length, 'target_bond_length_A': target})
    if seen != expected:
        raise InvalidArtifact('not every protein hydrogen has exactly one known attachment')
    return new, moves


def load_bonds(audit):
    verify(audit['forcefield'])
    with verify(audit['all_bonds']).open() as f:
        return list(csv.DictReader(f, delimiter='\t'))


def prepare(medium, large, bond_audit, agreement, output):
    from mace_hybrid import dry_run, check_atoms, write_xyz
    audit = read_json(bond_audit)
    original = read_json(verify(audit['source_state']))
    physical, moves = repair(original['physical_atoms'], load_bonds(audit))
    out = Path(output).resolve()
    out.mkdir(parents=True, exist_ok=False)
    prep = {'recipe': RECIPE, 'source_state': audit['source_state'], 'bond_audit': record(bond_audit),
            'forcefield': audit['forcefield'], 'physical_atoms': physical, 'moves': moves,
            'chemical_state_changed': False, 'heavy_coordinates_changed': False,
            'protein_H_changed': len(moves), 'quantum_core_reference': None}
    write_new(out / 'preparation.json', prep)
    xyzs, states = {}, []
    for metal in ('La', 'Ca'):
        coords = [(metal if a['element'] == 'M' else a['element'], *a['xyz_A']) for a in physical]
        p = out / f'1h4i_full_{metal}.xyz'
        write_xyz(p, coords)
        xyzs[metal] = record(p)
        state = {'metal': metal, 'source': original['source'], 'prepared_geometry': record(out/'preparation.json'),
                 'source_geometry_note': 'source protein H coordinates superseded by explicit preparation record',
                 'physical_atoms': physical, 'assembly': original['assembly'],
                 'microstate': original['microstate'], 'explicit_waters': original['explicit_waters']}
        write_new(out / f'{metal}_state.json', state)
        states.append(record(out / f'{metal}_state.json'))
    write_new(out / 'unavailable_DFT.json', {'rows': [], 'status': 'unavailable',
              'reason': 'no matched DFT endpoints exist for the new hydrogen coordinates'})
    write_new(out / 'physical_mapping.json', {'mapping': [{'physical_index':i, 'source_id':a['id']} for i,a in enumerate(physical)],
                                             'core_mapping': None})
    manifests = {}
    for variant, collection in (('medium', medium), ('large', large)):
        source = read_json(collection)
        sm = read_json(verify(source['manifest']))
        dry_run(verify(source['manifest']))
        if source['status'] != 'complete' or sm['schema_version'] != 'alquemia.mace_analytic.v1':
            raise InvalidArtifact('complete analytic parent required')
        if sm['source_states'][0] != audit['source_state']:
            raise InvalidArtifact('bond audit uses another physical state')
        d = out / variant
        d.mkdir()
        impl = d / 'implementation'
        impl.mkdir()
        pins = {}
        for name in (*sm['implementation'], 'mace_response_trace.py', 'mace_hydrogen.py'):
            shutil.copyfile(Path(__file__).with_name(name), impl / name)
            pins[name] = record(impl / name)
        m = copy.deepcopy(sm)
        tasks = []
        for metal in ('La', 'Ca'):
            t = copy.deepcopy(next(t for t in sm['tasks'] if t['kind'] == 'full' and t['metal'] == metal and t['variant'] == 'primary'))
            t['xyz'] = xyzs[metal]
            t['state'] = check_atoms(xyz(verify(t['xyz'])), t['charge'])
            tasks.append(t)
        m.update(schema_version=SCHEMA, protocol_id=f'mace_polar_1{variant[0]}_analytic_vacuum_protein_H_v1',
                 agreement=record(agreement), hydrogen_reference=record(collection),
                 hydrogen_preparation=record(out/'preparation.json'),
                 source_states=states, archived_endpoints=record(out/'unavailable_DFT.json'),
                 atom_mappings=record(out/'physical_mapping.json'), implementation=pins, tasks=tasks,
                 run_inventory={'distinct_mace_energy_force_calls': 2, 'new_DFT_endpoints': 0})
        m['historical_source_endpoint_manifest'] = m.pop('source_endpoint_manifest')
        m['model']['response_trace'] = TRACE_ID
        for t in tasks:
            t.pop('cache_key')
            t['cache_key'] = cache_key({'task': t, 'model': m['model'], 'software': m['software'], 'implementation': pins})
        write_new(d/'manifest.json', m)
        manifests[variant] = dry_run(d/'manifest.json')
    return manifests


def validate(m):
    from mace_hybrid import dry_run
    parent = read_json(verify(m['hydrogen_reference']))
    sm = read_json(verify(parent['manifest']))
    dry_run(verify(parent['manifest']))
    if not parent['analytic_comparison']['core_gate'] or not parent['checks'] or not all(x['pass'] for x in parent['checks']):
        raise InvalidArtifact('parent numerical checks are not valid')
    expected = copy.deepcopy(sm['model'])
    expected['response_trace'] = TRACE_ID
    if m['model'] != expected or m['software'] != sm['software'] or m['tolerances'] != sm['tolerances']:
        raise InvalidArtifact('hydrogen pilot changes electronic model/software/tolerances')
    prep = read_json(verify(m['hydrogen_preparation']))
    audit = read_json(verify(prep['bond_audit']))
    original = read_json(verify(prep['source_state']))
    predicted, moves = repair(original['physical_atoms'], load_bonds(audit))
    if prep['recipe'] != RECIPE or len(prep['physical_atoms']) != len(predicted) or len(prep['moves']) != len(moves):
        raise InvalidArtifact('hydrogen preparation is not the declared geometry operation')
    # Recomputing a three-vector norm on another CPU may change the last bit.
    # Stored inputs remain hash-exact; recipe replay uses the same 1e-12 A
    # tolerance as the declared geometry covariance checks.
    moving = {x['physical_index'] for x in moves}
    for i, (stored, expected_atom) in enumerate(zip(prep['physical_atoms'], predicted)):
        if {k:v for k,v in stored.items() if k!='xyz_A'} != {k:v for k,v in expected_atom.items() if k!='xyz_A'}:
            raise InvalidArtifact('hydrogen preparation changes atom identity')
        if i in moving:
            if not np.allclose(stored['xyz_A'],expected_atom['xyz_A'],atol=1e-12,rtol=0):
                raise InvalidArtifact('hydrogen preparation changes its prescribed position')
        elif stored['xyz_A'] != expected_atom['xyz_A']:
            raise InvalidArtifact('hydrogen preparation changes a fixed atom')
    for stored, expected_move in zip(prep['moves'], moves):
        numeric = ('after_A','old_bond_length_A')
        if {k:v for k,v in stored.items() if k not in numeric} != {k:v for k,v in expected_move.items() if k not in numeric}:
            raise InvalidArtifact('hydrogen move mapping or target differs')
        if any(not np.allclose(stored[k],expected_move[k],atol=1e-12,rtol=0) for k in numeric):
            raise InvalidArtifact('hydrogen move coordinates differ')
    if [t['task_id'] for t in m['tasks']] != [f'1h4i_full_{metal}_primary' for metal in ('La','Ca')]:
        raise InvalidArtifact('hydrogen task inventory changed')
    for t, state_ref in zip(m['tasks'], m['source_states']):
        old = next(x for x in sm['tasks'] if x['task_id'] == t['task_id'])
        if {k:v for k,v in t.items() if k not in ('cache_key','xyz')} != {k:v for k,v in old.items() if k not in ('cache_key','xyz')}:
            raise InvalidArtifact('hydrogen task changes chemical state')
        s = read_json(verify(state_ref))
        if s['physical_atoms'] != prep['physical_atoms'] or s['prepared_geometry'] != m['hydrogen_preparation']:
            raise InvalidArtifact('prepared state/source mapping mismatch')
        expected_rows = [(t['metal'] if a['element'] == 'M' else a['element'], *a['xyz_A']) for a in prep['physical_atoms']]
        rows = xyz(verify(t['xyz']))
        if [tuple(r) for r in rows] != expected_rows:
            raise InvalidArtifact('prepared full coordinates mismatch')
    if read_json(verify(m['archived_endpoints']))['rows']:
        raise InvalidArtifact('unmatched DFT endpoint supplied to hydrogen pilot')
    return {'status': 'pass'}


def collect_hydrogen(manifest):
    from mace_hybrid import accepted_attempt, EV_TO_KCAL
    mp = Path(manifest).resolve()
    m = read_json(mp)
    parent = read_json(verify(m['hydrogen_reference']))
    rows, attempts, changes = {}, [], []
    for t in m['tasks']:
        valid = []
        for a in sorted((mp.parent/'execution'/t['task_id']).glob('attempt_*')):
            r = accepted_attempt(a, t, mp)
            attempts.append({'task_id': t['task_id'], 'path': str(a), 'accepted': r is not None,
                             'receipt': record(a/'receipt.json') if (a/'receipt.json').exists() else None})
            if r is not None:
                valid.append(r)
        row = valid[-1] if valid else {'status': 'unavailable', 'energy_eV': None}
        rows[t['task_id']] = row
        if valid:
            old = parent['rows'][t['task_id']]
            changes.append({'task_id':t['task_id'], 'H_prepared_minus_original_energy_kcal_mol': (row['energy_eV'] - old['energy_eV'])*EV_TO_KCAL,
                            'forces_max_change_eV_A': float(np.max(np.abs(np.load(verify(row['forces']))-np.load(verify(old['forces'])))))})
    done = all(r['status']=='computed' for r in rows.values())
    result = {'status':'complete' if done else 'incomplete', 'protocol_id':m['protocol_id'],
              'manifest': record(mp), 'collection_implementation':record(__file__), 'rows':rows, 'attempts':attempts,
              'reference':None, 'S_kcal_mol':None, 'calibrated_class':None, 'direct':None, 'hybrid':None,
              'geometry_only_pilot':True, 'hybrid_unavailable_reason':'new core geometry lacks matching DFT endpoints',
              'geometry_changes':changes, 'hydrogen_preparation':m['hydrogen_preparation']}
    if done:
        la, ca = (rows[f'1h4i_full_{metal}_primary'] for metal in ('La','Ca'))
        diff = ca['energy_eV'] - la['energy_eV']
        result['direct'] = {'R_eV':diff, 'R_kcal_mol':diff*EV_TO_KCAL,
                            'change_from_original_R_kcal_mol':diff*EV_TO_KCAL-parent['direct']['R_kcal_mol']}
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for key in ('medium-collection','large-collection','bond-audit','agreement','output'):
        p.add_argument('--'+key,required=True)
    a=p.parse_args()
    print(json.dumps(prepare(a.medium_collection,a.large_collection,a.bond_audit,a.agreement,a.output),indent=2))
